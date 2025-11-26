from argparse import ArgumentParser
from collections import defaultdict
import logging
from pathlib import Path
import tomllib

from astroalign import estimate_transform, apply_transform
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.io.fits as pyfits
from astropy import units
from astropy.wcs import WCS
import duckdb
import numpy as np
import pandas as pd
import sep

from .constants import ANGLES, MATRIX
from .frames import Frame, PolFrame, StokesFrame, StokesParameter

type PolFrames = dict[int, PolFrame]
type StokesFrames = dict[StokesParameter, StokesFrame]


LOGLEVELS = {
    "warn": logging.WARNING,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}


debug = None

logger = logging.getLogger(__package__)


def demodulate(frames: PolFrames) -> StokesFrames:
    """Demodulate a set of polarization frames to a set of IQU Stokes frames

    Arguments
    ---------

    - `frames` is dict of `PolFrame`s with the polarizer angle as the
      key

    The list of frames should contain a set of unique angles, which
    should match one of {0, 60, 120} or {0, 45, 90, 135}. Integer
    angles are assumed.

    The function returns a dict of I, Q and U StokesFrame.

    """

    nframes = len(frames)
    try:
        matrix = MATRIX[nframes]
        angles = ANGLES[nframes]
    except KeyError:
        raise ValueError(f"incorrect number of frames to demodulate: {nframes}")

    # Set up the initial Stokes parameter frames
    frame = frames[0]
    if not isinstance(frame, PolFrame):
        raise ValueError("list of frames does not contain polarization frames")
    shape = frame.image.shape

    newframes = {}
    # Perform the matrix - vector multiplication,
    # assigning to the relevant stokes parameters
    parameters = [StokesParameter.I, StokesParameter.Q, StokesParameter.U]
    for row, parameter in zip(matrix, parameters):
        image = np.zeros(shape, dtype=float)  # single "Stokes parameter" image
        for factor, angle in zip(row, angles):
            frame = frames[angle]
            image += frame.image * factor
        # Copy and update the header from the polarization frame to the Stokes frame
        header = frame.header.copy()
        del header["polangle"]
        header["stokes"] = str(parameter)
        newframes[parameter] = StokesFrame(
            image=image, header=header, parameter=parameter
        )

    return newframes


def match_sources(sources1, sources2, maxsep=1):
    """Filter two input source lists to their matching sources"""

    if len({"ra", "dec"} & set(sources1.columns) & set(sources2.columns)) != 2:
        raise ValueError("missing WCS information")

    coords1 = SkyCoord(ra=sources1["ra"], dec=sources1["dec"], unit="degree")
    coords2 = SkyCoord(ra=sources2["ra"], dec=sources2["dec"], unit="degree")
    # Find indices into coords1 that matches with coords2
    logger.info("Matching coordinates")
    idx, d2d, _ = match_coordinates_sky(coords2, coords1)
    sel = d2d < maxsep * units.arcsec
    sources2 = sources2[sel]
    sources1 = sources1.iloc[idx[sel]]

    return sources1, sources2


def filter_by_catalog(
    sources: pd.DataFrame, catalog: str | Path, maxsep: float = 1
) -> pd.DataFrame:
    """Filter out sources in the dataframe that are not found in the catalog

    catalog is a DuckDB file

    maxsep is the maximum matching distance in arcseconds between sources and catalog sources

    """

    minra = sources["ra"].min()
    maxra = sources["ra"].max()
    mindec = sources["dec"].min()
    maxdec = sources["dec"].max()
    with duckdb.connect(catalog, read_only=True) as conn:
        # Provide a simple search for the chip area
        conn.execute(
            """\
select ra, dec from catalog
where ra between ? and ? and "dec" between ? and ?
        """,
            (
                minra,
                maxra,
                mindec,
                maxdec,
            ),
        )
        catsources = pd.DataFrame(conn.fetchall(), columns=["ra", "dec"])

    logger.info("Filtering sources using %s", catalog)
    sources, catsources = match_sources(sources, catsources, maxsep=maxsep)

    return sources


def detect_sources(frame: Frame, good_flags=None) -> pd.DataFrame:
    """Detect sources in a frame

    This function detects sources in a frame using SExtractor. It
    returns a dataframe with 'x' and 'y' columns for the pixel
    positions, and optinally 'ra' and 'dec' columns for WCS positions
    if the header of the input frame contains WCS information

    - good_flags: list of integer flags that are considered acceptable
    flags; only these sources are kept. The default is to keep all
    sources.

    """

    image = frame.image
    bkg = sep.Background(image)
    image = image - bkg

    logger.info("finding sources in image %s", frame.name)
    sources = sep.extract(image, 2, err=bkg.globalrms)
    sources = pd.DataFrame.from_records(sources)
    if good_flags:
        # Keep only sources with appropriate flags
        keep = np.zeros(len(sources), dtype=bool)
        for flag in good_flags:
            keep |= sources["flag"] == flag
        sources = sources[keep]

    # Transform x, y coordinates to RA & Dec if possible
    wcs = WCS(frame.header)
    if wcs.is_celestial:
        # 1-based origin for FITS image
        ra, dec = wcs.all_pix2world(sources[["x", "y"]].to_numpy(), 1).T
        sources["ra"], sources["dec"] = ra, dec

    return sources


def calc_misalignment(
    sources1: PolFrame | pd.DataFrame,
    sources2: PolFrame | pd.DataFrame,
    maxsep: float = 1,
) -> tuple[float, float]:
    """Determine the misalignment from sources1 to sources2

    - catalog: DuckDB file with catalog sources to use

    - maxsep is in arcsec

    Returns the median misalignment values dx and dy
    """

    # If input(s) are frames, obtain the source list
    if isinstance(sources1, PolFrame):
        sources1 = detect_sources(sources1)
    if isinstance(sources2, PolFrame):
        sources2 = detect_sources(sources2)

    if len({"ra", "dec"} & set(sources1.columns) & set(sources2.columns)) != 2:
        raise ValueError("missing WCS information")

    coords1 = SkyCoord(ra=sources1["ra"], dec=sources1["dec"], unit="degree")
    coords2 = SkyCoord(ra=sources2["ra"], dec=sources2["dec"], unit="degree")
    # Find indices into coords1 that matches with coords2
    idx, d2d, _ = match_coordinates_sky(coords2, coords1)
    # Avoid bad matches
    sel = d2d < maxsep * units.arcsec
    # Match the source lists
    sources2 = sources2[sel]
    sources1 = sources1.iloc[idx[sel]]

    dx = np.median(sources2["x"].to_numpy() - sources1["x"].to_numpy())
    dy = np.median(sources2["y"].to_numpy() - sources1["y"].to_numpy())

    return dx, dy


def calc_gradient(frame: StokesFrame, mask=None):
    """Calculate the gradient, but filter on detected & catalog sources

    The frame should be the Stokes I frame

    Only the areas within a boundingbox square around detected sources
    are used, to prevent sky background noise from dominating the
    gradient.

    If catalog is given, detected sources are further filtered on
    sources from this catalog (a database query); this can select
    sources that are deemed unpolarized.

    Set boundingbox to zero or a negative value to use the full image.

    """

    if mask is None:
        dx = np.gradient(frame.image, axis=1)
        dy = np.gradient(frame.image, axis=0)
    else:
        image = np.gradient(frame.image, axis=1)
        dx = image[mask]
        image = np.gradient(frame.image, axis=0)
        dy = image[mask]

    return dx, dy


def calc_mask(
    frame: Frame, catalog: str = "", boundingbox: int = 10, maxsep: float = 1
) -> np.ndarray | None:
    """Return a pixel mask for `frame` based on selected stars and the frame

    Leave catalog empty or set boundingbox to a negative value to not
    create a selection; this will return None instead.

    """

    if catalog and boundingbox > 0:
        sources = detect_sources(frame)
        sources = filter_by_catalog(sources, catalog, maxsep)

        halfwidth = boundingbox / 2

        height, width = frame.image.shape
        yy, xx = np.ogrid[:height, :width]
        mask = np.zeros(frame.image.shape, dtype=bool)
        for _, row in sources.iterrows():
            mask |= (np.abs(xx - row["x"]) <= halfwidth) & (
                np.abs(yy - row["y"]) <= halfwidth
            )
    else:
        mask = None

    return mask


def calc_induced_polarization_from_shift(
    frame1: PolFrame,
    frame2: PolFrame,
    frame_i: StokesFrame,
) -> tuple[StokesFrame, StokesFrame, np.ndarray | None]:
    """Calculate the induced Q and U polarization, given two
    polarization frames with a shift in between, and a Stokes I frame

    Returns a tuple of of (induced) Q and U frames

    """

    hx, hy = calc_misalignment(frame1, frame2)

    dx, dy = calc_gradient(frame_i)
    qx = hx / 3 * dx
    qy = hy / 3 * dy
    ux = hx / np.sqrt(3) * dx
    uy = hy / np.sqrt(3) * dy

    q = qx + qy
    u = ux + uy

    header = frame_i.header.copy()
    header["stokes"] = str(StokesParameter.Q)
    frame_q = StokesFrame(
        header, image=q, name="induced-Q", parameter=StokesParameter.Q
    )
    header["stokes"] = str(StokesParameter.U)
    frame_u = StokesFrame(
        header, image=u, name="induced-U", parameter=StokesParameter.U
    )

    return frame_q, frame_u


def fit_for_induced_polarization(
    frames: StokesFrames, mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Fit for alpha & beta between gradient of Stokes I and Stokes Q or U"""

    iframe = frames[StokesParameter.I]
    qframe = frames[StokesParameter.Q]
    uframe = frames[StokesParameter.U]
    if mask is None:
        dxfull = np.gradient(iframe.image, axis=1)
        dyfull = np.gradient(iframe.image, axis=0)
        dx = dxfull.flatten()
        dy = dyfull.flatten()
        q = qframe.image.flatten()
        u = uframe.image.flatten()
    else:
        dxfull = np.gradient(iframe.image, axis=1)
        dx = dxfull[mask]
        dyfull = np.gradient(iframe.image, axis=0)
        dy = dyfull[mask]
        q = qframe.image[mask]
        u = uframe.image[mask]

    matrix = np.column_stack([dx, dy])
    qpars, *_ = np.linalg.lstsq(matrix, q)
    upars, *_ = np.linalg.lstsq(matrix, u)
    logger.debug("qpars = %s, upars = %s", qpars, upars)

    dx = dxfull
    dy = dyfull

    q_ind = qpars[0] * dx + qpars[1] * dy
    u_ind = upars[0] * dx + upars[1] * dy

    if debug:
        pyfits.PrimaryHDU(data=(dx + dy), header=iframe.header).writeto(
            debug["gradient-file"].replace("<#>", debug["extno"]), overwrite=True
        )
        pyfits.PrimaryHDU(data=q_ind, header=qframe.header).writeto(
            debug["u-induced-file"].replace("<#>", debug["extno"]), overwrite=True
        )
        pyfits.PrimaryHDU(data=u_ind, header=uframe.header).writeto(
            debug["q-induced-file"].replace("<#>", debug["extno"]), overwrite=True
        )

    return q_ind, u_ind


def correct_induced_polarization(
    frames: StokesFrames,
    config: dict,
    polframes: PolFrames | None = None,
    inplace: bool = False,
) -> None | StokesFrames:
    """Correct induced polarization

    `frames` should be a dict of Stokes parameter (keys) and
    corresponding (demodulated) frames.

    `config` is a dict of correction parameters, supplied normally
    through the [correction] section in the configuration file.

    - If `polframes` is given, the induced polarization is calculated
      from the gradient of the Stokes I frame and that of the shift
      between the various `polframes`: the 0 and 60 degree frame, and
      the 0 and 120 degree frame. (At the moment, there is no correction
      implemented for a set of 0/45/90/135 frames.)

      The mask is calculated for the Stokes I frame, but also applied
      to the Stokes Q and U frames.

    If `inplace` is True, the input StokesFrames dict will be
    corrected for its Q and U values. Otherwise a copy is made before
    correction, to be returned.

    The dict of corrected StokesParameters is returned.

    """

    iframe = frames[StokesParameter.I]

    mask = calc_mask(
        iframe, config["catalog"], config["boundingbox"], maxsep=config["maxsep"]
    )

    if polframes is None:
        # Fit directly for alpha & beta,
        # from the gradient-I and Q & U frames
        q_ind, u_ind = fit_for_induced_polarization(frames, mask)

    else:
        # Calculate the shift, the matching alpha & beta, and
        # calculate the induced polarization
        frame_q60, frame_u60 = calc_induced_polarization_from_shift(
            polframes[0], polframes[60], frames[StokesParameter.I]
        )
        frame_q120, frame_u120 = calc_induced_polarization_from_shift(
            polframes[0], polframes[120], frames[StokesParameter.I]
        )
        # To do: verify that these two equations correct
        q_ind = frame_q60.image + frame_q120.image
        u_ind = frame_u60.image - frame_u120.image

    if not inplace:
        # Create a dict of new images, instead of using a deepcopy
        frames = {
            parameter: StokesFrame(
                image=frame.image.copy(),
                header=frame.header.copy(),
                name=frame.name,
                parameter=frame.parameter,
            )
            for parameter, frame in frames.items()
        }

    # Subtract induced polarization
    frames[StokesParameter.Q].image -= q_ind
    frames[StokesParameter.U].image -= u_ind

    return frames


def register_frames(frames: PolFrames, inplace=False) -> None | PolFrames:
    """This function uses AstroAlign (DOI 10.1016/j.ascom.2020.100384)
    functionality (which builds on scikit-image). It may not be fully
    flux conservative; an estimate from the authors is that it may
    introduce flux errors at less than half a percent of the flux in
    the worst cases.

    """

    if not inplace:
        newframes = {
            0: PolFrame(
                image=frames[0].image, header=frames[0].header, name="image-0", angle=0
            )
        }
    else:
        newframes = frames

    # Since we transform all images to the 0-angle one, we can set the WCS
    # of the transformed images to that WCS
    wcs = WCS(frames[0].header)
    # We can let AstroAlign do its full work, but we actually obtain
    # the source lists ourselves to calculate the transformation
    sources0 = detect_sources(frames[0])
    for angle, frame in frames.items():
        if angle == 0:
            continue
        sources = detect_sources(frame)
        target, sources = match_sources(sources0, sources)
        target = target[["x", "y"]].to_numpy()
        source = sources[["x", "y"]].to_numpy()
        logger.info(
            "estimating transformation between %s and %s", frames[0].name, frame.name
        )
        transform = estimate_transform("euclidian", source, target)
        logger.info("transforming %d angle frame to 0 angle frame")
        dest, _ = apply_transform(transform, frame.image, frames[0].image)
        header = frame.header
        header.update(wcs.to_header())
        newframes[angle] = PolFrame(
            image=dest, header=header, name=f"image-{angle}", angle=angle
        )

    return newframes


def run(frames: PolFrames, corr_config: dict) -> StokesFrames:
    """
    Calculate an induced correction and correct for it

    `method` should be one of "induced", "shift-induced", "interpolation"

    """

    method = corr_config["method"].lower()
    if method not in ["induced", "shift-induced", "interpolation", "none"]:
        raise ValueError(
            '`method` should be one of "induced", "shift-induced", "interpolation" or "none"'
        )

    if method == "interpolation":
        register_frames(frames, inplace=True)

    stokes_frames = demodulate(frames)

    polframes = frames if method == "shift-induced" else None
    if method in ("induced", "shift-induced"):
        correct_induced_polarization(
            stokes_frames, corr_config, polframes, inplace=True
        )

    return stokes_frames


def read_files(paths, headerkeys):
    """Read a set of input multi-HDU FITS files into a dict of dict,
    keys being the polarizer angle and the frame number

    """

    shapes = {}
    mosaics = defaultdict(dict)
    for path in paths:
        curangle = None
        with pyfits.open(path) as hdulist:
            for i, hdu in enumerate(hdulist):
                if (
                    isinstance(hdu, (pyfits.PrimaryHDU, pyfits.ImageHDU))
                    and isinstance(hdu.data, np.ndarray)
                    and hdu.data.ndim == 2
                ):
                    anglekey = headerkeys["polangle"]
                    try:
                        angle = hdu.header[anglekey]
                    except KeyError:
                        logging.error(
                            "missing polarization angle keyword '%s'", anglekey
                        )
                        raise KeyError(
                            f"missing polarization angle keyword '{anglekey}'"
                        )

                    if curangle is None:
                        curangle = angle
                    else:
                        if curangle != angle:
                            logging.error(
                                "polarization angle mismatch within FITS file"
                            )
                            raise ValueError(
                                "polarization angle mismatch within FITS file"
                            )
                    # Verify that matching HDUs between mosaics have the same image size
                    shape = shapes.get(i, hdu.data.shape)
                    if hdu.data.shape != shape:
                        logging.error("incorrect shape between mosaic HDUs")
                        raise ValueError("incorrect shape between mosaic HDUs")
                    shapes[i] = shape
                    frame = PolFrame(
                        image=hdu.data,
                        header=hdu.header,
                        name=f"image-{angle}",
                        path=path,
                        angle=angle,
                    )
                    mosaics[angle][i] = frame

    # Verify that each mosaic has the same number of frames
    if len({len(mosaic) for mosaic in mosaics.values()}) != 1:
        logging.error("input files have a non-matching nubmer of HDUs")
        raise ValueError("input files have a non-matching number of HDUs")

    # Invert the nested dict
    allframes = defaultdict(dict)
    for angle, mosaic in mosaics.items():
        for number, frame in mosaic.items():
            # Necessary transformation for images read with Astropy's io.fits module
            if not frame.image.dtype.isnative:
                frame.image = frame.image.view(
                    frame.image.dtype.newbyteorder("=")
                ).byteswap()
            allframes[number][angle] = frame

    return allframes


def save_stokes_frames_as_mosaic(
    frames: list[StokesFrames],
    path: str | Path = "polaris-stokes<stokes>.fits",
    overwrite=False,
):
    """Write a list of StokesFrames dict to 3 FITS files, each
    containing a single mosaic

    basepath is the first part of the output file. Appended to that is
    ext, with <parameter> the Stokes vector.

    """

    hdus = {}
    for stokes_frames in frames:
        for param, frame in stokes_frames.items():
            if param not in hdus:
                hdus[param] = [pyfits.PrimaryHDU(data=None, header=frame.header)]
            hdu = pyfits.ImageHDU(data=frame.image, header=frame.header)
            hdus[param].append(hdu)

    for param, hdus in hdus.items():
        outpath = Path(path.replace("<stokes>", str(param)))
        hdulist = pyfits.HDUList(hdus)
        hdulist.writeto(outpath, overwrite=overwrite)


def read_config():
    parser = ArgumentParser()
    parser.add_argument("config", help="Polaris configuration file")
    args = parser.parse_args()
    with open(args.config, "rb") as file:
        config = tomllib.load(file)
    return config


def setup_logging(level: str = "warning"):
    """Set up some default logging configuration"""

    level = LOGLEVELS[level.lower()]
    fmt = "%(asctime)s  [%(levelname)-5s] - %(module)s.%(funcName)s():%(lineno)d: %(message)s"
    formatter = logging.Formatter(fmt, datefmt="%y-%m-%d %H:%M:%S")
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    handler.setLevel(level)
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)


def main():
    config = read_config()

    if config["debug"]["on"]:  # Save intermediate files for later inspection
        global debug
        debug = config["debug"]

    setup_logging(level=config["logging"]["level"])

    paths = [Path(fname) for fname in config["files"]["input"]]
    nframes = len(paths)
    if nframes not in [3, 4]:  # this is already checked at argparse level
        logging.error("3 or 4 input frames are required")
        raise ValueError("3 or 4 input frames are required")

    allframes = read_files(paths, headerkeys=config["headerkeys"])

    stokes_frames = []  # list of dicts of sets of Stokes frames
    for i, frames in allframes.items():
        if debug:
            debug["extno"] = str(i)
        stokes_frames.append(run(frames, corr_config=config["correction"]))

    save_stokes_frames_as_mosaic(
        stokes_frames,
        path=config["files"]["output"],
        overwrite=config["files"]["overwrite"],
    )


if __name__ == "__main__":
    main()
