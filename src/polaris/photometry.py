import sys

from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
from astropy import units
from astropy.wcs import WCS
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
import sep

from .constants import ANGLES, MATRIX
from .frames import verify_hdus_equal, PolFrameSet, IQUFrameSet, StokesParameter
from .utils import add_logging_level, LOGLEVELS


# Stokes parameters of interest; circular polarization is not measured
STOKES_PARAMETERS = [StokesParameter.I, StokesParameter.Q, StokesParameter.U]
REL_STOKES_PARAMETERS = [StokesParameter.relQ, StokesParameter.relU]


logger = logging.getLogger(__package__)


add_logging_level(LOGLEVELS["notice"], "NOTICE")


def plot_image(image, sources=None):
    fig, ax = plt.subplots(figsize=(6, 12))
    mean, std = np.mean(image), np.std(image)
    ax.imshow(
        image,
        interpolation="nearest",
        cmap="gray",
        vmin=mean - std,
        vmax=mean + std,
        origin="lower",
    )

    if sources is not None:
        for i in range(len(sources)):
            e = Ellipse(
                xy=(sources["x"][i], sources["y"][i]),
                width=12 * sources["a"][i],
                height=12 * sources["b"][i],
                angle=sources["theta"][i] * 180.0 / np.pi,
            )
            e.set_facecolor("none")
            e.set_edgecolor("lightgreen")
            ax.add_artist(e)

    plt.tight_layout()
    plt.show()


def add_wcs(sources, header: pyfits.Header):
    """Add RA & declination information from the (x, y) pixel coordinates

    Note that we ignore the positional errors

    """

    wcs = WCS(header)
    ra, dec = wcs.all_pix2world(
        sources[["x", "y"]].to_numpy(), 1
    ).T  # 1-based origin from FITS image
    sources["ra"] = ra
    sources["dec"] = dec

    return sources


def find_sources(frameset: PolFrameSet | IQUFrameSet, aper=2.5):
    background = {}
    sources_ = {}
    subimage = {}
    for key, frame in frameset.frames.items():
        image = frame.image
        if not image.dtype.isnative:
            image = image.view(image.dtype.newbyteorder("=")).byteswap()
        bkg = sep.Background(image)
        background[key] = bkg.globalrms
        image = image - bkg
        sources = sep.extract(image, 2, err=bkg.globalrms)
        sources = pd.DataFrame.from_records(sources)
        width = (np.median(sources["a"]) + np.median(sources["b"])) / 2
        aper = aper * width
        flux, fluxerr, flag = sep.sum_circle(
            image, sources["x"], sources["y"], aper, err=bkg.globalrms, gain=1.0
        )
        sources["flux"] = flux
        sources["fluxerr"] = fluxerr
        sources["fluxflag"] = flag

        sources = add_wcs(sources, frame.header)

        # Drop sources that are flagged
        sources = sources[(sources["fluxflag"] == 0) & (sources["flag"] == 0)]
        sources_[key] = sources
        subimage[key] = image
        pyfits.PrimaryHDU(data=image, header=frame.header).writeto(
            f"frame{key}.fits", overwrite=True
        )
    return sources_, subimage, background


# to do: verify that there are no double matches, e.g.  a single
# source in one catalogue is matched with 2 sources in another
# catalogue, both within `radius` arcseconds
def match_sources(sources: dict[int, pd.DataFrame], radius=0.1):
    """Match the sources by their (x, y) pixel position.

    Sources not detected are discarded for all relevant frames.
    This does mean that a source with very strong polarization, so that
    it happens to be only visible at one or two polarizer angles, is
    also discarded

    radius is in arcseconds

    """

    radius = radius * units.arcsecond
    angles = list(sources)
    angle0 = angles[0]
    ra, dec = sources[angle0]["ra"], sources[angle0]["dec"]
    cat_ref = SkyCoord(ra=ra, dec=dec, unit="degree")
    # Find the intersection that matches between the first
    # and other source lists
    idrefs = []
    cats = {}
    for angle in angles[1:]:
        ra, dec = sources[angle]["ra"], sources[angle]["dec"]
        cat = SkyCoord(ra=ra, dec=dec, unit="degree")
        idref, d2d, _ = cat.match_to_catalog_sky(cat_ref)
        selection = d2d < radius
        idrefs.append(set(idref[selection]))
        cats[angle] = cat
    # Reduce the first source list to only matches
    # with all the other lists
    idref = list(set.intersection(*idrefs))
    sources[angle0] = sources[angle0].iloc[idref]
    sources[angle0].reset_index(inplace=True)
    cat_ref = cat_ref[idref]

    # Now reduce the other source lists, by matching again
    for angle in angles[1:]:
        cat = cats[angle]
        idx, d2d, _ = cat_ref.match_to_catalog_sky(cat)
        selection = d2d < radius
        sources[angle] = sources[angle].iloc[idx[selection]]
        sources[angle].reset_index(inplace=True)

    n = len(sources[angle0])
    for angle in angles[1:]:
        assert len(sources[angle]) == n

    return sources


def demodulate(sources):
    nframes = len(sources)
    try:
        matrix = MATRIX[nframes]
        angles = ANGLES[nframes]
    except KeyError:
        raise ValueError(f"incorrect number of frames to demodulate: {nframes}")

    dfcomb = pd.concat(list(sources.values()))
    sources0 = dfcomb.groupby(dfcomb.index).mean()
    sources0["flux"] = 0
    sources0["fluxerr"] = 0
    for row, parameter in zip(matrix, STOKES_PARAMETERS):
        # Average across the 3/4 dataframes
        # This is fine for coordinate information (they should be nearly the same);
        sources[parameter] = sources0.copy()
        # the flux is calculated below
        for factor, angle in zip(row, angles):
            sources[parameter]["flux"] += sources[angle]["flux"] * factor
            # weighted standard error
            sources[parameter]["fluxerr"] += sources[angle]["fluxerr"] ** 2 * factor**2
        sources[parameter]["fluxerr"] = np.sqrt(sources[parameter]["fluxerr"])

    return sources


def calc_relative(sources):
    """Calculate the relative Stokes Parameters"""

    iflux = sources[StokesParameter.I]["flux"]
    irelfluxerr = sources[StokesParameter.I]["relfluxerr"]
    for parameter in REL_STOKES_PARAMETERS:
        param = StokesParameter[
            parameter.name[3:]
        ]  # get the 'absolute' Stokes parameter
        sources[parameter] = sources[param].copy()
        sources[parameter]["flux"] = sources[param]["flux"] / iflux
        sources[parameter]["relfluxerr"] = np.sqrt(
            sources[param]["relfluxerr"] ** 2 + irelfluxerr**2
        )
        sources[parameter]["fluxerr"] = (
            sources[parameter]["relfluxerr"] * sources[parameter]["flux"]
        )
    return sources


def main():
    pd.set_option(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.width",
        None,
        "display.max_colwidth",
        None,
    )

    files = sys.argv[1:]

    indices = verify_hdus_equal(files)
    # `sources` is a dict of {HDU-number: {frameset-key: dataframe}}
    all_sources = {}
    for index in indices:
        hdus = []
        for file in files:
            with pyfits.open(file) as hdulist:
                hdus.append(hdulist[index].copy())
        frameset = PolFrameSet.from_hdus(hdus)
        angles = frameset.frames.keys()

        sources, *_ = find_sources(frameset, aper=5)

        sources = match_sources(sources)

        sources = demodulate(sources)

        for parameter in STOKES_PARAMETERS:
            sources[parameter]["relfluxerr"] = (
                sources[parameter]["fluxerr"] / sources[parameter]["flux"].abs()
            )
            sources[parameter]["s2n"] = (
                sources[parameter]["flux"].abs() / sources[parameter]["fluxerr"]
            )

        sources = calc_relative(sources)

        kappa = 10
        s = sources[StokesParameter.Q]
        selQ = s["s2n"] > kappa
        s = sources[StokesParameter.U]
        selU = s["s2n"] > kappa
        if selQ.sum() or selU.sumU():
            sel = selQ | selU
            print()
            print(f"HDU = {index}")
            for angle in angles:
                print(f"Angle {angle}:")
                print(sources[angle].loc[sel, ["x", "y", "ra", "dec", "flux"]])
            for parameter in STOKES_PARAMETERS:
                print(f"Stokes {parameter}:")
                print(
                    sources[parameter].loc[
                        sel, ["x", "y", "ra", "dec", "flux", "fluxerr", "s2n"]
                    ]
                )
        if selQ.sum():
            print("Relative Stokes Q:")
            print(
                sources[StokesParameter.relQ].loc[
                    selQ,
                    ["x", "y", "ra", "dec", "flux", "fluxerr", "relfluxerr", "s2n"],
                ],
            )

        if selU.sum():
            print("Relative Stokes U:")
            print(
                sources[StokesParameter.relU].loc[
                    selU,
                    ["x", "y", "ra", "dec", "flux", "fluxerr", "relfluxerr", "s2n"],
                ],
            )

        all_sources[index] = sources


if __name__ == "__main__":
    main()
