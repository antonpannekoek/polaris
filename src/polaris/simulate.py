"""Simulate a set of images at a given coordinate in the sky

A local GAIA database is used to obtain star positions.


View a (4, 3) array with DS9 as follows:

ds9 -tile grid layout 3 4 'polaris.fits[3]' 'polaris.fits[2]' 'polaris.fits[1]' 'polaris.fits[6]' 'polaris.fits[5]' 'polaris.fits[4]' 'polaris.fits[9]' 'polaris.fits[8]' 'polaris.fits[7]' 'polaris.fits[12]' 'polaris.fits[11]' 'polaris.fits[10]'

"""

from argparse import ArgumentParser
import logging
from pathlib import Path

from astropy.coordinates import SkyCoord
import astropy.io.fits as pyfits
from astropy.table import Table, vstack
from astropy import units
from astropy.wcs import WCS
import duckdb
import numpy as np
from photutils.datasets import make_model_image, apply_poisson_noise
from photutils.psf import MoffatPSF

from .config import read_config
from .utils import add_logging_level, LOGLEVELS


logger = logging.getLogger(__package__)


add_logging_level(LOGLEVELS["notice"], "NOTICE")


MATRIX = {
    # Unnormalized modulation matrix
    # inner lists are the rows
    3: [[1, 1, 0], [1, -0.5, -np.sqrt(3) / 2], [1, -0.5, np.sqrt(3) / 2]],
    # 4: [[1, 1, 1, 1], [2, 0, -2, 0], [0, -2, 0, -2]],
}


def mod_vector(angle: float):
    """Return a module vector for a specific angle"""

    angle = np.deg2rad(angle)
    return [1, np.cos(2 * angle), -np.sin(2 * angle)]


class CCDImage:
    def __init__(self, shape: tuple[int], pixelsize: float = 0.2, ron: float = 0.0):
        self.shape = shape
        self.pixelsize = pixelsize * units.arcsec
        self.ron = ron

    def run(self, table: Table, wcs: WCS, beta: float = 3.5):
        """Simulate a single CCD with a given WCS for stars given in `table`"""

        coords = SkyCoord(table["ra"], table["dec"], unit="degree")
        xcoords, ycoords = coords.to_pixel(wcs)

        table["x_0"] = xcoords
        table["y_0"] = ycoords
        sel = (xcoords > 0) & (xcoords < self.shape[0])
        sel &= (ycoords > 0) & (ycoords < self.shape[1])
        seltable = table[sel]
        psf_model = MoffatPSF(bbox_factor=150.0 / beta)
        image = make_model_image(tuple(self.shape[::-1]), psf_model, seltable)

        return image


def query_catalog(name, coord: SkyCoord, fov: dict[str, units.Quantity]):
    # radius = np.sqrt((fov['width'] / 2)**2 + (fov['height'] / 2)**2)
    conn = duckdb.connect(name, read_only=True)
    ra = coord.ra.degree
    decl = coord.dec.degree
    width = fov["width"].to(units.degree).value
    height = fov["height"].to(units.degree).value
    ra_range = [ra - width / 2, ra + width / 2]
    dec_range = [decl - height / 2, decl + height / 2]

    query = conn.execute(
        """
select ra, "dec", phot_g_mean_mag from gaia where
ra between ? and ? and "dec" between ? and ?
    """,
        (*ra_range, *dec_range),
    )
    df = query.df()
    df.rename(columns={"phot_g_mean_mag": "mag"}, inplace=True)
    catalog = Table.from_pandas(df)
    conn.close()

    return catalog


def create_wcs(center: SkyCoord, ccd: CCDImage):
    ra = center.ra.degree
    decl = center.dec.degree
    nx, ny = ccd.shape

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [nx / 2, ny / 2]
    wcs.wcs.cdelt = np.array([-ccd.pixelsize.value, ccd.pixelsize.value]) / 3600
    wcs.wcs.crval = [ra, decl]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = ["deg", "deg"]

    return wcs


class StarSimulation:
    def __init__(
        self,
        magrange: tuple[float],
        psfparams: dict[str, float],
        catalog: str | None = None,
    ):
        """

        - magminmax: magnitude range, [minimum, maximum]
             Note that the lowest value is always interpreted as the minimum magnitude

        - fwhm: full width half maximum in pixels

        - fwhm_std: spread in fwhm, in pixels

        - beta: "wings" parameter of the Moffat PSF function
            Alpha is determined from the fwhm

        - catalog: catalog name of the DuckDB database file

        """

        self.magrange = magrange
        self.fwhm = psfparams.get("fwhm", 2.5)
        self.fwhm_std = psfparams.get("fwhm_std", 0.0)
        self.beta = psfparams.get("beta", 3.5)
        self.catalog = catalog

    def copy(self):
        return StarSimulation(
            magrange=self.magrange[:],
            psfparams={"fwhm": self.fwhm, "fwhm": self.fwhm_std, "beta": self.beta},
            catalog=self.catalog,
        )

    def run(
        self,
        center: SkyCoord,
        fov: dict[str, units.Quantity],
        t_exp: float,
        nstars: int = 1000,
        zeropoint: float = 30.0,
        rng=None,
    ):
        """Create a list of random stars, and optionally add stars from an existing catalog to it"""

        rng = np.random.default_rng(rng)

        stars = Table()
        low = (center.ra - fov["width"] / 2).value
        high = (center.ra + fov["width"] / 2).value
        stars["ra"] = rng.uniform(low, high, nstars) * units.degree
        low = (center.dec - fov["height"] / 2).value
        high = (center.dec + fov["height"] / 2).value
        stars["dec"] = rng.uniform(low, high, nstars) * units.degree
        unirand = np.random.uniform(0, 1, nstars)
        mmin, mmax = sorted(self.magrange)
        stars["mag"] = (1 / 0.6) * np.log10(
            unirand * (10 ** (0.6 * mmax) - 10 ** (0.6 * mmin)) + 10 ** (0.6 * mmin)
        )
        if self.catalog:
            catalog = query_catalog(self.catalog, center, fov)
            stars = vstack([stars, catalog])

        stars["flux"] = t_exp * 10.0 ** (-0.4 * (stars["mag"] - zeropoint))
        stars["fwhm"] = rng.normal(self.fwhm, self.fwhm_std, len(stars))
        return stars


class Mosaic:
    def __init__(
        self, ccd: CCDImage, layout: tuple[int] = (1, 1), gaps: tuple[float] = (0, 0)
    ):
        self.layout = layout
        self.gaps = gaps
        self.ccd = ccd

        self.fov = {
            "width": ccd.pixelsize
            * (ccd.shape[0] * layout[0] + gaps[0] * (layout[0] - 1)),
            "height": ccd.pixelsize
            * (ccd.shape[1] * layout[1] + gaps[1] * (layout[1] - 1)),
        }

    def simulate(
        self,
        center: SkyCoord,
        starsim: StarSimulation,
        t_exp: float,
        stars: int | Table = 1000,
        skylevel: float = 1,
        zeropoint: float = 30.0,
        outfile: str | Path = "polaris.fits",
        angle: float = None,
        polstars=None,
        poisson_noise: bool = True,
        rng=None,
    ):
        if isinstance(stars, int):
            table = starsim.run(
                center, self.fov, t_exp, nstars, zeropoint=zeropoint, rng=rng
            )
        else:
            table = stars

        centers = self._calc_centers(center)

        cards = [
            ("magzp", zeropoint, ""),
            ("beta", starsim.beta, ""),
            ("bkg", skylevel, ""),
        ]
        if angle is not None:
            cards.append(("polangle", angle, "polarizer angle"))

        hdus = [pyfits.PrimaryHDU()]
        hdus[0].header.extend(cards)
        for key, center in centers.items():
            logger.notice("Creating single CCD image")
            wcs = create_wcs(center, self.ccd)
            image = self.ccd.run(table, wcs, starsim.beta)
            if polstars:
                # Add polarised stars
                image = image + self.ccd.run(polstars, wcs, starsim.beta)
            image = image + skylevel * t_exp
            if poisson_noise:
                image = apply_poisson_noise(image, seed=rng)
            # Add read-out noise
            image = image + rng.normal(0.0, self.ccd.ron, size=image.shape)

            hdu = pyfits.ImageHDU(image.astype(np.float64), header=wcs.to_header())
            hdu.header.extend(cards)
            hdus.append(hdu)
        logger.notice("Writing mosaic to file %s", outfile)
        pyfits.HDUList(hdus).writeto(outfile, overwrite=True)

    def _calc_centers(self, center: SkyCoord):
        # Global x and y center
        gxc = self.layout[0] / 2
        gyc = self.layout[1] / 2
        offsets = {}  # offsets in arcseconds from the FoV center
        for x in range(self.layout[0]):
            # local chip center
            xc = x - (gxc - 0.5)
            xc = self.ccd.pixelsize * (self.ccd.shape[0] + self.gaps[0]) * xc
            for y in range(self.layout[1]):
                yc = y - (gyc - 0.5)
                yc = self.ccd.pixelsize * (self.ccd.shape[1] + self.gaps[1]) * yc
                offsets[(x, y)] = (xc, yc)

        centers = {
            key: center.spherical_offsets_by(*offset) for key, offset in offsets.items()
        }

        return centers


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("config", help="configuration file")
    args = parser.parse_args()
    return args


def run(params):
    rng = np.random.default_rng(params["seed"])
    center = SkyCoord(params["sky"]["ra"], params["sky"]["dec"], unit="degree")
    psfparams = params["psf"]
    starsim = StarSimulation(
        magrange=params["magrange"], psfparams=psfparams, catalog=params["catalog"]
    )
    ccdparams = params["ccd"]
    ccd = CCDImage(
        shape=ccdparams["shape"], pixelsize=ccdparams["pixelsize"], ron=ccdparams["ron"]
    )
    mosaic = Mosaic(ccd, params["mosaic"]["layout"], params["mosaic"]["gaps"])
    stars = starsim.run(
        center,
        mosaic.fov,
        params["t_exp"],
        params["nstars"],
        zeropoint=params["zeropoint"],
        rng=rng,
    )

    if "polarisation" in params:
        polpars = params["polarisation"]
        starsim.magrange = polpars["magrange"]
        starsim.catalog = None
        nstars = polpars["nstars"]
        qrange = polpars["qrange"]
        urange = polpars["urange"]
        I = starsim.run(
            center,
            mosaic.fov,
            params["t_exp"],
            polpars["nstars"],
            zeropoint=params["zeropoint"],
            rng=rng,
        )
        Q = I.copy()
        U = I.copy()
        qfractions = rng.uniform(*qrange, size=nstars)
        Q["flux"] *= qfractions
        ufractions = rng.uniform(*urange, size=nstars)
        U["flux"] *= ufractions
        with open("ds9.reg", "w") as file:
            file.write("# Region file format: DS9 version 4.1\n")
            file.write(
                'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
            )
            file.write("fk5\n")
            for x, y, q, u in zip(I["ra"], I["dec"], qfractions, ufractions):
                file.write(f'circle({x},{y},3.0") # text={{{q:.3f}, {u:.3f}}}\n')
        for angle in polpars["angles"]:
            polstars = I.copy()
            modvec = mod_vector(angle)
            polstars["flux"] = (
                I["flux"] * modvec[0] + Q["flux"] * modvec[1] + U["flux"] * modvec[2]
            ) / 2
            logger.notice("Creating mosaic for polarization angle %.1f", angle)
            outfile = f"polaris-deg{angle:.1f}".replace(".", "_") + ".fits"
            mosaic.simulate(
                center,
                starsim,
                t_exp=params["t_exp"],
                stars=stars,
                skylevel=params["skylevel"],
                outfile=outfile,
                angle=angle,
                polstars=polstars,
                poisson_noise=params["poisson_noise"],
                rng=rng,
            )

    else:
        logger.notice("Creating mosaic")
        mosaic.simulate(
            center,
            starsim,
            t_exp=params["t_exp"],
            stars=stars,
            skylevel=params["skylevel"],
            rng=rng,
        )


def main():
    args = parse_args()
    config = read_config(args.config)
    params = config["simulation"]
    params["seed"] = None if params["random_seed"] < 0 else params["random_seed"]

    run(params)


if __name__ == "__main__":
    main()
