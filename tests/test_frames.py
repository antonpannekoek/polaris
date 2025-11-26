from polaris.frames import verify_hdus_equal, NonMatching, MisMatch

import astropy.io.fits as pyfits
from astropy.io.fits import Header
import numpy as np
import pytest


def test_verify_hdus_equal(tmp_path):
    file1 = tmp_path / "file1.fits"
    file2 = tmp_path / "file2.fits"
    file3 = tmp_path / "file3.fits"
    files = [file1, file2, file3]

    card = ("exptime", 300, "seconds")
    cards = [card, ("polangle", 0, "polarizer angle")]
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((12, 8)), header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
        ]
    ).writeto(file1, overwrite=True)
    cards = [card, ("polangle", 60, "polarizer angle")]
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((12, 8)), header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
        ]
    ).writeto(file2, overwrite=True)
    cards = [card, ("polangle", 120, "polarizer angle")]
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((12, 8)), header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
        ]
    ).writeto(file3, overwrite=True)

    verify_hdus_equal(files, keys=["polangle"])

    # Missing required card
    cards = [card]
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((12, 8)), header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
        ]
    ).writeto(file3, overwrite=True)

    with pytest.raises(NonMatching) as excinfo:
        verify_hdus_equal(files, keys=["polangle"])
        assert excinfo.value.code == MisMatch.CARD

    # Missing one HDU
    cards = [card, ("polangle", 0, "polarizer angle")]
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((12, 8)), header=Header(cards)),
        ]
    ).writeto(file3, overwrite=True)

    with pytest.raises(NonMatching) as excinfo:
        verify_hdus_equal(files, keys=["polangle"])
        assert excinfo.value.code == MisMatch.NHDUS

    # First image has an extra HDU with an image
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((12, 8)), header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
        ]
    ).writeto(file1, overwrite=True)
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((12, 8)), header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
        ]
    ).writeto(file2, overwrite=True)

    with pytest.raises(NonMatching) as excinfo:
        verify_hdus_equal(files, keys=["polangle"])
        assert excinfo.value.code == MisMatch.NHDUS

    # Shape in a single HDU differs
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((12, 8)), header=Header(cards)),
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
        ]
    ).writeto(file1, overwrite=True)
    pyfits.HDUList(
        [
            pyfits.PrimaryHDU(header=Header(cards)),
            pyfits.ImageHDU(
                data=np.ones((24, 16)), header=Header(cards)
            ),  # incorrect shape
            pyfits.ImageHDU(data=np.ones((24, 16)), header=Header(cards)),
        ]
    ).writeto(file2, overwrite=True)

    with pytest.raises(NonMatching) as excinfo:
        verify_hdus_equal(files, keys=["polangle"])
        assert excinfo.value.code == MisMatch.SHAPE
