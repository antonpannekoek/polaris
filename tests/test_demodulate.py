from polaris.frames import PolFrame, StokesParameter
from polaris.pipeline import demodulate

import astropy.io.fits as pyfits
import numpy as np


def test_demodulate_unpolarized():
    frames = {}
    # simple 2-D array with zero background,
    # and a single "point source"
    image = np.zeros((8, 8), dtype=float)
    image[2:6, 2:6] = 1
    image[3:5, 3:5] = 3
    hdu = pyfits.PrimaryHDU(image)
    hdu.header["polangle"] = 0
    frames[0] = PolFrame(header=hdu.header, image=hdu.data)
    hdu.header["polangle"] = 60
    frames[60] = PolFrame(header=hdu.header, image=hdu.data)
    hdu.header["polangle"] = 120
    frames[120] = PolFrame(header=hdu.header, image=hdu.data)

    stokes_frames = demodulate(frames)

    frame = stokes_frames[StokesParameter.I]
    assert frame.header.get("stokes") == "I"
    np.testing.assert_allclose(frame.image, image, strict=True)

    zeros = np.zeros((8, 8), dtype=float)
    frame = stokes_frames[StokesParameter.Q]
    assert frame.header.get("stokes") == "Q"
    np.testing.assert_allclose(frame.image, zeros, strict=True)

    frame = stokes_frames[StokesParameter.U]
    assert frame.header.get("stokes") == "U"
    np.testing.assert_allclose(frame.image, zeros, strict=True)


def test_demodulate_polarized():
    frames = {}
    # simple 2-D array with zero background,
    # and a single "point source"
    image0 = np.zeros((8, 8), dtype=float)
    image0[2:6, 2:6] = 1.0
    image0[3:5, 3:5] = 3.0
    hdu = pyfits.PrimaryHDU(image0)
    hdu.header["polangle"] = 0
    frames[0] = PolFrame(header=hdu.header, image=hdu.data)
    hdu.header["polangle"] = 60
    # Adjust the intensity of the 60 degrees frame only
    image60 = np.zeros_like(image0)
    image60[2:6, 2:6] = 0.5
    image60[3:5, 3:5] = 1.5
    frames[60] = PolFrame(header=hdu.header, image=image60)
    hdu.header["polangle"] = 120
    frames[120] = PolFrame(header=hdu.header, image=hdu.data)

    stokes_frames = demodulate(frames)

    frame = stokes_frames[StokesParameter.I]
    imageI = np.zeros_like(image0)
    imageI[2:6, 2:6] = 5 / 6
    imageI[3:5, 3:5] = 2.5
    np.testing.assert_allclose(frame.image, imageI, strict=True)

    frame = stokes_frames[StokesParameter.Q]
    imageQ = np.zeros_like(image0)
    imageQ[2:6, 2:6] = 1 / 6
    imageQ[3:5, 3:5] = 0.5
    np.testing.assert_allclose(frame.image, imageQ, strict=True)

    frame = stokes_frames[StokesParameter.U]
    imageU = np.zeros_like(image0)
    imageU[2:6, 2:6] = 0.5 / np.sqrt(3)
    imageU[3:5, 3:5] = 1.5 / np.sqrt(3)
    np.testing.assert_allclose(frame.image, imageU, strict=True)

    # Increase the intensity in the 120 degree frame as well
    frames[120] = PolFrame(header=hdu.header, image=image60)

    stokes_frames = demodulate(frames)

    frame = stokes_frames[StokesParameter.I]
    imageI = np.zeros_like(image0)
    imageI[2:6, 2:6] = 2 / 3
    imageI[3:5, 3:5] = 2.0
    np.testing.assert_allclose(frame.image, imageI, strict=True)

    frame = stokes_frames[StokesParameter.Q]
    imageQ = np.zeros_like(image0)
    imageQ[2:6, 2:6] = 1 / 3
    imageQ[3:5, 3:5] = 1.0
    np.testing.assert_allclose(frame.image, imageQ, strict=True)

    frame = stokes_frames[StokesParameter.U]
    imageU = np.zeros_like(image0)
    np.testing.assert_allclose(frame.image, imageU, strict=True)

    # Reset the intensity in the 60 degree frame to that of the 0 degrees frame.
    # Effectively, this flips the U parameter around its zeropoint.
    frames[60] = PolFrame(header=hdu.header, image=image0)

    stokes_frames = demodulate(frames)

    print(stokes_frames)

    frame = stokes_frames[StokesParameter.I]
    imageI = np.zeros_like(image0)
    imageI[2:6, 2:6] = 5 / 6
    imageI[3:5, 3:5] = 2.5
    np.testing.assert_allclose(frame.image, imageI, strict=True)

    frame = stokes_frames[StokesParameter.Q]
    imageQ = np.zeros_like(image0)
    imageQ[2:6, 2:6] = 1 / 6
    imageQ[3:5, 3:5] = 0.5
    np.testing.assert_allclose(frame.image, imageQ, strict=True)

    frame = stokes_frames[StokesParameter.U]
    imageU = np.zeros_like(image0)
    imageU[2:6, 2:6] = -0.5 / np.sqrt(3)
    imageU[3:5, 3:5] = -1.5 / np.sqrt(3)
    np.testing.assert_allclose(frame.image, imageU, strict=True)
