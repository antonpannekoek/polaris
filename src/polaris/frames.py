from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
import enum
from functools import partial
import logging
from multiprocessing import Pool
from pathlib import Path
import sys

import astropy.io.fits as pyfits
import numpy as np

from .config import config as cfg
from .constants import ANGLES, MATRIX
from .utils import NonePath
from .version import __version__


logger = logging.getLogger(__package__)


# # # Helper enum and exceptions


class StokesParameter(enum.Enum):
    I = enum.auto()
    U = enum.auto()
    Q = enum.auto()
    V = enum.auto()
    # relative parameters
    relU = enum.auto()
    relQ = enum.auto()
    relV = enum.auto()

    def __str__(self):
        return self.name


class MisMatch(enum.Enum):
    NHDUS = enum.auto()
    SHAPE = enum.auto()
    CARD = enum.auto()


class NonMatching(Exception):
    def __init__(self, message: str, code: MisMatch):
        super().__init__(message)
        self.message = message
        self.code = code


class UnAvailable(Exception):
    """Exception to indicate the resources needed are not available"""


# # # Two helper functions


def verify_polangles(
    hdus: list[pyfits.PrimaryHDU | pyfits.ImageHDU], req_angles: set[int]
):
    """Verify that the set of `req_angles` exists in the `hdus`"""

    test_angles = {hdu.header["polangle"] for hdu in hdus}
    if test_angles != req_angles:
        raise ValueError("required polarizer angles missing")


def verify_hdus_equal(files, keys=None) -> list[int]:
    """Verify that the individual HDUs of a set of files are the same size

    If one or more files contain a HDU list, then each image is
    separately matched against the relevant HDU in other images. Any
    missing HDU between images, or a mismatch in shape, results in a
    NonMatching exception.

    If `keys` is passed, each HDU is checked for the existince of a
    header keyword `key`; the value is ignored.

    HDUs that are not images or don't contain relevant (2D) image data
    are ignored.

    The list of matching (correct) HDU indices is returned.

    """

    keys = keys or []

    shapes = {}
    valid_hdus = []  # ensure the number (and type) of HDUs match
    valid_keys = []
    for file in files:
        with pyfits.open(file) as hdulist:
            test_hdus = []
            for i, hdu in enumerate(hdulist):
                if (
                    isinstance(hdu, (pyfits.PrimaryHDU, pyfits.ImageHDU))
                    and isinstance(hdu.data, np.ndarray)
                    and hdu.data.ndim == 2
                ):
                    if shape := shapes.get(i):
                        if hdu.data.shape != shape:
                            raise NonMatching("incorrect shape", MisMatch.SHAPE)
                        test_hdus.append(i)
                    else:
                        # first file; sets the base for comparison
                        shapes[i] = hdu.data.shape
                        valid_hdus.append(i)
                        test_hdus.append(i)
                    for key in keys:
                        # Are the given keys present? Their value doesn't matter
                        if key not in hdu.header:
                            raise NonMatching(
                                f"missing header card: {key}", MisMatch.CARD
                            )

        if valid_hdus != test_hdus:
            raise NonMatching("missing HDU", MisMatch.NHDUS)

    return valid_hdus


# # # Frame and FrameSet definitions


@dataclass  # (frozen=True)
class Frame:
    """Basically a simple alternative for an ImageHDU"""

    header: pyfits.Header
    image: np.ndarray
    path: Path

    @classmethod
    def from_file(cls, path: Path, hdu: int = 0):
        """Loads the relevant HDU from `path` into  memory and instantiates a frame"""
        with pyfits.open(path) as hdulist:
            hdu = hdulist[hdu].copy()
        return cls(hdu.header, hdu.data, path)

    @classmethod
    def from_hdu(cls, hdu: pyfits.ImageHDU, path: Path = NonePath):
        return cls(hdu.header, hdu.data, path)

    def __str__(self):
        return self.path


@dataclass  # (frozen=True)
class PolFrame(Frame):
    angle: int | None = None

    @classmethod
    def from_frame(cls, frame: Frame):
        return cls(frame.header, frame.image, frame.path)

    def __post_init__(self):
        if self.angle is None:
            self.angle = self.header.get("polangle", 0)


@dataclass
class StokesFrame(Frame):
    parameter: StokesParameter

    @classmethod
    def from_file(cls, path: Path, parameter: StokesParameter | str, hdu: int = 0):
        """Loads the relevant HDU from `path` into  memory and instantiates a frame"""
        with pyfits.open(path) as hdulist:
            hdu = hdulist[0].copy()
        if isinstance(parameter, str):
            parameter = StokesParameter[parameter.upper()]
        return cls(hdu.header, hdu.data, path, parameter)

    def __str__(self):
        return f"Stokes {self.parameter}"

    def demodulate(self, frames: "PolFrameSet"):
        """Calcualte a single Stokes parameter from a polarizer frameset

        This is usually not efficient: reducing a polarization
        frameset to a set of Stokes parameter is more efficient. See
        IQUFrameSet for this

        """
        nframes = len(frames)
        matrix = MATRIX[nframes]
        angles = ANGLES[nframes]

        # Set up the initial Stokes parameter frames

        icol = "IQU".index(self.parameter)
        row = matrix[icol]
        self.image = sum(
            frames[angle].image * factor for angle, factor in zip(angles, row)
        )

        # Set the header to that of an input frame
        # Update some fields, and the history
        self.header = frames[0].header
        self.header["stokes"] = self.parameter
        del self.header[cfg.header.polangle]
        self.header["polaris processed"] = str(datetime.now())
        self.header["polaris version"] = __version__
        self.header["polaris process"] = f"demodulate, {nframes} frames"

    def reduce(self, frames: "PolFrameSet", subsky: True):
        """Reduce a frameset to a single Stokes parameter

        This is usually not efficient: reducing a polarization
        frameset to a set of Stokes parameter is more efficient. See
        IQUFrameSet for this

        """
        if not self._verify(frames):
            logger.error("reduction cancelled")
            return
        if len(frames) in [3, 4]:
            self.demodulate(frames)
        else:
            logger.error("unable to demodulate frames")


@dataclass
class RelStokesFrame(StokesFrame):
    def __post_init__(self):
        valid = {
            StokesParameter.I,
            StokesParameter.relQ,
            StokesParameter.relU,
            StokesParameter.relV,
        }
        if self.parameter not in valid:
            raise ValueError(f"incorrect Stokes Parameter {frame.parameter}")


@dataclass
class FrameSet:
    frames: set[Frame]


@dataclass
class PolFrameSet(FrameSet):
    frames: dict[int, PolFrame]
    req_angles: frozenset[int] = frozenset()

    @classmethod
    def from_frames(cls, frames: list[Frame]):
        if not all(isinstance(frame, PolFrame) for frame in frames):
            raise ValueError("incorrect frame type")
        return cls(frames={frame.angle: PolFrame.from_frame(frame) for frame in frames})

    @classmethod
    def from_hdus(cls, hdus: list[pyfits.PrimaryHDU | pyfits.ImageHDU], path=NonePath):
        frames = {}
        for hdu in hdus:
            angle = hdu.header["polangle"]
            frames[angle] = PolFrame(hdu.header, hdu.data, path=path, angle=angle)
        return cls(frames=frames)

    # @classmethod
    def __post_init__(self):
        if not all(isinstance(frame, PolFrame) for frame in self.frames.values()):
            raise ValueError("incorrect frames")
        if not all(key == frame.angle for key, frame in self.frames.items()):
            raise ValueError("keys do not match angles")
        if self.req_angles:
            if self.req_angles != {angle for angle in self.frames}:
                raise ValueError("incorrect angles")


@dataclass
class Pol3FrameSet(PolFrameSet):
    req_angles: frozenset[int] = frozenset([0, 60, 120])


@dataclass
class Pol4FrameSet(PolFrameSet):
    req_angles: frozenset[int] = frozenset({0, 45, 90, 135})


@dataclass
class Pol3FrameSet(PolFrameSet):
    req_angles: set[int] = frozenset({0, 60, 120})


@dataclass
class IQUFrameSet(FrameSet):
    frames: dict[str, StokesFrame]

    @classmethod
    def from_polframeset(cls, frameset: PolFrameSet, subsky: bool = False):
        iquframeset = cls({})
        iquframeset.reduce(frameset, subsky=subsky)

        return cls(iquframeset.frames)

    @classmethod
    def from_hdus(
        cls,
        hdus: list[pyfits.PrimaryHDU | pyfits.ImageHDU],
        path: Path = NonePath,
        subsky: bool = False,
    ):
        """Create an IQUFrameSet from a set of input image HDUs

        The intermediate set of 3 or 4 polarizer frames is created first

        subsky: perform sky subtraction
        """

        nfiles = len(hdus)
        if nfiles == 4:
            polframeset = Pol4FrameSet.from_hdus(hdus, path)
        elif nfiles == 3:
            polframeset = Pol3FrameSet.from_hdus(hdus, path)
        else:
            raise UnAvailable("incorrect number of input files")
        return cls.from_polframeset(polframeset, subsky=subsky)

    @classmethod
    def from_files(cls, files: list[str | Path], hdu: int, subsky: bool = False):
        hdus = []
        for file in files:
            hdus.append(pyfits.open(file)[hdu])
        return cls.from_hdus(hdus, subsky=subsky)

    # @staticmethod
    def reduce(self, frameset: PolFrameSet, subsky: bool = False):
        if len(frameset.frames) in [3, 4]:
            return self.demodulate(frameset)
        else:
            logger.error("unable to demodulate frames")
            return set()

    # @staticmethod
    def demodulate(self, frameset: PolFrameSet):
        frames = frameset.frames
        nframes = len(frames)
        try:
            matrix = MATRIX[nframes]
            angles = ANGLES[nframes]
        except KeyError:
            raise ValueError(f"incorrect number of frames to demodulate: {nframes}")

        # Set up the initial Stokes parameter frames
        frame = next(iter(frames.values()))
        if isinstance(frame, Frame):
            shape = frame.image.shape
        else:
            raise ValueError("frameset does not contain Frames or FrameLists")

        # Perform the matrix - vector multiplication,
        # assigning to the relevant stokes parameters
        # First, turn the frameset into a dict with the angle as key
        # frames = {frame.header['polangle']: frame for frame in frameset.frames}
        parameters = [StokesParameter.I, StokesParameter.Q, StokesParameter.U]
        for row, parameter in zip(matrix, parameters):
            image = np.zeros(shape, dtype=float)  # single "Stokes parameter" image
            for factor, angle in zip(row, angles):
                frame = frames[angle]
                # factor /= nframes
                image += frame.image * factor
            header = frame.header.copy()
            del header["polangle"]
            header["stokes"] = str(parameter)
            path = Path(f"polaris-stokes{parameter}")
            self.frames[parameter] = StokesFrame(
                header=header, image=image, path=path, parameter=parameter
            )

    def to_hdus(self) -> dict:
        return {
            hdus[key]: pyfits.ImageHDU(frame.image, frame.header)
            for key, frame in self.frames.items()
        }


@dataclass
class RelIQUFrameSet(IQUFrameSet):
    frames: dict[str, RelStokesFrame]

    @classmethod
    def from_files(cls, files: list[str | Path], hdu: int, subsky: bool = False):
        hdus = []
        for file in files:
            hdus.append(pyfits.open(file)[hdu])
        return cls.from_hdus(hdus, subsky=subsky)

    @classmethod
    def from_hdus(
        cls,
        hdus: list[pyfits.PrimaryHDU | pyfits.ImageHDU],
        path: Path = NonePath,
        subsky: bool = False,
    ):
        """Create an IQUFrameSet from a set of input image HDUs

        The intermediate set of 3 or 4 polarizer frames is created first.
        Then an IQUFrameSet is created.
        Finally, the RelIQUFrameSet is created.

        subsky: perform sky subtraction
        """

        iquframeset = IQUFrameSet.from_hdus(hdus, path=path, subsky=subsky)
        frames = {}
        iframe = iquframeset.frames[StokesParameter.I]
        qframe = iquframeset.frames[StokesParameter.Q]
        uframe = iquframeset.frames[StokesParameter.U]

        # Since sky-subtracted frames are used in the denominator, we need to account for potential divisions by (near) zero.
        # We set NaNs and very large (absolute) values to zero in the resulting image.
        # Alternatively, we could have added 1 to the I image, so the median sky is 1, but this does add
        # an artificial flux to the stars as well.
        parameter = StokesParameter.I
        frames[parameter] = RelStokesFrame(
            header=iframe.header,
            image=iframe.image,
            path=iframe.path,
            parameter=parameter,
        )
        parameter = StokesParameter.relQ
        image = qframe.image / iframe.image
        image = np.nan_to_num(image)
        image[(image < -1e10) | (image > 1e0)] = 0.0
        header = qframe.header
        header["stokes"] = str(parameter)
        frames[parameter] = RelStokesFrame(
            header=header, image=image, path=qframe.path, parameter=parameter
        )
        parameter = StokesParameter.relU
        image = uframe.image / iframe.image
        image = np.nan_to_num(image)
        image[(image < -1e10) | (image > 1e0)] = 0.0
        header = uframe.header
        header["stokes"] = str(parameter)
        frames[parameter] = RelStokesFrame(
            header=header, image=image, path=uframe.path, parameter=parameter
        )
        return cls(frames)


class IQUMosaic:
    """Container for multiple IQUFrameSets, for multiple-HDU FITS files"""

    @classmethod
    def from_files(cls, files, subsky: bool = False, nproc=1):
        framesets = []
        ihdus = verify_hdus_equal(files, keys=None)
        hdus = defaultdict(list)
        paths = defaultdict(list)
        for ihdu in ihdus:
            for file in files:
                with pyfits.open(file) as hdulist:
                    hdus[ihdu].append(hdulist[ihdu].copy())

        if nproc > 1:
            with Pool(nproc) as pool:
                paths = [NonePath] * len(hdus)
                subskies = [subsky] * len(hdus)
                framesets = pool.starmap(
                    IQUFrameSet.create, zip(hdus.values(), paths, subskies)
                )
        else:
            framesets = [
                IQUFrameSet.create(hduset, NonePath, subsky) for hduset in hdus.values()
            ]

        return cls(framesets)

    def __init__(self, framesets: list[IQUFrameSet], *args, **kwargs):
        self.framesets = framesets

    def __str__(self):
        return str([str(frameset) for frameset in self.framesets])

    def write_fits(self, path: str | Path, overwrite=False):
        # Invert the framesets: from a list of
        # IQU frames, to a set of I, Q and U lists of frames
        hdulists = {}
        keys = []  # convenience list
        # Initialize dict with an empty PrimaryHDU
        for frameset in self.framesets[0]:
            header = frameset.header
            key = frameset.parameter
            hdulists[key] = pyfits.HDUList(pyfits.PrimaryHDU(header=header))
            keys.append(key)
        # Append the actual ImageHDUs
        for frameset in self.framesets:
            for frame in frameset:
                key = frame.parameter
                hdulists[key].append(pyfits.ImageHDU(frame.image, frame.header))
        if not path:
            for key in keys:
                path = Path(f"polaris-stokes{key}-mosaic.fits")
                hdulists[key].writeto(path, overwrite=overwrite)


def write_fits(
    framesets: list[IQUFrameSet],
    prefix="polaris-stokes",
    postfix="-mosaic",
    overwrite=False,
):
    # Invert the framesets: from a list of
    # IQU frames, to a set of I, Q and U lists of frames
    hdulists = {}
    keys = []  # convenience list
    # Initialize dict with an empty PrimaryHDU
    for key, frame in framesets[0].frames.items():
        header = frame.header
        key = frame.parameter
        hdulists[key] = pyfits.HDUList(pyfits.PrimaryHDU(header=header))
        keys.append(key)
    # Append the actual ImageHDUs
    for frameset in framesets:
        for key, frame in frameset.frames.items():
            hdulists[key].append(pyfits.ImageHDU(frame.image, frame.header))
    for key in keys:
        path = Path(f"{prefix}{key}{postfix}.fits")
        print(f"{key = }; {key}, {path = }")
        hdulists[key].writeto(path, overwrite=overwrite)


def main():
    print(1)

    frame = Frame.from_file("polaris-deg60_0.fits", hdu=2)
    print(2)
    frames = [
        Frame.from_file("polaris-deg0_0.fits", hdu=2),
        Frame.from_file("polaris-deg60_0.fits", hdu=2),
        Frame.from_file("polaris-deg120_0.fits", hdu=2),
    ]
    print(3)
    frameset = FrameSet(frames)
    print(4)
    try:
        polframeset = PolFrameSet.from_frames(frames)
        print("ERROR")
    except ValueError as exc:
        assert str(exc) == "incorrect frame type"
        print(5)
    frames = [
        PolFrame.from_file("polaris-deg0_0.fits", hdu=2),
        PolFrame.from_file("polaris-deg60_0.fits", hdu=2),
    ]
    polframeset = PolFrameSet.from_frames(frames)
    print(6)
    try:
        polframeset = Pol3FrameSet.from_frames(frames)
        print("ERROR")
    except ValueError as exc:
        assert str(exc) == "incorrect angles"
        print(7)
    frames = [
        PolFrame.from_file("polaris-deg0_0.fits", hdu=2),
        PolFrame.from_file("polaris-deg60_0.fits", hdu=2),
        PolFrame.from_file("polaris-deg120_0.fits", hdu=2),
    ]
    pol3frameset = Pol3FrameSet.from_frames(frames)
    print(8)
    stokesI = StokesFrame.from_file("polaris-stokesI.fits", parameter="i", hdu=2)

    files = ["polaris-deg0_0.fits", "polaris-deg60_0.fits", "polaris-deg120_0.fits"]
    framesetIQU = IQUFrameSet.from_polframeset(pol3frameset)
    print(9)

    indices = verify_hdus_equal(files)
    iquframesets = []

    import time

    start = time.time()
    for index in indices:
        iquframesets.append(IQUFrameSet.from_files(files, hdu=index))
    delta = time.time() - start
    print(10, delta)
    # Multiprocessing is slow as long as no shared memory is used.
    # In this particular case, it is the copying from the reduced frame
    # back to the main thread, into the `uqyframesets`.
    if len(sys.argv) > 1:
        nproc = int(sys.argv[1])
        func = partial(IQUFrameSet.from_files, files)
        iquframesets = []
        start = time.time()
        with Pool(nproc) as pool:
            iquframesets = pool.map(func, indices)
        delta = time.time() - start
        print(11, delta)
    write_fits(iquframesets, overwrite=True)

    reliquframesets = []
    for index in indices:
        reliquframesets.append(RelIQUFrameSet.from_files(files, hdu=index))
    print(12)
    write_fits(reliquframesets, overwrite=True)


if __name__ == "__main__":
    main()
