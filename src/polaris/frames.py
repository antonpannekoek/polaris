from dataclasses import dataclass
import enum
import logging

import astropy.io.fits as pyfits
import numpy as np


logger = logging.getLogger(__package__)


# # # Helper enum and exceptions


class StokesParameter(enum.Enum):
    I = enum.auto()  # noqa: E741
    Q = enum.auto()
    U = enum.auto()
    V = enum.auto()
    # relative parameters
    relU = enum.auto()
    relQ = enum.auto()
    relV = enum.auto()

    def __str__(self):
        return self.name


# Frame and FrameSet definitions


@dataclass
class Frame:
    """Basically a simple alternative for an ImageHDU"""

    header: pyfits.Header
    image: np.ndarray
    name: str = ""


@dataclass
class PolFrame(Frame):
    angle: int | None = None


@dataclass
class StokesFrame(Frame):
    parameter: StokesParameter | None = None


# Additional types

type PolFrames = dict[int, PolFrame]
type StokesFrames = dict[StokesParameter, StokesFrame]
