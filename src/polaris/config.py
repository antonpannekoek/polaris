from dataclasses import dataclass
from pathlib import Path
import tomllib


@dataclass
class Header:
    """Set of header keywords"""

    polangle: str = "polangle"


@dataclass
class Config:
    header: Header  # header keywords


def read_config(path: str | Path) -> dict:
    with open(path, "rb") as fp:
        config = tomllib.load(fp)
    return config


config = Config(
    Header(
        polangle="polangle",
    )
)
