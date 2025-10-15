from pathlib import Path
import tomllib


def read_config(path: str|Path) -> dict:
    with open(path, "rb") as fp:
        config = tomllib.load(fp)
    return config
