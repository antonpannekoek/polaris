from functools import partial, partialmethod
import logging
from pathlib import Path


LOGLEVELS = {
    "warn": logging.WARNING,
    "warning": logging.WARNING,
    "notice": (logging.WARNING + logging.INFO) // 2,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}


logger = logging.getLogger(__package__)


def add_logging_level(level, name):
    """Add an extra logging level

    See https://stackoverflow.com/questions/2183233/how-to-add-a-custom-loglevel-to-pythons-logging-facility
    for details

    """

    logging.addLevelName(level, name)
    setattr(logging, name, level)
    name = name.lower()
    logclass = logging.getLoggerClass()
    setattr(logclass, name, partialmethod(logclass.log, level))
    setattr(logger, name, partial(logger.log, level))


add_logging_level(LOGLEVELS["notice"], "NOTICE")


def set_logging(level: str = "NOTICE", logfile: str = "", file_level: str = ""):
    """Set up some default logging configuration.

    Note: this doesn't use a dict- or file-config; it is felt that the
    logging setup should still be relatively simple, with only options
    for the logging level and whether or not to (also) log to file.

    """

    if level is None:
        level = cfg.log.loglevel

    level = LOGLEVELS[level.lower()]
    fmt = "%(asctime)s  [%(levelname)-5s] - %(module)s.%(funcName)s():%(lineno)d: %(message)s"
    formatter = logging.Formatter(fmt, datefmt="%y-%m-%d %H:%M:%S")
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    handler.setLevel(level)
    logger.addHandler(handler)
    if logfile:
        file_level = LOGLEVELS[file_level.lower()]
        handler = logging.FileHandler(logfile)
        handler.setFormatter(formatter)
        handler.setLevel(file_level)
        logger.addHandler(handler)
    logger.setLevel(logging.DEBUG)


class NonePath(Path):
    def exists(self, **_):
        return False

    def is_file(self, **_):
        return False

    def is_dir(self, **_):
        return False

    def joinpath(self, *pathsegments):
        return self

    def with_segments(self, *pathsegments):
        return self

    def absolute(self):
        return self
