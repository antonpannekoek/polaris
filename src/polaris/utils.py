from functools import partial, partialmethod
import logging
from pathlib import Path


LOGLEVELS = {
    "warn": logging.WARNING,
    "warning": logging.WARNING,
    "notice": (logging.WARNING + logging.INFO),
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
