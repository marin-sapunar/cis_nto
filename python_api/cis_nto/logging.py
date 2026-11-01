"""Logging utilities for cis_nto Python API."""
import logging
import sys

_PACKAGE_LOGGER = "cis_nto"
_DEFAULT_FORMAT = "[%(name)s] %(levelname)s: %(message)s"


def get_logger(name=None, level=None):
    """Get or create a cis_nto logger."""
    full_name = f"{_PACKAGE_LOGGER}.{name}" if name else _PACKAGE_LOGGER
    logger = logging.getLogger(full_name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(_DEFAULT_FORMAT))
        logger.addHandler(handler)
    if level is not None:
        logger.setLevel(level if isinstance(level, int) else getattr(logging, level.upper()))
    elif logger.level == logging.NOTSET:
        logger.setLevel(logging.INFO)
    return logger


def set_verbosity(level):
    """Set global cis_nto verbosity. Accepts 'DEBUG', 'INFO', 'WARNING', 'ERROR'."""
    if isinstance(level, str):
        level = getattr(logging, level.upper())
    logging.getLogger(_PACKAGE_LOGGER).setLevel(level)
