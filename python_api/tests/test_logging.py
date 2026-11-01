"""Tests for logging."""
import logging
from cis_nto.logging import get_logger, set_verbosity


def test_get_logger_default():
    logger = get_logger("test")
    assert logger.name == "cis_nto.test"
    assert logger.level in (logging.NOTSET, logging.INFO)


def test_set_verbosity_debug():
    set_verbosity("DEBUG")
    logger = get_logger()
    assert logger.level == logging.DEBUG


def test_set_verbosity_int():
    set_verbosity(logging.WARNING)
    logger = get_logger()
    assert logger.level == logging.WARNING
