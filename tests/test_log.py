import pytest
import logging

import gromacs
import gromacs.log

import logging

# This is almost certainly not thread/parallel safe.


@pytest.fixture
def logger_with_file(tmp_path):
    logfile = tmp_path / "gmx.log"
    logger = gromacs.log.create("GMX", logfile=str(logfile))
    return logfile, logger


def test_create(logger_with_file):
    logfile, logger = logger_with_file
    logger.info("Jabberwock")
    logger.debug("Cheshire Cat")
    txt = logfile.read_text()
    assert "Jabberwock" in txt
    assert "Cheshire Cat" in txt


def test_clear_handlers(logger_with_file):
    logfile, logger = logger_with_file
    gromacs.log.clear_handlers(logger)
    logger.warning("Dodo")
    txt = logfile.read_text()
    assert "Dodo" not in txt


def test_NullHandler():
    h = gromacs.log.NullHandler()
    logger = logging.getLogger("GMX")
    logger.addHandler(h)
    logger.warning("screaming in silence")
    assert True  # not sure what to test here


@pytest.fixture
def gromacs_logger(tmp_path):
    logfile = tmp_path / "gromacs.log"
    gromacs.start_logging(logfile=str(logfile))
    logger = logging.getLogger("gromacs")
    logger.info("Running a test for logging")
    gromacs.stop_logging()
    return logfile


def _assert_msg_in_log(logfile, msg):
    output = logfile.read_text()
    assert msg in output


def test_start_logger(gromacs_logger):
    _assert_msg_in_log(
        gromacs_logger, "GromacsWrapper {} STARTED".format(gromacs.__version__)
    )


def test_using_logger(gromacs_logger):
    _assert_msg_in_log(gromacs_logger, "Running a test for logging")


def test_stop_logger(gromacs_logger):
    _assert_msg_in_log(
        gromacs_logger, "GromacsWrapper {} STOPPED".format(gromacs.__version__)
    )
