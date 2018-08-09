import pytest
import logging

import gromacs.log

# This is almost certainly not thread/parallel safe.

@pytest.fixture
def logger_with_file(tmpdir):
    logfile = tmpdir.join('gmx.log')
    logger = gromacs.log.create("GMX", logfile=str(logfile))
    return logfile, logger

def test_create(logger_with_file):
    logfile, logger = logger_with_file
    logger.info("Jabberwock")
    logger.debug("Cheshire Cat")
    txt = logfile.read('rb').decode('utf-8')
    assert "Jabberwock" in txt
    assert "Cheshire Cat" in txt

def test_clear_handlers(logger_with_file):
    logfile, logger = logger_with_file
    gromacs.log.clear_handlers(logger)
    logger.warn("Dodo")
    txt = logfile.read('rb').decode('utf-8')
    assert "Dodo" not in txt

def test_NullHandler():
    h = gromacs.log.NullHandler()
    logger = logging.getLogger("GMX")
    logger.addHandler(h)
    logger.warn("screaming in silence")
    assert True   # not sure what to test here
