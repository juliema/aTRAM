"""Common logging functions."""

from os.path import basename, splitext
import sys
import logging
import tempfile
import subprocess
from . import db

LOGGER = None  # Global logger so we can switch between queries & blast DBs
FORMATTER = logging.Formatter('%(asctime)s %(levelname)s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
NAME = 'atram_logger'


def setup(log_file, log_level, blast_db, query_file=''):
    """Logger setup."""
    log_file = file_name(log_file, blast_db, query_file)
    _setup(log_file, log_level)


def stitcher_setup(log_file, log_level):
    """Build a logger for the stitcher."""
    _setup(log_file, log_level)


def _setup(log_file, log_level):
    global LOGGER  # pylint: disable=global-statement

    if not LOGGER:
        handler = logging.FileHandler(log_file)
        handler.setFormatter(FORMATTER)
        handler.setLevel(logging.DEBUG)

        stream = logging.StreamHandler()
        stream.setFormatter(FORMATTER)
        stream.setLevel(logging.INFO)

        LOGGER = logging.getLogger(log_file)

        log_level = getattr(logging, log_level.upper())
        LOGGER.setLevel(log_level)

        LOGGER.addHandler(handler)
        LOGGER.addHandler(stream)

        info('#' * 80)
        info('aTRAM version: {}'.format(db.ATRAM_VERSION))
        info('Python version: {}'.format(' '.join(sys.version.split())))
        info(' '.join(sys.argv[:]))


def file_name(log_file, blast_db, query_file=''):
    """
    Create the log file name for each run.

    Honor user's argument if given.
    """
    if log_file:
        return log_file

    program = splitext(basename(sys.argv[0]))[0]

    if query_file:
        query_file = splitext(basename(query_file))[0]
        return '{}.{}.{}.log'.format(blast_db, query_file, program)

    return '{}.{}.log'.format(blast_db, program)


def subcommand(cmd, temp_dir, timeout=None):
    """
    Call a subprocess and log the output.

    Note: stdout=PIPE is blocking and large logs cause a hang.
    So we don't use it.
    """
    LOGGER.debug(cmd)

    with tempfile.NamedTemporaryFile(mode='w', dir=temp_dir) as log_output:
        try:
            subprocess.check_call(
                cmd,
                shell=True,
                timeout=timeout,
                stdout=log_output,
                stderr=log_output)
        except Exception as err:  # pylint: disable=broad-except
            error('Exception: {}'.format(err))
        finally:
            with open(log_output.name) as log_input:
                for line in log_input:
                    line = line.strip()
                    if line:
                        LOGGER.debug(line)


def info(msg):
    """Log an info message."""
    LOGGER.info(msg)


def error(msg):
    """Log an error message."""
    LOGGER.error(msg)


def fatal(msg):
    """Log an error message and exit."""
    error(msg)
    sys.exit(1)
