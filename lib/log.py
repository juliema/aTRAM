"""Common logging functions."""

from os.path import basename, splitext
import sys
import logging
import tempfile
import subprocess

LOGGER = None  # Global logger so we can switch between queries & blast DBs
FORMATTER = logging.Formatter('%(asctime)s %(levelname)s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')


def setup(log_file, blast_db, query):
    """Standard logger setup."""
    global LOGGER

    log_file = file_name(log_file, blast_db, query)

    handler = logging.FileHandler(log_file)
    handler.setFormatter(FORMATTER)
    handler.setLevel(logging.DEBUG)

    stream = logging.StreamHandler()
    stream.setFormatter(FORMATTER)
    stream.setLevel(logging.INFO)

    LOGGER = logging.getLogger(log_file)
    LOGGER.setLevel(logging.DEBUG)
    LOGGER.addHandler(handler)
    LOGGER.addHandler(stream)
    # logging.basicConfig(
    #     filename=log_file,
    #     level=logging.DEBUG,
    #     format='%(asctime)s %(levelname)s: %(message)s',
    #     datefmt='%Y-%m-%d %H:%M:%S')
    # logging.info(' '.join(sys.argv))


def file_name(log_file, blast_db, query_file=''):
    """Setup default log file name for each run.

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
    """Handle subprocess calls and log their output.

    Note: stdout=PIPE is blocking and large logs cause a hang.
    So we don't use it.
    """
    LOGGER.debug(cmd)

    with tempfile.NamedTemporaryFile(mode='w', dir=temp_dir) as log_output:
        try:
            subprocess.check_call(cmd,
                                  shell=True,
                                  timeout=timeout,
                                  stdout=log_output,
                                  stderr=log_output)
        except (subprocess.CalledProcessError, TimeoutError):
            raise
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
