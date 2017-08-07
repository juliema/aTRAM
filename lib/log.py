"""Common logging functions."""

from os.path import basename
import sys
import logging
import tempfile
import subprocess


def setup(name):
    """Standard logger setup."""

    logging.basicConfig(
        filename=name,
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))


def file_name(args, blast_db, query_file):
    """Setup default log file name for each run."""

    if args.log_file:
        return args.log_file

    program = basename(sys.argv[0][:-3])
    return '{}.{}.{}.log'.format(blast_db, query_file, program)


def subcommand(cmd, temp_dir, timeout=None):
    """Handle subprocess calls and log their output."""

    logging.info(cmd)

    with tempfile.NamedTemporaryFile(mode='w', dir=temp_dir) as log_output:

        # Note: stdout=PIPE is blocking and large logs cause a hang
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
                        logging.log(logging.INFO, line)


def info(msg, line_break=None):
    """Log and display an info message."""

    if line_break is not None:
        print(line_break)
    print(msg)

    if line_break is not None:
        logging.info(line_break)
    logging.info(msg)


def error(msg):
    """Log and display an error message."""

    print(msg)
    logging.error(msg)


def fatal(msg):
    """Log an error message and exit."""

    error(msg)
    sys.exit(1)
