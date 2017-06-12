"""Common logging functions."""

import sys
import logging
import tempfile
import subprocess


def setup(args):
    """Standard logger setup."""

    logging.basicConfig(
        filename=args.log_file,
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))


def subcommand(cmd, temp_dir):
    """Handle subprocess calls and log their output."""

    logging.info(cmd)

    with tempfile.NamedTemporaryFile(mode='w', dir=temp_dir) as log_output:

        # Note: stdout=PIPE is blocking and large logs cause a hang
        subprocess.check_call(
            cmd, shell=True, stdout=log_output, stderr=log_output)

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
