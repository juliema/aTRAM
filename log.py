"""
Functions that handle common aTRAM logging functions.
"""

import sys
import logging


def setup(args):
    """
    Set up the logs for a common format. We need the prefix of the output log file name and
    the command-line arguments for the starting message. Both are gotten from the user input.
    """
    logging.basicConfig(
        filename='{}{}.log'.format(args.out, sys.argv[0][:-3]),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))
