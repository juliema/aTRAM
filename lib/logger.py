"""Logger functions."""

import os
import sys
import logging


def setup(work_dir, blast_db):
    """Set up the logs for a common format. We need the prefix of the
    output log file name and the command-line arguments for the starting
    message. Both are gotten from the user input.
    """

    file_name = '{}.{}.log'.format(blast_db, sys.argv[0][:-3])

    logging.basicConfig(
        filename=os.path.join(work_dir, file_name),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))
