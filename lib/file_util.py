"""Utility functions."""

import os
import sys


def temp_root_dir(args, temp_dir):
    """Make a temporary root directory. Verify that it is clean."""

    if args.temp_dir:
        args.temp_dir = os.path.abspath(args.temp_dir)
        os.makedirs(args.temp_dir, exist_ok=True)
        if os.listdir(args.temp_dir):
            sys.exit('"{}" directory is not empty.'.format(args.temp_dir))
    else:
        args.temp_dir = temp_dir


def temp_iter_dir(temp_dir, iteration):
    """Make a temp directory for the iteration underneath the root temporary
    directory. Make sure that the iteration directory exists.
    """

    iter_name = 'iteration_{:02d}'.format(iteration)
    iter_dir = os.path.join(temp_dir, iter_name)
    os.makedirs(iter_dir, exist_ok=True)
    return iter_dir


def temp_iter_file(temp_dir, iteration, file_name):
    """Make a temporary file name that is inside of the iteration temporary
    directory.
    """

    iter_dir = temp_iter_dir(temp_dir, iteration)
    return os.path.join(iter_dir, file_name)


def output_file(args, file_suffix):
    """Build the output file name."""

    return '{}.{}'.format(args.output, file_suffix)
