"""Misc. utilities."""

import os
from os.path import exists
import sys


def update_temp_dir(temp_dir, args):
    """Handle the new temporary directory name."""
    args['temp_root'] = args['temp_dir']
    args['temp_dir'] = temp_dir
    os.environ['SQLITE_TMPDIR'] = temp_dir


def set_blast_batch_size(batch_size):
    """Use this setting to control blast memory usage & query concatenation."""
    if batch_size:
        os.environ['BATCH_SIZE'] = str(batch_size)


def temp_dir_exists(temp_dir, debug_dir):
    """Make sure the temporary directory exits."""
    if temp_dir and not exists(temp_dir):
        sys.exit('The temporary directory must exist.')
    if debug_dir and not exists(debug_dir):
        sys.exit('The temporary debug directory must exist.')
