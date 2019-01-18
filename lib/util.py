"""Misc. utilities."""

import os
import sys


def iter_dir(temp_dir, blast_db, query_name, iteration):
    """
    Get the work directory for the current iteration.

    We need to call this function in child processes so it cannot be in an
    object.
    """
    name = '{}_{}_{:02d}'.format(
        os.path.basename(blast_db), os.path.basename(query_name), iteration)

    return os.path.join(temp_dir, name)


def update_temp_dir(temp_dir, args):
    """Handle the new temporary directory name."""
    args['temp_dir'] = temp_dir
    os.environ['SQLITE_TMPDIR'] = temp_dir


def set_blast_batch_size(batch_size):
    """Use this setting to control blast memory usage & query concatenation."""
    if batch_size:
        os.environ['BATCH_SIZE'] = str(batch_size)


def temp_dir_exists(temp_dir):
    """Make sure the temporary directory exits."""
    if temp_dir and not os.path.exists(temp_dir):
        sys.exit('The temporary directory must exist.')
