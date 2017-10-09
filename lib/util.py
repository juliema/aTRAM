"""Misc. utilities."""

from os.path import basename, join


def iter_dir(temp_dir, blast_db, query_name, iteration):
    """
    Get the work directory for the current iteration.

    We need to call this function in child processes so it cannot be in an
    object.
    """
    name = '{}_{}_{:02d}'.format(
        basename(blast_db), basename(query_name), iteration)

    return join(temp_dir, name)
