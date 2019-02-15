"""Misc. utilities."""

import os
from os.path import exists
import sys
from shutil import rmtree
import gzip
import bz2
from contextlib import contextmanager
from tempfile import mkdtemp


def set_blast_batch_size(batch_size):
    """Use this setting to control blast memory usage & query concatenation."""
    if batch_size:
        os.environ['BATCH_SIZE'] = str(batch_size)


def write_fasta_record(out_file, seq_name, seq, seq_end=None):
    """Write a fasta record to the file."""
    out_file.write('>')
    out_file.write(seq_name)
    if seq_end:
        out_file.write('/')
        out_file.write(seq_end)
    out_file.write('\n')

    out_file.write(seq)
    out_file.write('\n')


def temp_dir_exists(temp_dir, debug_dir=None):
    """Make sure the temporary directory exits."""
    if temp_dir and not exists(temp_dir):
        sys.exit('The temporary directory must exist.')
    if debug_dir and not exists(debug_dir):
        sys.exit('The temporary debug directory must exist.')


def update_temp_dir(temp_dir, args):
    """Handle the new temporary directory name."""
    args['temp_root'] = args['temp_dir']
    args['temp_dir'] = temp_dir
    os.environ['SQLITE_TMPDIR'] = temp_dir


@contextmanager
def make_temp_dir(where=None, prefix=None, keep=False):
    """Handle creation and deletion of temporary directory."""
    temp_dir = mkdtemp(prefix=prefix, dir=where)
    try:
        yield temp_dir
    finally:
        if not keep or not where:
            rmtree(temp_dir)


@contextmanager
def open_file(args, file_name):
    """Handle creation and deletion of temporary directory."""
    if args.get('gzip'):
        stream = gzip.open(file_name, 'rt')
    elif args.get('bzip'):
        stream = bz2.open(file_name, 'rt')
    else:
        stream = open(file_name)

    try:
        yield stream
    finally:
        stream.close()
