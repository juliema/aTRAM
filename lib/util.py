"""Misc. utilities."""

import os
from os.path import exists, getsize, join, split
import io
import re
import sys
from shutil import rmtree
import gzip
import bz2
from contextlib import contextmanager
from tempfile import mkdtemp
from Bio.SeqIO.FastaIO import SimpleFastaParser


def shorten(text):
    """Collapse whitespace in a string."""
    return ' '.join(text.split())


def set_blast_batch_size(batch_size):
    """Use this to control blast memory usage & query concatenation."""
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
    args['temp_dir'] = str(temp_dir)
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


def clean_name(name):
    """Replace problem characters in file names."""
    return re.sub(r'[^\w.]+', '_', name.strip())


def as_word(number):
    """Convert a number in a word.

    If this gets complex we will add the inflect module instead.
    """
    ordinal = {
        1: 'First',
        2: 'Second',
        3: 'Third'}
    return ordinal.get(number, '{}th'.format(number))


def fasta_file_is_empty(fasta_path):
    """Check if a fasta file is either empty or does not have a sequence."""
    if os.stat(fasta_path).st_size == 0:
        return True

    with open(fasta_path) as fasta_file:
        _, seq = next(SimpleFastaParser(fasta_file))

    if not seq:
        return True

    return False


def is_fastq_file(args, file_name):
    """Check if this a FASTQ file."""
    if args.get('fasta'):
        return False
    if args.get('fastq'):
        return True

    parts = file_name.lower().split('.')
    index = -2 if re.search(r'[zp2]$', parts[-1]) and len(parts) > 2 else -1
    return parts[index].startswith('f') and parts[index].endswith('q')


def shard_file_size(args, file_name):
    """Calculate shard file size for FASTA/Q files in raw or zipped format."""
    file_size = getsize(file_name)

    if args.get('gzip'):
        with gzip.open(file_name, 'rb') as zippy:
            file_size = zippy.seek(0, io.SEEK_END)
    elif args.get('bzip'):
        with bz2.open(file_name, 'rb') as zippy:
            file_size = zippy.seek(0, io.SEEK_END)

    if is_fastq_file(args, file_name):
        file_size /= 2  # Guessing that fastq files ~2x fasta files

    return file_size


def prefix_file(prefix, name):
    """Calculate the output path."""
    dir_, file_ = split(prefix)
    file_ += '.' if file_ and file_[-1] != '.' else ''
    return join(dir_, file_ + name)
