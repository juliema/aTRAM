#!/usr/bin/env python3
"""
Start the atram exon framer.

This wrapper module parses the input arguments and passes them to the module
that does the actual framing (core_framer.py).
"""

from os.path import join
from datetime import date
import argparse
import textwrap
import lib.db as db
import lib.util as util
import lib.core_framer as framer


def parse_command_line():
    """Process command-line arguments."""
    description = """
        This program will align contigs to a reference sequence and put them
        into the correct reading frame.
        """

    parser = argparse.ArgumentParser(
        fromfile_prefix_chars='@',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(db.ATRAM_VERSION))

    parser.add_argument(
        '-T', '--taxa', metavar='TAXA', required=True,
        help="""A text file of all of your taxon names.""")

    parser.add_argument(
        '-r', '--reference-genes', '--refs', metavar='FASTA', required=True,
        help="""Reference amino acid sequences in a FASTA file.""")

    parser.add_argument(
        '-a', '--assemblies-dir', metavar='PATH', required=True,
        help="""The path to the DNA contigs.""")

    parser.add_argument(
        '-m', '--min-length', metavar='LENGTH', default=100, type=int,
        help="""Remove contigs that are less than this length. The default is
            100.""")

    parser.add_argument(
        '-t', '--temp-dir', metavar='DIR',
        help="""Place temporary files in this directory. All files will be
            deleted after aTRAM completes. The directory must exist.""")

    parser.add_argument(
        '--keep-temp-dir', action='store_true',
        help="""This flag will keep the temporary files in the --temp-dir
        around for debugging.""")

    parser.add_argument(
        '-l', '--log-file',
        help="""Log file (full path). The default is
            "atram_framer_<date>.log".""")
    parser.add_argument(
        '--log-level', choices=['debug', 'info', 'error'], default='info',
        help="""Log messages of the given level (or above). 'debug' shows the
            most messages and 'error' shows the least. The default is
            'info'""")

    parser.add_argument(
        '-o', '--output-prefix',
        help="""This is the prefix of all of the output files. So you can
            identify different stitcher output file sets. You may include a
            directory as part of the prefix. The stitcher will add suffixes to
            differentiate output files.""")

    parser.add_argument(
        '-f', '--file-filter', default='*.fasta',
        help="""Use this to filter files in the assemblies directory. For
            example '*filtered*.fasta' will select all fasta files in the
            assemblies directory with the word filtered in them. The default
            is to select all fasta files in the assemblies directory
            '*.fasta'.""")

    parser.add_argument(
        '--reference-name', action='store_true',
        help="""Prepend the reference name to the final assembled gene name?
            if false the gene name in the reference file with just be the
            <taxon-name> if you select this then the assembled gene name
            will be <reference-name>.<taxon-name>.""")

    parser.add_argument(
        '--long-contig', type=float, default=0.7,
        help="""A long contig is considered to be this fraction [0-1] of the
            longest contig assembled by exonerate. The default is 0.7.""")

    args = parser.parse_args()

    util.temp_dir_exists(args.temp_dir)

    if not args.output_prefix:
        args.output_prefix = join(
            '.', 'atram_framer_' + date.today().isoformat())

    if not args.log_file and args.output_prefix[-1] == '/':
        args.log_file = join(
            args.output_prefix,
            'atram_framer_' + date.today().isoformat() + '.log')
    else:
        args.log_file = args.output_prefix + '.log'

    return args


if __name__ == '__main__':
    ARGS = parse_command_line()
    framer.frame(ARGS)
