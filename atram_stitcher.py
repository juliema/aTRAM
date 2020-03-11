#!/usr/bin/env python3
"""
Start the atram exon stitcher.

This wrapper module parses the input arguments and passes them to the module
that does the actual stitching (core_stitcher.py).
"""

from os.path import join
from datetime import date
import argparse
import textwrap
import lib.db as db
import lib.log as log
import lib.util as util
import lib.core_stitcher as stitcher


def parse_command_line():
    """Process command-line arguments."""
    description = """
        This program will find and stitch together exons from targeted
        assemblies using amino acid targets and DNA assemblies.
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
        '-O', '--overlap', type=int, default=10,
        help="""Contigs must overlap by this many codons before it is
            considered a real overlap.""")

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
            "atram_stitcher_<date>.log".""")
    parser.add_argument(
        '--log-level', choices=['debug', 'info', 'error'], default='info',
        help="""Log messages of the given level (or above). 'debug' shows the
            most messages and 'error' shows the least. The default is
            'info'""")

    parser.add_argument(
        '-i', '--iterations', type=int, default=2, metavar='N',
        help="""The number of times to run the main stitcher loop. This
            must be either 1 or 2, the default is 2.""")

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

    args = parser.parse_args()

    util.temp_dir_exists(args.temp_dir)

    if not args.output_prefix:
        args.output_prefix = join(
            '.', 'atram_stitcher_' + date.today().isoformat())

    if not args.log_file and args.output_prefix[-1] == '/':
        args.log_file = join(
            args.output_prefix,
            'atram_stitcher_' + date.today().isoformat() + '.log')
    else:
        args.log_file = args.output_prefix + '.log'

    if 1 > args.iterations > 2:
        log.fatal('The iterations must be either 1 or 2.')

    return args


if __name__ == '__main__':
    ARGS = parse_command_line()
    stitcher.stitch(ARGS)
