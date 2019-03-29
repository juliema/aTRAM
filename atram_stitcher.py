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
import lib.util as util
from lib.core_stitcher import Sticher


def parse_command_line():
    """Process command-line arguments."""
    description = """
        This program will find and stitch together exons from targeted
        assemblies using amino acid targets and DNA assemblies.
        """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(db.ATRAM_VERSION))

    parser.add_argument(
        '-T', '--taxa-list', metavar='TAXA', required=True,
        help="""A text file of all your taxon names.""")

    parser.add_argument(
        '-r', '--reference-genes', metavar='FASTA', required=True,
        help="""Reference amino acid sequences in a FASTA file.""")

    parser.add_argument(
        '-a', '--assemblies-dir', metavar='PATH', required=True,
        help="""The path to the target assemblies directory.""")

    parser.add_argument(
        '-O', '--overlap', default=10,
        help=""".""")

    parser.add_argument(
        '-t', '--temp-dir', metavar='DIR',
        help="""Place temporary files in this directory. All files will be
            deleted after aTRAM completes. The directory must exist.""")

    parser.add_argument(
        '--keep-temp-dir', action='store_true',
        help="""This flag will keep the temporary files in the --temp-dir
        around for debugging.""")

    log_file = join('.', 'atram_stitcher_' + date.today().isoformat() + '.log')
    parser.add_argument(
        '-l', '--log-file', default=log_file,
        help="""Log file (full path). The default is
            "atram_stitcher_<date>.log".""")

    parser.add_argument(
        '-o', '--output-prefix', required=True,
        help="""This is the prefix of all of the output files. So you can
            identify different stitcher output file sets. You may include a
            directory as part of the prefix. The stitcher will add suffixes to
            differentiate output files.""")

    args = vars(parser.parse_args())

    util.temp_dir_exists(args['temp_dir'])

    return args


if __name__ == '__main__':
    ARGS = parse_command_line()
    Sticher(ARGS).stitch()
