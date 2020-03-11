#!/usr/bin/env python3
"""
Start the atram preprocessor.

This wrapper module parses the input arguments and passes them to the module
that does the actual preprocessing (core_preprocessor.py).
"""

import os
from os.path import join
from glob import glob
from itertools import chain
import argparse
import textwrap
from datetime import date
import lib.db as db
import lib.util as util
import lib.blast as blast
from lib.core_preprocessor import preprocess


def parse_command_line():
    """Process command-line arguments."""
    description = """
        This script prepares data for use by the atram.py
        script. It takes fasta or fastq files of paired-end (or
        single-end) sequence reads and creates a set of atram
        databases.

        You need to prepare the sequence read archive files so that the
        header lines contain only a sequence ID with the optional
        paired-end suffix at the end of the header line. The separator
        for the optional trailing paired-end suffix may be a space,
        a slash "/", a dot ".", or an underscore "_".

        For example:

            >DBRHHJN1:427:H9YYAADXX:1:1101:10001:77019/1
            GATTAA...
            >DBRHHJN1:427:H9YYAADXX:1:1101:10001:77019/2
            ATAGCC...
            >DBRHHJN1:427:H9YYAADXX:1:1101:10006:63769/2
            CGAAAA...
        """

    parser = argparse.ArgumentParser(
        fromfile_prefix_chars='@',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(db.ATRAM_VERSION))

    parser.add_argument(
        '--end-1', '-1', metavar='FASTA/Q', action='append',
        help="""Sequence read archive files that have only end 1 sequences. The
            sequence names do not need an end suffix, we will assume the suffix
            is always 1. The files are in fasta or fastq format. You may
            repeat this argument or use wildcards.
            """)

    parser.add_argument(
        '--end-2', '-2', metavar='FASTA/Q', action='append',
        help="""Sequence read archive files that have only end 2 sequences.
            The sequence names do not need an end suffix, we will assume the
            suffix is always 2. The files are in fasta or fastq format. You
            may repeat this argument or use wildcards.
            """)

    parser.add_argument(
        '--mixed-ends', '-m', metavar='FASTA/Q', action='append',
        help="""Sequence read archive files that have a mix of both end 1 and
            end 2 sequences (or single ends). The files are in fasta or fastq
            format. You may repeat this argument or use wildcards.
            """)

    parser.add_argument(
        '--single-ends', '-0', metavar='FASTA/Q', action='append',
        help="""Sequence read archive files that have only unpaired sequences.
            Any sequence suffix will be ignored. The files are in fasta or
            fastq format. You may repeat this argument or use wildcards.
            """)

    group = parser.add_argument_group('preprocessor arguments')

    blast_db = join('.', 'atram_' + date.today().isoformat())
    group.add_argument(
        '-b', '--blast-db', '--db', default=blast_db, metavar='DB',
        help="""This is the prefix of all of the blast database files. So you
            can identify different blast database sets. You may include a
            directory as part of the prefix. The default is "{}".
            """.format(blast_db))

    cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
    group.add_argument(
        '--cpus', '--processes', '--max-processes', type=int, default=cpus,
        help="""Number of CPU threads to use. On this machine the default is
        ("{}")""".format(cpus))

    group.add_argument(
        '-t', '--temp-dir', metavar='DIR',
        help="""Place temporary files in this directory. All files will be
            deleted after aTRAM completes. The directory must exist.""")

    group.add_argument(
        '--keep-temp-dir', action='store_true',
        help="""This flag will keep the temporary files in the --temp-dir
        around for debugging.""")

    group.add_argument(
        '-l', '--log-file',
        help="""Log file (full path). The default is to use the DB and program
            name to come up with a name like "<DB>_atram_preprocessor.log".""")
    group.add_argument(
        '--log-level', choices=['debug', 'info', 'error'], default='info',
        help="""Log messages of the given level (or above). 'debug' shows the
            most messages and 'error' shows the least. The default is
            'info'""")

    group.add_argument(
        '-s', '--shards', '--number', type=int, metavar='SHARDS',
        dest='shard_count',
        help="""Number of blast DB shards to create. The default is to have
            each shard contain roughly 250MB of sequence data.""")

    group.add_argument(
        '--path',
        help="""If makeblastdb is not in your $PATH then use this to prepend
            directories to your path.""")

    group.add_argument(
        '--fasta', action='store_true',
        help="""Are these fasta files? If you do not specify either --fasta or
            --fastq then aTRAM will guess the file type by looking at the last
            character of the file name.""")

    group.add_argument(
        '--fastq', action='store_true',
        help="""Are these fastq files? If you do not specify either --fasta or
            --fastq then aTRAM will guess the file type by looking at the last
            character of the file name.""")

    group.add_argument(
        '--gzip', action='store_true',
        help="""Are these gzip files?""")

    group.add_argument(
        '--bzip', action='store_true',
        help="""Are these bzip files?""")

    args = vars(parser.parse_args())

    # Prepend to PATH environment variable if requested
    if args['path']:
        os.environ['PATH'] = '{}:{}'.format(args['path'], os.environ['PATH'])

    all_files = []
    for ends in ['mixed_ends', 'end_1', 'end_2', 'single_ends']:
        if args.get(ends):
            end_files = [glob(p) for p in args[ends]]
            end_files = sorted(list(chain.from_iterable(end_files)))
            args[ends] = end_files
            all_files.extend(end_files)

    args['shard_count'] = blast.default_shard_count(args, all_files)

    blast.make_blast_output_dir(args['blast_db'])

    blast.find_program('makeblastdb')

    util.temp_dir_exists(args['temp_dir'])

    return args


if __name__ == '__main__':
    ARGS = parse_command_line()
    preprocess(ARGS)
