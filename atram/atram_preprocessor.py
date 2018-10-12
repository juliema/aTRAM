#!/bin/env python
"""
Start the atram preprocessor.

This wrapper module parses the input arguments and passes them to the module
that does the actual preprocessing (core_preprocessor.py).
"""

import os
import argparse
import textwrap
import tempfile
from datetime import date
import lib.core_preprocessor as core_preprocessor
import lib.db as db
import lib.blast as blast


def parse_command_line(temp_dir_default):
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
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(db.ATRAM_VERSION))

    parser.add_argument('--end-1', '-1', metavar='FASTA_or_FASTQ', nargs='+',
                        help='''Sequence read archive files that have only
                             end 1 sequences. The sequence names do not need an
                             end suffix, we will assume the suffix is always 1.
                             The files are in fasta or fastq format. You may
                             enter more than one file or you may use wildcards.
                             ''')

    parser.add_argument('--end-2', '-2', metavar='FASTA_or_FASTQ', nargs='+',
                        help='''Sequence read archive files that have only
                             end 2 sequences. The sequence names do not need an
                             end suffix, we will assume the suffix is always 2.
                             The files are in fasta or fastq format. You may
                             enter more than one file or you may use wildcards.
                             ''')

    parser.add_argument('--mixed-ends', '-m', metavar='FASTA_or_FASTQ',
                        nargs='+',
                        help='''Sequence read archive files that have a mix of
                             both end 1 and end 2 sequences (or single ends).
                             The files are in fasta or fastq format. You may
                             enter more than one file or you may use wildcards.
                             ''')

    parser.add_argument('--single-ends', '-0', metavar='FASTA_or_FASTQ',
                        nargs='+',
                        help='''Sequence read archive files that have only
                             unpaired sequences. Any sequence suffix will be
                             ignored. The files are in fasta or fastq format.
                             You may enter more than one file or you may use
                             wildcards.''')

    group = parser.add_argument_group('preprocessor arguments')

    blast_db = os.path.join('.', 'atram_' + date.today().isoformat())
    group.add_argument('-b', '--blast-db', '--output', '--db',
                       default=blast_db, metavar='DB',
                       help='''This is the prefix of all of the blast
                            database files. So you can identify
                            different blast database sets. You may include
                            a directory as part of the prefix. The default
                            is "{}".'''.format(blast_db))

    cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
    group.add_argument('--cpus', '--processes', '--max-processes',
                       type=int, default=cpus,
                       help='''Number of CPU threads to use. On this
                            machine the default is ("{}")'''.format(cpus))

    group.add_argument('-t', '--temp-dir', metavar='DIR',
                       help='''You may save intermediate files for debugging
                            in this directory. The directory must be empty.''')

    group.add_argument('-l', '--log-file',
                       help='''Log file (full path). The default is to use the
                            DB and program name to come up with a name like
                            "<DB>_atram_preprocessor.log"''')

    group.add_argument('-s', '--shards', '--number',
                       type=int, metavar='SHARDS',
                       dest='shard_count',
                       help='''Number of blast DB shards to create.
                            The default is to have each shard contain
                            roughly 250MB of sequence data.''')

    group.add_argument('--path',
                       help='''If blast or makeblastdb is not in your $PATH
                            then use this to prepend directories to your
                            path.''')

    group.add_argument('--sqlite-temp-dir', metavar='DIR',
                       help='''Use this directory to save temporary SQLITE3
                            files. This is a possible fix for "database or
                            disk is full" errors.''')

    args = vars(parser.parse_args())

    # Prepend to PATH environment variable if requested
    if args['path']:
        os.environ['PATH'] = '{}:{}'.format(args['path'], os.environ['PATH'])

    # Add an sqlite3 temporary directory if requested
    if args['sqlite_temp_dir']:
        os.environ['SQLITE_TMPDIR'] = args['sqlite_temp_dir']
    elif args['temp_dir']:
        os.environ['SQLITE_TMPDIR'] = args['temp_dir']

    # Setup temp dir
    if not args['temp_dir']:
        args['temp_dir'] = temp_dir_default
    else:
        os.makedirs(args['temp_dir'], exist_ok=True)

    all_files = []
    for ends in ['mixed_ends', 'end_1', 'end_2', 'single_ends']:
        if args.get(ends):
            all_files.extend([i for i in args[ends]])

    args['shard_count'] = blast.default_shard_count(
        args['shard_count'], all_files)

    blast.make_blast_output_dir(args['blast_db'])

    blast.find_program('makeblastdb')

    return args


if __name__ == '__main__':

    with tempfile.TemporaryDirectory(prefix='atram_') as TEMP_DIR_DEFAULT:
        ARGS = parse_command_line(TEMP_DIR_DEFAULT)
        core_preprocessor.preprocess(ARGS)
