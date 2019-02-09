#!/usr/bin/env python3
"""Convert data in an atram sqlite database into fasta files."""

import sys
from os.path import exists, splitext
import argparse
import textwrap
import lib.db as db
import lib.blast as blast
import lib.util as util


def create_fasta_files(args):
    """Convert data in an atram sqlite database into fasta files."""
    if not exists(db.get_db_name(args['blast_db'])):
        sys.exit('Could not find the database.')

    with db.connect(args['blast_db'], check_version=True) as cxn:
        try:
            files = open_fasta_files(args, cxn)

            for rec in db.get_all_sequences(cxn):
                util.write_fasta_record(files[rec[1]], rec[0], rec[2], rec[1])

        finally:
            close_fasta_files(files)


def open_fasta_files(args, cxn):
    """Open one fasta file for each sequence end."""
    files = {}
    for end in [e[0] for e in db.get_sequence_ends(cxn)]:
        name = '{}{}.{}'.format(args['fasta_root'], end, args['fasta_ext'])
        files[end] = open(name, 'w')
    return files


def close_fasta_files(files):
    """Close all fasta files."""
    for file in files.values():
        file.close()


def parse_command_line():
    """Process command-line arguments."""
    description = """
        This will read through the aTRAM SQLite database and create fasta
        files. One for each end. So end 1 will be named <file-name>1.fasta
        etc. If there is no end in the DB (i.e. the DB was built with the
        --single-ends option) then the file name will be <file-name>.fasta.
        """
    parser = argparse.ArgumentParser(
        fromfile_prefix_chars='@',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(db.ATRAM_VERSION))

    parser.add_argument('-b', '--blast-db', '--sra', '--db', '--database',
                        required=True, metavar='DB',
                        help="""This needs to match the DB prefix you
                             entered for atram_preprocessor.py.""")

    parser.add_argument('-f', '--fasta', required=True,
                        help="""What to name the output fasta files without
                            then end indicator.""")

    args = vars(parser.parse_args())

    args['blast_db'] = blast.touchup_blast_db_names([args['blast_db']])[0]

    (args['fasta_root'], args['fasta_ext']) = splitext(args['fasta'])
    args['fasta_ext'] = args['fasta_ext'] if args['fasta_ext'] else 'fasta'

    return args


if __name__ == '__main__':
    ARGS = parse_command_line()
    create_fasta_files(ARGS)
