"""
Format the data so that atram can use it later in atram itself.

It takes sequence read archive (SRA) files and converts them into coordinated
blast and sqlite3 databases.
"""

import os
import sys
import argparse
import textwrap
import tempfile
import multiprocessing
from shutil import which
from datetime import date
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
import lib.db as db
import lib.log as log
import lib.blast as blast

__all__ = ('preprocess', )


def preprocess(args):
    """Run the preprocessor."""
    log.setup(args['log_file'], args['blast_db'])

    with db.connect(args['blast_db'], clean=True) as db_conn:
        db.create_metadata_table(db_conn)

        db.create_sequences_table(db_conn)
        load_seqs(args, db_conn)

        log.info('Creating an index for the sequence table')
        db.create_sequences_index(db_conn)

        shard_list = assign_seqs_to_shards(db_conn, args['shard_count'])

    create_all_blast_shards(args, shard_list)


def load_seqs(args, db_conn):
    """Load sequences from a fasta/fastq files into the atram database."""
    # We have to clamp the end suffix depending on the file type.
    for (arg, clamp) in [('mixed_ends', ''), ('end_1', '1'),
                         ('end_2', '2'), ('single_ends', '')]:
        if args.get(arg):
            for file_name in args[arg]:
                load_one_file(db_conn, file_name, arg, clamp)


def load_one_file(db_conn, file_name, arg, seq_end_clamp=''):
    """Load sequences from a fasta/fastq file into the atram database."""
    log.info('Loading "{}" into sqlite database'.format(file_name))

    is_fastq = file_name.lower().endswith('q')
    parser = FastqGeneralIterator if is_fastq else SimpleFastaParser

    with open(file_name) as sra_file:
        batch = []

        for rec in parser(sra_file):
            title = rec[0].strip()
            seq = rec[1]

            match = blast.PARSE_HEADER.match(title)
            if match.group(2):
                seq_name = match.group(1)
                if arg == 'mixed_ends':
                    seq_end = match.group(2)
                else:
                    seq_end = seq_end_clamp
            else:
                seq_name = title
                seq_end = seq_end_clamp

            batch.append((seq_name, seq_end, seq))

            if len(batch) >= db.BATCH_SIZE:
                db.insert_sequences_batch(db_conn, batch)
                batch = []

        db.insert_sequences_batch(db_conn, batch)


def assign_seqs_to_shards(db_conn, shard_count):
    """Assign sequences to blast DB shards."""
    log.info('Assigning sequences to shards')

    total = db.get_sequence_count(db_conn)
    offsets = np.linspace(0, total - 1, dtype=int, num=shard_count + 1)
    cuts = [db.get_shard_cut(db_conn, offset) for offset in offsets]

    # Make sure the last sequence gets included
    cuts[-1] = cuts[-1] + 'z'

    # Now organize the list into pairs of sequence names
    pairs = [(cuts[i - 1], cuts[i]) for i in range(1, len(cuts))]

    return pairs


def create_all_blast_shards(args, shard_list):
    """
    Assign processes to make the blast DBs.

    One process for each blast DB shard.
    """
    log.info('Making blast DBs')

    with multiprocessing.Pool(processes=args['cpus']) as pool:
        results = []
        for idx, shard_params in enumerate(shard_list, 1):
            results.append(pool.apply_async(
                create_one_blast_shard,
                (args, shard_params, idx)))

        _ = [result.get() for result in results]  # noqa

    log.info('Finished making blast DBs')


def create_one_blast_shard(args, shard_params, shard_index):
    """
    Create a blast DB from the shard.

    We fill a fasta file with the appropriate sequences and hand things off
    to the makeblastdb program.
    """
    shard = '{}.{:03d}.blast'.format(args['blast_db'], shard_index)
    fasta_name = '{}_{:03d}.fasta'.format(os.path.basename(sys.argv[0][:-3]),
                                          shard_index)
    fasta_path = os.path.join(args['temp_dir'], fasta_name)

    fill_blast_fasta(args['blast_db'], fasta_path, shard_params)

    blast.create_db(args['temp_dir'], fasta_path, shard)


def fill_blast_fasta(blast_db, fasta_path, shard_params):
    """
    Fill the fasta file used as input into blast.

    Use sequences from the sqlite3 DB. We use the shard partitions passed in to
    determine which sequences to get for this shard.
    """
    with db.connect(blast_db) as db_conn:
        limit, offset = shard_params

        with open(fasta_path, 'w') as fasta_file:
            for row in db.get_sequences_in_shard(db_conn, limit, offset):
                seq_end = '/{}'.format(row[1]) if row[1] else ''
                fasta_file.write('>{}{}\n'.format(row[0], seq_end))
                fasta_file.write('{}\n'.format(row[2]))


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

    args = vars(parser.parse_args())

    # Prepend to PATH environment variable if requested
    if args['path']:
        os.environ['PATH'] = '{}:{}'.format(args['path'], os.environ['PATH'])

    # Setup temp dir
    if not args['temp_dir']:
        args['temp_dir'] = temp_dir_default
    else:
        os.makedirs(args['temp_dir'], exist_ok=True)

    all_files = []
    for arg in ['mixed_ends', 'end_1', 'end_2', 'single_ends']:
        if args.get(arg):
            all_files.extend([i for i in args[arg]])

    args['shard_count'] = blast.default_shard_count(
        args['shard_count'], all_files)

    blast.make_blast_output_dir(args['blast_db'])

    find_programs()

    return args


def find_programs():
    """Make sure we can find the programs needed by the assembler and blast."""
    if not which('makeblastdb'):
        err = ('We could not find the programs "makeblastdb". You either need '
               'to install it or you need adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        sys.exit(err)


if __name__ == '__main__':

    with tempfile.TemporaryDirectory(prefix='atram_') as TEMP_DIR_DEFAULT:
        ARGS = parse_command_line(TEMP_DIR_DEFAULT)
        preprocess(ARGS)
