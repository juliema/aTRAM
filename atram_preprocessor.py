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
from datetime import date
import numpy as np
import lib.db as db
import lib.log as log
import lib.blast as blast

__all__ = ('preprocess', )


def preprocess(args):
    """Run the preprocessor."""
    log.setup(args['log_file'], args['blast_db'])

    with db.connect(args['blast_db'], bulk_mode=True) as db_conn:
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
    for (arg, clamp) in [('mixed_ends', None), ('end_1', '1'),
                         ('end_2', '2'), ('single_ends', '')]:
        if args.get(arg):
            for file_name in args[arg]:
                load_one_file(db_conn, file_name, clamp)


def load_one_file(db_conn, file_name, seq_end_clamp=None):
    """Load sequences from a fasta/fastq file into the atram database."""
    # A hand rolled version of "Bio.SeqIO". It's faster because we can take
    # shortcuts due to its limited functionality.

    # We're using a very simple state machine on lines to do the parsing.
    #     1) header      (exactly 1 line)  Starts with a '>' or an '@'
    #     2) sequence    (1 or more lines) Starts with a letter
    #     3) fastq stuff (0 or more lines) Starts with a '+' on the 1st line
    #     4) Either go back to 1 or end
    log.info('Loading "{}" into sqlite database'.format(file_name))

    with open(file_name) as sra_file:
        batch = []      # The batch of records to insert

        seq = []        # It will become the sequence string. A DB field
        seq_end = ''    # Which end? 1 or 2. A DB field
        seq_name = ''   # Name from the fasta file. A DB field
        is_seq = True   # The state machine's only flag for state

        for line in sra_file:

            # Skip blank lines
            if not line:
                pass

            # Handle a header line
            elif line[0] in ['>', '@']:
                if seq:  # Append last record to the batch?
                    batch.append((seq_name, seq_end, ''.join(seq)))

                seq = []        # Reset the sequence
                is_seq = True   # Reset the state flag

                # Get data from the header
                match = blast.PARSE_HEADER.match(line)
                seq_name = match.group(1) if match else line[1:].strip()

                # Sequence end either equals its clamp value or we get it
                # from the header line
                end = match.group(2)
                if seq_end_clamp is None:
                    seq_end = end if end else ''
                else:
                    seq_end = end if end else seq_end_clamp
                    # seq_end = seq_end_clamp

            # Handle fastq stuff, so skip these lines
            elif line[0] == '+':
                is_seq = False

            # Handle sequence lines
            elif line[0].isalpha() and is_seq:
                seq.append(line.rstrip())

            if len(batch) >= db.BATCH_SIZE:
                db.insert_sequences_batch(db_conn, batch)
                batch = []

        if seq:
            batch.append((seq_name, seq_end, ''.join(seq)))

        db.insert_sequences_batch(db_conn, batch)


def assign_seqs_to_shards(db_conn, shard_count):
    """
    Assign sequences to blast DB shards.

    We are dividing all of the input sequences into shard_count buckets of
    sequences. If there are two ends of a sequence we have to make sure that
    both ends (1 & 2) wind up in the same shard.

    What we do is look at two sequence names at the shard boundary. If they are
    the same then the chard boundary is fine where it is. But if they are
    different then the second sequence is either the start of a sequence pair
    or a singleton so we can safely start the shard at the second sequence.

    Note: This will only work for sequence pairs. Which is all we care about.
    """
    log.info('Assigning sequences to shards')

    total = db.get_sequence_count(db_conn)

    # This creates a list of roughly equal partition indexes
    offsets = np.linspace(0, total, dtype=int, num=shard_count + 1)

    # Checking to make sure we don't split up the ends of a sequence
    for i in range(1, len(offsets) - 1):

        # Get the first two sequences of a possible partition
        first, second = db.get_shard_cut_pair(db_conn, offsets[i])

        # If both have the same name then both sequences will be in the
        # same partition and we're done. If they have different names then
        # we know that second sequence either starts a new pair or is a
        # singleton, so we can safely start the partition at the second one
        if first != second:
            offsets[i] += 1

    # Get the length of each partition
    limits = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]

    return list(zip(limits, offsets[:-1]))


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

    parser.add_argument('--mixed-ends', '-m', metavar='FASTA', nargs='+',
                        help='''Sequence read archive files that have a mix of
                             both end 1 and end 2 sequences. The sequence names
                             MUST have an end suffix like "/1" or "_2". The
                             files are in fasta or fastq format. You may enter
                             more than one file or you may use wildcards.''')

    parser.add_argument('--end-1', '-1', metavar='FASTA', nargs='+',
                        help='''Sequence read archive files that have only
                             end 1 sequences. The sequence names do not need an
                             end suffix, we will assume the suffix is 1 if it
                             is missing. The files are in fasta or fastq
                             format. You may enter more than one file or you
                             may use wildcards.''')

    parser.add_argument('--end-2', '-2', metavar='FASTA', nargs='+',
                        help='''Sequence read archive files that have only
                             end 2 sequences. The sequence names do not need an
                             end suffix, we will assume the suffix is 2 if it
                             is missing. The files are in fasta or fastq
                             format. You may enter more than one file or you
                             may use wildcards.''')

    # parser.add_argument('--single-ends', '-S', metavar='FASTA', nargs='+',
    #                     help='''Sequence read archive files that have only
    #                          unpaired sequences. Any sequence end suffixes
    #                          will
    #                          be ignored. The files are in fasta or fastq
    #                          format. You may enter more than one file or you
    #                          may use wildcards.''')

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(db.ATRAM_VERSION))

    group = parser.add_argument_group('preprocessor arguments')

    blast_db = os.path.join('.', 'atram_' + date.today().isoformat())
    group.add_argument('-b', '--blast-db', '--output', '--db',
                       default=blast_db, metavar='DB',
                       help='''This is the prefix of all of the blast
                            database files. So you can identify
                            different blast database sets. You may include
                            a directory as part of the prefix. The default
                            is "{}".'''.format(blast_db))

    cpus = os.cpu_count() - 4 if os.cpu_count() > 4 else 1
    group.add_argument('-c', '--cpus', '--processes', '--max-processes',
                       type=int, default=cpus,
                       help='''Number of cpus to use. Defaults to: Total CPUS
                            - 4 = "{}"'''.format(cpus))

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

    args = vars(parser.parse_args())

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

    return args


if __name__ == '__main__':

    with tempfile.TemporaryDirectory(prefix='atram_') as TEMP_DIR_DEFAULT:
        ARGS = parse_command_line(TEMP_DIR_DEFAULT)
        preprocess(ARGS)
