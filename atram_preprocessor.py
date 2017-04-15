"""Format the data so that atram can use it later. It takes sequence read
archive files and converts them into blast and sqlite3 databases.
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


def run(args):
    """Run the preprocessor."""

    log.setup(args)

    db_conn = db.connect(args.blast_db)
    db.create_sequences_table(db_conn)
    load_seqs(args, db_conn)
    db.create_sequences_index(db_conn)

    shard_list = assign_seqs_to_shards(args, db_conn)
    create_blast_dbs(args, shard_list)

    db_conn.close()


def load_seqs(args, db_conn):
    """A hand rolled version of "Bio.SeqIO". It's faster because we can
    take shortcuts due to its limited use.

    We're using a very simple state machine on lines to do the parsing.
        1) header      (exactly 1 line)  Starts with a '>' or an '@'
        2) sequence    (1 or more lines) Starts with a letter
        3) fastq stuff (0 or more lines) Starts with a '+', Look for header
        4) Either go back to 1 or end
    """

    for file_name in args.sra_files:

        log.info('Loading "%s" into sqlite database' % file_name)

        with open(file_name) as sra_file:
            batch = []      # The batch of records to insert

            seq = ''        # The sequence string. A DB field
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
                        batch.append((seq_name, seq_end, seq))

                    seq = ''        # Reset the sequence
                    is_seq = True   # Reset the state flag

                    # Get data from the header
                    match = blast.PARSE_HEADER.match(line)
                    seq_name = match.group(1) if match else ''
                    seq_end = match.group(2) if match else ''

                # Handle fastq stuff, so skip these lines
                elif line[0] == '+':
                    is_seq = False

                # Handle sequence lines
                elif line[0].isalpha() and is_seq:
                    seq += line.rstrip()  # ''.join([strs])

                if len(batch) >= db.BATCH_SIZE:
                    db.insert_sequences_batch(db_conn, batch)
                    batch = []

            if seq:
                batch.append((seq_name, seq_end, seq))

            db.insert_sequences_batch(db_conn, batch)


def assign_seqs_to_shards(args, db_conn):
    """Put the sequences into the DB shards. What we doing is dividing all
    of the input sequences into shard_count bucket of sequences. If there
    are two ends of a sequence we have to make sure that both ends (1 & 2)
    wind up in the same shard. These shards will then be turned into blast
    databases.

    This will build up an array of "LIMIT len OFFSET start" parameters for
    SQL SELECT statements that are used for building the shard fasta files
    that get input into the makeblastdb statements.
    """

    log.info('Assigning sequences to shards')

    total = db.get_sequence_count(db_conn)

    # This creates a list of roughly equal partition indexes
    offsets = np.linspace(0, total, dtype=int, num=args.shard_count + 1)

    # Checking to make sure we don't split up the ends of a sequence
    for i in range(1, len(offsets) - 1):

        # Get the first two sequences of a partition
        first, second = db.get_two_sequences(db_conn, offsets[i])

        # If both have the same name then both sequences will be in the
        # same partition and we're done. If they have different names then
        # we know that second sequence either starts a new pair or is a
        # singleton, so we can safely start the partition at the second one
        if first != second:
            offsets[i] += 1

    # Get the length of each partition
    limits = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]

    return list(zip(limits, offsets[:-1]))


def create_blast_dbs(args, shard_list):
    """Assign processes to make the blast DBs. One process for each blast
    DB shard.
    """

    log.info('Making blast DBs')

    with multiprocessing.Pool(processes=args.cpus) as pool:
        results = []
        for idx, params in enumerate(shard_list, 1):
            results.append(pool.apply_async(
                create_blast_db, (args.blast_db, args.temp_dir, params, idx)))

        _ = [result.get() for result in results]

    log.info('Finished making blast DBs')


def create_blast_db(blast_db, temp_dir, shard_params, shard_index):
    """Create a blast DB from the shard. We fill a fasta file with the
    appropriate sequences and hand things off to the makeblastdb program.
    """
    # NOTE: Because this is called in a child process, the address space is not
    # shared with the parent (caller) hence we cannot share object variables.

    shard_path = blast.shard_path(blast_db, shard_index)
    fasta_name = '{}_{:03d}.fasta'.format(os.path.basename(sys.argv[0][:-3]),
                                          shard_index)
    fasta_path = os.path.join(temp_dir, fasta_name)

    with open(fasta_path, 'w') as fasta_file:
        fill_blast_fasta(blast_db, fasta_file, shard_params)
        blast.create_db(temp_dir, fasta_path, shard_path)


def fill_blast_fasta(blast_db, fasta_file, shard_params):
    """Fill the fasta file used as input into blast with shard sequences from
    the sqlite3 DB. We use the shard partitions passed in to determine
    which sequences to get for this shard.
    """
    # NOTE: Because this is called in a child process, the address space is not
    # shared with the parent (caller) hence we cannot share object variables.

    db_conn = db.connect(blast_db)

    limit, offset = shard_params

    for row in db.get_sequences_in_shard(db_conn, limit, offset):
        fasta_file.write('>{}{}\n'.format(row[0], row[1]))
        fasta_file.write('{}\n'.format(row[2]))

    db_conn.close()


def parse_command_line(temp_dir):
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

    parser.add_argument(dest='sra_files', metavar='FASTA_FILE', nargs='+',
                        help='Sequence read archive files in fasta or '
                             'fastq format. You may enter more than '
                             'one file and you may use wildcards.')

    parser.add_argument('--version', action='version', version='%(prog)s 2.0')

    group = parser.add_argument_group('preprocessor arguments')

    blast_db = os.path.join('.', 'atram_' + date.today().isoformat())
    group.add_argument('-b', '--blast-db', '--output', '--db',
                       default=blast_db, metavar='DB',
                       help=('This is the prefix of all of the blast '
                             'database files. So you can identify '
                             'different blast database sets. You may include '
                             'a directory as part of the prefix. The default '
                             'is "{}".').format(blast_db))

    cpus = os.cpu_count() - 4 if os.cpu_count() > 4 else 1
    group.add_argument('-c', '--cpus', '--processes', '--max-processes',
                       type=int, default=cpus,
                       help=('Number of cpus to use. Defaults to: Total CPUS '
                             ' - 4 = "{}"').format(cpus))

    group.add_argument('-t', '--temp-dir', metavar='DIR',
                       help='You may save intermediate files for '
                            'debugging in this directory.')

    group.add_argument('-l', '--log-file',
                       help='Log file (full path). The default is to use the '
                            'DB and program name to come up with a name like '
                            '"DB_atram_preprocessor.log"')

    group.add_argument('-s', '--shards', '--number',
                       type=int, metavar='SHARDS',
                       dest='shard_count',
                       help='Number of blast DB shards to create. '
                            'The default is to have each shard contain '
                            'roughly 250MB of sequence data.')

    args = parser.parse_args()

    # Set default --shard-count
    if not args.shard_count:
        total_fasta_size = 0
        for file_name in args.sra_files:
            file_size = os.path.getsize(file_name)
            if file_name.lower().endswith('q'):
                file_size /= 2  # Guessing that fastq files ~2x fasta files
            total_fasta_size += file_size
        args.shard_count = int(total_fasta_size / 2.5e8)
        args.shard_count = args.shard_count if args.shard_count else 1

    # Make output directory
    output_dir = os.path.dirname(args.blast_db)
    if output_dir and output_dir not in ['.']:
        os.makedirs(output_dir, exist_ok=True)

    # Set default log file name
    if not args.log_file:
        args.log_file = '{}.{}.log'.format(
            args.blast_db, os.path.basename(sys.argv[0][:-3]))

    # Make temp directory
    if args.temp_dir:
        os.makedirs(args.temp_dir, exist_ok=True)
    else:
        args.temp_dir = temp_dir

    return args


if __name__ == '__main__':

    with tempfile.TemporaryDirectory(prefix='atram_') as TEMP_DIR:
        ARGS = parse_command_line(TEMP_DIR)
        run(ARGS)
