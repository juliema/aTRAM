"""Format the data so that atram can use it later. It takes given sequence
read archive files and converts them.
"""

import logging
import subprocess
import multiprocessing
import numpy as np
import lib.db as db
import lib.bio as bio
from lib.filer import Filer
from lib.configure import Configure


class AtramPreprocessor:
    """Class to handle getting the data into a format that atram can use."""

    def __init__(self):
        self.config = None
        self.db_conn = None
        self.shard_list = None

    def run(self):
        """The main method."""

        self.config = Configure().parse_command_line(
            description="""This script prepares data for use by the atram.py
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
                """,
            args='sra_files db_prefix work_dir shard_count cpus')

        filer = Filer(self.config.work_dir, self.config.db_prefix)
        filer.log_setup()

        self.db_conn = db.connect(filer)
        db.create_sequences_table(self.db_conn)
        self.load_seqs()
        db.create_sequences_index(self.db_conn)

        self.shard_list = self.assign_seqs_to_shards()
        self.create_blast_dbs()

        self.db_conn.close()

    def load_seqs(self):
        """A hand rolled version of "Bio.SeqIO". It's faster because we can
        take shortcuts due to its limited use.

        We're using a very simple state machine on lines to do the parsing.
            1) header      (exactly 1 line)  Starts with a '>' or an '@'
            2) sequence    (1 or more lines) Starts with a letter
            3) fastq stuff (0 or more lines) Starts with a '+', Wait for header
            4) Either go back to 1 or end
        """

        for file_name in self.config.sra_files:

            logging.info('Loading "%s" into sqlite database', file_name)

            with open(file_name) as sra_file:
                sequence_batch = []  # The batch of records to insert

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
                            sequence_batch.append((seq_name, seq_end, seq))

                        seq = ''        # Reset the sequence
                        is_seq = True   # Reset the state flag

                        # Get data from the header
                        match = bio.PARSE_HEADER.match(line)
                        seq_name = match.group(1) if match else ''
                        seq_end = match.group(2) if match else ''

                    # Handle fastq stuff, so skip these lines
                    elif line[0] == '+':
                        is_seq = False

                    # Handle sequence lines
                    elif line[0].isalpha() and is_seq:
                        seq += line.rstrip()  # ''.join([strs])

                    if len(sequence_batch) >= db.BATCH_SIZE:
                        db.insert_sequences_batch(self.db_conn, sequence_batch)
                        sequence_batch = []

                if seq:
                    sequence_batch.append((seq_name, seq_end, seq))

                db.insert_sequences_batch(self.db_conn, sequence_batch)

    def assign_seqs_to_shards(self):
        """Put the sequences into the DB shards. What we doing is dividing all
        of the input sequences into shard_count bucket of sequences. If there
        are two ends of a sequence we have to make sure that both ends (1 & 2)
        wind up in the same shard. These shards will then be turned into blast
        databases.

        This will build up an array of "LIMIT len OFFSET start" parameters for
        SQL SELECT statements that are used for building the shard fasta files
        that get input into the makeblastdb statements.
        """

        logging.info('Assigning sequences to shards')

        total = db.get_sequence_count(self.db_conn)

        # This creates a list of roughly equal partition indexes
        offsets = np.linspace(0, total, dtype=int,
                              num=self.config.shard_count + 1)

        # Checking to make sure we don't split up the ends of a sequence
        for i in range(1, len(offsets) - 1):

            # Get the first two sequences of a partition
            first, second = db.get_two_sequences(self.db_conn, offsets[i])

            # If both have the same name then both sequences will be in the
            # same partition and we're done. If they have different names then
            # we know that second sequence either starts a new pair or is a
            # singleton, so we can safely start the partition at the second one
            if first != second:
                offsets[i] += 1

        # Get the length of each partition
        limits = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]

        return list(zip(limits, offsets[:-1]))

    def create_blast_dbs(self):
        """Assign processes to make the blast DBs. One process for each blast
        DB shard.
        """

        logging.info('Making blast DBs')

        with multiprocessing.Pool(processes=self.config.cpus) as pool:
            results = []
            for shard_index, shard_params in enumerate(self.shard_list):
                results.append(pool.apply_async(
                    create_blast_db,
                    (self.config.work_dir, self.config.db_prefix,
                     shard_params, shard_index)))

            _ = [result.get() for result in results]

        logging.info('Finished making blast DBs')


def create_blast_db(work_dir, db_prefix, shard_params, shard_index):
    """Create a blast DB from the shard. We fill a fasta file with the
    appropriate sequences and hand things off to the makeblastdb program.

    Because this is called in a child process, the address space is not
    shared with the parent (caller) hence we cannot use instance variables.
    """

    filer = Filer(work_dir=work_dir, db_prefix=db_prefix)
    blast_db = filer.blast_shard_name(shard_index)

    with filer.temp_file() as fasta_file:
        fill_blast_fasta(fasta_file, filer, shard_params)

        cmd = 'makeblastdb -dbtype nucl -in {} -out {}'
        cmd = cmd.format(fasta_file.name, blast_db)
        subprocess.check_call(cmd, shell=True)


def fill_blast_fasta(fasta_file, filer, shard_params):
    """Fill the fasta file used as input into blast with shard sequences from
    the sqlite3 DB. We use the shard partitions passed in to determine
    which sequences to get for this shard.

    Because this is called in a child process, the address space is not
    shared with the parent (caller) hence we cannot use instance variables.
    """

    db_conn = db.connect(filer)
    result = db.get_sequences_in_shard(
        db_conn, shard_params[0], shard_params[1])
    for row in result:
        fasta_file.write('>{}{}\n'.format(row[0], row[1]))
        fasta_file.write('{}\n'.format(row[2]))
    db_conn.close()


if __name__ == '__main__':

    AtramPreprocessor().run()
