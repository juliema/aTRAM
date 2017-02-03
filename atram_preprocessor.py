"""
Format the data so that atram can use it later. It takes given sequence read
archive files and converts them.
"""

import re
import sqlite3
import logging
import subprocess
import multiprocessing
import numpy as np
from lib.filer import Filer
from lib.configure import Configure


class AtramPreprocessor:
    """Class to handle getting the data into a format that atram can use."""

    batch_size = 1e7  # How many sequence records to insert at a time

    # Try to get the sequence name and which end it is from the fasta header
    parse_header = re.compile(r'^ [>@] \s* ( .* ) ( [\s\/_] [12] )',
                              re.VERBOSE)

    def __init__(self):
        self.config = None
        self.db_conn = None
        self.shard_list = None

    def preprocess(self):
        """The main method."""

        self.config = Configure().parse_command_line(
            description="""
                This script prepares data for use by the atram.py script.
                It takes fasta or fastq files of paired-end (or single-end)
                sequence reads and creates a set of atram databases.
                """,
            args='sra_files file_prefix work_dir shard_count cpus')

        filer = Filer(self.config.work_dir, self.config.file_prefix)
        filer.log_setup()

        self.db_conn = self.connect_db(filer)
        self.create_table()
        self.load_seqs()
        self.create_index()

        self.shard_list = self.assign_seqs_to_shards()
        self.create_blast_dbs()

        self.db_conn.close()

    @staticmethod
    def connect_db(filer):
        """
        Setup the DB for our processing needs and return a DB connection.

        Because this is called in a child process, the address space is not
        shared with the parent (caller) and we cannot use instance variables.
        """

        db_path = filer.db_file_name()
        db_conn = sqlite3.connect(db_path)
        db_conn.execute("PRAGMA page_size = {}".format(2**16))
        db_conn.execute("PRAGMA journal_mode = 'off'")
        db_conn.execute("PRAGMA synchronous = 'off'")
        return db_conn

    def create_table(self):
        """Reset the DB. Delete the tables and recreate them."""

        logging.info('Creating sqlite tables')
        self.db_conn.execute('''DROP INDEX IF EXISTS seq_names''')
        self.db_conn.execute('''DROP TABLE IF EXISTS sequences''')
        sql = '''
            CREATE TABLE IF NOT EXISTS sequences
            (seq_name TEXT, seq_end TEXT, seq TEXT)
            '''
        self.db_conn.execute(sql)

    def create_index(self):
        """
        Create indexes after we build the table. This speeds up the program
        significantly.
        """

        logging.info('Creating sqlite indices')
        self.db_conn.execute('''
            CREATE INDEX IF NOT EXISTS seq_names ON sequences (seq_name)
            ''')

    def load_seqs(self):
        """
        A hand rolled version of "Bio.SeqIO". It's faster because we can take
        shortcuts due to its limited use.

        We're using a very simple state machine on lines to do the parsing.
            1) header      (exactly 1 line)  Starts with a '>' or an '@'
            2) sequence    (1 or more lines) Starts with a letter
            3) fastq stuff (0 or more lines) Starts with a '+', Wait for header
            4) Either go back to 1 or end
        """

        for file_name in self.config.sra_files:

            logging.info('Loading "%s" into sqlite database', file_name)

            with open(file_name) as sra_file:
                record_batch = []  # The batch of records to insert

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
                            record_batch.append((seq_name, seq_end, seq))

                        seq = ''        # Reset the sequence
                        is_seq = True   # Reset the state flag

                        # Get data from the header
                        match = AtramPreprocessor.parse_header.match(line)
                        seq_name = match.group(1) if match else ''
                        seq_end = match.group(2) if match else ''

                    # Handle fastq stuff, so skip these lines
                    elif line[0] == '+':
                        is_seq = False

                    # Handle sequence lines
                    elif line[0].isalpha() and is_seq:
                        seq += line.rstrip()  # ''.join([strs])

                    if len(record_batch) >= AtramPreprocessor.batch_size:
                        self.insert_record_batch(record_batch)
                        record_batch = []

                if seq:
                    record_batch.append((seq_name, seq_end, seq))

                self.insert_record_batch(record_batch)

    def insert_record_batch(self, record_batch):
        """Insert a batch of sequence records into the sqlite database."""

        if record_batch:
            sql = '''
                INSERT INTO sequences (seq_name, seq_end, seq) VALUES (?, ?, ?)
                '''
            self.db_conn.executemany(sql, record_batch)
            self.db_conn.commit()

    def assign_seqs_to_shards(self):
        """
        Put the sequences into the DB shards. What we doing is dividing all of
        the input sequences into shard_count bucket of sequences. If there are
        two ends of a sequence we have to make sure that both ends (1 & 2) wind
        up in the same shard. These shards will then be turned into blast
        databases.

        This will build up an array of "LIMIT len OFFSET start" parameters for
        SQL SELECT statements that are used for building the shard fasta files
        that get input into the makeblastdb statements.
        """

        logging.info('Assigning sequences to shards')

        result = self.db_conn.execute('SELECT COUNT(*) FROM sequences')
        total = result.fetchone()[0]

        # This creates a list of roughly equal partition indexes
        offsets = np.linspace(0, total, dtype=int,
                              num=self.config.shard_count + 1)

        # Checking to make sure we don't split up the ends of a sequence
        for i in range(1, len(offsets) - 1):

            # Get the first two sequences of a partition
            sql = '''
                SELECT seq_name FROM sequences
                ORDER BY seq_name LIMIT 2 OFFSET {}
                '''
            result = self.db_conn.execute(sql.format(offsets[i]))
            first = result.fetchone()[0]
            second = result.fetchone()[0]

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
        """
        Assign processes to make the blast DBs. One process for each blast DB
        shard.
        """

        logging.info('Making blast DBs')

        config = self.config

        with multiprocessing.Pool(processes=self.config.cpus) as pool:
            results = []
            for shard_index, shard_params in enumerate(self.shard_list):
                results.append(pool.apply_async(
                    self.create_blast_db, (config.work_dir, config.file_prefix,
                                           shard_params, shard_index)))

            _ = [result.get() for result in results]

        logging.info('Finished making blast DBs')

    @staticmethod
    def create_blast_db(work_dir, file_prefix, shard_params, shard_index):
        """
        Create a blast DB from the shard. We fill a fasta file with the
        appropriate sequences and hand things off to the makeblastdb program.

        Because this is called in a child process, the address space is not
        shared with the parent (caller) and we cannot use instance variables.
        """

        filer = Filer(work_dir=work_dir, file_prefix=file_prefix)
        blast_db = filer.shard_db_name(shard_index)

        with filer.temp_file() as fasta_file:
            AtramPreprocessor.fill_blast_fasta(fasta_file, filer, shard_params)

            cmd = 'makeblastdb -dbtype nucl -in {} -out {}'
            cmd = cmd.format(fasta_file.name, blast_db)
            subprocess.check_call(cmd, shell=True)

    @staticmethod
    def fill_blast_fasta(fasta_file, filer, shard_params):
        """
        Fill the fasta file used as input into blast with shard sequences from
        the sqlite3 DB. We use the shard partitions passed in to determine
        which sequences to get for this shard.

        Because this is called in a child process, the address space is not
        shared with the parent (caller) and we cannot use instance variables.
        """

        db_conn = AtramPreprocessor.connect_db(filer)
        sql = '''
            SELECT seq_name, seq_end, seq FROM sequences
            ORDER BY seq_name LIMIT {} OFFSET {}
            '''
        sql = sql.format(shard_params[0], shard_params[1])
        result = db_conn.execute(sql)
        for row in result:
            fasta_file.write('>{}{}\n'.format(row[0], row[1]))
            fasta_file.write('{}\n'.format(row[2]))
        db_conn.close()


if __name__ == '__main__':

    AtramPreprocessor().preprocess()
