"""Format the data so that atram can use it."""

import re
import sqlite3
import logging
import tempfile
import subprocess
import multiprocessing
import numpy as np
from filer import Filer
from configure import Configure


class AtramPreprocessor:
    """Class to handle getting the data into a form that atram can use."""

    batch_size = 1e7  # How many sequence records to insert at a time

    # Try to get the sequence name and which end it is from the fasta header
    parse_header = re.compile(r'^ [>@] \s* ( .* ) ( [\s\/_] [12] )',
                              re.VERBOSE)

    def __init__(self):
        self.config = None
        self.db_conn = None
        self.shard_list = None

    def preprocess(self):
        """The main program."""

        self.config = Configure().parse_command_line(
            description="""
                This script prepares data for use by the atram.py script.
                It takes fasta or fastq files of paired-end (or single-end)
                sequence reads and creates a set of aTRAM databases.
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
        Setup the DB for our processing needs.

        This is called in a child process so the address space is not shared
        and we cannot use instance variables.
        """

        db_path = filer.db_file()
        db_conn = sqlite3.connect(db_path)
        db_conn.execute("PRAGMA page_size = {}".format(2**16))
        db_conn.execute("PRAGMA journal_mode = 'off'")
        db_conn.execute("PRAGMA synchronous = 'off'")
        return db_conn

    def create_table(self):
        """Reset the DB. Delete the table and read it."""

        logging.info('Creating sqlite tables')
        self.db_conn.execute('''DROP INDEX IF EXISTS seq_names''')
        self.db_conn.execute('''DROP TABLE IF EXISTS sequences''')
        sql = '''
            CREATE TABLE IF NOT EXISTS sequences
            (seq_name TEXT, seq_end TEXT, seq TEXT)
            '''
        self.db_conn.execute(sql)

    def create_index(self):
        """Create the index after we build the table."""

        logging.info('Creating sqlite indices')
        self.db_conn.execute('''
            CREATE INDEX IF NOT EXISTS seq_names ON sequences (seq_name)
            ''')

    def load_seqs(self):
        """
        A hand rolled version of "Bio.SeqIO". It's faster because we can take
        shortcuts due to its limited use.
        """

        for file_name in self.config.sra_files:
            logging.info('Loading "%s" into sqlite database', file_name)
            with open(file_name) as sra_file:
                recs, seq, seq_end, seq_name, is_seq = [], '', '', 0, True
                for line in sra_file:
                    if not line:
                        pass
                    elif line[0] in ['>', '@']:
                        if seq:
                            recs.append((seq_name, seq_end, seq))
                        seq = ''
                        match = AtramPreprocessor.parse_header.match(line)
                        is_seq = True
                        if match:
                            seq_name = match.group(1)
                            seq_end = match.group(2)
                        else:
                            seq_name = line[1:]
                            seq_end = ''
                    elif line[0] == '+':
                        is_seq = False
                    elif line[0].isalpha() and is_seq:
                        seq += line.rstrip()  # ''.join([strs])

                    if len(recs) >= AtramPreprocessor.batch_size:
                        self.bulk_insert(recs)
                        recs = []

                if seq:
                    recs.append((seq_name, seq_end, seq))
                self.bulk_insert(recs)

    def bulk_insert(self, recs):
        """Insert a batch of sequence records into the sqlite database."""

        if recs:
            sql = '''
                INSERT INTO sequences (seq_name, seq_end, seq) VALUES (?, ?, ?)
                '''
            self.db_conn.executemany(sql, recs)
            self.db_conn.commit()

    def assign_seqs_to_shards(self):
        """Put the sequences into the DB shards."""

        logging.info('Assigning sequences to shards')
        result = self.db_conn.execute('SELECT COUNT(*) FROM sequences')
        total = result.fetchone()[0]
        offsets = np.linspace(0, total,
                              num=self.config.shard_count + 1, dtype=int)
        for i in range(1, len(offsets) - 1):
            sql = '''
                SELECT seq_name FROM sequences
                ORDER BY seq_name LIMIT 2 OFFSET {}
                '''
            result = self.db_conn.execute(sql.format(offsets[i]))
            first = result.fetchone()[0]
            second = result.fetchone()[0]
            if first != second:
                offsets[i] += 1
        limits = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]
        return list(zip(limits, offsets[:-1]))

    def create_blast_dbs(self):
        """Assign processes to make the blast DBs."""

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
        Create a blast DB from the shard.

        This is called in a child process so the address space is not shared
        and we cannot use instance variables.
        """

        filer = Filer(work_dir=work_dir, file_prefix=file_prefix)
        blast_db = filer.shard_db_name(shard_index)

        with tempfile.NamedTemporaryFile(mode='w') as fasta_file:
            AtramPreprocessor.fill_blast_fasta(fasta_file, filer, shard_params)

            cmd = 'makeblastdb -dbtype nucl -in {} -out {}'
            cmd = cmd.format(fasta_file.name, blast_db)
            subprocess.check_call(cmd, shell=True)

    @staticmethod
    def fill_blast_fasta(fasta_file, filer, shard_params):
        """
        Fill the fasta file used as input into blast with shard sequences from
        the DB.

        This is called in a child process so the address space is not shared
        and we cannot use instance variables.
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
