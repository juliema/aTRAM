"""Create the initial blast databases that will be used by aTRAM."""

import re
import sqlite3
import logging
import tempfile
import subprocess
import multiprocessing
import numpy as np
from filer import Filer
from configure import Configure


DEFAULT_BATCH_SIZE = 1e7
FRAGMENT = re.compile(r'^ [>@] \s* ( .* ) ( [\s\/_] [12] )', re.VERBOSE)


def preprocessor():
    """The main program."""

    config = Configure().parse_command_line(
        description="""
            This script prepares data for use by the atram.py script. It takes
            fasta or fastq files of paired-end (or single-end) sequence reads
            and creates a set of aTRAM databases.
            """,
        args='sra_files file_prefix work_dir shard_count cpus')

    filer = Filer(config.work_dir, config.file_prefix)
    filer.log_setup()

    db_conn = connect_db(filer)
    create_table(db_conn)
    load_seqs(db_conn, config.sra_files)
    create_index(db_conn)

    shard_list = assign_seqs_to_shards(db_conn, config.shard_count)
    create_blast_dbs(config.cpus, config.work_dir, config.file_prefix,
                     shard_list)

    db_conn.close()


def bulk_insert(db_conn, recs):
    """Insert a batch of sequence records into the sqlite database."""

    if recs:
        sql = '''INSERT INTO frags (frag, frag_end, seq) VALUES (?, ?, ?)'''
        db_conn.executemany(sql, recs)
        db_conn.commit()


def load_seqs(db_conn, sra_files):
    """
    A hand rolled version of "Bio.SeqIO". It's faster because we can take
    shortcuts due to its limited use.
    """

    for file_name in sra_files:
        logging.info('Loading "%s" into sqlite database', file_name)
        with open(file_name) as sra_file:
            recs, seq, frag_end, frag, is_seq = [], '', '', 0, True
            for line in sra_file:
                if not line:
                    pass
                elif line[0] in ['>', '@']:
                    if seq:
                        recs.append((frag, frag_end, seq))
                    seq = ''
                    match = FRAGMENT.match(line)
                    is_seq = True
                    if match:
                        frag = match.group(1)
                        frag_end = match.group(2)
                    else:
                        frag = line[1:]
                        frag_end = ''
                elif line[0] == '+':
                    is_seq = False
                elif line[0].isalpha() and is_seq:
                    seq += line.rstrip()  # ''.join([strs])

                if len(recs) >= DEFAULT_BATCH_SIZE:
                    bulk_insert(db_conn, recs)
                    recs = []

            if seq:
                recs.append((frag, frag_end, seq))
            bulk_insert(db_conn, recs)


def create_table(db_conn):
    """Reset the DB. Delete the table and read it."""

    logging.info('Creating sqlite tables')
    db_conn.execute('''DROP INDEX IF EXISTS frag''')
    db_conn.execute('''DROP TABLE IF EXISTS frags''')
    sql = '''
        CREATE TABLE IF NOT EXISTS frags (frag TEXT, frag_end TEXT, seq TEXT)
    '''
    db_conn.execute(sql)


def create_index(db_conn):
    """Create the index after we build the table."""

    logging.info('Creating sqlite indices')
    db_conn.execute('''CREATE INDEX IF NOT EXISTS frag ON frags (frag)''')


def connect_db(filer):
    """Setup the DB for our processing needs."""

    db_path = filer.db_file()
    db_conn = sqlite3.connect(db_path)
    db_conn.execute("PRAGMA page_size = {}".format(2**16))
    db_conn.execute("PRAGMA journal_mode = 'off'")
    db_conn.execute("PRAGMA synchronous = 'off'")
    return db_conn


def assign_seqs_to_shards(db_conn, shard_count):
    """Put the sequences into the DB shards."""

    logging.info('Assigning sequences to shards')
    result = db_conn.execute('SELECT COUNT(*) FROM frags')
    total = result.fetchone()[0]
    offsets = np.linspace(0, total, num=shard_count + 1, dtype=int)
    for i in range(1, len(offsets) - 1):
        sql = 'SELECT frag FROM frags ORDER BY frag LIMIT 2 OFFSET {}'
        result = db_conn.execute(sql.format(offsets[i]))
        first = result.fetchone()[0]
        second = result.fetchone()[0]
        if first != second:
            offsets[i] += 1
    limits = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]
    return list(zip(limits, offsets[:-1]))


def fill_blast_fasta(fasta_file, filer, shard_params):
    """
    Fill the fasta file used as input into blast with shard sequences from
    the DB.
    """

    db_conn = connect_db(filer)
    sql = '''
        SELECT frag, frag_end, seq FROM frags ORDER BY frag LIMIT {} OFFSET {}
    '''
    sql = sql.format(shard_params[0], shard_params[1])
    result = db_conn.execute(sql)
    for row in result:
        fasta_file.write('>{}{}\n'.format(row[0], row[1]))
        fasta_file.write('{}\n'.format(row[2]))
    db_conn.close()


def create_blast_db(work_dir, file_prefix, shard_params, shard_index):
    """Create a blast DB from the shard."""

    filer = Filer(work_dir=work_dir, file_prefix=file_prefix)
    blast_db = filer.shard_db_name(shard_index)

    with tempfile.NamedTemporaryFile(mode='w') as fasta_file:
        fill_blast_fasta(fasta_file, filer, shard_params)

        cmd = 'makeblastdb -dbtype nucl -in {} -out {}'
        cmd = cmd.format(fasta_file.name, blast_db)
        subprocess.check_call(cmd, shell=True)


def create_blast_dbs(cpus, work_dir, file_prefix, shard_list):
    """Assign processes to make the blast DBs."""

    logging.info('Making blast DBs')
    with multiprocessing.Pool(processes=cpus) as pool:
        results = [pool.apply_async(
            create_blast_db,
            (work_dir, file_prefix, shard_params, shard_index))
                   for shard_index, shard_params in enumerate(shard_list)]
        _ = [result.get() for result in results]
    logging.info('Finished making blast DBs')


if __name__ == '__main__':

    preprocessor()
