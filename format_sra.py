#!/usr/bin/python3

"""Create the initial blast databases that will be used by aTRAM."""

import re
import os
import sqlite3
import logging
import argparse
import subprocess
import multiprocessing
import numpy as np
import configure
import log
# Bio.SeqIO  # 2-3x slower than load_seqs() :(

FRAGMENT = re.compile(r'^ [>@] \s* ( .* ) ( [\s\/_] [12] )', re.VERBOSE)
DEFAULT_BATCH_SIZE = 1e7


def bulk_insert(db_conn, recs):
    """Insert a batch of sequence records into the sqlite database."""
    if recs:
        db_conn.executemany('''INSERT INTO frags (frag, frag_end, seq) VALUES (?, ?, ?)''', recs)
        db_conn.commit()


def load_seqs(db_conn, config):
    """A faster version of "Bio.SeqIO. We can take shortcuts because of the limited use."""
    for file_name in config.sra_files:
        logging.info('Loading "%s" into sqlite database', file_name)
        with open(file_name, 'r') as sra_file:
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
                    seq += line.rstrip()  # Almost always singletons so ''.join([strs]) is no help

                if len(recs) >= DEFAULT_BATCH_SIZE:
                    bulk_insert(db_conn, recs)
                    recs = []

            if seq:
                recs.append((frag, frag_end, seq))
            bulk_insert(db_conn, recs)


def create_table(db_conn):
    """Reset the DB. Delete the table and readd it."""
    logging.info('Creating sqlite tables')
    db_conn.execute('''DROP INDEX IF EXISTS frag''')
    db_conn.execute('''DROP TABLE IF EXISTS frags''')
    db_conn.execute('''CREATE TABLE IF NOT EXISTS frags (frag TEXT, frag_end TEXT, seq TEXT)''')


def create_index(db_conn):
    """Create the index after we build the table."""
    logging.info('Creating sqlite indices')
    db_conn.execute('''CREATE INDEX IF NOT EXISTS frag ON frags (frag)''')


def connect_db(config):
    """Setup the DB for our processing needs."""
    db_path = '{}sqlite.db'.format(config.out)
    db_conn = sqlite3.connect(db_path)
    db_conn.execute("PRAGMA page_size = {}".format(2**16))
    db_conn.execute("PRAGMA journal_mode = 'off'")
    db_conn.execute("PRAGMA synchronous = 'off'")
    return db_conn


def assign_seqs_to_shards(db_conn, config):
    """Put the sequences into the DB shards."""
    logging.info('Assigning sequences to shards')
    connection = db_conn.execute('SELECT COUNT(*) FROM frags')
    total = connection.fetchone()[0]
    offsets = np.linspace(0, total, num=config.shards + 1, dtype=int)
    for i in range(1, len(offsets) - 1):
        connection = db_conn.execute(
            'SELECT frag FROM frags ORDER BY frag LIMIT 2 OFFSET {}'.format(offsets[i]))
        first = connection.fetchone()[0]
        second = connection.fetchone()[0]
        if first != second:
            offsets[i] += 1
    limits = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]
    connection.close()
    return list(zip(limits, offsets[:-1]))


def create_blast_db(config, shard, i):
    """Create a blast DB from the shard."""
    db_conn = connect_db(config)
    sql = ('SELECT frag, frag_end, seq FROM frags '
           'ORDER BY frag LIMIT {} OFFSET {}').format(shard[0], shard[1])
    connection = db_conn.execute(sql)
    fasta_file = '{}temp_seqs_{}.format'.format(config.out, i)
    blast_prefix = '{}blast_{}'.format(config.out, i)
    with open(fasta_file, 'w') as out_file:
        for row in connection:
            out_file.write('>{}{}\n'.format(row[0], row[1]))
            out_file.write('{}\n'.format(row[2]))
    db_conn.close()
    subprocess.check_call(('makeblastdb -dbtype nucl -in {} -out {}').format(
        fasta_file, blast_prefix), shell=True)
    os.remove(fasta_file)


def create_blast_dbs(config, shards):
    """Assign processes to make the blast DBs."""
    logging.info('Making blast DBs')
    with multiprocessing.Pool(processes=config.processes) as pool:
        for i, shard in enumerate(shards):
            proc = pool.Process(target=create_blast_db, args=(config, shard, i))
            proc.start()
    pool.join()


def parse_args():
    """Parse the input arguments and assign defaults."""
    parser = argparse.ArgumentParser(
        description=('Takes fasta or fastq files of paired-end (or single-end) '
                     'short reads and creates an aTRAM database.'))
    parser.add_argument('sra_files', nargs='+',
                        help='short read archives in fasta or fastq format. May contain wildcards.')
    parser.add_argument('-o', '--out',
                        help=('output aTRAM files with this prefix. '
                              'May include a directory in the prefix.'))
    parser.add_argument('-s', '--shards', type=int, help='number of shards to create')
    parser.add_argument('-p', '--processes', type=int, help='number of processes to create',
                        default=configure.default_process_count())
    config = parser.parse_args()

    if not config.shards:
        config.shards = configure.default_shard_count(config.sra_files)

    return config


if __name__ == '__main__':
    ARGS = parse_args()
    log.setup(ARGS)

    DB = connect_db(ARGS)
    create_table(DB)
    load_seqs(DB, ARGS)
    create_index(DB)

    SHARDS = assign_seqs_to_shards(DB, ARGS)
    create_blast_dbs(ARGS, SHARDS)

    DB.close()
