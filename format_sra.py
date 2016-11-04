import os
import re
import sys
import glob
import time
import sqlite3
import argparse
import subprocess
from multiprocessing import Pool
import numpy as np

BATCH_SIZE = 750000
FRAGMENT = re.compile(r'[>@] \s* ( .*? ) ( [\s\/_] [12] )', re.VERBOSE)

DEFAULT_SHARDS = 8  # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DEFAULT_PROCESSES = DEFAULT_SHARDS  # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!
DB_PATH = 'data/temp_sqlite.db'  # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print('%s function took %0.3f s' % (f.__name__, (time2 - time1)))
        return ret
    return wrap


def log(s):
    print(s)  # TODO


def create_temp_file(shard):
    return 1


def bulk_insert(db, recs, count):
    if recs:
        db.executemany('''INSERT INTO ends (frag, seq, seq_end) VALUES (?, ?, ?)''', recs)
        db.commit()


def append_rec(recs, frag, seq, end, count):
    if seq:
        recs.append((frag, seq, end))
        count += 1
    return count


# @timing
def load_seqs(db, sra_files):
    count = 0
    for file_name in sra_files:
        log('Loading "{}" into sqlite database'.format(file_name))
        with open(file_name, 'r') as sra_file:
            recs = []
            seq, end, frag = '', '', 0
            for line in sra_file:
                if not line:
                    pass
                elif line[0] in ['>', '@']:
                    count = append_rec(recs, frag, seq, end, count)
                    seq = ''
                    match = FRAGMENT.match(line)
                    if match:
                        frag = match.group(1)
                        end = match.group(2)
                    else:
                        frag = line[1:]
                        end = ''
                elif line[0].isalpha():
                    seq += line.rstrip()
                if len(recs) >= BATCH_SIZE:
                    bulk_insert(db, recs, count)
                    recs = []
            else:
                count = append_rec(recs, frag, seq, end, count)
                bulk_insert(db, recs, count)


def create_table(db):
    log('Creating sqlite tables')
    db.execute('''CREATE TABLE IF NOT EXISTS ends (frag TEXT, seq TEXT, seq_end TEXT)''')


# @timing
def create_index(db):
    log('Creating sqlite indices')
    db.execute('''CREATE INDEX IF NOT EXISTS frag ON ends (frag)''')


def connect_db(db_path):
    return sqlite3.connect(db_path)


def assign_seqs_to_shards(db, shard_count):
    log('Assigning sequences to shards')
    connection = db.execute('SELECT COUNT(*) FROM ends')
    total = connection.fetchone()[0]
    offsets = np.linspace(0, total, num=shard_count + 1, dtype=int)
    for i in range(1, len(offsets) - 1):
        connection = db.execute('SELECT frag FROM ends ORDER BY frag LIMIT 2 OFFSET {}'.format(offsets[i]))
        first = connection.fetchone()[0]
        second = connection.fetchone()[0]
        if first != second:
            offsets[i] += 1
    limits = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]
    connection.close()
    return list(zip(limits, offsets[:-1]))


def create_blast_db(shard, i):
    db = connect_db(DB_PATH)
    sql = 'SELECT frag, seq_end, seq FROM ends ORDER BY frag LIMIT {} OFFSET {}'.format(shard[0], shard[1])
    connection = db.execute(sql)
    fasta_file = 'data/temp_seqs_{}.fasta'.format(i)
    blast_prefix = 'data/blast_{}'.format(i)
    with open(fasta_file, 'w') as out_file:  # TODO
        for row in connection:
            out_file.write('>{}{}\n'.format(row[0], row[1]))
            out_file.write('{}\n'.format(row[2]))
    db.close()
    result = subprocess.check_call(('makeblastdb -dbtype nucl -in {} -out {}').format(
        fasta_file, blast_prefix), shell=True)
    os.remove(fasta_file)


def create_blast_dbs(process_count, shards):
    log('Making blast DBs')
    with Pool(processes=process_count) as pool:
        for i, shard in enumerate(shards):
            proc = pool.Process(target=create_blast_db, args=(shard, i))
            proc.start()
    pool.join()


def delete_temp_files():
    log('Deleting temporary files')
    try:
        for f in glob.glob(DB_PATH + '*'):  # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!
            os.remove(f)
    except:
        log('Did not find the temporary sqlite database files')

    try:
        for f in glob.glob('data/blast_*'):  # TODO !!!!!!!!!!!!!!!!!!!!!!!!!!
            os.remove(f)
    except:
        log('Did not find the temporary blast files')


def parse_args():
    parser = argparse.ArgumentParser(description='''
        Takes fasta or fastq files of paired-end (or single-end) short reads and creates an aTRAM database.
    ''')
    parser.add_argument('sra_files', nargs='+', help='short read archive in fasta or fastq format.')
    parser.add_argument('-o', '--out', help='output aTRAM files with this prefix')
    parser.add_argument('-l', '--log', help='save log information to this file')
    parser.add_argument('-n', '--shards', type=int, help='number of shards to create')
    parser.add_argument('-p', '--processes', type=int, help='number of processes to create')
    args = parser.parse_args()

    args.shards = args.shards if args.shards else DEFAULT_SHARDS
    args.processes = args.processes if args.processes else DEFAULT_PROCESSES

    # TODO sra_files must exist and be in fasta/q format
    # TODO sra_files should be globs

    return args


if __name__ == "__main__":
    args = parse_args()

    db = connect_db(DB_PATH)  # TODO: make based on args

    create_table(db)
    load_seqs(db, args.sra_files)
    create_index(db)

    shards = assign_seqs_to_shards(db, args.shards)
    create_blast_dbs(args.processes, shards)
    delete_temp_files()

    db.close()
