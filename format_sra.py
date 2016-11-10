import re
import os
import sys
import time
import sqlite3
import argparse
import datetime
import functools
import subprocess
import multiprocessing
import numpy as np
# from Bio import SeqIO  # 2-3x slower

FRAGMENT = re.compile(r'^ [>@] \s* ( .* ) ( [\s\/_] [12] )', re.VERBOSE)

DEFAULT_BATCH_SIZE = 1e7
LOG_FILE = None


def create_log_file(args):
    global LOG_FILE
    log_file = '{}log.txt'.format(args.out)
    LOG_FILE = open(log_file, 'w')


def log(s):
    global LOG_FILE
    LOG_FILE.write('{}: {}\n'.format(datetime.datetime.now().isoformat(' '), s))


def bulk_insert(db, recs):
    if recs:
        db.executemany('''INSERT INTO frags (frag, frag_end, seq) VALUES (?, ?, ?)''', recs)
        db.commit()


def load_seqs(db, args):
    for file_name in args.sra_files:
        log('Loading "{}" into sqlite database'.format(file_name))
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
                    seq += line.rstrip()
                if len(recs) >= DEFAULT_BATCH_SIZE:
                    bulk_insert(db, recs)
                    recs = []
            else:
                if seq:
                    recs.append((frag, frag_end, seq))
                bulk_insert(db, recs)


def create_table(db):
    log('Creating sqlite tables')
    db.execute('''DROP INDEX IF EXISTS frag''')
    db.execute('''DROP TABLE IF EXISTS frags''')
    db.execute('''CREATE TABLE IF NOT EXISTS frags (frag TEXT, frag_end TEXT, seq TEXT)''')


def create_index(db):
    log('Creating sqlite indices')
    db.execute('''CREATE INDEX IF NOT EXISTS frag ON frags (frag)''')


def connect_db(args):
    db_path = '{}sqlite.db'.format(args.out)
    db = sqlite3.connect(db_path)
    db.execute("PRAGMA page_size = {}".format(2**16))
    db.execute("PRAGMA journal_mode = 'off'")
    db.execute("PRAGMA synchronous = 'off'")
    return db


def assign_seqs_to_shards(db, args):
    log('Assigning sequences to shards')
    connection = db.execute('SELECT COUNT(*) FROM frags')
    total = connection.fetchone()[0]
    offsets = np.linspace(0, total, num=args.shards + 1, dtype=int)
    for i in range(1, len(offsets) - 1):
        connection = db.execute('SELECT frag FROM frags ORDER BY frag LIMIT 2 OFFSET {}'.format(offsets[i]))
        first = connection.fetchone()[0]
        second = connection.fetchone()[0]
        if first != second:
            offsets[i] += 1
    limits = [offsets[i + 1] - offsets[i] for i in range(len(offsets) - 1)]
    connection.close()
    return list(zip(limits, offsets[:-1]))


def create_blast_db(args, shard, i):
    db = connect_db(args)
    sql = 'SELECT frag, frag_end, seq FROM frags ORDER BY frag LIMIT {} OFFSET {}'.format(shard[0], shard[1])
    connection = db.execute(sql)
    fasta_file = '{}temp_seqs_{}.format'.format(args.out, i)
    blast_prefix = '{}blast_{}'.format(args.out, i)
    with open(fasta_file, 'w') as out_file:
        for row in connection:
            out_file.write('>{}{}\n'.format(row[0], row[1]))
            out_file.write('{}\n'.format(row[2]))
    db.close()
    result = subprocess.check_call(('makeblastdb -dbtype nucl -in {} -out {}').format(
        fasta_file, blast_prefix), shell=True)
    os.remove(fasta_file)


def create_blast_dbs(args, shards):
    log('Making blast DBs')
    with multiprocessing.Pool(processes=args.processes) as pool:
        for i, shard in enumerate(shards):
            proc = pool.Process(target=create_blast_db, args=(args, shard, i))
            proc.start()
    pool.join()


def parse_args():
    parser = argparse.ArgumentParser(description='''
        Takes fasta or fastq files of paired-end (or single-end) short reads and creates an aTRAM database.
    ''')
    parser.add_argument('sra_files', nargs='+', help='short read archives in fasta or fastq format. May contain wildcards.')
    parser.add_argument('-o', '--out', help='output aTRAM files with this prefix. May include a directory in the prefix.')
    parser.add_argument('-s', '--shards', type=int, help='number of shards to create')
    parser.add_argument('-p', '--processes', type=int, help='number of processes to create')
    args = parser.parse_args()

    size = 0
    for f in args.sra_files:
        size += os.path.getsize(f) / 2 if f.lower().endswith('.fastq') else os.path.getsize(f)
    size = int(size / 2.5e8)
    default_shards = size if size else 1
    args.shards = args.shards if args.shards else default_shards

    default_processes = os.cpu_count() - 2 if os.cpu_count() > 2 else 1
    args.processes = args.processes if args.processes else default_processes

    return args


if __name__ == '__main__':
    args = parse_args()

    create_log_file(args)
    log(' '.join(sys.argv))

    db = connect_db(args)

    create_table(db)
    load_seqs(db, args)
    create_index(db)

    shards = assign_seqs_to_shards(db, args)
    create_blast_dbs(args, shards)

    db.close()
    LOG_FILE.close()
