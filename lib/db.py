"""Handle database functions."""

import sqlite3
import sys
import os
from os.path import basename, join, exists


ATRAM_VERSION = '2.0.alpha.5'
DB_VERSION = '2.0'

BATCH_SIZE = 1e6  # How many sequence records to insert at a time


def connect(blast_db, check_version=False, clean=False):
    """Create DB connection."""
    db_name = '{}.sqlite.db'.format(blast_db)

    if clean and exists(db_name):
        os.remove(db_name)

    db_conn = sqlite3.connect(db_name)

    db_conn.execute("PRAGMA page_size = {}".format(2**16))
    db_conn.execute("PRAGMA busy_timeout = 10000")
    db_conn.execute("PRAGMA journal_mode = OFF")
    db_conn.execute("PRAGMA synchronous = OFF")

    if check_version:
        check_versions(db_conn)

    return db_conn


def aux_db(db_conn, temp_dir, blast_db, query_name):
    """Create & attach an temporary database to the current DB connection."""
    db_dir = join(temp_dir, 'db')
    os.makedirs(db_dir, exist_ok=True)

    db_name = '{}_{}_temp.sqlite.db'.format(
        basename(blast_db), basename(query_name))
    db_name = join(db_dir, db_name)

    sql = """ATTACH DATABASE '{}' AS aux""".format(db_name)
    db_conn.execute(sql)


def aux_detach(db_conn):
    """Detach the temporary database."""
    db_conn.execute('DETACH DATABASE aux')


# ########################### misc functions #################################
# DB_VERSION != version in DB. Don't force DB changes until required. So
# this version will tend to lag ATRAM_VERSION.

def check_versions(db_conn):
    """Make sure the database version matches what we built it with."""
    version = get_version(db_conn)
    if version != DB_VERSION:
        err = ('The database was built with version {} but you are running '
               'version {}. You need to rebuild the atram database by '
               'running atram_preprocessor.py again.').format(
                   version, DB_VERSION)
        sys.exit(err)


# ########################## metadata table ##################################

def create_metadata_table(db_conn):
    """
    Create the metadata table.

    A single record used to tell if we are running atram.py against the
    schema version we built with atram_preprocessor.py.
    """
    db_conn.execute('''DROP TABLE IF EXISTS metadata''')
    sql = 'CREATE TABLE metadata (label TEXT, value TEXT)'
    db_conn.execute(sql)

    with db_conn:
        sql = '''INSERT INTO metadata (label, value) VALUES (?, ?)'''
        db_conn.execute(sql, ('version', DB_VERSION))
        db_conn.commit()


def get_version(db_conn):
    """Get the current database version."""
    sql = '''SELECT value FROM metadata WHERE label = ?'''
    try:
        result = db_conn.execute(sql, ('version', ))
        return result.fetchone()[0]
    except sqlite3.OperationalError:
        return '1.0'


# ########################## sequences table ##################################

def create_sequences_table(db_conn):
    """Create the sequence table."""
    db_conn.execute('''DROP TABLE IF EXISTS sequences''')
    sql = 'CREATE TABLE sequences (seq_name TEXT, seq_end TEXT, seq TEXT)'
    db_conn.execute(sql)


def create_sequences_index(db_conn):
    """
    Create the sequences index after we build the table.

    This speeds up the program significantly.
    """
    sql = 'CREATE INDEX sequences_index ON sequences (seq_name, seq_end)'
    db_conn.execute(sql)


def insert_sequences_batch(db_conn, batch):
    """Insert a batch of sequence records into the database."""
    if batch:
        sql = '''INSERT INTO sequences (seq_name, seq_end, seq)
                      VALUES (?, ?, ?)
            '''
        db_conn.executemany(sql, batch)
        db_conn.commit()


def get_sequence_count(db_conn):
    """Get the number of sequences in the table."""
    result = db_conn.execute('SELECT COUNT(*) FROM sequences')
    return result.fetchone()[0]


def get_shard_cut(db_conn, offset):
    """Get the sequence name at the given offset."""
    sql = 'SELECT seq_name FROM sequences ORDER BY seq_name LIMIT 1 OFFSET {}'
    result = db_conn.execute(sql.format(offset))
    cut = result.fetchone()[0]
    return cut


def get_sequences_in_shard(db_conn, start, end):
    """Get all sequences in a shard."""
    sql = '''
        SELECT seq_name, seq_end, seq
          FROM sequences
         WHERE seq_name >= ?
           AND seq_name < ?
        '''
    return db_conn.execute(sql, (start, end))


# ######################## sra_blast_hits table ###############################

def create_sra_blast_hits_table(db_conn):
    """
    Reset the DB.

    Delete the tables and recreate them.
    """
    db_conn.execute('''DROP TABLE IF EXISTS aux.sra_blast_hits''')
    sql = '''
        CREATE TABLE aux.sra_blast_hits
                   (iteration INTEGER, seq_name TEXT, seq_end TEXT, shard TEXT)
        '''
    db_conn.execute(sql)

    sql = '''
        CREATE INDEX aux.sra_blast_hits_index
                  ON sra_blast_hits (iteration, seq_name, seq_end)
        '''
    db_conn.execute(sql)


def insert_blast_hit_batch(db_conn, batch):
    """Insert a batch of blast hit records into the database."""
    if batch:
        sql = '''
            INSERT INTO aux.sra_blast_hits
                        (iteration, seq_end, seq_name, shard)
                        VALUES (?, ?, ?, ?)
            '''
        db_conn.executemany(sql, batch)
        db_conn.commit()


def sra_blast_hits_count(db_conn, iteration):
    """Count the blast hist for select the iteration."""
    sql = '''
        SELECT COUNT(*) AS count
          FROM aux.sra_blast_hits
         WHERE iteration = ?
        '''

    result = db_conn.execute(sql, (iteration, ))
    return result.fetchone()[0]


def get_sra_blast_hits(db_conn, iteration):
    """Get all blast hits for the iteration."""
    sql = '''
        SELECT seq_name, seq_end, seq
          FROM sequences
         WHERE seq_name IN (SELECT DISTINCT seq_name
                              FROM aux.sra_blast_hits
                             WHERE iteration = ?)
      ORDER BY seq_name, seq_end
        '''

    db_conn.row_factory = sqlite3.Row
    return db_conn.execute(sql, (iteration, ))


def get_blast_hits_by_end_count(db_conn, iteration, end_count):
    """Get all blast hits for the iteration."""
    sql = '''
        SELECT seq_name, seq_end, seq
          FROM sequences
         WHERE seq_name IN (SELECT seq_name
                              FROM sequences
                             WHERE seq_name IN (SELECT DISTINCT seq_name
                                                  FROM aux.sra_blast_hits
                                                 WHERE iteration = ?)
                          GROUP BY seq_name
                            HAVING COUNT(*) = ?)
      ORDER BY seq_name, seq_end
        '''

    db_conn.row_factory = sqlite3.Row
    return db_conn.execute(sql, (iteration, end_count))


# ####################### contig_blast_hits table #############################

def create_contig_blast_hits_table(db_conn):
    """Reset the database. Delete the tables and recreate them."""
    db_conn.execute('''DROP TABLE IF EXISTS aux.contig_blast_hits''')
    sql = '''
        CREATE TABLE aux.contig_blast_hits
                     (iteration INTEGER, contig_id TEXT, description TEXT,
                      bit_score NUMERIC, len INTEGER,
                      query_from INTEGER, query_to INTEGER, query_strand TEXT,
                      hit_from INTEGER, hit_to INTEGER, hit_strand TEXT)
        '''
    db_conn.execute(sql)

    sql = '''
        CREATE INDEX aux.contig_blast_hits_index
                  ON contig_blast_hits (iteration, bit_score, len)
        '''
    db_conn.execute(sql)


def insert_contig_hit_batch(db_conn, batch):
    """Insert a batch of blast hit records into the database."""
    if batch:
        sql = '''
            INSERT INTO aux.contig_blast_hits
                        (iteration, contig_id, description, bit_score, len,
                         query_from, query_to, query_strand,
                         hit_from, hit_to, hit_strand)
                 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            '''
        db_conn.executemany(sql, batch)
        db_conn.commit()


def get_contig_blast_hits(db_conn, iteration):
    """Get all blast hits for the iteration."""
    sql = '''
        SELECT iteration, contig_id, description, bit_score, len,
               query_from, query_to, query_strand,
               hit_from, hit_to, hit_strand
          FROM aux.contig_blast_hits
         WHERE iteration = ?
        '''

    db_conn.row_factory = sqlite3.Row
    return db_conn.execute(sql, (iteration, ))


# ####################### assembled_contigs table #############################

def create_assembled_contigs_table(db_conn):
    """Reset the database. Delete the tables and recreate them."""
    db_conn.execute('''DROP TABLE IF EXISTS aux.assembled_contigs''')
    sql = '''
        CREATE TABLE aux.assembled_contigs
                     (iteration INTEGER, contig_id TEXT, seq TEXT,
                      description TEXT, bit_score NUMERIC, len INTEGER,
                      query_from INTEGER, query_to INTEGER, query_strand TEXT,
                      hit_from INTEGER, hit_to INTEGER, hit_strand TEXT)
        '''
    db_conn.execute(sql)

    sql = '''
        CREATE INDEX aux.assembled_contigs_index
                  ON assembled_contigs (iteration, contig_id)
        '''
    db_conn.execute(sql)


def assembled_contigs_count(db_conn, iteration, bit_score, length):
    """Count the blast hist for the iteration."""
    sql = '''
        SELECT COUNT(*) AS count
          FROM aux.assembled_contigs
         WHERE iteration = ?
           AND bit_score >= ?
           AND len >= ?
        '''

    result = db_conn.execute(sql, (iteration, bit_score, length))
    return result.fetchone()[0]


def iteration_overlap_count(db_conn, iteration, bit_score, length):
    """Count how many assembled contigs match what's in the last iteration."""
    sql = '''
        SELECT COUNT(*) AS overlap
          FROM aux.assembled_contigs AS curr_iter
          JOIN aux.assembled_contigs AS prev_iter
            ON (     curr_iter.contig_id = prev_iter.contig_id
                 AND curr_iter.iteration = prev_iter.iteration + 1)
         WHERE curr_iter.iteration = ?
           AND curr_iter.seq = prev_iter.seq
           AND curr_iter.bit_score >= ?
           AND prev_iter.bit_score >= ?
           AND curr_iter.len >= ?
        '''
    result = db_conn.execute(
        sql, (iteration, bit_score, bit_score, length))
    return result.fetchone()[0]


def insert_assembled_contigs_batch(db_conn, batch):
    """Insert a batch of blast hit records into the database."""
    if batch:
        sql = '''
            INSERT INTO aux.assembled_contigs
                        (iteration, contig_id, seq, description,
                         bit_score, len,
                         query_from, query_to, query_strand,
                         hit_from, hit_to, hit_strand)
                 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            '''
        db_conn.executemany(sql, batch)
        db_conn.commit()


def get_assembled_contigs(db_conn, iteration, bit_score, length):
    """
    Get all assembled contigs for the iteration.

    We will use them as the queries in the next atram iteration.
    """
    sql = '''
        SELECT contig_id, seq
          FROM aux.assembled_contigs
         WHERE iteration = ?
           AND bit_score >= ?
           AND len >= ?
         '''
    return db_conn.execute(sql, (iteration, bit_score, length))


def get_all_assembled_contigs(db_conn, bit_score=0, length=0):
    """Get all assembled contigs."""
    sql = '''
        SELECT iteration, contig_id, seq, description, bit_score, len,
               query_from, query_to, query_strand,
               hit_from, hit_to, hit_strand
          FROM aux.assembled_contigs
         WHERE bit_score >= ?
           AND len >= ?
      ORDER BY bit_score DESC, iteration
        '''

    db_conn.row_factory = sqlite3.Row
    return db_conn.execute(sql, (bit_score, length))
