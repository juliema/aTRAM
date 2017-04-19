"""Handle SQL functions."""

import sqlite3

BATCH_SIZE = 1e7  # How many sequence records to insert at a time


def connect(blast_db):
    """Setup the DB for our processing needs and return a DB connection."""

    db_name = '{}.sqlite.db'.format(blast_db)

    db_conn = sqlite3.connect(db_name)

    db_conn.execute("PRAGMA busy_timeout = 10000")
    db_conn.execute("PRAGMA page_size = {}".format(2**16))
    db_conn.execute("PRAGMA journal_mode = 'off'")
    db_conn.execute("PRAGMA synchronous = 'off'")
    return db_conn


# ########################## sequences table ##################################

def create_sequences_table(db_conn):
    """Create the sequence table."""

    db_conn.execute('''DROP INDEX IF EXISTS seq_names''')
    db_conn.execute('''DROP TABLE IF EXISTS sequences''')
    sql = 'CREATE TABLE sequences (seq_name TEXT, seq_end TEXT, seq TEXT)'
    db_conn.execute(sql)


def create_sequences_index(db_conn):
    """Create the sequences index after we build the table. This speeds up the
    program significantly.
    """

    sql = 'CREATE INDEX sequences_index ON sequences (seq_name, seq_end)'
    db_conn.execute(sql)


def insert_sequences_batch(db_conn, batch):
    """Insert a batch of sequence records into the sqlite database."""

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


def get_two_sequences(db_conn, offset):
    """Get two sequences at the given offset. This is used for calculating
    where to cut-off a shard in the table. We want the shard to contain both
    paired ends of a sequence.
    """

    sql = 'SELECT seq_name FROM sequences ORDER BY seq_name LIMIT 2 OFFSET {}'
    result = db_conn.execute(sql.format(offset))
    first = result.fetchone()[0]
    second = result.fetchone()[0]
    return first, second


def get_sequences_in_shard(db_conn, limit, offset):
    """Get all sequences in a shard."""

    sql = '''
        SELECT seq_name, seq_end, seq FROM sequences
        ORDER BY seq_name LIMIT {} OFFSET {}
        '''
    sql = sql.format(limit, offset)

    return db_conn.execute(sql)


# ########################## blast_hits table #################################

def create_blast_hits_table(db_conn):
    """Reset the DB. Delete the tables and recreate them."""

    db_conn.execute('''DROP INDEX IF EXISTS blast_hits_index''')
    db_conn.execute('''DROP TABLE IF EXISTS blast_hits''')
    sql = '''CREATE TABLE blast_hits
        (iteration INTEGER, seq_name TEXT, seq_end TEXT, shard TEXT)
        '''
    db_conn.execute(sql)

    sql = 'CREATE INDEX blast_hits_index ON blast_hits (iteration, seq_name)'
    db_conn.execute(sql)


def insert_blast_hit_batch(db_conn, batch):
    """Insert a batch of blast hit records into the sqlite database."""

    if batch:
        sql = '''INSERT INTO blast_hits
            (iteration, seq_end, seq_name, shard) VALUES (?, ?, ?, ?)
            '''
        db_conn.executemany(sql, batch)
        db_conn.commit()


def blast_hits_count(db_conn, iteration):
    """Count the blast hist for the iteration."""

    sql = '''SELECT COUNT(*) AS count FROM blast_hits WHERE iteration = ?'''

    result = db_conn.execute(sql, str(iteration))
    return result.fetchone()[0]


def get_blast_hits(db_conn, iteration):
    """Get all blast hits for the iteration."""

    sql = '''SELECT seq_name, seq_end, seq FROM sequences
        WHERE seq_name IN (SELECT DISTINCT seq_name
                             FROM blast_hits
                             WHERE iteration = ?)
        ORDER BY seq_name, seq_end
        '''

    db_conn.row_factory = sqlite3.Row
    return db_conn.execute(sql, str(iteration))


# ####################### assembled_contigs table #############################

def create_assembled_contigs_table(db_conn):
    """Reset the DB. Delete the tables and recreate them."""

    db_conn.execute('''DROP INDEX IF EXISTS assembled_contigs_index''')
    db_conn.execute('''DROP TABLE IF EXISTS assembled_contigs''')
    sql = '''CREATE TABLE assembled_contigs
        (iteration INTEGER, contig_id TEXT, seq TEXT, description TEXT,
         bit_score NUMERIC, len INTEGER,
         query_from INTEGER, query_to INTEGER, query_strand TEXT,
         hit_from INTEGER, hit_to INTEGER, hit_strand TEXT)
        '''
    db_conn.execute(sql)

    sql = '''CREATE INDEX assembled_contigs_index
        ON assembled_contigs (iteration, contig_id)
        '''
    db_conn.execute(sql)


def assembled_contigs_count(db_conn, iteration):
    """Count the blast hist for the iteration."""

    sql = '''SELECT COUNT(*) AS count
        FROM assembled_contigs WHERE iteration = ?
        '''

    result = db_conn.execute(sql, str(iteration))
    return result.fetchone()[0]


def iteration_overlap_count(db_conn, iteration):
    """Count how many assembled contigs match what was in the last iteration.
    """

    sql = '''SELECT COUNT(*) AS overlap
        FROM assembled_contigs AS curr_iter
        JOIN assembled_contigs AS prev_iter
          ON (    curr_iter.contig_id = prev_iter.contig_id
              AND curr_iter.iteration = prev_iter.iteration + 1)
        WHERE curr_iter.iteration = ?
          AND curr_iter.seq = prev_iter.seq
        '''
    result = db_conn.execute(sql, str(iteration))
    return result.fetchone()[0]


def insert_assembled_contigs_batch(db_conn, batch):
    """Insert a batch of blast hit records into the sqlite database."""

    if batch:
        sql = '''INSERT INTO assembled_contigs
            (iteration, contig_id, seq, description, bit_score, len,
             query_from, query_to, query_strand, hit_from, hit_to, hit_strand)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            '''
        db_conn.executemany(sql, batch)
        db_conn.commit()


def get_assembled_contigs(db_conn, iteration):
    """Get all assembled contigs for the iteration so that we can use them
    as the targets in the next atram iteration.
    """

    sql = 'SELECT contig_id, seq FROM assembled_contigs WHERE iteration = ?'
    return db_conn.execute(sql, str(iteration))


def get_all_assembled_contigs(db_conn):
    """Get all assembled contigs for the iteration so that we can use them
    as the targets in the next atram iteration.
    """

    sql = '''SELECT iteration, contig_id, seq, description, bit_score, len,
                    query_from, query_to, query_strand,
                    hit_from, hit_to, hit_strand
        FROM assembled_contigs ORDER BY bit_score DESC, iteration
        '''

    db_conn.row_factory = sqlite3.Row
    return db_conn.execute(sql)
