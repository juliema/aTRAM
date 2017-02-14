"""Handle SQL functions."""

import sqlite3

BATCH_SIZE = 1e7  # How many sequence records to insert at a time


def connect(filer):
    """Setup the DB for our processing needs and return a DB connection."""

    db_path = filer.db_file_name()
    db_conn = sqlite3.connect(db_path)
    db_conn.execute("PRAGMA page_size = {}".format(2**16))
    db_conn.execute("PRAGMA journal_mode = 'off'")
    db_conn.execute("PRAGMA synchronous = 'off'")
    return db_conn


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

    db_conn.execute('CREATE INDEX sequences_index ON sequences (seq_name)')


def create_blast_hits_table(db_conn):
    """Reset the DB. Delete the tables and recreate them."""

    db_conn.execute('''DROP INDEX IF EXISTS blast_hits_index''')
    db_conn.execute('''DROP TABLE IF EXISTS blast_hits''')
    sql = '''CREATE TABLE blast_hits
        (iteration INTEGER, seq_name TEXT, seq_end TEXT, shard_name TEXT)
        '''
    db_conn.execute(sql)


def create_blast_hits_index(db_conn):
    """Create the blast_hits_index after we build the table."""

    sql = 'CREATE INDEX blast_hits_index ON blast_hits (iteration, seq_name)'
    db_conn.execute(sql)


def insert_sequences_batch(db_conn, sequence_batch):
    """Insert a batch of sequence records into the sqlite database."""

    if sequence_batch:
        sql = '''INSERT INTO sequences (seq_name, seq_end, seq)
            VALUES (?, ?, ?)
            '''
        db_conn.executemany(sql, sequence_batch)
        db_conn.commit()


def insert_blast_hit_batch(db_conn, blast_hit_batch):
    """Insert a batch of blast hit records into the sqlite database."""

    if blast_hit_batch:
        sql = '''INSERT INTO blast_hits
            (iteration, seq_end, seq_name, shard_name) VALUES (?, ?, ?, ?)
            '''
        db_conn.executemany(sql, blast_hit_batch)
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


def get_blast_hits(db_conn, iteration):
    """Get all blast hits for the iteration."""

    print(iteration, type(iteration))
    sql = '''WITH hits AS
        (SELECT DISTINCT seq_name FROM blast_hits WHERE iteration = ?)
         SELECT seq_name, seq_end, seq FROM sequences WHERE seq_name IN hits
         ORDER BY seq_name, seq_end
        '''
    return db_conn.execute(sql, str(iteration))
