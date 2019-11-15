"""Database functions for the preprocessor."""

from lib.db import DB_VERSION


# ########################## metadata table ##################################

def create_metadata_table(cxn, args):
    """
    Create the metadata table.

    Information used to tell how aTRAM was set up.
    """
    cxn.executescript("""
        DROP TABLE IF EXISTS metadata;

        CREATE TABLE metadata (
            label TEXT,
            value TEXT);
        """)

    with cxn:
        sql = """INSERT INTO metadata (label, value) VALUES (?, ?)"""
        cxn.execute(sql, ('version', DB_VERSION))
        cxn.execute(sql, ('single_ends', bool(args.get('single_ends'))))


# ########################## sequences table ##################################

def create_sequences_table(cxn):
    """Create a table to hold the raw input sequences."""
    cxn.executescript("""
        DROP TABLE IF EXISTS sequences;

        CREATE TABLE sequences (
            seq_name TEXT,
            seq_end  TEXT,
            seq      TEXT);
        """)


def create_sequences_index(cxn):
    """
    Create the sequences index after we build the table.

    This speeds up the program significantly.
    """
    sql = 'CREATE INDEX sequences_index ON sequences (seq_name, seq_end)'
    cxn.execute(sql)


def insert_sequences_batch(cxn, batch):
    """Insert a batch of sequence records into the database."""
    sql = """INSERT INTO sequences (seq_name, seq_end, seq)
                VALUES (?, ?, ?)
        """
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


def get_sequence_count(cxn):
    """Get the number of sequences in the table."""
    result = cxn.execute('SELECT COUNT(*) FROM sequences')
    return result.fetchone()[0]


def get_shard_cut(cxn, offset):
    """Get the sequence name at the given offset."""
    sql = 'SELECT seq_name FROM sequences ORDER BY seq_name LIMIT 1 OFFSET {}'
    result = cxn.execute(sql.format(offset))
    cut = result.fetchone()[0]
    return cut


def get_sequences_in_shard(cxn, start, end):
    """Get all sequences in a shard."""
    sql = """
        SELECT seq_name, seq_end, seq
          FROM sequences
         WHERE seq_name >= ?
           AND seq_name < ?
        """
    return cxn.execute(sql, (start, end))
