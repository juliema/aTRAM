"""Database functions for the preprocessor."""

import os
from os.path import join

from .db import DB_VERSION


def aux_db(cxn, temp_dir):
    """Create & attach a temporary database to the current DB connection."""
    db_dir = join(temp_dir, 'db')
    os.makedirs(db_dir, exist_ok=True)

    db_name = join(db_dir, 'temp.sqlite.db')

    sql = """ATTACH DATABASE '{}' AS aux""".format(db_name)
    cxn.execute(sql)


def aux_detach(cxn):
    """Detach the temporary database."""
    cxn.execute('DETACH DATABASE aux')


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
        sql = """INSERT INTO metadata (label, value) VALUES (?, ?);"""
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
    cxn.executescript("""
        CREATE INDEX sequences_index ON sequences (seq_name, seq_end);
        """)


def insert_sequences_batch(cxn, batch):
    """Insert a batch of sequence records into the database."""
    sql = """INSERT INTO sequences (seq_name, seq_end, seq)
                VALUES (?, ?, ?);"""
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


# ########################## sequence names ##################################

def create_seq_names_table(cxn):
    """Create the sequence names table and index."""
    cxn.executescript("""
        CREATE TABLE aux.seq_names AS SELECT DISTINCT seq_name FROM sequences;
        CREATE INDEX aux.name_index ON seq_names (seq_name);
        """)


def get_sequences_in_shard(cxn, shard_count, shard_index):
    """Split the sequences by row ID to shuffle them into different shards."""
    sql = """
        SELECT seq_name, seq_end, seq
          FROM sequences
         WHERE seq_name IN (
               SELECT seq_name FROM aux.seq_names WHERE (rowid % ?) = ?);
        """
    return cxn.execute(sql, (shard_count, shard_index))
