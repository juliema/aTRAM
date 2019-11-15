"""Handle common database functions."""

import sqlite3
import sys
import os
from os.path import basename, join, exists


ATRAM_VERSION = 'v2.3.0'

# DB_VERSION != ATRAM_VERSION
# We don't force DB changes until required.
# Therefore DB_VERSION <= ATRAM_VERSION.
DB_VERSION = '2.0'

BATCH_SIZE = 1e6  # How many sequence records to insert at a time


def connect(blast_db, check_version=False, clean=False):
    """Create DB connection."""
    db_name = get_db_name(blast_db)

    if clean and exists(db_name):
        os.remove(db_name)

    if check_version and not exists(db_name):
        err = 'Could not find the database file "{}".'.format(db_name)
        sys.exit(err)

    if check_version:
        with db_setup(db_name) as cxn:
            check_versions(cxn)

    return db_setup(db_name)


def get_db_name(db_prefix):
    """Build the SQLite DB name from the prefix."""
    return '{}.sqlite.db'.format(db_prefix)


def aux_db(cxn, temp_dir, blast_db, query_name):
    """Create & attach a temporary database to the current DB connection."""
    db_dir = join(temp_dir, 'db')
    os.makedirs(db_dir, exist_ok=True)

    db_name = '{}_{}_temp.sqlite.db'.format(
        basename(blast_db), basename(query_name))
    db_name = join(db_dir, db_name)

    sql = """ATTACH DATABASE '{}' AS aux""".format(db_name)
    cxn.execute(sql)


def subprocess_db_name(temp_dir, blast_db, query_name, shard):
    """Build the subprocess database name."""
    db_name = '{}_{}_{}_temp.sqlite.db'.format(
        basename(blast_db), basename(query_name), basename(shard))
    db_name = join(temp_dir, 'db', db_name)
    return db_name


def subprocess_db(temp_dir, blast_db, query_name, shard):
    """Create a temporary database for a shard/subprocess."""
    db_name = subprocess_db_name(temp_dir, blast_db, query_name, shard)
    return db_setup(db_name)


def subprocess_attach(cxn, temp_dir, blast_db, query_name, shard):
    """Attach the subprocess database to the main database."""
    db_name = subprocess_db_name(temp_dir, blast_db, query_name, shard)
    sql = """ATTACH DATABASE '{}' AS subprocess""".format(db_name)
    cxn.execute(sql)


def subprocess_detach(cxn):
    """Detach the temporary database."""
    cxn.execute('DETACH DATABASE subprocess')


def aux_detach(cxn):
    """Detach the temporary database."""
    cxn.execute('DETACH DATABASE aux')


def temp_db(temp_dir, db_prefix):
    """Create a temporary database."""
    db_name = join(temp_dir, get_db_name(db_prefix))
    return db_setup(db_name)


def db_setup(db_name):
    """Database setup."""
    cxn = sqlite3.connect(db_name, timeout=30.0)
    cxn.execute("PRAGMA page_size = {}".format(2**16))
    cxn.execute("PRAGMA journal_mode = WAL")
    return cxn


# ########################### misc functions #################################

def check_versions(cxn):
    """Make sure the database version matches what we built it with."""
    version = get_version(cxn)
    if version != DB_VERSION:
        err = ('The database was built with version {} but you are running '
               'version {}. You need to rebuild the atram database by '
               'running atram_preprocessor.py again.').format(
                   version, DB_VERSION)
        sys.exit(err)


# ########################## metadata table ##################################

def get_metadata(cxn, key, default=''):
    """Get the current database version."""
    sql = """SELECT value FROM metadata WHERE label = ?"""
    try:
        result = cxn.execute(sql, (key, ))
        result = result.fetchone()
        return default if not result else result[0]
    except sqlite3.OperationalError:
        return default


def get_version(cxn):
    """Get the current database version."""
    return get_metadata(cxn, 'version', default='1.0')


def is_single_end(cxn):
    """Was the database build for single ends."""
    result = get_metadata(cxn, 'single_ends')
    return result != '0'


# ########################## sequences table ##################################

def get_sequence_ends(cxn):
    """Get a list of all seq_ends in the database."""
    return cxn.execute('SELECT DISTINCT seq_end FROM sequences')


def get_all_sequences(cxn):
    """Get a list of all seq_ends in the database."""
    return cxn.execute('SELECT * FROM sequences')
