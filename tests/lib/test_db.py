"""Testing functions in lib/db."""

import sqlite3
import lib.db as db
import lib.db_atram as db_atram
import lib.db_preprocessor as db_preprocessor


CXN = sqlite3.connect(':memory:')


def setUpModule():
    """Build the database for testing."""
    CXN.execute("""ATTACH DATABASE ':memory:' AS aux""")
    db_preprocessor.create_metadata_table(CXN, {})
    db_preprocessor.create_sequences_table(CXN)
    db_atram.create_sra_blast_hits_table(CXN)
    db_atram.create_contig_blast_hits_table(CXN)
    db_atram.create_assembled_contigs_table(CXN)


def test_get_db_name_01():
    """It prepends the blast_db name to the database name."""
    assert db.get_db_name('test_db') == 'test_db.sqlite.db'


def test_get_version_01():
    """It returns the current DB version."""
    assert db.get_version(CXN) == '2.0'


def test_get_version_02():
    """It returns a default version if there is no metadata table."""
    CXN.execute("""DROP TABLE IF EXISTS metadata""")
    assert db.get_version(CXN) == '1.0'
