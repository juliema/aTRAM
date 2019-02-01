"""Testing functions in lib/db."""

# pylint: disable=invalid-name

import sqlite3
import lib.db as db


cxn = sqlite3.connect(':memory:')


def setUpModule():
    """Setup the database for testing."""
    cxn.execute("""ATTACH DATABASE ':memory:' AS aux""")
    db.create_metadata_table(cxn)
    db.create_sequences_table(cxn)
    db.create_sra_blast_hits_table(cxn)
    db.create_contig_blast_hits_table(cxn)
    db.create_assembled_contigs_table(cxn)

def test_get_db_name_01():
    """It prepends the blast_db name to the database name."""
    assert db.get_db_name('test_db') == 'test_db.sqlite.db'

def test_get_version_01():
    """It returns the current DB version."""
    assert db.get_version(cxn) == '2.0'

def test_get_version_02():
    """It returns a default version if there is no metadata table."""
    cxn.execute("""DROP TABLE IF EXISTS metadata""")
    assert db.get_version(cxn) == '1.0'
