"""Testing functions in lib/log."""

# pylint: disable=invalid-name


from os.path import basename, splitext
import lib.log as log


log_file = 'my_log_file'
blast_db = 'my_blast_db'
query_file = 'my_query_file'


def test_file_name_01():
    """It returns the given log file name."""
    actual = log.file_name(log_file, blast_db, query_file=query_file)
    assert log_file == actual


def test_file_name_02():
    """It returns a default log file name."""
    actual = log.file_name('', blast_db)
    expect = '{}.pytest.log'.format(blast_db)
    assert expect == actual


def test_file_name_03():
    """It adds the query file to the log file name."""
    actual = log.file_name('', blast_db, query_file=query_file)
    query = splitext(basename(query_file))[0]
    expect = '{}.{}.pytest.log'.format(blast_db, query)
    assert expect == actual
