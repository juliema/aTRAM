"""Testing functions in lib/log."""

# pylama: ignore=D103

import lib.log as log


def test_file_name_default():
    log_file = log.file_name('', 'blast_db', 'query_file')
    assert log_file == 'blast_db.query_file.pytest.log'


def test_file_name_arg():
    file_name = log.file_name('test.log', 'blast_db', 'query_file')
    assert file_name == 'test.log'


def test_info_default():
    pass
