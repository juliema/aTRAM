"""Testing functions in lib/log.py"""

# pylint: disable=missing-docstring

from argparse import Namespace
import lib.log as log


def test_file_name_default():
    args = Namespace(log_file=None)
    file_name = log.file_name(args, 'blast_db', 'query_file')
    assert file_name == 'blast_db.query_file.pytest.log'


def test_file_name_arg():
    args = Namespace(log_file='test.log')
    file_name = log.file_name(args, 'blast_db', 'query_file')
    assert file_name == 'test.log'
