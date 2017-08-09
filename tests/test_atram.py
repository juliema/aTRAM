"""Testing functions in atram.py"""

# pylint: disable=missing-docstring
# pylint: disable=unused-argument
# pylint: disable=global-statement
# pylint: disable=too-few-public-methods

import atram


# ============================================================================
# Mock functions in the db module

# Mocking the sequences table. It will hold the tuples that get inserted into
# the table.

WRITE_SEQ_DATA = []
WRITE_SEQ_COUNT = 0


def mock_write_seq(file_name, seq_id, seq):
    global WRITE_SEQ_DATA, WRITE_SEQ_COUNT
    WRITE_SEQ_COUNT += 1
    WRITE_SEQ_DATA.append({'file': file_name, 'id': seq_id, 'seq': seq})


# ============================================================================
# Tests start here


def test_split_queries(monkeypatch):
    global WRITE_SEQ_DATA, WRITE_SEQ_COUNT

    monkeypatch.setattr(atram, 'write_seq', mock_write_seq)

    file_names = ['tests/data/split_queries1.txt']

    queries = list(atram.split_queries(
        {'temp_file': 'tests', 'query': file_names}))

    assert WRITE_SEQ_COUNT == 3

    for query in queries:
        assert query == 'tests/sequence_01.fasta'

    for i, datum in enumerate(WRITE_SEQ_DATA, 1):
        assert datum['file'] == 'tests/sequence_01.fasta'
        assert datum['id'] == 'seq{}'.format(i)
        assert datum['seq'] == ['', 'A', 'C', 'G'][i] * (i * 10)
