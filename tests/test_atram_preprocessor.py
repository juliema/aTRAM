"""Testing functions in atram_preprocessor.py"""

# pylint: disable=missing-docstring
# pylint: disable=unused-argument
# pylint: disable=global-statement

# import builtins
from os.path import join
import lib.db as db
import atram_preprocessor


# ============================================================================
# Mock functions in the db module

# Mocking the sequences table. It will hold the tuples that get inserted into
# the table.

INSERT_SEQUENCES_BATCH = []
INSERT_SEQUENCES_BATCH_COUNT = 0


def mock_insert_sequences_batch(db_conn, batch):
    global INSERT_SEQUENCES_BATCH, INSERT_SEQUENCES_BATCH_COUNT
    INSERT_SEQUENCES_BATCH_COUNT += 1
    INSERT_SEQUENCES_BATCH += batch


def mock_get_sequence_count(db_conn):
    return 100


# Mocking get_two_sequences function it returns two consecutive sequence names.

GET_TWO_SEQUENCES = {
    33: ('seq1', 'seq2'),
    66: ('seq3', 'seq3')}


def mock_get_two_sequences(db_conn, offset):
    return GET_TWO_SEQUENCES[offset]


# ============================================================================
# Tests start here

def test_load_seq(monkeypatch):
    global INSERT_SEQUENCES_BATCH, INSERT_SEQUENCES_BATCH_COUNT
    INSERT_SEQUENCES_BATCH = []

    db.BATCH_SIZE = 5

    monkeypatch.setattr(
        db, 'insert_sequences_batch', mock_insert_sequences_batch)

    file_1 = join('tests', 'data', 'load_seq1.txt')
    file_2 = join('tests', 'data', 'load_seq2.txt')
    atram_preprocessor.load_seqs(True, [file_1, file_2])

    assert INSERT_SEQUENCES_BATCH_COUNT == 3
    assert INSERT_SEQUENCES_BATCH == [
        ('seq1', '1', 'AAAAAAAAAA'),
        ('seq2', '1', 'AAAAAAAAAAGGGGGGGGGG'),
        ('seq3', '1', 'AAAAAAAAAA'),
        ('seq4', '1', 'AAAAAAAAAA'),
        ('seq5/3', '', 'AAAAAAAAAAGGGGGGGGGG'),
        ('seq1', '2', 'AAAAAAAAAA'),
        ('seq2', '2', 'AAAAAAAAAAGGGGGGGGGG'),
        ('seq3', '2', 'AAAAAAAAAA'),
        ('seq4', '2', 'AAAAAAAAAAGGGGGGGGGG'),
        ('seq6', '1', 'TTTTTTTTTT'),
        ('seq7', '1', 'TTTTTTTTTTCCCCCCCCCC'),
        ('seq8/a.1 suffix', '', 'TTTTTTTTTT'),
        ('seq8', '2', 'TTTTTTTTTTCCCCCCCCCC')]


def test_assign_seqs_to_shards(monkeypatch):
    monkeypatch.setattr(db, 'get_sequence_count', mock_get_sequence_count)
    monkeypatch.setattr(db, 'get_two_sequences', mock_get_two_sequences)

    shard_count = 3
    shard_list = atram_preprocessor.assign_seqs_to_shards(True, shard_count)

    # A list of pairs of LIMIT and OFFSET for queries that build shards
    assert shard_list == [(34, 0), (32, 34), (34, 66)]
