"""Testing functions in atram."""

# pylama: ignore=D103

import atram


# ============================================================================
# Mock functions in atram itself

# Mocking write_seq so that we don't need to write to files

WRITE_QUERY_SEQ_DATA = []
WRITE_QUERY_SEQ_COUNT = 0


def mock_write_query_seq(file_name, seq_id, seq):
    global WRITE_QUERY_SEQ_DATA, WRITE_QUERY_SEQ_COUNT
    WRITE_QUERY_SEQ_COUNT += 1
    WRITE_QUERY_SEQ_DATA.append({'file': file_name, 'id': seq_id, 'seq': seq})


# ============================================================================
# Tests start here


def test_split_queries(monkeypatch):
    global WRITE_QUERY_SEQ_DATA, WRITE_QUERY_SEQ_COUNT

    monkeypatch.setattr(atram, 'write_query_seq', mock_write_query_seq)

    file_names = ['tests/data/split_queries1.txt']

    queries = atram.split_queries(
        {'temp_dir': 'temp_dir', 'query': file_names})

    fasta = [
        'temp_dir/queries/split_queries1_seq1_1_1.fasta',
        'temp_dir/queries/split_queries1_seq2_2_2_2.fasta',
        'temp_dir/queries/split_queries1_seq3_3.fasta',
        'temp_dir/queries/split_queries1_seq1_1_4.fasta']

    assert WRITE_QUERY_SEQ_COUNT == 4

    assert queries == fasta

    assert WRITE_QUERY_SEQ_DATA == [
        {'file': fasta[0], 'id': 'seq1/1', 'seq': 'A' * 10},
        {'file': fasta[1], 'id': 'seq2:2/2', 'seq': 'C' * 20},
        {'file': fasta[2], 'id': 'seq3', 'seq': 'G' * 30},
        {'file': fasta[3], 'id': 'seq1+1', 'seq': 'T' * 10}]
