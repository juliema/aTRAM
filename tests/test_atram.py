"""Testing functions in atram."""

# pylama: ignore=D103

import atram
import tests.mock as mock


def test_split_queries():
    mock.it(atram, 'write_query_seq')

    file_names = ['tests/data/split_queries1.txt']

    queries = atram.split_queries(
        {'temp_dir': 'temp_dir', 'query': file_names})

    fasta = [
        'temp_dir/queries/split_queries1_seq1_1_1.fasta',
        'temp_dir/queries/split_queries1_seq2_2_2_2.fasta',
        'temp_dir/queries/split_queries1_seq3_3.fasta',
        'temp_dir/queries/split_queries1_seq1_1_4.fasta']

    assert queries == fasta

    assert mock.history == [
        {'module': 'atram', 'func': 'write_query_seq',
            'file_name': fasta[0], 'seq_id': 'seq1/1', 'seq': 'A' * 10},
        {'module': 'atram', 'func': 'write_query_seq',
            'file_name': fasta[1], 'seq_id': 'seq2:2/2', 'seq': 'C' * 20},
        {'module': 'atram', 'func': 'write_query_seq',
            'file_name': fasta[2], 'seq_id': 'seq3', 'seq': 'G' * 30},
        {'module': 'atram', 'func': 'write_query_seq',
            'file_name': fasta[3], 'seq_id': 'seq1+1', 'seq': 'T' * 10}]
