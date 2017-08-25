"""Testing functions in atram."""

# pylama: ignore=D103,D101,D102

from os.path import join
import tempfile
import atram
import lib.db as db
import lib.log as log
import lib.assembler as assembly
import tests.mock as mock


class MockAssembler:
    def write_final_output(self, *args):
        hist = {a: args[i] for i, a in enumerate(['blast_db', 'query'])}
        hist['module'] = 'MockAssembler'
        hist['func'] = 'write_final_output'
        mock.history.append(hist)


def test_assemble():
    args = {
        'split_queries': False,
        'query': ['my_query1', 'my_query2'],
        'blast_db': ['my_blast_db1', 'my_blast_db2'],
        'log_file': 'my_log_file'}
    mock.context(db, 'connect', 'my_connection')
    mock.it(atram, 'clean_database')
    mock.it(log, 'file_name', 'my_new_log_file')
    mock.it(log, 'setup')
    mock.it(assembly, 'factory', MockAssembler())
    mock.it(atram, 'assembly_loop')

    atram.assemble(args)

    history = mock.filter('atram', 'clean_database')
    assert history == [{'db_conn': 'my_connection'}] * 4

    history = mock.filter('lib.log', 'file_name')
    assert history == [
        {'blast_db': 'my_blast_db1',
         'log_file': 'my_log_file',
         'query_file': 'my_query1'},
        {'blast_db': 'my_blast_db1',
         'log_file': 'my_log_file',
         'query_file': 'my_query2'},
        {'blast_db': 'my_blast_db2',
         'log_file': 'my_log_file',
         'query_file': 'my_query1'},
        {'blast_db': 'my_blast_db2',
         'log_file': 'my_log_file',
         'query_file': 'my_query2'}]

    history = mock.filter('lib.log', 'setup')
    assert history == [{'log_file': 'my_new_log_file'}] * 4

    history = mock.filter('lib.assembler', 'factory')
    assert history == [{'args': args, 'db_conn': 'my_connection'}] * 4

    history = mock.filter('MockAssembler', 'write_final_output')
    assert history == [
        {'blast_db': 'my_blast_db1', 'query': 'my_query1'},
        {'blast_db': 'my_blast_db1', 'query': 'my_query2'},
        {'blast_db': 'my_blast_db2', 'query': 'my_query1'},
        {'blast_db': 'my_blast_db2', 'query': 'my_query2'}]


def test_split_queries():
    mock.it(atram, 'write_query_seq')

    file_names = ['tests/data/split_queries1.txt']

    fasta = [
        'temp_dir/queries/split_queries1_seq1_1_1.fasta',
        'temp_dir/queries/split_queries1_seq2_2_2_2.fasta',
        'temp_dir/queries/split_queries1_seq3_3.fasta',
        'temp_dir/queries/split_queries1_seq1_1_4.fasta']

    queries = atram.split_queries(
        {'temp_dir': 'temp_dir', 'query': file_names})

    history = mock.filter('atram', 'write_query_seq')

    assert queries == fasta

    assert history == [
        {'file_name': fasta[0], 'seq_id': 'seq1/1', 'seq': 'A' * 10},
        {'file_name': fasta[1], 'seq_id': 'seq2:2/2', 'seq': 'C' * 20},
        {'file_name': fasta[2], 'seq_id': 'seq3', 'seq': 'G' * 30},
        {'file_name': fasta[3], 'seq_id': 'seq1+1', 'seq': 'T' * 10}]


def test_write_query_seq():
    with tempfile.TemporaryDirectory(prefix='test_') as temp_dir:
        path = join(temp_dir, 'test_query.fasta')

        atram.write_query_seq(path, 'my sequence name', 'aaaacccgggtt')

        with open(path) as test_file:
            assert test_file.read() == (
                '>my sequence name\n'
                'aaaacccgggtt\n')


def test_clean_database():
    mock.it(db, 'create_sra_blast_hits_table')
    mock.it(db, 'create_contig_blast_hits_table')
    mock.it(db, 'create_assembled_contigs_table')

    atram.clean_database('my_db_conn')

    history = mock.filter('lib.db')
    assert history == [
        {'func': 'create_sra_blast_hits_table', 'db_conn': 'my_db_conn'},
        {'func': 'create_contig_blast_hits_table', 'db_conn': 'my_db_conn'},
        {'func': 'create_assembled_contigs_table', 'db_conn': 'my_db_conn'}]
