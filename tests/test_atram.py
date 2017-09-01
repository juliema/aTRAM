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
        hist = {'blast_db': args[0], 'query': args[1]}
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
    mock.it(log, 'setup')
    mock.it(assembly, 'factory', MockAssembler())
    mock.it(atram, 'assembly_loop')

    atram.assemble(args)

    expect = [{'db_conn': 'my_connection'}] * 4
    assert expect == mock.filter('atram', 'clean_database')

    expect = [
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
    assert expect == mock.filter('lib.log', 'setup')

    expect = [{'args': args, 'db_conn': 'my_connection'}] * 4
    assert expect == mock.filter('lib.assembler', 'factory')

    expect = [
        {'blast_db': 'my_blast_db1', 'query': 'my_query1'},
        {'blast_db': 'my_blast_db1', 'query': 'my_query2'},
        {'blast_db': 'my_blast_db2', 'query': 'my_query1'},
        {'blast_db': 'my_blast_db2', 'query': 'my_query2'}]
    assert expect == mock.filter('MockAssembler', 'write_final_output')


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

    assert fasta == queries

    expect = [
        {'file_name': fasta[0], 'seq_id': 'seq1/1', 'seq': 'A' * 10},
        {'file_name': fasta[1], 'seq_id': 'seq2:2/2', 'seq': 'C' * 20},
        {'file_name': fasta[2], 'seq_id': 'seq3', 'seq': 'G' * 30},
        {'file_name': fasta[3], 'seq_id': 'seq1+1', 'seq': 'T' * 10}]
    assert expect == mock.filter('atram', 'write_query_seq')


def test_write_query_seq():
    with tempfile.TemporaryDirectory(prefix='test_') as temp_dir:
        path = join(temp_dir, 'test_query.fasta')

        atram.write_query_seq(path, 'my sequence name', 'aaaacccgggtt')

        with open(path) as test_file:
            expect = (
                '>my sequence name\n'
                'aaaacccgggtt\n')
            assert expect == test_file.read()


def test_clean_database():
    mock.it(db, 'create_sra_blast_hits_table')
    mock.it(db, 'create_contig_blast_hits_table')
    mock.it(db, 'create_assembled_contigs_table')

    atram.clean_database('my_db_conn')

    expect = [
        {'func': 'create_sra_blast_hits_table', 'db_conn': 'my_db_conn'},
        {'func': 'create_contig_blast_hits_table', 'db_conn': 'my_db_conn'},
        {'func': 'create_assembled_contigs_table', 'db_conn': 'my_db_conn'}]
    assert expect == mock.filter('lib.db')


# def test_assembly_loop_one_iteration():
#     pass
#
#
# def test_assembly_loop_multiple_iterations():
#     pass
#
#
# def test_assembly_loop_blast_only():
#     pass
#
#
# def test_assembly_loop_no_blast_hits():
#     pass
#
#
# def test_assembly_loop_nothing_assembled():
#     pass
#
#
# def test_assembly_loop_no_assembled_contigs_count():
#     pass
#
#
# def test_assembly_loop_no_new_contigs():
#     pass
#
#
# def test_shard_fraction():
#     pass
#
#
# def test_blast_query_against_one_shard():
#     pass
#
#
# def test_filter_contigs():
#     pass
#
#
# def test_save_blast_against_contigs():
#     pass
#
#
# def test_save_contigs():
#     pass
#
#
# def test_create_query_from_contigs():
#     pass
