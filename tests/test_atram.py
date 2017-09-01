"""Testing functions in atram."""

# pylama: ignore=D103,D101,D102

from os.path import join
import tempfile
import atram
import lib.db as db
import lib.log as log
import lib.blast as blast
import lib.file_util as file_util
import lib.assembler as assembly
from lib.assemblers.base import BaseAssembler
import tests.mock as mock


def test_assemble():
    args = {
        'split_queries': False,
        'query': ['my_query1', 'my_query2'],
        'blast_db': ['my_blast_db1', 'my_blast_db2'],
        'log_file': 'my_log_file'}
    assembler = BaseAssembler(args, 'my_connection')
    mock.context(db, 'connect', 'my_connection')
    mock.it(atram, 'clean_database')
    mock.it(log, 'setup')
    mock.it(assembly, 'factory', assembler)
    mock.it(assembler, 'write_final_output')
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

    from pprint import pprint
    pprint(mock.history)
    expect = [
        {'blast_db': 'my_blast_db1', 'query': 'my_query1'},
        {'blast_db': 'my_blast_db1', 'query': 'my_query2'},
        {'blast_db': 'my_blast_db2', 'query': 'my_query1'},
        {'blast_db': 'my_blast_db2', 'query': 'my_query2'}]
    assert expect == mock.filter('BaseAssembler', 'write_final_output')


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


def test_assembly_loop_one_iter():
    args = {'temp_dir': 'my_temp_dir', 'iterations': 1}
    assembler = BaseAssembler(args, 'my_connection')
    assembler.blast_only = False

    mock.it(log, 'info')
    mock.it(file_util, 'temp_iter_dir')
    mock.it(assembler, 'initialize_iteration')
    mock.it(atram, 'blast_query_against_all_shards')
    mock.it(assembler, 'no_blast_hits', False)
    mock.it(assembler, 'write_input_files')
    mock.it(assembler, 'run')
    mock.it(assembler, 'nothing_assembled', False)
    mock.it(atram, 'filter_contigs', 99)
    mock.it(assembler, 'assembled_contigs_count', 11)
    mock.it(assembler, 'no_new_contigs', False)
    mock.it(atram, 'create_query_from_contigs', 'my_query')

    atram.assembly_loop(args, 'a_blast_db', 'a_query', 'a_db_conn', assembler)

    expect = [{'msg': ('aTRAM blast DB = "a_blast_db", '
                       'query = "a_query", iteration 1')},
              {'msg': 'All iterations completed'}]
    assert expect == mock.filter('lib.log', 'info')

    expect = [{'temp_dir': 'my_temp_dir',
               'iteration': 1,
               'query_file': 'a_query',
               'blast_db': 'a_blast_db'}]
    assert expect == mock.filter('lib.file_util', 'temp_iter_dir')

    expect = [
        {'blast_db': 'a_blast_db', 'iteration': 1, 'query_file': 'a_query'}]
    assert expect == mock.filter('BaseAssembler', 'initialize_iteration')

    expect = [{'args': {'temp_dir': 'my_temp_dir', 'iterations': 1},
               'blast_db': 'a_blast_db',
               'query': 'a_query',
               'iteration': 1}]
    assert expect == mock.filter('atram', 'blast_query_against_all_shards')

    assert [{}] == mock.filter('BaseAssembler', 'no_blast_hits')
    assert [{}] == mock.filter('BaseAssembler', 'write_input_files')
    assert [{}] == mock.filter('BaseAssembler', 'run')
    assert [{}] == mock.filter('BaseAssembler', 'nothing_assembled')
    assert [{'count': 11}] == mock.filter('BaseAssembler', 'no_new_contigs')

    expect = [{'args': {'temp_dir': 'my_temp_dir', 'iterations': 1},
               'db_conn': 'a_db_conn',
               'assembler': assembler,
               'iteration': 1}]
    assert expect == mock.filter('atram', 'create_query_from_contigs')


# def test_assembly_loop_multiple_iterations():


# def test_assembly_loop_blast_only():


# def test_assembly_loop_no_blast_hits():


# def test_assembly_loop_nothing_assembled():


# def test_assembly_loop_no_assembled_contigs_count():


# def test_assembly_loop_no_new_contigs():


def test_shard_fraction_one():
    args = {'fraction': 1.0}
    returns = ['1st', '2nd', '3rd', '4th']
    mock.it(blast, 'all_shard_paths', returns=[returns])

    shards = atram.shard_fraction(args, 'my_blast_db')

    assert returns == shards

    expect = [{'blast_db': 'my_blast_db'}]
    assert expect == mock.filter('lib.blast', 'all_shard_paths')


def test_shard_fraction_half():
    args = {'fraction': 0.5}
    mock.it(blast, 'all_shard_paths', returns=[['1st', '2nd', '3rd', '4th']])

    shards = atram.shard_fraction(args, 'my_blast_db')

    assert ['1st', '2nd'] == shards

    expect = [{'blast_db': 'my_blast_db'}]
    assert expect == mock.filter('lib.blast', 'all_shard_paths')


# def test_blast_query_against_one_shard():


# def test_filter_contigs():


# def test_save_blast_against_contigs():


# def test_save_contigs():


# def test_create_query_from_contigs():
