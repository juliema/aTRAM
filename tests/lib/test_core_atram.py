"""Testing functions in core_atram."""

# pylint: disable=too-many-arguments,unused-variable

from os.path import join
from unittest.mock import patch, MagicMock, call
import tempfile
import lib.core_atram as core_atram
from lib.assemblers.base import BaseAssembler


def set_up():
    """Setup a generic assembler."""
    cxn = 'cxn'
    args = {
        'query': ['query_file_1', 'query_file_2'],
        'blast_db': ['blast_db_1', 'blast_db_2'],
        'iterations': 1,
        'log_file': 'log_file_1',
        'temp_dir': 'temp_dir_1'}
    assembler = BaseAssembler(args, cxn)
    return args, cxn, assembler

@patch('lib.assembler.factory')
@patch('lib.db.connect')
@patch('lib.db.aux_db')
@patch('lib.db.aux_detach')
@patch('lib.log.setup')
@patch('lib.core_atram.split_queries')
@patch('lib.core_atram.clean_database')
@patch('lib.core_atram.assembly_loop')
def test_assemble_01(assembly_loop, clean_database, split_queries,
                     setup, aux_detach, aux_db, connect, factory):
    """It runs the assemble function."""
    args, cxn, assembler = set_up()
    connect.return_value.__enter__ = lambda x: 'my_db'
    assembler.write_final_output = MagicMock()
    split_queries.return_value = args['query']
    factory.return_value = assembler

    core_atram.assemble(args)

    split_queries.assert_called_once_with(args)

    calls = [
        call(args['blast_db'][0], check_version=True),
        call().__exit__(None, None, None),
        call(args['blast_db'][1], check_version=True),
        call().__exit__(None, None, None)]
    connect.assert_has_calls(calls)

    calls = [call('my_db'), call('my_db'), call('my_db'), call('my_db')]
    clean_database.assert_has_calls(calls)

    log_file = args['log_file']
    blast_db = args['blast_db']
    query = args['query']

    calls = [
        call('my_db', args['temp_dir'], blast_db[0], query[0]),
        call('my_db', args['temp_dir'], blast_db[0], query[1]),
        call('my_db', args['temp_dir'], blast_db[1], query[0]),
        call('my_db', args['temp_dir'], blast_db[1], query[1])]
    aux_db.assert_has_calls(calls)

    calls = [call('my_db'), call('my_db'), call('my_db'), call('my_db')]
    aux_detach.assert_has_calls(calls)

    calls = [
        call(log_file, blast_db[0], query[0]),
        call(log_file, blast_db[0], query[1]),
        call(log_file, blast_db[1], query[0]),
        call(log_file, blast_db[1], query[1])]
    setup.assert_has_calls(calls)

    calls = [
        call(assembler, blast_db[0], query[0]),
        call(assembler, blast_db[0], query[1]),
        call(assembler, blast_db[1], query[0]),
        call(assembler, blast_db[1], query[1])]
    assembly_loop.assert_has_calls(calls)

    calls = [
        call(blast_db[0], query[0]),
        call(blast_db[0], query[1]),
        call(blast_db[1], query[0]),
        call(blast_db[1], query[1])]
    assembler.write_final_output.assert_has_calls(calls)

@patch('lib.core_atram.write_query_seq')
def test_split_queries_01(write_query_seq):
    """Test split queries where there are no fasta files to split."""
    args, cxn, _ = set_up()
    args['query_split'] = []

    queries = core_atram.split_queries(args)

    write_query_seq.assert_not_called()

    assert args['query'] == queries

@patch('lib.core_atram.write_query_seq')
def test_split_queries_02(write_query_seq):
    """Test split queries where there are fasta files to split."""
    args, cxn, assembler = set_up()
    args['query_split'] = ['tests/data/split_queries1.txt']

    queries = core_atram.split_queries(args)

    split_files = [
        'temp_dir_1/queries/split_queries1_seq1_1_1.fasta',
        'temp_dir_1/queries/split_queries1_seq2_2_2_2.fasta',
        'temp_dir_1/queries/split_queries1_seq3_3.fasta',
        'temp_dir_1/queries/split_queries1_seq1_1_4.fasta']

    calls = [
        call(split_files[0], 'seq1/1', 'A' * 10),
        call(split_files[1], 'seq2:2/2', 'C' * 20),
        call(split_files[2], 'seq3', 'G' * 30),
        call(split_files[3], 'seq1+1', 'T' * 10)]
    write_query_seq.assert_has_calls(calls)

    assert split_files == queries

def test_write_query_seq_01():
    """It writes a sequence to a fasta file."""
    args, cxn, assembler = set_up()
    with tempfile.TemporaryDirectory(prefix='test_') as temp_dir:
        path = join(temp_dir, 'test_query.fasta')

        core_atram.write_query_seq(
            path,
            'my sequence name',
            'aaaacccgggtt')

        with open(path) as test_file:
            expect = (
                '>my sequence name\n'
                'aaaacccgggtt\n')
            assert expect == test_file.read()

@patch('lib.db.create_sra_blast_hits_table')
@patch('lib.db.create_contig_blast_hits_table')
@patch('lib.db.create_assembled_contigs_table')
def test_clean_database_01(
        create_assembled_contigs_table,
        create_contig_blast_hits_table,
        create_sra_blast_hits_table):
    """It runs the clean_database function."""
    args, cxn, assembler = set_up()
    dbh = 'my_db'
    core_atram.clean_database(dbh)

    create_assembled_contigs_table.assert_called_once_with(dbh)
    create_contig_blast_hits_table.assert_called_once_with(dbh)
    create_sra_blast_hits_table.assert_called_once_with(dbh)

@patch('lib.log.info')
@patch('os.makedirs')
@patch('lib.core_atram.blast_query_against_all_shards')
@patch('lib.core_atram.create_query_from_contigs')
@patch('lib.core_atram.filter_contigs')
def test_assembly_loop_01(
        filter_contigs,
        create_query_from_contigs,
        blast_query_against_all_shards,
        b_makedirs,
        info):
    """It iterates over the assembly processes."""
    args, _, assembler = set_up()

    iter_dir = 'my_iter_dir'

    assembler.blast_only = False
    assembler.state['query_file'] = args['query'][0]
    assembler.state['blast_db'] = args['blast_db'][0]

    assembler.initialize_iteration = MagicMock()
    assembler.iter_dir = MagicMock(return_value=iter_dir)
    assembler.count_blast_hits = MagicMock(return_value=1)
    assembler.write_input_files = MagicMock()
    assembler.run = MagicMock()
    assembler.nothing_assembled = MagicMock(return_value=False)
    assembler.assembled_contigs_count = MagicMock(return_value=11)
    assembler.no_new_contigs = MagicMock(return_value=False)

    core_atram.assembly_loop(
        assembler, args['blast_db'][0], args['query'][0])

    calls = [
        call('aTRAM blast DB = "{}", query = "{}", iteration {}'.format(
            args['blast_db'][0], args['query'][0], 1)),
        call('All iterations completed')]
    info.assert_has_calls(calls)

    assembler.initialize_iteration.assert_called_once_with(
        args['blast_db'][0], args['query'][0], 1)
    assert assembler.iter_dir.call_count == 1
    b_makedirs.assert_called_once_with(iter_dir, exist_ok=True)
    blast_query_against_all_shards.assert_called_once_with(assembler)
    assert assembler.count_blast_hits.call_count == 1
    assembler.no_new_contigs.assert_called_once_with(11)

    create_query_from_contigs.create_query_from_contigs(assembler)
    filter_contigs.create_query_from_contigs(assembler)

# def test_assembly_loop_multiple_iterations():

# def test_assembly_loop_blast_only():

# def test_assembly_loop_no_blast_hits():

# def test_assembly_loop_nothing_assembled():

# def test_assembly_loop_no_assembled_contigs_count():

# def test_assembly_loop_no_new_contigs():

@patch('lib.blast.all_shard_paths')
def test_shard_fraction_01(all_shard_paths):
    """It gets the shards we are using when there is no split."""
    args, cxn, assembler = set_up()
    returns = ['1st', '2nd', '3rd', '4th']
    assembler.state['blast_db'] = args['blast_db'][0]
    assembler.args['fraction'] = 1.0

    all_shard_paths.return_value = returns

    shards = core_atram.shard_fraction(assembler)

    assert returns == shards

    all_shard_paths.assert_called_once_with(args['blast_db'][0])

@patch('lib.blast.all_shard_paths')
def test_shard_fraction_02(all_shard_paths):
    """It gets the shards we are using when there is a split."""
    args, cxn, assembler = set_up()
    assembler.args['fraction'] = 0.5
    assembler.state['blast_db'] = args['blast_db'][0]
    returns = ['1st', '2nd', '3rd', '4th']

    all_shard_paths.return_value = returns

    shards = core_atram.shard_fraction(assembler)

    assert ['1st', '2nd'] == shards

    all_shard_paths.assert_called_once_with(args['blast_db'][0])

# def test_blast_query_against_one_shard():

# def test_filter_contigs():

# def test_save_blast_against_contigs():

# def test_save_contigs():


# def test_create_query_from_contigs(self):
