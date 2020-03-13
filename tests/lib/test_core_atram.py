"""Testing functions in core_atram."""

# pylint: disable=too-many-arguments,unused-variable

from os.path import join
from unittest.mock import patch, MagicMock, call
import tempfile
import lib.core_atram as core_atram
from lib.assemblers.base import BaseAssembler


def set_up():
    """Build a generic assembler."""
    cxn = 'cxn'
    args = {
        'query': ['query_file_1', 'query_file_2'],
        'blast_db': ['blast_db_1', 'blast_db_2'],
        'iterations': 1,
        'log_file': 'log_file_1',
        'log_level': 'info',
        'temp_dir': 'temp_dir_1'}
    assembler = BaseAssembler(args, cxn)
    return args, cxn, assembler


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
    args['protein'] = True

    with tempfile.TemporaryDirectory(prefix='test_') as temp_dir:
        args['temp_dir'] = temp_dir
        queries = core_atram.split_queries(args)

    split_files = [
        join(temp_dir, 'queries', 'split_queries1_seq1_1_1.fasta'),
        join(temp_dir, 'queries', 'split_queries1_seq2_2_2_2.fasta'),
        join(temp_dir, 'queries', 'split_queries1_seq3_3.fasta'),
        join(temp_dir, 'queries', 'split_queries1_seq1_1_4.fasta')]

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


@patch('lib.db_atram.create_sra_blast_hits_table')
@patch('lib.db_atram.create_contig_blast_hits_table')
@patch('lib.db_atram.create_assembled_contigs_table')
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


@patch('lib.core_atram.blast_query_against_all_shards')
@patch('lib.core_atram.create_query_from_contigs')
@patch('lib.core_atram.filter_contigs')
def test_assembly_loop_iteration_01(
        filter_contigs,
        create_query_from_contigs,
        blast_query_against_all_shards):
    """It iterates over the assembly processes."""
    args, _, assembler = set_up()

    temp_dir = 'my_temp_dir'

    assembler.blast_only = False
    assembler.state['query_file'] = args['query'][0]
    assembler.state['blast_db'] = args['blast_db'][0]
    assembler.state['iter_dir'] = 'my_iter_dir'

    assembler.init_iteration = MagicMock()
    assembler.count_blast_hits = MagicMock(return_value=1)
    assembler.write_input_files = MagicMock()
    assembler.run = MagicMock()
    assembler.nothing_assembled = MagicMock(return_value=False)
    assembler.assembled_contigs_count = MagicMock(return_value=11)
    assembler.no_new_contigs = MagicMock(return_value=False)

    core_atram.assembly_loop_iteration(args, assembler)

    blast_query_against_all_shards.assert_called_once_with(assembler)
    assert assembler.count_blast_hits.call_count == 1
    assembler.no_new_contigs.assert_called_once_with(11)

    create_query_from_contigs.create_query_from_contigs(assembler)
    filter_contigs.create_query_from_contigs(assembler)


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
