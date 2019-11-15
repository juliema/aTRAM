"""Testing functions in core_preprocessor."""

# pylint: disable=too-many-arguments,too-many-locals,unused-argument

from os.path import join
import tempfile
from unittest.mock import patch, call
import lib.db as db
import lib.db_preprocessor as db_preprocessor
import lib.core_preprocessor as core_preprocessor


def set_up():
    """Build structures for the tests."""
    args = {
        'blast_db': 'blast_db_1',
        'log_file': 'log_file_1',
        'sra_files': 'sra_files_1',
        'shard_count': 4,
        'temp_dir': 'temp_dir_1',
        'keep_temp_dir': False}
    return 'cxn', args


@patch('lib.util.make_temp_dir')
@patch('lib.util.update_temp_dir')
@patch('lib.log.info')
@patch('lib.log.setup')
@patch('lib.db.connect')
@patch('lib.db_preprocessor.create_metadata_table')
@patch('lib.db_preprocessor.create_sequences_table')
@patch('lib.db_preprocessor.create_sequences_index')
@patch('lib.core_preprocessor.load_seqs')
@patch('lib.core_preprocessor.assign_seqs_to_shards')
@patch('lib.core_preprocessor.create_all_blast_shards')
def test_preprocess_01(
        create_all_blast_shards, assign_seqs_to_shards, load_seqs,
        create_sequences_index, create_sequences_table,
        create_metadata_table, connect, setup, info, update_temp_dir,
        make_temp_dir):
    """It runs the function to build the databases required by atram."""
    _, args = set_up()
    shard_list = ['shard_1', 'shard_4', 'shard_3', 'shard_4']
    dbh = 'my_db'
    assign_seqs_to_shards.return_value = shard_list
    connect.return_value.__enter__ = lambda x: dbh

    core_preprocessor.preprocess(args)

    setup.assert_called_once_with(
        args['log_file'], args['blast_db'])

    calls = [
        call(args['blast_db'], clean=True),
        call().__exit__(None, None, None)]
    connect.assert_has_calls(calls)

    create_metadata_table.assert_called_once_with(dbh, args)
    create_sequences_table.assert_called_once_with(dbh)
    load_seqs.assert_called_once_with(args, dbh)
    info.assert_called_once_with(
        'Creating an index for the sequence table')
    create_sequences_index.assert_called_once_with(dbh)
    assign_seqs_to_shards.assert_called_once_with(
        dbh, args['shard_count'])
    create_all_blast_shards.assert_called_once_with(args, shard_list)


@patch('lib.log.info')
@patch('lib.db_preprocessor.insert_sequences_batch')
def test_load_one_file_01(insert_sequences_batch, info):
    """Load mixed sequences from a file into the atram database."""
    cxn, args = set_up()
    db.BATCH_SIZE = 5

    file_1 = join('tests', 'data', 'load_seq1.txt')
    core_preprocessor.load_one_file(args, cxn, file_1, 'mixed_ends')

    msg = 'Loading "{}" into sqlite database'.format(file_1)
    info.assert_called_once_with(msg)

    calls = [
        call(cxn, [
            ('seq1', '1', 'AAAAAAAAAA'),
            ('seq2', '1', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3', '1', 'AAAAAAAAAA'),
            ('seq4', '1', 'AAAAAAAAAA'),
            ('seq5/3', '', 'AAAAAAAAAAGGGGGGGGGG')]),
        call(cxn, [
            ('seq1', '2', 'AAAAAAAAAA'),
            ('seq2', '2', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3', '2', 'AAAAAAAAAA'),
            ('seq4', '2', 'AAAAAAAAAAGGGGGGGGGG')])]
    insert_sequences_batch.assert_has_calls(calls)


@patch('lib.log.info')
@patch('lib.db_preprocessor.insert_sequences_batch')
def test_load_one_file_02(insert_sequences_batch, info):
    """Load end 1 sequences from a file into the atram database."""
    cxn, args = set_up()
    db.BATCH_SIZE = 5

    file_1 = join('tests', 'data', 'load_seq1.txt')
    core_preprocessor.load_one_file(args, cxn, file_1, 'end_1', '1')

    msg = 'Loading "{}" into sqlite database'.format(file_1)
    info.assert_called_once_with(msg)

    calls = [
        call(cxn, [
            ('seq1', '1', 'AAAAAAAAAA'),
            ('seq2', '1', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3', '1', 'AAAAAAAAAA'),
            ('seq4', '1', 'AAAAAAAAAA'),
            ('seq5/3', '1', 'AAAAAAAAAAGGGGGGGGGG')]),
        call(cxn, [
            ('seq1', '1', 'AAAAAAAAAA'),
            ('seq2', '1', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3', '1', 'AAAAAAAAAA'),
            ('seq4', '1', 'AAAAAAAAAAGGGGGGGGGG')])]
    insert_sequences_batch.assert_has_calls(calls)


@patch('lib.log.info')
@patch('lib.db_preprocessor.insert_sequences_batch')
def test_load_one_file_03(insert_sequences_batch, info):
    """Load end 2 sequences from a file into the atram database."""
    cxn, args = set_up()
    db.BATCH_SIZE = 5

    file_1 = join('tests', 'data', 'load_seq1.txt')
    core_preprocessor.load_one_file(args, cxn, file_1, 'end_2', '2')

    msg = 'Loading "{}" into sqlite database'.format(file_1)
    info.assert_called_once_with(msg)

    calls = [
        call(cxn, [
            ('seq1', '2', 'AAAAAAAAAA'),
            ('seq2', '2', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3', '2', 'AAAAAAAAAA'),
            ('seq4', '2', 'AAAAAAAAAA'),
            ('seq5/3', '2', 'AAAAAAAAAAGGGGGGGGGG')]),
        call(cxn, [
            ('seq1', '2', 'AAAAAAAAAA'),
            ('seq2', '2', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3', '2', 'AAAAAAAAAA'),
            ('seq4', '2', 'AAAAAAAAAAGGGGGGGGGG')])]
    insert_sequences_batch.assert_has_calls(calls)


@patch('lib.log.info')
@patch('lib.db_preprocessor.insert_sequences_batch')
def test_load_one_file_04(insert_sequences_batch, info):
    """Load single end sequences from a file into the atram database."""
    cxn, args = set_up()
    db.BATCH_SIZE = 5

    file_1 = join('tests', 'data', 'load_seq1.txt')
    core_preprocessor.load_one_file(
        args, cxn, file_1, 'single_ends', '')

    msg = 'Loading "{}" into sqlite database'.format(file_1)
    info.assert_called_once_with(msg)

    calls = [
        call(cxn, [
            ('seq1/1', '', 'AAAAAAAAAA'),
            ('seq2 1', '', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3 1', '', 'AAAAAAAAAA'),
            ('seq4_1', '', 'AAAAAAAAAA'),
            ('seq5/3', '', 'AAAAAAAAAAGGGGGGGGGG')]),
        call(cxn, [
            ('seq1/2', '', 'AAAAAAAAAA'),
            ('seq2 2', '', 'AAAAAAAAAAGGGGGGGGGG'),
            ('seq3 2', '', 'AAAAAAAAAA'),
            ('seq4_2', '', 'AAAAAAAAAAGGGGGGGGGG')])]
    insert_sequences_batch.assert_has_calls(calls)


@patch('lib.log.info')
@patch('lib.db_preprocessor.insert_sequences_batch')
def test_load_one_file_05(insert_sequences_batch, info):
    """Load single end sequences from a file into the atram database."""
    cxn, args = set_up()
    db.BATCH_SIZE = 5

    file_1 = join('tests', 'data', 'load_seq2.txt')
    core_preprocessor.load_one_file(args, cxn, file_1, 'mixed_ends')

    msg = 'Loading "{}" into sqlite database'.format(file_1)
    info.assert_called_once_with(msg)

    calls = [
        call(cxn, [
            ('seq6', '1', 'TTTTTTTTTT'),
            ('seq7', '1', 'TTTTTTTTTTCCCCCCCCCC'),
            ('seq8/a', '1', 'TTTTTTTTTT'),
            ('seq8', '2', 'TTTTTTTTTTCCCCCCCCCC')])]
    insert_sequences_batch.assert_has_calls(calls)


@patch('lib.log.info')
@patch('lib.db_preprocessor.get_sequence_count')
@patch('lib.db_preprocessor.get_shard_cut')
def test_assign_seqs_to_shards_01(
        get_shard_cut, get_sequence_count, info):
    """It assigns sequences to blast DB shards."""
    get_sequence_count.return_value = 100
    get_shard_cut.side_effect = ['seq1', 'seq2', 'seq3', 'seq4']

    shard_list = core_preprocessor.assign_seqs_to_shards(True, 3)

    assert [
        ('seq1', 'seq2'),
        ('seq2', 'seq3'),
        ('seq3', 'seq4z')] == shard_list

    msg = 'Assigning sequences to shards'
    info.assert_called_once_with(msg)


@patch('lib.blast.create_db')
@patch('lib.core_preprocessor.fill_blast_fasta')
def test_create_one_blast_shard_01(fill_blast_fasta, create_db):
    """It creates a blast DB from the shard."""
    _, args = set_up()
    shard_params = ['limit', 'offset']

    core_preprocessor.create_one_blast_shard(args, shard_params, 11)

    path = join(args['temp_dir'], 'pytest_011.fasta')

    fill_blast_fasta.assert_called_once_with(
        args['blast_db'], path, shard_params)

    shard = '{}.{:03d}.blast'.format(args['blast_db'], 11)
    create_db.assert_called_once_with(args['temp_dir'], path, shard)


@patch('lib.db.connect')
@patch('lib.db_preprocessor.get_sequences_in_shard')
def test_fill_blast_fasta_01(get_sequences_in_shard, connect):
    """It fills the fasta file used as blast input."""
    dbh = 'my_db'

    connect.return_value.__enter__ = lambda x: dbh
    get_sequences_in_shard.return_value = [
        ('seq1', '1', 'AAAAAAAAAA'),
        ('seq1', '2', 'CCCCCCCCCC'),
        ('seq2', '', 'GGGGGGGGGG'),
        ('seq3', '1', 'TTTTTTTTTT')]

    with tempfile.TemporaryDirectory(prefix='test_') as temp_dir:
        blast_db = 'test_blast_db'
        fasta_path = join(temp_dir, 'test_output.fasta')
        shard_params = (100, 200)  # limit and offset

        core_preprocessor.fill_blast_fasta(
            blast_db, fasta_path, shard_params)

        with open(fasta_path) as test_file:
            expect = (
                '>seq1/1\n'
                'AAAAAAAAAA\n'
                '>seq1/2\n'
                'CCCCCCCCCC\n'
                '>seq2\n'
                'GGGGGGGGGG\n'
                '>seq3/1\n'
                'TTTTTTTTTT\n')
            assert expect == test_file.read()

    calls = [
        call(blast_db),
        call().__exit__(None, None, None)]
    connect.assert_has_calls(calls)

    get_sequences_in_shard.assert_called_once_with(dbh, 100, 200)
