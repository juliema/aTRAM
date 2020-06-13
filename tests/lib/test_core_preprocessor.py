"""Testing functions in core_preprocessor."""

# pylint: disable=too-many-arguments,too-many-locals,unused-argument

from os.path import join
from unittest.mock import call, patch

import lib.core_preprocessor as core_preprocessor
import lib.db as db


def set_up():
    """Build structures for the tests."""
    args = {
        'blast_db': 'blast_db_1',
        'log_file': 'log_file_1',
        'log_level': 'info',
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
@patch('lib.core_preprocessor.create_all_blast_shards')
def test_preprocess_01(
        create_all_blast_shards, load_seqs,
        create_sequences_index, create_sequences_table,
        create_metadata_table, connect, setup, info, update_temp_dir,
        make_temp_dir):
    """It runs the function to build the databases required by atram."""
    _, args = set_up()
    shard_list = ['shard_1', 'shard_4', 'shard_3', 'shard_4']
    dbh = 'my_db'
    connect.return_value.__enter__ = lambda x: dbh

    core_preprocessor.preprocess(args)

    setup.assert_called_once_with(
        args['log_file'], args['log_level'], args['blast_db'])

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
    create_all_blast_shards.assert_called_once_with(args, dbh, len(shard_list))


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


@patch('lib.blast.create_db')
def test_create_one_blast_shard_01(create_db):
    """It creates a blast DB from the shard."""
    _, args = set_up()
    path = join(args['temp_dir'], 'pytest_011.fasta')

    core_preprocessor.create_one_blast_shard(args, path, 11)

    shard = '{}.{:03d}.blast'.format(args['blast_db'], 12)
    create_db.assert_called_once_with(args['temp_dir'], path, shard)
