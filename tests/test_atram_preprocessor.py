"""Testing functions in atram_preprocessor."""

from os.path import join
import tempfile
import unittest
from unittest.mock import patch, call
import lib.db as db

import atram_preprocessor


class TestAtramPreprocessor(unittest.TestCase):

    def setUp(self):
        self.db_conn = 'db_conn'
        self.args = {
            'blast_db': 'blast_db_1',
            'log_file': 'log_file_1',
            'sra_files': 'sra_files_1',
            'shard_count': 4,
            'temp_dir': 'temp_dir_1'}

    @patch('lib.log.info')
    @patch('lib.log.setup')
    @patch('lib.db.connect')
    @patch('lib.db.create_metadata_table')
    @patch('lib.db.create_sequences_table')
    @patch('lib.db.create_sequences_index')
    @patch('atram_preprocessor.load_seqs')
    @patch('atram_preprocessor.assign_seqs_to_shards')
    @patch('atram_preprocessor.create_all_blast_shards')
    def test_preprocess(
            self, create_all_blast_shards, assign_seqs_to_shards, load_seqs,
            create_sequences_index, create_sequences_table,
            create_metadata_table, connect, setup, info):
        shard_list = ['shard_1', 'shard_4', 'shard_3', 'shard_4']
        dbh = 'my_db'
        assign_seqs_to_shards.return_value = shard_list
        connect.return_value.__enter__ = lambda x: dbh

        atram_preprocessor.preprocess(self.args)

        setup.assert_called_once_with(
            self.args['log_file'], self.args['blast_db'])

        calls = [
            call(self.args['blast_db'], clean=True),
            call().__exit__(None, None, None)]
        connect.assert_has_calls(calls)

        create_metadata_table.assert_called_once_with(dbh)
        create_sequences_table.assert_called_once_with(dbh)
        load_seqs.assert_called_once_with(self.args, dbh)
        info.assert_called_once_with(
            'Creating an index for the sequence table')
        create_sequences_index.assert_called_once_with(dbh)
        assign_seqs_to_shards.assert_called_once_with(
            dbh, self.args['shard_count'])
        create_all_blast_shards.assert_called_once_with(self.args, shard_list)

    @patch('lib.log.info')
    @patch('lib.db.insert_sequences_batch')
    def test_load_one_file_mixed(self, insert_sequences_batch, info):
        db.BATCH_SIZE = 5

        file_1 = join('tests', 'data', 'load_seq1.txt')
        atram_preprocessor.load_one_file(self.db_conn, file_1, 'mixed_ends')

        msg = 'Loading "{}" into sqlite database'.format(file_1)
        info.assert_called_once_with(msg)

        calls = [
            call(self.db_conn, [
                ('seq1', '1', 'AAAAAAAAAA'),
                ('seq2', '1', 'AAAAAAAAAAGGGGGGGGGG'),
                ('seq3', '1', 'AAAAAAAAAA'),
                ('seq4', '1', 'AAAAAAAAAA'),
                ('seq5/3', '', 'AAAAAAAAAAGGGGGGGGGG')]),
            call(self.db_conn, [
                ('seq1', '2', 'AAAAAAAAAA'),
                ('seq2', '2', 'AAAAAAAAAAGGGGGGGGGG'),
                ('seq3', '2', 'AAAAAAAAAA'),
                ('seq4', '2', 'AAAAAAAAAAGGGGGGGGGG')])]
        insert_sequences_batch.assert_has_calls(calls)

    @patch('lib.log.info')
    @patch('lib.db.insert_sequences_batch')
    def test_load_one_file_end1(self, insert_sequences_batch, info):
        db.BATCH_SIZE = 5

        file_1 = join('tests', 'data', 'load_seq1.txt')
        atram_preprocessor.load_one_file(self.db_conn, file_1, 'end_1', '1')

        msg = 'Loading "{}" into sqlite database'.format(file_1)
        info.assert_called_once_with(msg)

        calls = [
            call(self.db_conn, [
                ('seq1', '1', 'AAAAAAAAAA'),
                ('seq2', '1', 'AAAAAAAAAAGGGGGGGGGG'),
                ('seq3', '1', 'AAAAAAAAAA'),
                ('seq4', '1', 'AAAAAAAAAA'),
                ('seq5/3', '1', 'AAAAAAAAAAGGGGGGGGGG')]),
            call(self.db_conn, [
                ('seq1', '1', 'AAAAAAAAAA'),
                ('seq2', '1', 'AAAAAAAAAAGGGGGGGGGG'),
                ('seq3', '1', 'AAAAAAAAAA'),
                ('seq4', '1', 'AAAAAAAAAAGGGGGGGGGG')])]
        insert_sequences_batch.assert_has_calls(calls)

    @patch('lib.log.info')
    @patch('lib.db.insert_sequences_batch')
    def test_load_one_file_end2(self, insert_sequences_batch, info):
        db.BATCH_SIZE = 5

        file_1 = join('tests', 'data', 'load_seq1.txt')
        atram_preprocessor.load_one_file(self.db_conn, file_1, 'end_2', '2')

        msg = 'Loading "{}" into sqlite database'.format(file_1)
        info.assert_called_once_with(msg)

        calls = [
            call(self.db_conn, [
                ('seq1', '2', 'AAAAAAAAAA'),
                ('seq2', '2', 'AAAAAAAAAAGGGGGGGGGG'),
                ('seq3', '2', 'AAAAAAAAAA'),
                ('seq4', '2', 'AAAAAAAAAA'),
                ('seq5/3', '2', 'AAAAAAAAAAGGGGGGGGGG')]),
            call(self.db_conn, [
                ('seq1', '2', 'AAAAAAAAAA'),
                ('seq2', '2', 'AAAAAAAAAAGGGGGGGGGG'),
                ('seq3', '2', 'AAAAAAAAAA'),
                ('seq4', '2', 'AAAAAAAAAAGGGGGGGGGG')])]
        insert_sequences_batch.assert_has_calls(calls)

    @patch('lib.log.info')
    @patch('lib.db.insert_sequences_batch')
    def test_load_one_file_single(self, insert_sequences_batch, info):
        db.BATCH_SIZE = 5

        file_1 = join('tests', 'data', 'load_seq1.txt')
        atram_preprocessor.load_one_file(
            self.db_conn, file_1, 'single_ends', '')

        msg = 'Loading "{}" into sqlite database'.format(file_1)
        info.assert_called_once_with(msg)

        calls = [
            call(self.db_conn, [
                ('seq1', '', 'AAAAAAAAAA'),
                ('seq2', '', 'AAAAAAAAAAGGGGGGGGGG'),
                ('seq3', '', 'AAAAAAAAAA'),
                ('seq4', '', 'AAAAAAAAAA'),
                ('seq5/3', '', 'AAAAAAAAAAGGGGGGGGGG')]),
            call(self.db_conn, [
                ('seq1', '', 'AAAAAAAAAA'),
                ('seq2', '', 'AAAAAAAAAAGGGGGGGGGG'),
                ('seq3', '', 'AAAAAAAAAA'),
                ('seq4', '', 'AAAAAAAAAAGGGGGGGGGG')])]
        insert_sequences_batch.assert_has_calls(calls)

    @patch('lib.log.info')
    @patch('lib.db.insert_sequences_batch')
    def test_load_one_file_2(self, insert_sequences_batch, info):
        db.BATCH_SIZE = 5

        file_1 = join('tests', 'data', 'load_seq2.txt')
        atram_preprocessor.load_one_file(self.db_conn, file_1, 'mixed_ends')

        msg = 'Loading "{}" into sqlite database'.format(file_1)
        info.assert_called_once_with(msg)

        calls = [
            call(self.db_conn, [
                ('seq6', '1', 'TTTTTTTTTT'),
                ('seq7', '1', 'TTTTTTTTTTCCCCCCCCCC'),
                ('seq8/a.1 suffix', '', 'TTTTTTTTTT'),
                ('seq8', '2', 'TTTTTTTTTTCCCCCCCCCC')])]
        insert_sequences_batch.assert_has_calls(calls)

    @patch('lib.log.info')
    @patch('lib.db.get_sequence_count')
    @patch('lib.db.get_shard_cut')
    def test_assign_seqs_to_shards(
            self, get_shard_cut, get_sequence_count, info):

        get_sequence_count.return_value = 100
        get_shard_cut.side_effect = ['seq1', 'seq2', 'seq3', 'seq4']

        shard_list = atram_preprocessor.assign_seqs_to_shards(True, 3)

        assert [
            ('seq1', 'seq2'),
            ('seq2', 'seq3'),
            ('seq3', 'seq4z')] == shard_list

        msg = 'Assigning sequences to shards'
        info.assert_called_once_with(msg)

    @patch('lib.blast.create_db')
    @patch('atram_preprocessor.fill_blast_fasta')
    def test_create_one_blast_shard(self, fill_blast_fasta, create_db):
        shard_params = ['limit', 'offset']

        atram_preprocessor.create_one_blast_shard(self.args, shard_params, 11)

        path = join(self.args['temp_dir'], 'pyt_011.fasta')

        fill_blast_fasta.assert_called_once_with(
            self.args['blast_db'], path, shard_params)

        shard = '{}.{:03d}.blast'.format(self.args['blast_db'], 11)
        create_db.assert_called_once_with(self.args['temp_dir'], path, shard)

    @patch('lib.db.connect')
    @patch('lib.db.get_sequences_in_shard')
    def test_fill_blast_fasta(self, get_sequences_in_shard, connect):
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

            atram_preprocessor.fill_blast_fasta(
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
