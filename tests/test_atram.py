"""Testing functions in atram."""

# pylint: disable=too-many-arguments
# pylint: disable=missing-docstring

import unittest
from unittest.mock import patch, MagicMock, call
from lib.assemblers.base import BaseAssembler
# import lib.db
import atram


class TestAtram(unittest.TestCase):
    """Testing functions in atram."""

    def setUp(self):
        self.db_conn = 'db_conn'
        self.args = {
            'query': ['query_file_1', 'query_file_2'],
            'blast_db': ['blast_db_1', 'blast_db_2'],
            'log_file': 'log_file_1',
            'temp_dir': 'temp_dir_1'}
        self.assembler = BaseAssembler(self.args, self.db_conn)

    @patch('lib.assembler.factory')
    @patch('lib.db.connect')
    @patch('lib.log.setup')
    @patch('atram.split_queries')
    @patch('atram.clean_database')
    @patch('atram.assembly_loop')
    def test_assemble(self, assembly_loop, clean_database, split_queries,
                      setup, connect, factory):
        connect.return_value.__enter__ = lambda x: 'my_db'
        self.assembler.write_final_output = MagicMock()
        split_queries.return_value = self.args['query']
        factory.return_value = self.assembler

        atram.assemble(self.args)

        split_queries.assert_called_once_with(self.args)

        calls = [
            call(self.args['blast_db'][0], check_version=True),
            call().__exit__(None, None, None),
            call(self.args['blast_db'][1], check_version=True),
            call().__exit__(None, None, None)]
        connect.assert_has_calls(calls)

        calls = [call('my_db'), call('my_db'), call('my_db'), call('my_db')]
        clean_database.assert_has_calls(calls)

        log_file = self.args['log_file']
        blast_db = self.args['blast_db']
        query = self.args['query']

        calls = [
            call(log_file, blast_db[0], query[0]),
            call(log_file, blast_db[0], query[1]),
            call(log_file, blast_db[1], query[0]),
            call(log_file, blast_db[1], query[1])]
        setup.assert_has_calls(calls)

        calls = [
            call(self.assembler, blast_db[0], query[0]),
            call(self.assembler, blast_db[0], query[1]),
            call(self.assembler, blast_db[1], query[0]),
            call(self.assembler, blast_db[1], query[1])]
        assembly_loop.assert_has_calls(calls)

        calls = [
            call(blast_db[0], query[0]),
            call(blast_db[0], query[1]),
            call(blast_db[1], query[0]),
            call(blast_db[1], query[1])]
        self.assembler.write_final_output.assert_has_calls(calls)

    @patch('atram.write_query_seq')
    def test_split_queries_none(self, write_query_seq):
        """Test split queries where there are no fasta files to split."""

        self.args['query_split'] = []

        queries = atram.split_queries(self.args)

        write_query_seq.assert_not_called()

        assert self.args['query'] == queries

    @patch('atram.write_query_seq')
    def test_split_queries_some(self, write_query_seq):
        """Test split queries where there are fasta files to split."""

        self.args['query_split'] = ['tests/data/split_queries1.txt']

        queries = atram.split_queries(self.args)

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

        expected = self.args['query'] + split_files
        assert expected == queries
