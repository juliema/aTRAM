"""Testing functions in atram."""

# pylint: disable=missing-docstring,too-many-arguments,no-self-use

from os.path import join
import tempfile
import unittest
from unittest.mock import patch, MagicMock, call
from lib.assemblers.base import BaseAssembler
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

    def test_write_query_seq(self):
        with tempfile.TemporaryDirectory(prefix='test_') as temp_dir:
            path = join(temp_dir, 'test_query.fasta')

            atram.write_query_seq(path, 'my sequence name', 'aaaacccgggtt')

            with open(path) as test_file:
                expect = (
                    '>my sequence name\n'
                    'aaaacccgggtt\n')
                assert expect == test_file.read()

    @patch('lib.db.create_sra_blast_hits_table')
    @patch('lib.db.create_contig_blast_hits_table')
    @patch('lib.db.create_assembled_contigs_table')
    def test_clean_database(
            self,
            create_assembled_contigs_table,
            create_contig_blast_hits_table,
            create_sra_blast_hits_table):

        dbh = 'my_db'
        atram.clean_database(dbh)

        create_assembled_contigs_table.assert_called_once_with(dbh)
        create_contig_blast_hits_table.assert_called_once_with(dbh)
        create_sra_blast_hits_table.assert_called_once_with(dbh)

    # def test_assembly_loop_one_iter(self):
    #     db_conn = st.text()
    #     blast_db = st.text()
    #     query = st.text()
    #     iter_dir = st.text()
    #     args = {'temp_dir': st.text(), 'iterations': 1}
    #
    #     assembler = BaseAssembler(args, db_conn)
    #     assembler.blast_only = False
    #     assembler.state['query_file'] = query
    #     assembler.state['blast_db'] = blast_db
    #
    #     mock.it(log, 'info')
    #     mock.it(assembler, 'initialize_iteration')
    #     mock.it(assembler, 'iter_dir', iter_dir)
    #     mock.it(os, 'makedirs')
    #     mock.it(atram, 'blast_query_against_all_shards')
    #     mock.it(assembler, 'no_blast_hits', False)
    #     mock.it(assembler, 'write_input_files')
    #     mock.it(assembler, 'run')
    #     mock.it(assembler, 'nothing_assembled', False)
    #     mock.it(atram, 'filter_contigs', 99)
    #     mock.it(assembler, 'assembled_contigs_count', 11)
    #     mock.it(assembler, 'no_new_contigs', False)
    #     mock.it(atram, 'create_query_from_contigs', 'my_query')
    #
    #     atram.assembly_loop(assembler, blast_db, query)
    #
    #     expect = [{'msg': 'aTRAM blast DB = "{}", query = "{}", '
    #                       'iteration {}'.format(blast_db, query, 1)},
    #               {'msg': 'All iterations completed'}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #     expect = [
    #         {'blast_db': blast_db, 'iteration': 1, 'query_file': query}]
    #     assert expect == mock.filter('BaseAssembler', 'initialize_iteration')
    #
    #     assert [{}] == mock.filter('BaseAssembler', 'iter_dir')
    #
    #     expect = [{'name': iter_dir, 'exist_ok': True}]
    #     assert expect == mock.filter('os', 'makedirs')
    #
    #     expect = [{'assembler': assembler}]
    #     assert expect == mock.filter('atram', 'blast_query_against_all_shards')
    #
    #     assert [{}] == mock.filter('BaseAssembler', 'no_blast_hits')
    #     assert [{}] == mock.filter('BaseAssembler', 'write_input_files')
    #     assert [{}] == mock.filter('BaseAssembler', 'run')
    #     assert [{}] == mock.filter('BaseAssembler', 'nothing_assembled')
    #     assert [{'count': 11}] == mock.filter('BaseAssembler', 'no_new_contigs')
    #
    #     expect = [{'assembler': assembler}]
    #     assert expect == mock.filter('atram', 'create_query_from_contigs')
    #
    # # def test_assembly_loop_multiple_iterations(self):
    #
    # # def test_assembly_loop_blast_only(self):
    #
    # # def test_assembly_loop_no_blast_hits(self):
    #
    # # def test_assembly_loop_nothing_assembled(self):
    #
    # # def test_assembly_loop_no_assembled_contigs_count(self):
    #
    # # def test_assembly_loop_no_new_contigs(self):
    #
    # def test_shard_fraction_one():
    #     args = {'fraction': 1.0}
    #     returns = ['1st', '2nd', '3rd', '4th']
    #     assembler = BaseAssembler(args, 'db_conn')
    #     assembler.state['blast_db'] = 'my_blast_db'
    #
    #     mock.it(blast, 'all_shard_paths', returns=[returns])
    #
    #     shards = atram.shard_fraction(assembler)
    #
    #     assert returns == shards
    #
    #     expect = [{'blast_db': 'my_blast_db'}]
    #     assert expect == mock.filter('lib.blast', 'all_shard_paths')
    #
    # def test_shard_fraction_half(self):
    #     args = {'fraction': 0.5}
    #     assembler = BaseAssembler(args, 'db_conn')
    #     assembler.state['blast_db'] = 'my_blast_db'
    #
    #     mock.it(blast, 'all_shard_paths', returns=[['1st', '2nd', '3rd', '4th']])
    #
    #     shards = atram.shard_fraction(assembler)
    #
    #     assert ['1st', '2nd'] == shards
    #
    #     expect = [{'blast_db': 'my_blast_db'}]
    #     assert expect == mock.filter('lib.blast', 'all_shard_paths')
    #
    # # def test_blast_query_against_one_shard(self):
    #
    # # def test_filter_contigs(self):
    #
    # # def test_save_blast_against_contigs(self):
    #
    # # def test_save_contigs(self):
    #
    #
    # # def test_create_query_from_contigs(self):
