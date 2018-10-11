"""Testing functions in atram."""

from os.path import join
import tempfile
import unittest
from unittest.mock import patch, MagicMock, call
from lib.assemblers.base import BaseAssembler
import atram


class TestAtram(unittest.TestCase):
    """Testing functions in atram."""

    def setUp(self):
        self.cxn = 'cxn'
        self.args = {
            'query': ['query_file_1', 'query_file_2'],
            'blast_db': ['blast_db_1', 'blast_db_2'],
            'iterations': 1,
            'log_file': 'log_file_1',
            'temp_dir': 'temp_dir_1'}
        self.assembler = BaseAssembler(self.args, self.cxn)

    @patch('lib.assembler.factory')
    @patch('lib.db.connect')
    @patch('lib.db.aux_db')
    @patch('lib.db.aux_detach')
    @patch('lib.log.setup')
    @patch('atram.split_queries')
    @patch('atram.clean_database')
    @patch('atram.assembly_loop')
    def test_assemble(self, assembly_loop, clean_database, split_queries,
                      setup, aux_detach, aux_db, connect, factory):
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
            call('my_db', self.args['temp_dir'], blast_db[0], query[0]),
            call('my_db', self.args['temp_dir'], blast_db[0], query[1]),
            call('my_db', self.args['temp_dir'], blast_db[1], query[0]),
            call('my_db', self.args['temp_dir'], blast_db[1], query[1])]
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

        assert split_files == queries

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

    @patch('lib.log.info')
    @patch('os.makedirs')
    @patch('atram.blast_query_against_all_shards')
    @patch('atram.create_query_from_contigs')
    @patch('atram.filter_contigs')
    def test_assembly_loop_one_iter(
            self,
            filter_contigs,
            create_query_from_contigs,
            blast_query_against_all_shards,
            b_makedirs,
            info):
        iter_dir = 'my_iter_dir'

        self.assembler.blast_only = False
        self.assembler.state['query_file'] = self.args['query'][0]
        self.assembler.state['blast_db'] = self.args['blast_db'][0]
        self.assembler.run = MagicMock()
        self.assembler.write_input_files = MagicMock()
        self.assembler.initialize_iteration = MagicMock()
        self.assembler.iter_dir = MagicMock(return_value=iter_dir)
        self.assembler.no_blast_hits = MagicMock(return_value=False)
        self.assembler.no_new_contigs = MagicMock(return_value=False)
        self.assembler.nothing_assembled = MagicMock(return_value=False)
        self.assembler.assembled_contigs_count = MagicMock(return_value=11)

        atram.assembly_loop(
            self.assembler, self.args['blast_db'][0], self.args['query'][0])

        calls = [
            call('aTRAM blast DB = "{}", query = "{}", iteration {}'.format(
                self.args['blast_db'][0], self.args['query'][0], 1)),
            call('All iterations completed')]
        info.assert_has_calls(calls)

        self.assembler.initialize_iteration.assert_called_once_with(
            self.args['blast_db'][0], self.args['query'][0], 1)

        assert self.assembler.iter_dir.call_count == 1

        b_makedirs.assert_called_once_with(iter_dir, exist_ok=True)

        blast_query_against_all_shards.assert_called_once_with(self.assembler)

        assert self.assembler.no_blast_hits.call_count == 1
        self.assembler.no_blast_hits.write_input_files()
        self.assembler.no_blast_hits.run()
        self.assembler.no_blast_hits.nothing_assembled()

        self.assembler.no_new_contigs.assert_called_once_with(11)

        create_query_from_contigs.create_query_from_contigs(self.assembler)
        filter_contigs.create_query_from_contigs(self.assembler)

    # def test_assembly_loop_multiple_iterations(self):

    # def test_assembly_loop_blast_only(self):

    # def test_assembly_loop_no_blast_hits(self):

    # def test_assembly_loop_nothing_assembled(self):

    # def test_assembly_loop_no_assembled_contigs_count(self):

    # def test_assembly_loop_no_new_contigs(self):

    @patch('lib.blast.all_shard_paths')
    def test_shard_fraction_one(self, all_shard_paths):
        returns = ['1st', '2nd', '3rd', '4th']
        self.assembler.state['blast_db'] = self.args['blast_db'][0]
        self.assembler.args['fraction'] = 1.0

        all_shard_paths.return_value = returns

        shards = atram.shard_fraction(self.assembler)

        assert returns == shards

        all_shard_paths.assert_called_once_with(self.args['blast_db'][0])

    @patch('lib.blast.all_shard_paths')
    def test_shard_fraction_half(self, all_shard_paths):
        self.assembler.args['fraction'] = 0.5
        self.assembler.state['blast_db'] = self.args['blast_db'][0]
        returns = ['1st', '2nd', '3rd', '4th']

        all_shard_paths.return_value = returns

        shards = atram.shard_fraction(self.assembler)

        assert ['1st', '2nd'] == shards

        all_shard_paths.assert_called_once_with(self.args['blast_db'][0])

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
