"""Test object in lib/assemblers/base."""

# pylint: disable=missing-docstring,too-many-arguments,no-self-use

import re
from os.path import basename, join
import subprocess
import pytest
import unittest
from unittest.mock import patch, MagicMock, call
from hypothesis import given, assume
import hypothesis.strategies as st
from lib.assemblers.base import BaseAssembler


PATH_PATTERN = r'(?:[\w.-]+/)[\w.-]+'


class TestAssemblersBase(unittest.TestCase):

    def build_assembler(self):
        args = {'bit_score': 44, 'contig_length': 55}
        state = {'blast_db': 'my_blast_db',
                 'query_file': 'my_query_file',
                 'query_target': 'my_query_name',
                 'iteration': 99}
        assembler = BaseAssembler(args, 'my_db_conn')
        assembler.set_state(
            state['blast_db'],
            state['query_file'],
            state['iteration'])
        return assembler

    @given(args=st.text(), db_conn=st.text())
    def test_init(self, args, db_conn):
        assembler = BaseAssembler(args, db_conn)

        assert args == assembler.args
        assert assembler.blast_only is False
        assert [] == assembler.steps
        assert {} == assembler.file

        expected = {
            'iteration': 0,
            'query_target': '',
            'query_file': '',
            'blast_db': '',
            'db_conn': db_conn}
        assert expected == assembler.state

    @given(
        blast_db=st.from_regex(PATH_PATTERN),
        query_file=st.from_regex(PATH_PATTERN),
        iteration=st.integers())
    def test_initialize_iteration(self, blast_db, query_file, iteration):
        output_files = ['prefix/dir/' + f for f
                        in ['output.fasta',
                            'paired_1.fasta',
                            'paired_2.fasta',
                            'single_1.fasta',
                            'single_2.fasta',
                            'single_any.fasta']]
        assembler = BaseAssembler('args', 'db_conn')
        assembler.iter_file = MagicMock(side_effect=output_files)

        assembler.initialize_iteration(blast_db, query_file, iteration)

        assert blast_db == assembler.state['blast_db']
        assert query_file == assembler.state['query_file']
        assert iteration == assembler.state['iteration']
        assert assembler.file['long_reads'] == ''
        assert output_files[0] == assembler.file['output']
        assert output_files[1] == assembler.file['paired_1']
        assert output_files[2] == assembler.file['paired_2']
        assert output_files[3] == assembler.file['single_1']
        assert output_files[4] == assembler.file['single_2']
        assert output_files[5] == assembler.file['single_any']
        assert assembler.file['paired_count'] == 0
        assert assembler.file['single_1_count'] == 0
        assert assembler.file['single_2_count'] == 0
        assert assembler.file['single_any_count'] == 0

    @given(
        args=st.fixed_dictionaries({'temp_dir': st.from_regex(PATH_PATTERN)}),
        blast_db=st.from_regex(PATH_PATTERN),
        query_file=st.from_regex(PATH_PATTERN),
        iteration=st.integers())
    def test_iter_dir(self, args, blast_db, query_file, iteration):
        assume(re.search(r'\w', query_file))
        assembler = BaseAssembler(args, 'db_conn')
        assembler.set_state(blast_db, query_file, iteration)

        base_blast_db = basename(blast_db)
        base_query_file = basename(query_file) if iteration == 1 else ''

        dir_name = '{}_{}_{:02d}'.format(
            base_blast_db, base_query_file, iteration)

        expect = join(args['temp_dir'], dir_name)
        assert expect == assembler.iter_dir()

    @given(
        args=st.fixed_dictionaries({'temp_dir': st.from_regex(PATH_PATTERN)}),
        blast_db=st.from_regex(PATH_PATTERN),
        query_file=st.from_regex(PATH_PATTERN),
        iteration=st.integers(),
        file_name=st.from_regex(r'\w[\w.-]+'))
    def test_iter_file(
            self, args, blast_db, query_file, iteration, file_name):
        assume(re.search(r'\w', query_file))
        assembler = BaseAssembler(args, 'db_conn')
        assembler.set_state(blast_db, query_file, iteration)

        base_blast_db = basename(blast_db)
        base_query_file = basename(query_file) if iteration == 1 else ''

        dir_name = '{}_{}_{:02d}'.format(
            base_blast_db, base_query_file, iteration)

        expect = join(args['temp_dir'], dir_name, file_name)
        assert expect == assembler.iter_file(file_name)

    def test_work_path(self):
        iter_dir = 'whatever'
        assembler = BaseAssembler('args', 'db_conn')
        assembler.iter_dir = MagicMock(return_value=iter_dir)

        assert iter_dir == assembler.work_path()

    @patch('lib.log.info')
    def test_run_ok(self, info):
        args = {'assembler': 'my_assembler', 'timeout': 10}

        assembler = BaseAssembler(args, 'db_conn')
        assembler.set_state('blast_db', 'query_file', 99)

        assembler.assemble = MagicMock()

        assembler.run()

        info.assert_called_once_with('Assembling shards with {}: '
                                     'iteration {}'.format(
                                         args['assembler'],
                                         assembler.state['iteration']))

    @patch('lib.log.info')
    @patch('lib.log.error')
    def test_run_timeout(self, error, info):
        args = {'assembler': 'my_assembler', 'timeout': 10}

        assembler = BaseAssembler(args, 'db_conn')
        assembler.set_state('blast_db', 'query_file', 99)

        assembler.assemble = MagicMock(side_effect=TimeoutError())

        with pytest.raises(TimeoutError) as timeout_error:
            assembler.run()

        error_msg = str(timeout_error.value)
        if error_msg[-1] != '.':
            error_msg += '.'
        assert error_msg == (
            'Time ran out for the assembler after 0:00:10 (HH:MM:SS).')

        expect = 'Assembling shards with {}: iteration {}'.format(
            args['assembler'], assembler.state['iteration'])
        info.assert_called_once_with(expect)

        # Python 3.6 formats exceptions differently so we need to do this
        assert error.call_count == 1
        regex = re.compile(
            r'Time ran out for the assembler after 0:00:10 \(HH:MM:SS\)')
        assert regex.match(error.call_args[0][0])

    @patch('lib.log.info')
    @patch('lib.log.error')
    def test_run_called_process_error(self, error, info):
        args = {'assembler': 'my_assembler', 'timeout': 10}
        error_code = 88
        cmd = 'my command'
        error = subprocess.CalledProcessError(error_code, cmd)

        assembler = BaseAssembler(args, 'db_conn')
        assembler.set_state('blast_db', 'query_file', 99)

        assembler.assemble = MagicMock(side_effect=error)

        with pytest.raises(RuntimeError) as runtime_error:
            assembler.run()

        error_msg = str(runtime_error.value)
        if error_msg[-1] != '.':
            error_msg += '.'
        assert error_msg == (
            "The assembler failed with error: Command 'my command' "
            "returned non-zero exit status 88.")

        expect = 'Assembling shards with {}: iteration {}'.format(
            args['assembler'], assembler.state['iteration'])
        info.assert_called_once_with(expect)

    @patch('lib.log.info')
    @patch('lib.db.sra_blast_hits_count')
    def test_no_blast_hits_true(self, sra_blast_hits_count, info):
        assembler = BaseAssembler('args', 'db_conn')
        assembler.set_state('blast_db', 'query_file', 99)

        sra_blast_hits_count.return_value = 0

        assert assembler.no_blast_hits()

        info.assert_called_once_with('No blast hits in iteration {}'.format(
            assembler.state['iteration']))

    @patch('lib.log.info')
    @patch('lib.db.sra_blast_hits_count')
    def test_no_blast_hits_false(self, sra_blast_hits_count, info):
        assembler = BaseAssembler('args', 'db_conn')
        assembler.set_state('blast_db', 'query_file', 99)

        sra_blast_hits_count.return_value = 1

        assert not assembler.no_blast_hits()

        info.assert_not_called()

    @patch('lib.log.info')
    def test_nothing_assembled_missing(self, info):
        assembler = BaseAssembler('args', 'db_conn')
        assembler.set_state('blast_db', 'query_file', 99)
        assembler.file['output'] = 'tests/data/missing_file.txt'

        assert assembler.nothing_assembled()

        info.assert_called_once_with(
            'No new assemblies in iteration {}'.format(
                assembler.state['iteration']))

    @patch('lib.log.info')
    def test_nothing_assembled_empty(self, info):
        assembler = BaseAssembler('args', 'db_conn')
        assembler.set_state('blast_db', 'query_file', 99)
        assembler.file['output'] = 'tests/data/empty_file.txt'

        assert assembler.nothing_assembled()

        expect = 'No new assemblies in iteration {}'.format(
            assembler.state['iteration'])
        info.assert_called_once_with(expect)

    @patch('lib.log.info')
    def test_nothing_assembled_false(self, info):
        assembler = BaseAssembler('args', 'db_conn')
        assembler.set_state('blast_db', 'query_file', 99)
        assembler.file['output'] = 'tests/data/load_seq1.txt'

        assert not assembler.nothing_assembled()

        info.assert_not_called()

    @patch('lib.log.info')
    @patch('lib.db.assembled_contigs_count')
    def test_assembled_contigs_count_0(self, assembled_contigs_count, info):
        high_score = 5
        assembler = self.build_assembler()

        assembled_contigs_count.return_value = 0

        assert assembler.assembled_contigs_count(high_score) == 0

        assembled_contigs_count.assert_called_once_with(
            assembler.state['db_conn'],
            assembler.state['iteration'],
            assembler.args['bit_score'],
            assembler.args['contig_length'])

        expect = ('No contigs had a bit score greater than {} and are at '
                  'least {} long in iteration {}. The highest score for '
                  'this iteration is {}').format(
                      assembler.args['bit_score'],
                      assembler.args['contig_length'],
                      assembler.state['iteration'],
                      high_score)
        info.assert_called_once_with(expect)

    @patch('lib.log.info')
    @patch('lib.db.assembled_contigs_count')
    def test_assembled_contigs_count_1(self, assembled_contigs_count, info):
        high_score = 5
        count = 1
        assembler = self.build_assembler()

        assembled_contigs_count.return_value = count

        assert assembler.assembled_contigs_count(high_score) == count

        assembled_contigs_count.assert_called_once_with(
            assembler.state['db_conn'],
            assembler.state['iteration'],
            assembler.args['bit_score'],
            assembler.args['contig_length'])

        info.assert_not_called()

    @patch('lib.log.info')
    @patch('lib.db.iteration_overlap_count')
    def test_no_new_contigs_ne(self, iteration_overlap_count, info):
        count = 1
        assembler = self.build_assembler()

        iteration_overlap_count.return_value = count + 1

        assert not assembler.no_new_contigs(count)

        iteration_overlap_count.assert_called_once_with(
            assembler.state['db_conn'],
            assembler.state['iteration'],
            assembler.args['bit_score'],
            assembler.args['contig_length'])

        info.assert_not_called()

    @patch('lib.log.info')
    @patch('lib.db.iteration_overlap_count')
    def test_no_new_contigs_eq(self, iteration_overlap_count, info):
        count = 1
        assembler = self.build_assembler()

        iteration_overlap_count.return_value = count

        assert assembler.no_new_contigs(count)

        iteration_overlap_count.assert_called_once_with(
            assembler.state['db_conn'],
            assembler.state['iteration'],
            assembler.args['bit_score'],
            assembler.args['contig_length'])

        expect = 'No new contigs were found in iteration {}'.format(
            assembler.state['iteration'])
        info.assert_called_once_with(expect)

    @patch('lib.log.subcommand')
    def test_assemble(self, subcommand):
        assembler = self.build_assembler()
        assembler.args['temp_dir'] = 'my_temp_dir'
        assembler.args['timeout'] = 333

        assembler.post_assembly = MagicMock()
        assembler.step1 = MagicMock(return_value='step1')
        assembler.step2 = MagicMock(return_value='step2')

        assembler.steps = [assembler.step1, assembler.step2]

        assembler.assemble()

        calls = [
            call('step1',
                 assembler.args['temp_dir'],
                 assembler.args['timeout']),
            call('step2',
                 assembler.args['temp_dir'],
                 assembler.args['timeout'])]
        subcommand.assert_has_calls(calls)
