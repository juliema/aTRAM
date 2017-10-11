"""Test object in lib/assemblers/base."""

# pylint: disable=missing-docstring,too-many-arguments

import re
from os.path import basename, join
import unittest
from unittest.mock import patch, call
from subprocess import CalledProcessError
from hypothesis import given, assume
import hypothesis.strategies as st
import lib.db as db
import lib.log as log
from lib.assemblers.base import BaseAssembler


PATH_PATTERN = r'(?:[\w.-]+/)[\w.-]+'

class TestAtramPreprocessor(unittest.TestCase):
    pass

    # def build_assembler():
    #     args = {'bit_score': 44, 'contig_length': 55}
    #     state = {'blast_db': 'my_blast_db',
    #              'query_file': 'my_query_file',
    #              'query_target': 'my_query_name',
    #              'iteration': 99}
    #     assembler = BaseAssembler(args, 'my_db_conn')
    #     assembler.set_state(
    #         state['blast_db'],
    #         state['query_file'],
    #         state['iteration'])
    #     return assembler
    #
    #
    # @given(args=st.text(), db_conn=st.text())
    # def test_init(args, db_conn):
    #     assembler = BaseAssembler(args, db_conn)
    #
    #     assert args == assembler.args
    #     assert assembler.blast_only is False
    #     assert [] == assembler.steps
    #     assert {} == assembler.file
    #
    #     expected = {
    #         'iteration': 0,
    #         'query_target': '',
    #         'query_file': '',
    #         'blast_db': '',
    #         'db_conn': db_conn}
    #     assert expected == assembler.state
    #
    #
    # @given(
    #     blast_db=st.from_regex(PATH_PATTERN),
    #     query_file=st.from_regex(PATH_PATTERN),
    #     iteration=st.integers())
    # def test_initialize_iteration(blast_db, query_file, iteration):
    #     output_files = ['prefix/dir/' + f for f
    #                     in ['output.fasta',
    #                         'paired_1.fasta',
    #                         'paired_2.fasta',
    #                         'single_1.fasta',
    #                         'single_2.fasta',
    #                         'single_any.fasta']]
    #     assembler = BaseAssembler('args', 'db_conn')
    #
    #     mock.it(assembler, 'iter_file', output_files)
    #
    #     assembler.initialize_iteration(blast_db, query_file, iteration)
    #
    #     assert blast_db == assembler.state['blast_db']
    #     assert query_file == assembler.state['query_file']
    #     assert iteration == assembler.state['iteration']
    #     assert assembler.file['long_reads'] == ''
    #     assert output_files[0] == assembler.file['output']
    #     assert output_files[1] == assembler.file['paired_1']
    #     assert output_files[2] == assembler.file['paired_2']
    #     assert output_files[3] == assembler.file['single_1']
    #     assert output_files[4] == assembler.file['single_2']
    #     assert output_files[5] == assembler.file['single_any']
    #     assert assembler.file['paired_count'] == 0
    #     assert assembler.file['single_1_count'] == 0
    #     assert assembler.file['single_2_count'] == 0
    #     assert assembler.file['single_any_count'] == 0
    #
    #
    # @given(
    #     args=st.fixed_dictionaries({'temp_dir': st.from_regex(PATH_PATTERN)}),
    #     blast_db=st.from_regex(PATH_PATTERN),
    #     query_file=st.from_regex(PATH_PATTERN),
    #     iteration=st.integers())
    # def test_iter_dir(args, blast_db, query_file, iteration):
    #     assume(re.search(r'\w', query_file))
    #     assembler = BaseAssembler(args, 'db_conn')
    #     assembler.set_state(blast_db, query_file, iteration)
    #
    #     base_blast_db = basename(blast_db)
    #     base_query_file = basename(query_file) if iteration == 1 else ''
    #
    #     dir_name = '{}_{}_{:02d}'.format(
    #         base_blast_db, base_query_file, iteration)
    #
    #     expect = join(args['temp_dir'], dir_name)
    #     assert expect == assembler.iter_dir()
    #
    #
    # @given(
    #     args=st.fixed_dictionaries({'temp_dir': st.from_regex(PATH_PATTERN)}),
    #     blast_db=st.from_regex(PATH_PATTERN),
    #     query_file=st.from_regex(PATH_PATTERN),
    #     iteration=st.integers(),
    #     file_name=st.from_regex(r'\w[\w.-]+'))
    # def test_iter_file(args, blast_db, query_file, iteration, file_name):
    #     assume(re.search(r'\w', query_file))
    #     assembler = BaseAssembler(args, 'db_conn')
    #     assembler.set_state(blast_db, query_file, iteration)
    #
    #     base_blast_db = basename(blast_db)
    #     base_query_file = basename(query_file) if iteration == 1 else ''
    #
    #     dir_name = '{}_{}_{:02d}'.format(
    #         base_blast_db, base_query_file, iteration)
    #
    #     expect = join(args['temp_dir'], dir_name, file_name)
    #     assert expect == assembler.iter_file(file_name)
    #
    #
    # def test_work_path():
    #     iter_dir = 'whatever'
    #     assembler = BaseAssembler('args', 'db_conn')
    #     mock.it(assembler, 'iter_dir', iter_dir)
    #
    #     assert iter_dir == assembler.work_path()
    #
    #
    # def test_run_ok():
    #     args = {'assembler': 'my_assembler', 'timeout': 10}
    #
    #     assembler = BaseAssembler(args, 'db_conn')
    #     assembler.set_state('blast_db', 'query_file', 99)
    #
    #     mock.it(assembler, 'assemble')
    #     mock.it(log, 'info')
    #
    #     assembler.run()
    #
    #     expect = [{'msg': 'Assembling shards with {}: iteration {}'.format(
    #         args['assembler'], assembler.state['iteration'])}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #
    # def test_run_timeout():
    #     args = {'assembler': 'my_assembler', 'timeout': 10}
    #
    #     assembler = BaseAssembler(args, 'db_conn')
    #     assembler.set_state('blast_db', 'query_file', 99)
    #
    #     mock.exception(assembler, 'assemble', TimeoutError())
    #     mock.it(log, 'info')
    #     mock.it(log, 'fatal')
    #
    #     assembler.run()
    #
    #     expect = [{'msg': 'Assembling shards with {}: iteration {}'.format(
    #         args['assembler'], assembler.state['iteration'])}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #     fatal = mock.filter('lib.log', 'fatal')
    #     assert len(fatal) == 1
    #
    #     # Python 3.6 formats exceptions differently
    #     expect = 'Time ran out for the assembler after 0:00:{} (HH:MM:SS)'.format(
    #         args['timeout'])
    #     actual = re.sub(r'\.$', '', fatal[0]['msg'])
    #     assert expect == actual
    #
    #
    # def test_run_called_process_error():
    #     args = {'assembler': 'my_assembler', 'timeout': 10}
    #     error_code = 88
    #     cmd = 'my command'
    #     error = CalledProcessError(error_code, cmd)
    #
    #     assembler = BaseAssembler(args, 'db_conn')
    #     assembler.set_state('blast_db', 'query_file', 99)
    #
    #     mock.exception(assembler, 'assemble', error)
    #     mock.it(log, 'info')
    #     mock.it(log, 'fatal')
    #
    #     assembler.run()
    #
    #     expect = [{'msg': 'Assembling shards with {}: iteration {}'.format(
    #         args['assembler'], assembler.state['iteration'])}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #     fatal = mock.filter('lib.log', 'fatal')
    #     assert len(fatal) == 1
    #
    #     expect = ("The assembler failed with error: Command '{}' "
    #               'returned non-zero exit status {}').format(cmd, error_code)
    #     actual = re.sub(r'\.$', '', fatal[0]['msg'])
    #
    #     assert expect == actual
    #
    #
    # def test_no_blast_hits_true():
    #     assembler = BaseAssembler('args', 'db_conn')
    #     assembler.set_state('blast_db', 'query_file', 99)
    #
    #     mock.it(db, 'sra_blast_hits_count', 0)
    #     mock.it(log, 'info')
    #
    #     assert assembler.no_blast_hits()
    #
    #     expect = [{'msg': 'No blast hits in iteration {}'.format(
    #         assembler.state['iteration'])}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #
    # def test_no_blast_hits_false():
    #     assembler = BaseAssembler('args', 'db_conn')
    #     assembler.set_state('blast_db', 'query_file', 99)
    #
    #     mock.it(db, 'sra_blast_hits_count', 1)
    #     mock.it(log, 'info')
    #
    #     assert not assembler.no_blast_hits()
    #
    #     assert [] == mock.filter('lib.log', 'info')
    #
    #
    # def test_nothing_assembled_missing():
    #     assembler = BaseAssembler('args', 'db_conn')
    #     assembler.set_state('blast_db', 'query_file', 99)
    #     assembler.file['output'] = 'tests/data/missing_file.txt'
    #
    #     mock.it(log, 'info')
    #
    #     assert assembler.nothing_assembled()
    #
    #     expect = [{'msg': 'No new assemblies in iteration {}'.format(
    #         assembler.state['iteration'])}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #
    # def test_nothing_assembled_empty():
    #     assembler = BaseAssembler('args', 'db_conn')
    #     assembler.set_state('blast_db', 'query_file', 99)
    #     assembler.file['output'] = 'tests/data/empty_file.txt'
    #
    #     mock.it(log, 'info')
    #
    #     assert assembler.nothing_assembled()
    #
    #     expect = [{'msg': 'No new assemblies in iteration {}'.format(
    #         assembler.state['iteration'])}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #
    # def test_nothing_assembled_false():
    #     assembler = BaseAssembler('args', 'db_conn')
    #     assembler.set_state('blast_db', 'query_file', 99)
    #     assembler.file['output'] = 'tests/data/load_seq1.txt'
    #
    #     mock.it(log, 'info')
    #
    #     assert not assembler.nothing_assembled()
    #
    #     assert [] == mock.filter('lib.log', 'info')
    #
    #
    # def test_assembled_contigs_count_0():
    #     high_score = 5
    #     assembler = build_assembler()
    #
    #     mock.it(db, 'assembled_contigs_count', 0)
    #     mock.it(log, 'info')
    #
    #     assert assembler.assembled_contigs_count(high_score) == 0
    #
    #     expect = [{
    #         'db_conn': assembler.state['db_conn'],
    #         'iteration': assembler.state['iteration'],
    #         'bit_score': assembler.args['bit_score'],
    #         'length': assembler.args['contig_length']}]
    #     assert expect == mock.filter('lib.db', 'assembled_contigs_count')
    #
    #     expect = [{'msg': 'No contigs had a bit score greater than {} and are at '
    #                       'least {} long in iteration {}. The highest score for '
    #                       'this iteration is {}'.format(
    #                           assembler.args['bit_score'],
    #                           assembler.args['contig_length'],
    #                           assembler.state['iteration'],
    #                           high_score)}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #
    # def test_assembled_contigs_count_1():
    #     high_score = 5
    #     count = 1
    #     assembler = build_assembler()
    #
    #     mock.it(db, 'assembled_contigs_count', count)
    #     mock.it(log, 'info')
    #
    #     assert assembler.assembled_contigs_count(high_score) == count
    #
    #     expect = [{
    #         'db_conn': assembler.state['db_conn'],
    #         'iteration': assembler.state['iteration'],
    #         'bit_score': assembler.args['bit_score'],
    #         'length': assembler.args['contig_length'],
    #     }]
    #     assert expect == mock.filter('lib.db', 'assembled_contigs_count')
    #
    #     assert [] == mock.filter('lib.log', 'info')
    #
    #
    # def test_no_new_contigs_ne():
    #     count = 1
    #     assembler = build_assembler()
    #
    #     mock.it(db, 'iteration_overlap_count', count + 1)
    #     mock.it(log, 'info')
    #
    #     assert not assembler.no_new_contigs(count)
    #
    #     expect = [{
    #         'db_conn': assembler.state['db_conn'],
    #         'iteration': assembler.state['iteration'],
    #         'bit_score': assembler.args['bit_score'],
    #         'length': assembler.args['contig_length'],
    #     }]
    #     assert expect == mock.filter('lib.db', 'iteration_overlap_count')
    #
    #     assert [] == mock.filter('lib.log', 'info')
    #
    #
    # def test_no_new_contigs_eq():
    #     count = 1
    #     assembler = build_assembler()
    #
    #     mock.it(db, 'iteration_overlap_count', count)
    #     mock.it(log, 'info')
    #
    #     assert assembler.no_new_contigs(count)
    #
    #     expect = [{
    #         'db_conn': assembler.state['db_conn'],
    #         'iteration': assembler.state['iteration'],
    #         'bit_score': assembler.args['bit_score'],
    #         'length': assembler.args['contig_length'],
    #     }]
    #     assert expect == mock.filter('lib.db', 'iteration_overlap_count')
    #
    #     expect = [{'msg': 'No new contigs were found in iteration {}'.format(
    #         assembler.state['iteration'])}]
    #     assert expect == mock.filter('lib.log', 'info')
    #
    #
    # def test_assemble():
    #     assembler = build_assembler()
    #     assembler.args['temp_dir'] = 'my_temp_dir'
    #     assembler.args['timeout'] = 333
    #
    #     mock.it(log, 'subcommand')
    #     mock.it(assembler, 'post_assembly')
    #
    #     def step1():
    #         mock.history.append({'module': 'none', 'func': 'step1'})
    #         return 'step1'
    #
    #     def step2():
    #         mock.history.append({'module': 'none', 'func': 'step2'})
    #         return 'step2'
    #     assembler.steps = [step1, step2]
    #
    #     assembler.assemble()
    #
    #     assert mock.filter('none', 'step1') == [{}]
    #     assert mock.filter('none', 'step2') == [{}]
    #
    #     expect = [{
    #         'cmd': 'step1',
    #         'temp_dir': assembler.args['temp_dir'],
    #         'timeout': assembler.args['timeout']
    #     }, {
    #         'cmd': 'step2',
    #         'temp_dir': assembler.args['temp_dir'],
    #         'timeout': assembler.args['timeout']
    #     }]
    #     assert expect == mock.filter('lib.log', 'subcommand')
