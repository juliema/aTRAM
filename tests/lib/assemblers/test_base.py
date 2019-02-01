"""Test object in lib/assemblers/base."""

import re
import subprocess
from unittest.mock import MagicMock, patch, call
import pytest
from lib.assemblers.base import BaseAssembler


def build_assembler():
    """Build a generic assembler."""
    args = {'bit_score': 44, 'contig_length': 55}
    state = {'blast_db': 'my_blast_db',
             'query_file': 'my_query_file.seq',
             'query_target': 'my_query_name',
             'iteration': 99}
    assembler = BaseAssembler(args, 'my_db_conn')
    assembler.set_state(
        state['blast_db'],
        state['query_file'],
        state['iteration'])
    return assembler

def test_init_01():
    """It initializes an assembler object."""
    assembler = BaseAssembler('args', 'cxn')

    assert assembler.args == 'args'
    assert assembler.blast_only is False
    assert assembler.steps == []
    assert assembler.file == {}
    expected = {
        'iteration': 0,
        'query_target': '',
        'query_file': '',
        'blast_db': '',
        'cxn': 'cxn'}
    assert assembler.state == expected

def test_initialize_iteration_01():
    """It sets up assembler variables for an iteration."""
    blast_db = 'my_blast_db'
    query_file = 'my_query_file.seq'
    iteration = 99
    temp_dir = 'my/temp/dir/'
    args = {'temp_dir': temp_dir}

    assembler = BaseAssembler(args, 'cxn')
    assembler.set_state(blast_db, query_file, iteration)

    assembler.initialize_iteration(blast_db, query_file, iteration)
    expected_dir = '{}{}_{}_{:02d}/'.format(
        temp_dir, blast_db, assembler.state['query_target'], iteration)

    assert assembler.state['blast_db'] == blast_db
    assert assembler.state['query_file'] == query_file
    assert assembler.state['iteration'] == 99
    assert assembler.file['long_reads'] == ''
    assert assembler.file['output'] == expected_dir + 'output.fasta'
    assert assembler.file['paired_1'] == expected_dir + 'paired_1.fasta'
    assert assembler.file['paired_2'] == expected_dir + 'paired_2.fasta'
    assert assembler.file['single_1'] == expected_dir + 'single_1.fasta'
    assert assembler.file['single_2'] == expected_dir + 'single_2.fasta'
    assert assembler.file['single_any'] == expected_dir + 'single_any.fasta'
    assert assembler.file['paired_count'] == 0
    assert assembler.file['single_1_count'] == 0
    assert assembler.file['single_2_count'] == 0
    assert assembler.file['single_any_count'] == 0

def test_iter_dir_01():
    """It builds the directory for the iteration."""
    blast_db = 'my_blast_db'
    query_file = 'my_query_file.seq'
    iteration = 99
    temp_dir = 'my/temp/dir/'
    args = {'temp_dir': temp_dir}

    assembler = BaseAssembler(args, 'cxn')
    assembler.set_state(blast_db, query_file, iteration)

    expected_dir = '{}{}_{}_{:02d}'.format(
        temp_dir, blast_db, assembler.state['query_target'], iteration)

    assert assembler.iter_dir() == expected_dir

def test_iter_file_01():
    """It builds a file name with the iter dir."""
    blast_db = 'my_blast_db'
    query_file = 'my_query_file.seq'
    iteration = 99
    temp_dir = 'my/temp/dir/'
    args = {'temp_dir': temp_dir}
    file_name = 'my_file.txt'

    assembler = BaseAssembler(args, 'cxn')
    assembler.set_state(blast_db, query_file, iteration)

    expected_dir = '{}{}_{}_{:02d}/'.format(
        temp_dir, blast_db, assembler.state['query_target'], iteration)

    assert assembler.iter_file(file_name) == expected_dir + file_name

def test_work_path():
    """It builds the work directory."""
    blast_db = 'my_blast_db'
    query_file = 'my_query_file.seq'
    iteration = 99
    temp_dir = 'my/temp/dir/'
    args = {'temp_dir': temp_dir}

    assembler = BaseAssembler(args, 'cxn')
    assembler.set_state(blast_db, query_file, iteration)

    expected_dir = '{}{}_{}_{:02d}'.format(
        temp_dir, blast_db, assembler.state['query_target'], iteration)

    assert assembler.work_path() == expected_dir

@patch('lib.log.info')
def test_run_01(info):
    """It performs a normal flow."""
    args = {'assembler': 'my_assembler', 'timeout': 10}

    assembler = BaseAssembler(args, 'cxn')
    assembler.set_state('blast_db', 'query_file', 99)

    assembler.assemble = MagicMock()

    assembler.run()

    info.assert_called_once_with('Assembling shards with {}: '
                                 'iteration {}'.format(
                                     args['assembler'],
                                     assembler.state['iteration']))

@patch('lib.log.info')
@patch('lib.log.error')
def test_run_02(error, info):
    """It handles a timeout error."""
    args = {'assembler': 'my_assembler', 'timeout': 10}

    assembler = BaseAssembler(args, 'cxn')
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
def test_run_03(error, info):
    """It handles a subprocess error."""
    args = {'assembler': 'my_assembler', 'timeout': 10}
    error_code = 88
    cmd = 'my command'
    error = subprocess.CalledProcessError(error_code, cmd)

    assembler = BaseAssembler(args, 'cxn')
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
def test_count_blast_hits_01(sra_blast_hits_count, info):
    """It handles no blast hits."""
    assembler = BaseAssembler('args', 'cxn')
    assembler.set_state('blast_db', 'query_file', 99)

    sra_blast_hits_count.return_value = 0

    assert assembler.count_blast_hits() == 0

    info.assert_called_once_with('0 blast hits in iteration {}'.format(
        assembler.state['iteration']))

@patch('lib.log.info')
@patch('lib.db.sra_blast_hits_count')
def test_count_blast_hits_02(sra_blast_hits_count, info):
    """It handles one blast hits."""
    assembler = BaseAssembler('args', 'cxn')
    assembler.set_state('blast_db', 'query_file', 99)

    sra_blast_hits_count.return_value = 1

    assert assembler.count_blast_hits() == 1

    info.assert_called_once_with('1 blast hits in iteration {}'.format(
        assembler.state['iteration']))

@patch('lib.log.info')
def test_nothing_assembled_01(info):
    """It handles when nothing is assembled."""
    assembler = BaseAssembler('args', 'cxn')
    assembler.set_state('blast_db', 'query_file', 99)
    assembler.file['output'] = 'tests/data/missing_file.txt'

    assert assembler.nothing_assembled()

    info.assert_called_once_with(
        'No new assemblies in iteration {}'.format(
            assembler.state['iteration']))

@patch('lib.log.info')
def test_nothing_assembled_02(info):
    """It handles an empty assembly."""
    assembler = BaseAssembler('args', 'cxn')
    assembler.set_state('blast_db', 'query_file', 99)
    assembler.file['output'] = 'tests/data/empty_file.txt'

    assert assembler.nothing_assembled()

    expect = 'No new assemblies in iteration {}'.format(
        assembler.state['iteration'])
    info.assert_called_once_with(expect)

@patch('lib.log.info')
def test_nothing_assembled_03(info):
    """It handles when something is assembled."""
    assembler = BaseAssembler('args', 'cxn')
    assembler.set_state('blast_db', 'query_file', 99)
    assembler.file['output'] = 'tests/data/load_seq1.txt'

    assert not assembler.nothing_assembled()

    info.assert_not_called()

@patch('lib.log.info')
@patch('lib.db.assembled_contigs_count')
def test_assembled_contigs_count_01(assembled_contigs_count, info):
    """Handle when here are no contigs."""
    high_score = 5
    assembler = build_assembler()

    assembled_contigs_count.return_value = 0

    assert assembler.assembled_contigs_count(high_score) == 0

    assembled_contigs_count.assert_called_once_with(
        assembler.state['cxn'],
        assembler.state['iteration'],
        assembler.args['bit_score'],
        assembler.args['contig_length'])

    expect = ('No contigs had a bit score greater than {} and are at '
              'least {} bp long in iteration {}. The highest score for '
              'this iteration is {}').format(
                  assembler.args['bit_score'],
                  assembler.args['contig_length'],
                  assembler.state['iteration'],
                  high_score)
    info.assert_called_once_with(expect)

@patch('lib.log.info')
@patch('lib.db.assembled_contigs_count')
def test_assembled_contigs_count_02(assembled_contigs_count, info):
    """Handle when here is one contig."""
    high_score = 5
    count = 1
    assembler = build_assembler()

    assembled_contigs_count.return_value = count

    assert assembler.assembled_contigs_count(high_score) == count

    assembled_contigs_count.assert_called_once_with(
        assembler.state['cxn'],
        assembler.state['iteration'],
        assembler.args['bit_score'],
        assembler.args['contig_length'])

    info.assert_not_called()

@patch('lib.log.info')
@patch('lib.db.iteration_overlap_count')
def test_no_new_contigs_01(iteration_overlap_count, info):
    """It handles when there are new contigs."""
    count = 1
    assembler = build_assembler()

    iteration_overlap_count.return_value = count + 1

    assert not assembler.no_new_contigs(count)

    iteration_overlap_count.assert_called_once_with(
        assembler.state['cxn'],
        assembler.state['iteration'],
        assembler.args['bit_score'],
        assembler.args['contig_length'])

    info.assert_not_called()

@patch('lib.log.info')
@patch('lib.db.iteration_overlap_count')
def test_no_new_contigs_02(iteration_overlap_count, info):
    """It handles when there are no new contigs."""
    count = 1
    assembler = build_assembler()

    iteration_overlap_count.return_value = count

    assert assembler.no_new_contigs(count)

    iteration_overlap_count.assert_called_once_with(
        assembler.state['cxn'],
        assembler.state['iteration'],
        assembler.args['bit_score'],
        assembler.args['contig_length'])

    expect = 'No new contigs were found in iteration {}'.format(
        assembler.state['iteration'])
    info.assert_called_once_with(expect)

@patch('lib.log.subcommand')
def test_assemble(subcommand):
    """It runs the assembler."""
    assembler = build_assembler()
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
