"""Test object in lib/assemblers/base."""

import re
from os.path import join
import subprocess
from unittest.mock import MagicMock, patch, call
import pytest
from lib.assemblers.base import BaseAssembler


def build_assembler():
    """Build a generic assembler."""
    args = {
        'bit_score': 44,
        'contig_length': 55,
        'temp_dir': 'my_temp_dir',
        'assembler': 'my_assembler'}
    asm = BaseAssembler(args, 'my_db_conn')
    asm.state = {
        'blast_db': 'my_blast_db',
        'query_file': 'my_query_file.seq',
        'query_target': 'my_query_name',
        'iter_dir': 'my_iter_dir',
        'iteration': 99,
        'cxn': 'my_cxn'}
    return asm


def test_init_01():
    """It initializes an assembler object."""
    asm = BaseAssembler('args', 'cxn')

    assert asm.args == 'args'
    assert asm.blast_only is False
    assert asm.steps == []
    assert asm.file == {}
    expected = {
        'iteration': 0,
        'query_target': '',
        'query_file': '',
        'blast_db': '',
        'iter_dir': '',
        'cxn': 'cxn'}
    assert asm.state == expected


def test_init_iteration_01():
    """It sets up assembler variables for an iteration."""
    asm = build_assembler()

    blast_db = 'another_blast_db'
    query_file = 'another_query_file'
    iteration = 99

    asm.init_iteration(blast_db, query_file, iteration)

    assert asm.state['blast_db'] == blast_db
    assert asm.state['query_file'] == query_file
    assert asm.state['iteration'] == iteration


def test_setup_files_01():
    """It sets up assembler file names for an iteration."""
    asm = build_assembler()

    iter_dir = 'my_iter_dir'

    asm.setup_files(iter_dir)

    assert asm.state['iter_dir'] == iter_dir
    assert asm.file['long_reads'] == ''
    assert asm.file['output'].endswith(join(iter_dir, 'output.fasta'))
    assert asm.file['paired_1'].endswith(join(iter_dir, 'paired_1.fasta'))
    assert asm.file['paired_2'].endswith(join(iter_dir, 'paired_2.fasta'))
    assert asm.file['single_1'].endswith(join(iter_dir, 'single_1.fasta'))
    assert asm.file['single_2'].endswith(join(iter_dir, 'single_2.fasta'))
    assert asm.file['single_any'].endswith(join(iter_dir, 'single_any.fasta'))
    assert asm.file['paired_count'] == 0
    assert asm.file['single_1_count'] == 0
    assert asm.file['single_2_count'] == 0
    assert asm.file['single_any_count'] == 0


def test_file_prefix_01():
    """It builds the directory for the iteration."""
    asm = build_assembler()

    expected = '{}_{}_{:02d}_'.format(
        asm.state['blast_db'],
        asm.state['query_target'],
        asm.state['iteration'])

    assert asm.file_prefix() == expected


def test_iter_file_01():
    """It builds a file name with the iter dir."""
    file_name = 'my_file.txt'
    asm = build_assembler()
    expected = join(asm.state['iter_dir'], file_name)
    assert asm.iter_file(file_name).endswith(expected)


def test_work_path_01():
    """It builds the work directory."""
    asm = build_assembler()
    assert asm.work_path() == asm.state['iter_dir']


@patch('lib.log.info')
def test_run_01(info):
    """It performs a normal flow."""
    asm = build_assembler()

    asm.assemble = MagicMock()

    asm.run()

    info.assert_called_once_with(
        'Assembling shards with {}: iteration {}'.format(
            asm.args['assembler'], asm.state['iteration']))


@patch('lib.log.info')
@patch('lib.log.error')
def test_run_02(error, info):
    """It handles a timeout error."""
    args = {'assembler': 'my_assembler', 'timeout': 10}

    asm = BaseAssembler(args, 'cxn')
    asm.assemble = MagicMock(side_effect=TimeoutError())

    with pytest.raises(TimeoutError) as timeout_error:
        asm.run()

    error_msg = str(timeout_error.value)
    if error_msg[-1] != '.':
        error_msg += '.'
    assert error_msg == (
        'Time ran out for the assembler after 0:00:10 (HH:MM:SS).')

    expect = 'Assembling shards with {}: iteration {}'.format(
        args['assembler'], asm.state['iteration'])
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
    error_code = 88
    cmd = 'my command'
    error = subprocess.CalledProcessError(error_code, cmd)
    asm = build_assembler()
    asm.assemble = MagicMock(side_effect=error)

    with pytest.raises(RuntimeError) as runtime_error:
        asm.run()

    error_msg = str(runtime_error.value)
    if error_msg[-1] != '.':
        error_msg += '.'
    assert error_msg == (
        "The assembler failed with error: Command 'my command' "
        "returned non-zero exit status 88.")

    expect = 'Assembling shards with {}: iteration {}'.format(
        asm.args['assembler'], asm.state['iteration'])
    info.assert_called_once_with(expect)


@patch('lib.log.info')
@patch('lib.db_atram.sra_blast_hits_count')
def test_count_blast_hits_01(sra_blast_hits_count, info):
    """It handles no blast hits."""
    asm = build_assembler()

    sra_blast_hits_count.return_value = 0

    assert asm.count_blast_hits() == 0

    info.assert_called_once_with('0 blast hits in iteration {}'.format(
        asm.state['iteration']))


@patch('lib.log.info')
@patch('lib.db_atram.sra_blast_hits_count')
def test_count_blast_hits_02(sra_blast_hits_count, info):
    """It handles one blast hits."""
    asm = build_assembler()

    sra_blast_hits_count.return_value = 1

    assert asm.count_blast_hits() == 1

    info.assert_called_once_with('1 blast hits in iteration {}'.format(
        asm.state['iteration']))


@patch('lib.log.info')
def test_nothing_assembled_01(info):
    """It handles when nothing is assembled."""
    asm = build_assembler()
    asm.file['output'] = 'tests/data/missing_file.txt'

    assert asm.nothing_assembled()

    info.assert_called_once_with(
        'No new assemblies in iteration {}'.format(
            asm.state['iteration']))


@patch('lib.log.info')
def test_nothing_assembled_02(info):
    """It handles an empty assembly."""
    asm = build_assembler()
    asm.file['output'] = 'tests/data/empty_file.txt'

    assert asm.nothing_assembled()

    expect = 'No new assemblies in iteration {}'.format(
        asm.state['iteration'])
    info.assert_called_once_with(expect)


@patch('lib.log.info')
def test_nothing_assembled_03(info):
    """It handles when something is assembled."""
    asm = build_assembler()
    asm.file['output'] = 'tests/data/load_seq1.txt'

    assert not asm.nothing_assembled()

    info.assert_not_called()


@patch('lib.log.info')
@patch('lib.db_atram.assembled_contigs_count')
def test_assembled_contigs_count_01(assembled_contigs_count, info):
    """Handle when here are no contigs."""
    high_score = 5
    asm = build_assembler()

    assembled_contigs_count.return_value = 0

    assert asm.assembled_contigs_count(high_score) == 0

    assembled_contigs_count.assert_called_once_with(
        asm.state['cxn'],
        asm.state['iteration'],
        asm.args['bit_score'],
        asm.args['contig_length'])

    expect = ('No contigs had a bit score greater than {} and are at '
              'least {} bp long in iteration {}. The highest score for '
              'this iteration is {}').format(
                  asm.args['bit_score'],
                  asm.args['contig_length'],
                  asm.state['iteration'],
                  high_score)
    info.assert_called_once_with(expect)


@patch('lib.log.info')
@patch('lib.db_atram.assembled_contigs_count')
def test_assembled_contigs_count_02(assembled_contigs_count, info):
    """Handle when here is one contig."""
    high_score = 5
    count = 1
    asm = build_assembler()

    assembled_contigs_count.return_value = count

    assert asm.assembled_contigs_count(high_score) == count

    assembled_contigs_count.assert_called_once_with(
        asm.state['cxn'],
        asm.state['iteration'],
        asm.args['bit_score'],
        asm.args['contig_length'])

    info.assert_not_called()


@patch('lib.log.info')
@patch('lib.db_atram.iteration_overlap_count')
def test_no_new_contigs_01(iteration_overlap_count, info):
    """It handles when there are new contigs."""
    count = 1
    asm = build_assembler()

    iteration_overlap_count.return_value = count + 1

    assert not asm.no_new_contigs(count)

    iteration_overlap_count.assert_called_once_with(
        asm.state['cxn'],
        asm.state['iteration'],
        asm.args['bit_score'],
        asm.args['contig_length'])

    info.assert_not_called()


@patch('lib.log.info')
@patch('lib.db_atram.iteration_overlap_count')
def test_no_new_contigs_02(iteration_overlap_count, info):
    """It handles when there are no new contigs."""
    count = 1
    asm = build_assembler()

    iteration_overlap_count.return_value = count

    assert asm.no_new_contigs(count)

    iteration_overlap_count.assert_called_once_with(
        asm.state['cxn'],
        asm.state['iteration'],
        asm.args['bit_score'],
        asm.args['contig_length'])

    expect = 'No new contigs were found in iteration {}'.format(
        asm.state['iteration'])
    info.assert_called_once_with(expect)


@patch('lib.log.subcommand')
def test_assemble(subcommand):
    """It runs the assembler."""
    asm = build_assembler()
    asm.args['timeout'] = 333

    asm.post_assembly = MagicMock()
    asm.step1 = MagicMock(return_value='step1')
    asm.step2 = MagicMock(return_value='step2')

    asm.steps = [asm.step1, asm.step2]

    asm.assemble()

    calls = [
        call('step1',
             asm.args['temp_dir'],
             asm.args['timeout']),
        call('step2',
             asm.args['temp_dir'],
             asm.args['timeout'])]
    subcommand.assert_has_calls(calls)
