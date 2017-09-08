"""Test object in lib/assemblers/base."""

# pylint: disable=missing-docstring

from os.path import basename, join
from hypothesis import given
import hypothesis.strategies as st
from lib.assemblers.base import BaseAssembler
import tests.mock as mock


PATH_PATTERN = r'(?:[\w.-]+/)*[\w.-]+'


@given(args=st.text(), db_conn=st.text())
def test_init(args, db_conn):
    assembler = BaseAssembler(args, db_conn)

    assert args == assembler.args
    assert assembler.blast_only is False
    assert [] == assembler.steps
    assert {} == assembler.file

    expected = {
        'iteration': 0,
        'query_file': '',
        'blast_db': '',
        'db_conn': db_conn}
    assert expected == assembler.state


@given(
    blast_db=st.from_regex(PATH_PATTERN),
    query_file=st.from_regex(PATH_PATTERN),
    iteration=st.integers())
def test_initialize_iteration(blast_db, query_file, iteration):
    output_files = ['prefix/dir/' + f for f
                    in ['output.fasta',
                        'paired_1.fasta',
                        'paired_2.fasta',
                        'single_1.fasta',
                        'single_2.fasta',
                        'single_any.fasta']]
    assembler = BaseAssembler('args', 'db_conn')

    mock.it(assembler, 'iter_file', output_files)

    assembler.initialize_iteration(blast_db, query_file, iteration)

    assert basename(blast_db) == assembler.state['blast_db']
    assert basename(query_file) == assembler.state['query_file']
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
def test_iter_dir(args, blast_db, query_file, iteration):
    assembler = BaseAssembler(args, 'db_conn')
    assembler.set_state(blast_db, query_file, iteration)

    base_blast_db = blast_db.split('/')[-1]
    base_query_file = query_file.split('/')[-1]

    file_name = '{}_{}_iteration_{:02d}'.format(
        base_blast_db, base_query_file, iteration)
    expect = join(args['temp_dir'], file_name)
    assert expect == assembler.iter_dir()
