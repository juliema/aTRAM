"""Testing functions in lib/log."""

# pylint: disable=missing-docstring, too-many-arguments, no-self-use


# from os.path import basename, splitext
# from hypothesis import given
# import hypothesis.strategies as st
# import lib.log as log
#
#
# @given(log_file=st.text(min_size=1), blast_db=st.text(), query_file=st.text())
# def test_file_name_default(log_file, blast_db, query_file):
#     actual = log.file_name(log_file, blast_db, query_file=query_file)
#     assert log_file == actual
#
#
# @given(log_file=st.just(''), blast_db=st.text())
# def test_file_name_no_query(log_file, blast_db):
#     actual = log.file_name(log_file, blast_db)
#     expect = '{}.pytest.log'.format(blast_db)
#     assert expect == actual
#
#
# @given(log_file=st.just(''), blast_db=st.text(), query=st.text(min_size=1))
# def test_file_name_with_query(log_file, blast_db, query):
#     actual = log.file_name(log_file, blast_db, query_file=query)
#     query = splitext(basename(query))[0]
#     expect = '{}.{}.pytest.log'.format(blast_db, query)
#     assert expect == actual
