"""Testing functions in lib/blast."""

# pylint: disable=missing-docstring, too-many-arguments, no-self-use

# from hypothesis import given
# import hypothesis.strategies as st
# import lib.blast as blast
#
#
# @given(
#     title=st.text(min_size=1),
#     ends=st.data.draw(sampled_from([])),
#     seq_end_clamp=st.data.draw(sampled_from([])))
# def test_parse_fasta_title_01():
#     expected_seq_name = ''
#     expected_seq_end = ''
#     actual_seq_name, actual_seq_end = blast.parse_fasta_title(
#         '', '', '')
#     assert expected_seq_name == actual_seq_name
#     assert expected_seq_end == actual_seq_end
