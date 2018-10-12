"""Testing functions in lib/bio."""

# pylint: disable=missing-docstring, too-many-arguments, no-self-use
# flake8: noqa


from hypothesis import given
import hypothesis.strategies as st
import lib.bio as bio

DNA = 'ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx'
REVERSIBLE = 'ACGTWSMKRYBDHVNXacgtwsmkrybdhvnx'  # 'U's removed
PROTEIN = 'EFILPQefilpq'
HAS_PROTEIN = '[{}]*[{}]+[{}]*'.format(DNA, PROTEIN, DNA)


@given(seq=st.text(alphabet=REVERSIBLE))
def test_reverse_complement_twice(seq):
    actual = bio.reverse_complement(bio.reverse_complement(seq))
    assert seq == actual


@given(seq=st.text(alphabet=DNA))
def test_is_protein_no(seq):
    assert not bio.is_protein(seq)


@given(seq=st.from_regex(HAS_PROTEIN))
def test_is_protein_yes(seq):
    assert bio.is_protein(seq)
