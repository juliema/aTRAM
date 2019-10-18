"""Testing functions in lib/bio."""

import lib.bio as bio


def test_reverse_complement_01():
    """It complements then reverses the string."""
    seq = 'ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx'
    actual = bio.reverse_complement(seq)
    assert actual == 'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx'[::-1]


def test_reverse_complement_02():
    """A reverse complement twice yields the original sequence, except U."""
    seq = 'ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx'
    actual = bio.reverse_complement(bio.reverse_complement(seq))
    assert actual == 'ACGTTWSMKRYBDHVNXacgttwsmkrybdhvnx'


def test_is_protein_no():
    """DNA, RNA, and wildcards are not protein."""
    seq = 'ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx'
    assert not bio.is_protein(seq)


def test_is_protein_yes():
    """Any protein character makes the whole sequence a protein."""
    seq = 'ACGTUWSMKRYBeDHVNXacgtuwsmkrybdhvnx'
    assert bio.is_protein(seq)
