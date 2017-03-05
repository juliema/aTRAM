"""Utilities for working with sequences."""

COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx')


def reverse_complement(seq):
    """Reverse complement a nucleotide squence. We added some wildcards."""
    return seq.translate(COMPLEMENT)[::-1]
