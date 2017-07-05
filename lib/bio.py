"""Utilities for working with sequences."""

import re

COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx')

IS_PROTEIN = re.compile(r'[EFILPQ]', re.IGNORECASE)


def reverse_complement(seq):
    """Reverse complement a nucleotide sequence. We added some wildcards."""

    return seq.translate(COMPLEMENT)[::-1]


def is_protein(seq):
    """A simple to test to see if the sequence is a protein."""

    return IS_PROTEIN.search(seq)
