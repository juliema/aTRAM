"""Utilities for working with sequences."""

import re
from Bio import SeqIO

COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx')

IS_PROTEIN = re.compile(r'[EFILPQ]', re.IGNORECASE)


def reverse_complement(seq):
    """Reverse complement a nucleotide sequence. We added some wildcards."""

    return seq.translate(COMPLEMENT)[::-1]


def is_protein(seq):
    """A simple to test to see if the sequence is a protein."""

    return IS_PROTEIN.search(seq)


def fasta_file_has_protein(has_protein, query_file):
    """Check if the fasta file has an obviously protein sequences in it. It the
    user has told us that we have a protein then return that."""

    if not has_protein:
        return True

    with open(query_file) as in_file:
        for query in SeqIO.parse(in_file, 'fasta'):
            if is_protein(str(query.seq)):
                return True

    return False
