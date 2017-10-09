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
    """Is the sequence a protein."""
    return IS_PROTEIN.search(seq)


def fasta_file_has_protein(query_files):
    """
    Search for protein characters in a fasta file.

    If the user has told us that we have a protein then return that.
    """
    for query_file in query_files:
        with open(query_file) as in_file:
            for query in SeqIO.parse(in_file, 'fasta'):
                if is_protein(str(query.seq)):
                    return True

    return False
