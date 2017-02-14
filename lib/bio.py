"""Utilities for filling in the gaps of BioPython or written for speed."""

import re


# Try to get the sequence name and which end it is from the fasta header
PARSE_HEADER = re.compile(r'^ [>@] \s* ( .* ) ( [\s/._] [12] ) \s* $',
                          re.VERBOSE)

# Parse blast hits file
PARSE_BLAST_RESULTS = re.compile(r'^ ( .* ) ( [\s\/_] [12] )', re.VERBOSE)


COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx')


def reverse_complement(seq):
    """Reverse complement a nucleotide squence. We added some wildcards."""
    return seq.translate(COMPLEMENT)[::-1]


# def makeblastdb(config, **kwargs):
#     """Build the command line for the makeblastdb program. BioPython seems
#     to be missing this.
#     """
#     pass
