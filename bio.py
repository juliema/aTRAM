"""Utilities for filling in the gaps of BioPython or written for speed."""


COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx')


def reverse_complement(seq):
    """Reverse complement a nucleotide squence. We added some wildcards."""
    return seq.translate(COMPLEMENT)[::-1]


def makeblastdb(config, **kwargs):
    """
    Build the command line for the makeblastdb program. BioPython seems to be missing this.
    """
    cmd = [config.get('makeblastdb')]
    for key, val in kwargs.items():
        cmd.append('-{} {}'.format(key, val))
    return ' '.join(cmd)
