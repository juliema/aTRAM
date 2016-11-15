import sys
import logging
import argparse
import util


COMPLEMENT = str.maketrans('ACGTUWSMKRYBDHVNXacgtuwsmkrybdhvnx',
                           'TGCAAWSKMYRVHDBNXtgcaawskmyrvhdbnx')


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_args():
    parser = argparse.ArgumentParser(
        description=''' ''')
    parser.add_argument('-o', '--out',
                        help='output aTRAM files with this prefix. May include a directory in the prefix.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    util.setup_log(args)
