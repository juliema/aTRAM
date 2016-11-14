import sys
import argparse
import util


def parse_args():
    parser = argparse.ArgumentParser(
        description=''' ''')
    parser.add_argument('-o', '--out',
                        help='output aTRAM files with this prefix. May include a directory in the prefix.')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    util.open_log_file(args.out)
    util.log(' '.join(sys.argv))

    util.close_log_file()
