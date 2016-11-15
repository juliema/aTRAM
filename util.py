import sys
import logging


def setup_log(args):
    logging.basicConfig(
        filename='{}{}.log'.format(args.out, sys.argv[0][:-3]),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))
