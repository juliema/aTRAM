"""Functions that handle common aTRAM chores."""

import sys
import logging
import glob


def log_setup(config):
    """
    Set up the logs for a common format. We need the prefix of the output log file name and
    the command-line arguments for the starting message. Both are gotten from the user input.
    """
    logging.basicConfig(
        filename='{}{}.log'.format(config['blast_db_prefix'], sys.argv[0][:-3]),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))


def blast_shard_name(config, iteration):
    """Standardize the BLAST shard DB names"""
    return '{}blast_{}'.format(config['blast_db_prefix'], str(iteration + 1).zfill(3))


def get_blast_shards(config):
    """Get all of the BLAST DB names built by format_sra."""
    pattern = '{}blast_*.nhr'.format(config['blast_db_prefix'])
    files = glob.glob(pattern)
    return sorted([f[:-4] for f in files])


def blast_result_file(shard, iteration):
    """Get the name of the blast result file. Can use wild"""
    return '{}_{}.txt'.format(shard, str(iteration).zfill(2))
