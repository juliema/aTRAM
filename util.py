"""Functions that handle common aTRAM chores."""

import os
import sys
import logging
import glob


def log_setup(config):
    """
    Set up the logs for a common format. We need the prefix of the output log file name and
    the command-line arguments for the starting message. Both are gotten from the user input.
    """
    file_name = '{}{}.log'.format(config['file_prefix'], sys.argv[0][:-3])
    logging.basicConfig(
        filename=os.path.join(config['work_dir'], file_name),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))


def db_name(config):
    """Create a standard DB name used across all programs."""
    file_name = '{}sqlite.db'.format(config['file_prefix'])
    return os.path.join(config['work_dir'], file_name)


def blast_shard_name(config, iteration):
    """Standardize the BLAST shard DB names."""
    file_name = '{}blast_{}'.format(config['file_prefix'], str(iteration + 1).zfill(3))
    return os.path.join(config['work_dir'], file_name)


def get_blast_shards(config):
    """Get all of the BLAST DB names built by format_sra."""
    file_name = '{}blast_*.nhr'.format(config['file_prefix'])
    pattern = os.path.join(config['work_dir'], file_name)
    files = glob.glob(pattern)
    return sorted([f[:-4] for f in files])


def blast_result_file(shard, iteration):
    """Get the name of the blast result file. Can use wildcards."""
    return '{}_{}.txt'.format(shard, str(iteration).zfill(2))


def paired_end_file(config, iteration, end):
    """Standardize the name of the paired file."""
    file_name = '{}matching_seqs_{}_{}.fasta'.format(config['file_prefix'], iteration, end)
    return os.path.join(config['work_dir'], file_name)
