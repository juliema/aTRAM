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
    print(config)
    file_name = '{}{}.log'.format(config['file_prefix'], sys.argv[0][:-3])
    logging.basicConfig(
        filename=os.path.join(config['work_dir'], file_name),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))


def get_shard_file_names(config):
    """Get all of the BLAST DB names built by format_sra."""
    file_name = '{}blast_*.nhr'.format(config['file_prefix'])
    pattern = os.path.join(config['work_dir'], file_name)
    files = glob.glob(pattern)
    return sorted([f[:-4] for f in files])


def path(config, file_name):
    """Standardize file names with a work directory and a file prefix."""
    file_name = config['file_prefix'] + file_name
    return os.path.join(config['work_dir'], file_name)


def db_file(config):
    """Create an SQL DB name."""
    file_name = 'sqlite.db'
    return path(config, file_name)


def blast_shard_file(config, iteration):
    """Standardize the BLAST shard DB names."""
    file_name = 'blast_{}'.format(str(iteration).zfill(2))
    return path(config, file_name)


def blast_contig_file(config, iteration):
    """Get the name of the blast result file."""
    file_name = 'contig_scores_{}.txt'.format(str(iteration).zfill(2))
    return path(config, file_name)


def contig_score_file(config, iteration):
    """Get the name of the blast DB for the assembled contigs."""
    file_name = 'blast_contigs_{}'.format(str(iteration).zfill(2))
    return path(config, file_name)


def paired_end_file(config, iteration, end):
    """Standardize the name of the paired file."""
    file_name = 'matching_seqs_{}_{}.fasta'.format(str(iteration).zfill(2), end)
    return path(config, file_name)


def raw_contig_file(config, iteration):
    """Standardize the contig file name from before it is filtered."""
    file_name = 'raw_contigs_{}.fasta'.format(str(iteration).zfill(2))
    return path(config, file_name)


def contig_file(config, iteration):
    """Standardize the contig file name from after it is filtered."""
    file_name = 'contigs_{}.fasta'.format(str(iteration).zfill(2))
    return path(config, file_name)


def blast_result_file(shard, iteration):
    """Get the name of the blast result file."""
    return '{}_{}.txt'.format(shard, str(iteration).zfill(2))
