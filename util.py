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


def path(file_name, config, iteration=None):
    """Standardize file names with a work directory and a file prefix."""
    if iteration:
        file_name = file_name.format(str(iteration).zfill(2))
    file_name = '{}{}'.format(config['file_prefix'], file_name)
    return os.path.join(config['work_dir'], file_name)


def db_file(config):
    """Create an SQLite3 DB name."""
    return path('sqlite.db', config)


def shard_db_names(config):
    """Get all of the BLAST DB names built by format_sra."""
    pattern = path('blast_*.nhr', config)
    files = glob.glob(pattern)
    return sorted([f[:-4] for f in files])


def shard_db_name(config, shard_index):
    """Create the BLAST shard DB names."""
    file_name = 'blast_{}'.format(str(shard_index + 1).zfill(3))
    return path(file_name, config)


def paired_end_file(config, iteration, end):
    """Create the file name of the paired end file."""
    file_name = 'matching_seqs_{}_{}.fasta'.format('{}', end)
    return path(file_name, config, iteration=iteration)


def contig_blast_file(config, iteration):
    """Create the file name of the blast DB for the assembled contigs."""
    return path('blast_contigs_{}', config, iteration=iteration)


def contig_score_file(config, iteration):
    """Create the contig blast result file name."""
    return path('contig_scores_{}.txt', config, iteration=iteration)


def contig_unfiltered_file(config, iteration):
    """Create the file name for the contigs before they are filtered."""
    return path('raw_contigs_{}.fasta', config, iteration=iteration)


def contig_filtered_file(config, iteration):
    """Create the file name for the contigs after they are filtered."""
    return path('contigs_{}.fasta', config, iteration=iteration)


def blast_result_file(shard_name, iteration):
    """Get the file name of the blast result file."""
    return '{}_{}.txt'.format(shard_name, str(iteration).zfill(2))
