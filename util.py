import os
import sys
import logging

DEFAULT_SHARD_SIZE = 2.5e8


def setup_log(args):
    logging.basicConfig(
        filename='{}{}.log'.format(args.out, sys.argv[0][:-3]),
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))


def default_shard_count(sra_files):
    total_fasta_size = 0
    for sra_file in sra_files:
        file_size = os.path.getsize(sra_file)
        if sra_file.lower().endswith('.fastq'):
            file_size /= 2  # Guessing that fastq files are about twice as big as fasta files
        total_fasta_size += file_size
    shard_count = int(total_fasta_size / DEFAULT_SHARD_SIZE)
    return shard_count if shard_count else 1  # We need at least one shard


def default_process_count():
    return os.cpu_count() - 2 if os.cpu_count() > 2 else 1
