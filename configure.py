"""Handle configuration options."""

import os
import configparser

DEFAULT = {
    'shard_size': 2.5e8,
}


class DictAttrs(dict):
    """Allow dot notation for accessing dict items."""
    __getattr__ = dict.__getitem__


def read_config_file(config, args):
    """Read in the config file and hoist the data into the config dict."""
    if not args.get('config_file', None):
        pass
    parser = configparser.ConfigParser()
    parser.read(args['config_file'])
    return config


def get_configs(args=None):
    """Setup program configuration based on defaults, a config file, and command-line args."""
    if not args:
        args = {}
    config = DEFAULT

    DictAttrs(config)


def default_shard_count(sra_files):
    """If we're not given an input shard count use the fasta file size / 250MB."""
    total_fasta_size = 0
    for sra_file in sra_files:
        file_size = os.path.getsize(sra_file)
        if sra_file.lower().endswith('.fastq'):
            file_size /= 2  # Guessing that fastq files are about twice as big as fasta files
        total_fasta_size += file_size
    shard_count = int(total_fasta_size / DEFAULT['shard_size'])
    return shard_count if shard_count else 1  # We need at least one shard


def default_process_count():
    """If we're not given a default process count use the number of CPUs minus 2."""
    return os.cpu_count() - 2 if os.cpu_count() > 2 else 1


# Read the config file
# Overwrite command-line args
