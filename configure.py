"""Handle configuration options and setup."""

# ??? Read the config file
# ??? Overwrite config file settings with command-line args

import os
import configparser
from dict_attrs import DictAttrs

DEFAULT = {
    'shard_size': 2.5e8,
}


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


def default_shard_count(config):
    """If we're not given an input shard count use the fasta file size / 250MB."""
    if not config.shards:
        total_fasta_size = 0
        for sra_file in config.sra_files:
            file_size = os.path.getsize(sra_file)
            if sra_file.lower().endswith('.fastq'):
                file_size /= 2  # Guessing that fastq files are about twice as big as fasta files
            total_fasta_size += file_size
        shard_count = int(total_fasta_size / DEFAULT['shard_size'])
        config.shards = shard_count if shard_count else 1  # We need at least one shard


def default_process_count():
    """If we're not given a default process count use the number of CPUs minus 2."""
    return os.cpu_count() - 2 if os.cpu_count() > 2 else 1


def setup_config_file():
    """Make our best guess for the configurations."""


def add_argument(parser, arg):
    """Add command-line arguments. We want to keep them consistent between programs."""
    if arg == 'out':
        parser.add_argument(
            '-o', '--out',
            help='output aTRAM files with this prefix. May include a directory in the prefix.')
    if arg == 'processes':
        parser.add_argument(
            '-p', '--processes', type=int, help='number of processes to create',
            default=default_process_count())
    if arg == 'shards':
        parser.add_argument('-s', '--shards', type=int, help='number of shards to create')
    if arg == 'sra_files':
        parser.add_argument(
            'sra_files', nargs='+',
            help='short read archives in fasta or fastq format. May contain wildcards.')


if __name__ == '__main__':
    setup_config_file()
