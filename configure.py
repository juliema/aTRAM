"""Handle configuration options and setup."""

# ??? Read the config file
# ??? Overwrite config file settings with command-line args

import os
import configparser
from dict_attrs import DictAttrs

DEFAULT = DictAttrs({
    'shard_size': 2.5e8,
    'evalue': 1e-9,
    'iterations': 5,
    'max_target_seqs': 100000000,
})


def read_config_file(config, args):
    """Read in the config file and hoist the data into the config dict."""
    if not args.get('config_file', None):
        pass
    parser = configparser.ConfigParser()
    parser.read(args['config_file'])
    return config


def default_shard_count(config):
    """If we're not given an input shard count use the fasta file size / 250MB."""
    total_fasta_size = 0
    for sra_file in config.sra_files:
        file_size = os.path.getsize(sra_file)
        if sra_file.lower().endswith('.fastq'):
            file_size /= 2  # Guessing that fastq files are about twice as big as fasta files
        total_fasta_size += file_size
    shard_count = int(total_fasta_size / DEFAULT.shard_size)
    return shard_count if shard_count else 1  # We need at least one shard


def default_process_count():
    """If we're not given a default process count use the number of CPUs minus 2."""
    return os.cpu_count() - 2 if os.cpu_count() > 2 else 1


def setup_config_file():
    """Make our best guess for the configurations."""


def parse_args(parser):
    """Parse the commandline arguments and return a dict attribute object."""
    args = parser.parse_args()

    # Default for shard count requires calulation after the args are parsed
    shards = args.get('shards', None)
    if shards is not None and shards <= 0:
        args['shards'] = default_shard_count(args)

    return DictAttrs(vars(args))


def add_argument(parser, arg):
    """Add command-line arguments. We want to keep them consistent between programs."""
    if arg == 'evalue':
        parser.add_argument(
            '-e', '--evalue', default=DEFAULT.evalue, type=float,
            help='The default evalue is {}.'.format(DEFAULT.evalue))

    elif arg == 'iterations':
        parser.add_argument(
            '-i', '--iterations', default=DEFAULT.iterations, type=int,
            help=('The number of pipline iterations. '
                  'The default is {}.').format(DEFAULT.iterations))

    elif arg == 'out':
        parser.add_argument(
            '-o', '--out',
            help='Output aTRAM files with this prefix. May include a directory in the prefix.')

    elif arg == 'max_target_seqs':
        parser.add_argument(
            '-M', '--max-target-seqs', type=int, default=DEFAULT.max_target_seqs,
            help='Maximum hit sequences per shard. Default is {}.'.format(DEFAULT.max_target_seqs))

    elif arg == 'process_count':
        parser.add_argument(
            '-p', '--processes-count', type=int, help='Number of processes to create.',
            default=default_process_count())

    if arg == 'protein':
        parser.add_argument(
            '-p', '--protein', nargs='?', const=True, default=False,
            help='Is the target sequence a protein?')

    elif arg == 'shard_count':
        parser.add_argument(
            '-s', '--shards-count', type=int, help='Number of SRA shards to create.', default=-1)

    elif arg == 'sra_files':
        parser.add_argument(
            'sra_files', nargs='+',
            help='Short read archives in fasta or fastq format. May contain wildcards.')


if __name__ == '__main__':
    setup_config_file()
