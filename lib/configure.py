"""Handle configuration options and setup."""

import os
import argparse
import textwrap
from datetime import date
import psutil
from lib.dict_attrs import DictAttrs


class Configure:
    """Class for handling configuration options."""

    default = {
        'abyss': 'abyss-pe',
        'bit_score': 70.0,
        'data_prefix': 'atram_db_' + date.today().isoformat(),
        'evalue': 1e-9,
        'genetic_code': 1,
        'iterations': 5,
        'kmer': 31,
        'max_memory': '50G',
        'max_target_seqs': 100000000,
        'shard_size': 2.5e8,
    }

    args = {
        'assembler': lambda p: p.add_argument(
            '-a', '--assembler', required=True,
            help=('Which assembler to use. Required. '
                  '(One of: Trinity, Abyss, Velvet)')),

        'bit_score': lambda p: p.add_argument(
            '-b', '--bit-score', type=float,
            default=Configure.default['bit_score'],
            help=('Remove contigs that have a value less than this. '
                  'The default is "{}"').format(
                    Configure.default['bit_score'])),

        'cpus': lambda p: p.add_argument(
            '-c', '--cpus', type=int, default=0,
            help='Number of cpus to use.'),

        'data_prefix_atram': lambda p: p.add_argument(
            '-d', '--db-prefix', default=Configure.default['data_prefix'],
            help=('This needs to match the --db-prefix you entered for '
                  'atram_preprocessor.py.')),

        'data_prefix_preprocessor': lambda p: p.add_argument(
            '-d', '--db-prefix', default=Configure.default['data_prefix'],
            help=('This will get prepended to all blast and database files '
                  'so you can identify different database sets. These files '
                  'will be placed into the --work-dir. '
                  'The default is "{}"').format(
                    Configure.default['data_prefix'])),

        'evalue': lambda p: p.add_argument(
            '-e', '--evalue', type=float, default=Configure.default['evalue'],
            help='The default evalue is "{}".'.format(
                Configure.default['evalue'])),

        'genetic_code': lambda p: p.add_argument(
            '-g', '--genetic-code', type=int,
            default=Configure.default['genetic_code'],
            help=('The genetic code to use during blast runs. '
                  'The default is "{}".').format(
                Configure.default['genetic_code'])),

        'iterations': lambda p: p.add_argument(
            '-i', '--iterations', type=int,
            default=Configure.default['iterations'],
            help=('The number of pipline iterations. '
                  'The default is "{}".').format(
                      Configure.default['iterations'])),

        'kmer': lambda p: p.add_argument(
            '-k', '--kmer', type=int, default=Configure.default['kmer'],
            help=('k-mer size for assembers that use it. '
                  'The default is "{}".').format(Configure.default['kmer'])),

        'max_memory': lambda p: p.add_argument(
            '-m', '--max_memory', default=0,
            help='Maximum amount of memory to use.'),

        'max_target_seqs': lambda p: p.add_argument(
            '-M', '--max-target-seqs',
            type=int, default=Configure.default['max_target_seqs'],
            help='Maximum hit sequences per shard. Default is "{}".'.format(
                Configure.default['max_target_seqs'])),

        'output_prefix': lambda p: p.add_argument(
            '-o', '--output-prefix', required=True,
            help='The give the aTRAM output files this prefix. Required.'),

        'protein': lambda p: p.add_argument(
            '-p', '--protein', nargs='?', const=True, default=False,
            help='Are the target sequences protein?'),

        'shard_count': lambda p: p.add_argument(
            '-s', '--shard-count', type=int, default=0,
            help='Number of blast DB shards to create.'),

        'sra_files': lambda p: p.add_argument(
            'sra_files', nargs='+',
            help=('Sequence read archives in fasta or fastq format. '
                  'You may use wildcards or list multiple files.')),

        'target': lambda p: p.add_argument(
            '-t', '--target', required=True,
            help=('The path to the fasta file with sequences of interest. '
                  'Required.')),

        'work_dir': lambda p: p.add_argument(
            '-w', '--work-dir', default='.',
            help=('Where to store temporary files and files needed by other '
                  'aTRAM programs and other temporary files. Defaults to '
                  'the current working directory.')),
    }

    def __init__(self):
        self.config = None

    def parse_command_line(self, description='', args=None):
        """Parse the commandline arguments and return a dict attribute object.
        """

        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=textwrap.dedent(description))
        self.add_arguments(parser, args)

        self.config = vars(parser.parse_args())

        self.set_environment_defaults()

        return DictAttrs(self.config)

    @staticmethod
    def add_arguments(parser, args):
        """Add command-line arguments. We're given a string of argument keys,
        and we split that and look up the command line parser action in the
        dispatch table above.
        """

        for arg in args.split():
            Configure.args[arg](parser)

    def set_environment_defaults(self):
        """Add defaults that require calulation base on the environment."""

        environment_defaults = {
            'cpus': self.default_cpu_count,
            'max_memory': self.default_max_memory,
            'shard_count': self.default_shard_count,
        }

        for arg, func in environment_defaults.items():
            if arg in self.config and self.config[arg] < 1:
                self.config[arg] = func()

    def default_shard_count(self):
        """Default the shard count to the total fasta file size / 250MB.
        """

        total_fasta_size = 0
        for sra_file in self.config['sra_files']:
            file_size = os.path.getsize(sra_file)
            if sra_file.lower().endswith('.fastq'):
                file_size /= 2  # Guessing that fastq files ~2x fasta files
            total_fasta_size += file_size

        shard_count = int(total_fasta_size / Configure.default['shard_size'])

        return shard_count if shard_count else 1  # We need at least one shard

    @staticmethod
    def default_cpu_count():
        """If we're not given a default process count use the number of CPUs
        minus 2.
        """

        return os.cpu_count() - 2 if os.cpu_count() > 2 else 1

    @staticmethod
    def default_max_memory():
        """Some assemblers want to know how much memory they can use.
        Default to available memory minus 2G.
        """
        return '{}G'.format(int(psutil.virtual_memory()[0] / (2**30)) - 2)
