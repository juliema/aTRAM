"""Wrapper for the Spades assembler."""

import os
import shutil
from os.path import join

import psutil

from .base import BaseAssembler


class SpadesAssembler(BaseAssembler):
    """Wrapper for the Spades assembler."""

    def __init__(self, args, cxn, log):
        """Build the assembler."""
        super().__init__(args, cxn, log)
        self.steps = [self.spades]

    def work_path(self):
        """
        Create output directory name.

        It has has unique requirements.
        """
        return join(self.state['iter_dir'], 'spades')

    def spades(self):
        """Build the command for assembly."""
        cmd = ['spades.py ',
               '--only-assembler',
               '--threads {}'.format(self.args['spades_threads']),
               '--memory {}'.format(self.args['spades_memory']),
               '--cov-cutoff {}'.format(self.args['spades_cov_cutoff']),
               '-o {}'.format(self.work_path())]

        if self.args['spades_careful']:
            cmd.append('--careful')

        if self.file['paired_count']:
            cmd.append("--pe1-1 '{}'".format(self.file['paired_1']))
            cmd.append("--pe1-2 '{}'".format(self.file['paired_2']))

        if self.file['single_1_count']:
            cmd.append("--s1 '{}'".format(self.file['single_1']))
        if self.file['single_2_count']:
            cmd.append("--s1 '{}'".format(self.file['single_2']))
        if self.file['single_any_count']:
            cmd.append("--s1 '{}'".format(self.file['single_any']))

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output."""
        src = join(self.work_path(), 'contigs.fasta')
        shutil.move(src, self.file['output'])

    @staticmethod
    def command_line_args(parser):
        """Add command-line arguments for this assembler."""
        group = parser.add_argument_group('optional Spades arguments')

        cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
        group.add_argument(
            '--spades-threads', type=int, default=cpus,
            help="""Number of threads to use.
                (default %(default)s)""".format(cpus))

        total_mem = psutil.virtual_memory().available >> 30
        max_mem = max(1.0, total_mem >> 1)
        group.add_argument(
            '--spades-memory', default=max_mem, metavar='MEMORY', type=int,
            help="""Maximum amount of memory to use in gigabytes. We will 
                use {} out of {} GB of free/unused memory. It passes
                --memory <int> argument to Spades.""".format(
                max_mem, total_mem))

        group.add_argument(
            '--spades-careful', action='store_true',
            help="""Tries to reduce the number of mismatches and short indels.
                It passes --careful.""")

        group.add_argument(
            '--spades-cov-cutoff', default='off',
            help="""Read coverage cutoff value. Must be a positive float value,
                or "auto", or "off". It passes  --cov-cutoff
                <keyword or int>. (default %(default)s)""")

    @staticmethod
    def validate_cov_cutoff(log, cov_cutoff):
        """Calculate default coverage cutoff argument."""
        if cov_cutoff in ['off', 'auto']:
            return cov_cutoff

        err = ('Read coverage cutoff value. Must be a positive '
               'float value, or "auto", or "off"')
        value = None
        try:
            value = float(cov_cutoff)
        except ValueError:
            log.fatal(err)

        if value < 0:
            log.fatal(err)

        return cov_cutoff
