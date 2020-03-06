"""Wrapper for the Spades assembler."""

from os.path import join
import shutil
from .base import BaseAssembler


class SpadesAssembler(BaseAssembler):
    """Wrapper for the Spades assembler."""

    def __init__(self, args, cxn):
        """Build the assembler."""
        super().__init__(args, cxn)
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
               '--threads {}'.format(self.args['cpus']),
               '--memory {}'.format(self.args['max_memory']),
               '--cov-cutoff {}'.format(self.args['cov_cutoff']),
               '-o {}'.format(self.work_path())]

        if self.args['careful']:
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
