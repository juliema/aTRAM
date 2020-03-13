"""Wrapper for the Trinity assembler."""

from os.path import join
from shutil import move
from .base import BaseAssembler


class TrinityAssembler(BaseAssembler):
    """Wrapper for the trinity assembler."""

    def __init__(self, args, cxn):
        """Build the assembler."""
        super().__init__(args, cxn)
        self.steps = [self.trinity]

    def work_path(self):
        """
        Create output directory name.

        It has has unique requirements.
        """
        return join(self.state['iter_dir'], 'trinity')

    def trinity(self):
        """Build the command for assembly."""
        cmd = ['Trinity',
               '--seqType fa',
               '--max_memory {}G'.format(self.args['max_memory']),
               '--CPU {}'.format(self.args['cpus']),
               "--output '{}'".format(self.work_path()),
               '--full_cleanup']

        if not self.args['bowtie2']:
            cmd.append('--no_bowtie')

        if self.file['paired_count']:
            cmd.append("--left '{}'".format(self.file['paired_1']))
            cmd.append("--right '{}'".format(self.file['paired_2']))
        else:
            single_ends = self.get_single_ends()
            if single_ends:
                cmd.append("--single '{}'".format(','.join(single_ends)))

        if self.file['long_reads'] and not self.args['no_long_reads']:
            cmd.append("--long_reads '{}'".format(self.file['long_reads']))

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output."""
        src = join(self.state['iter_dir'], 'trinity.Trinity.fasta')
        move(src, self.file['output'])
