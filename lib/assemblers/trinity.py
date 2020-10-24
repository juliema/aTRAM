"""Wrapper for the Trinity assembler."""

from os.path import join
from shutil import move

import psutil

from .base import BaseAssembler


class TrinityAssembler(BaseAssembler):
    """Wrapper for the trinity assembler."""

    def __init__(self, args, cxn, log):
        """Build the assembler."""
        super().__init__(args, cxn, log)
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
               '--max_memory {}G'.format(self.args['trinity_max_memory']),
               '--CPU {}'.format(self.args['cpus']),
               "--output '{}'".format(self.work_path()),
               '--full_cleanup']

        if not self.args['trinity_bowtie2']:
            cmd.append('--no_bowtie')

        if self.file['paired_count']:
            cmd.append("--left '{}'".format(self.file['paired_1']))
            cmd.append("--right '{}'".format(self.file['paired_2']))
        else:
            single_ends = self.get_single_ends()
            if single_ends:
                cmd.append("--single '{}'".format(','.join(single_ends)))

        if self.file['long_reads'] and not self.args['trinity_no_long_reads']:
            cmd.append("--long_reads '{}'".format(self.file['long_reads']))

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output."""
        src = join(self.state['iter_dir'], 'trinity.Trinity.fasta')
        move(src, self.file['output'])

    @staticmethod
    def command_line_args(parser):
        """Add command-line arguments for this assembler."""
        group = parser.add_argument_group('optional Trinity arguments')

        group.add_argument(
            '--trinity-no-long-reads', action='store_true',
            help="""Do not use long reads during assembly. for the assemblers
                is to use long reads. So this argument will stop the following:
                It passes --long_reads '<file>'.""")

        total_mem = psutil.virtual_memory().available >> 30
        max_mem = max(1.0, total_mem >> 1)
        group.add_argument(
            '--trinity-max-memory', default=max_mem, metavar='MEMORY',
            type=int,
            help="""Maximum amount of memory to use in gigabytes. We will 
                use {} out of {} GB of free/unused memory. It passes
                --max_memory <int>G.""".format(max_mem, total_mem))

        group.add_argument(
            '--trinity-bowtie2', action='store_true',
            help="""Use bowtie2 during assembly. This will prevent --no_bowtie
                from being passed to Trinity.""")
