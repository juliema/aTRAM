"""Wrapper for the Abyss assembler."""

from os.path import realpath
from shutil import copyfile
from .base import BaseAssembler


class AbyssAssembler(BaseAssembler):
    """Wrapper for the Abyss assembler."""

    def __init__(self, args, cxn):
        """Build the assembler."""
        super().__init__(args, cxn)
        self.steps = [self.abyss]

    def abyss(self):
        """Build the command for assembly."""
        cmd = ['abyss-pe',
               "-C '{}'".format(self.work_path()),
               'E=0',
               'k={}'.format(self.args['kmer'])]

        if self.args.get('abyss_p') is not None:
            cmd.append('p={}'.format(self.args['abyss_p']))

        cmd.append("name='{}'".format(self.file['output']))

        if self.args['mpi']:
            cmd.append('np={}'.format(self.args['cpus']))

        if self.args.get('abyss_paired_ends'):
            if self.file['paired_count']:
                cmd.append("in='{} {}'".format(
                    self.file['paired_1'], self.file['paired_2']))
            single_ends = self.get_single_ends()
            if single_ends:
                cmd.append("se='{}'".format(' '.join(single_ends)))
        else:
            in_files = []
            if self.file['paired_count']:
                in_files += [self.file['paired_1'], self.file['paired_2']]
            in_files += self.get_single_ends()
            cmd.append("se='{}'".format(' '.join(in_files)))

        if self.file['long_reads'] and not self.args['no_long_reads']:
            cmd.append("long='LONGREADS'")
            cmd.append("LONGREADS='{}'".format(self.file['long_reads']))

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output into the temp directory."""
        src = realpath(self.file['output'] + '-unitigs.fa')

        copyfile(src, self.file['output'])
