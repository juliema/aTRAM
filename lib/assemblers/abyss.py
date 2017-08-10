"""Wrapper for the Abyss assembler."""

import os
import shutil
from lib.assemblers.base import BaseAssembler


class AbyssAssembler(BaseAssembler):
    """Wrapper for the Abyss assembler."""

    def __init__(self, args, db_conn):
        """Build the assembler."""
        super().__init__(args)
        self.steps = [self.abyss]

    def abyss(self):
        """Build the command for assembly."""
        cmd = ['abyss-pe',
               "-C '{}'".format(self.work_path),
               'E=0',
               'k={}'.format(self.args.kmer),
               "name='{}'".format(self.file['output'])]

        if self.args.mpi:
            cmd.append('np={}'.format(self.args.cpus))

        if self.file['paired_count']:
            cmd.append("in='{} {}'".format(
                self.file['paired_1'], self.file['paired_2']))

        single_ends = []
        if self.file['single_1_count']:
            single_ends.append(self.file['single_1'])
        if self.file['single_2_count']:
            single_ends.append(self.file['single_2'])
        if self.file['single_any_count']:
            single_ends.append(self.file['single_any'])
        if single_ends:
            cmd.append("se='{}'".format(' '.join(single_ends)))

        if self.file['long_reads'] and not self.args.no_long_reads:
            cmd.append("long='LONGREADS'")
            cmd.append("LONGREADS='{}'".format(self.file['long_reads']))

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output into the temp directory."""
        src = os.path.realpath(self.file['output'] + '-unitigs.fa')

        shutil.copyfile(src, self.file['output'])
