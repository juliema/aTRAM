"""Wrapper for the Velvet assembler."""

import shutil
from lib.assemblers.base import BaseAssembler


class VelvetAssembler(BaseAssembler):
    """Wrapper for the Velvet assembler."""

    def __init__(self, args, db_conn):
        """Build the assembler."""
        super().__init__(args, db_conn)
        self.steps = [self.velveth, self.velvetg]

    @staticmethod
    def parse_contig_id(header):
        """Given a fasta header line return the contig ID."""
        return header

    def velveth(self):
        """Build the velveth for the first assembly step."""
        cmd = ['velveth',
               '{}'.format(self.work_path()),
               '{}'.format(self.args['kmer']),
               '-fasta']

        if self.file['paired_count']:
            cmd.append("-shortPaired '{}' '{}'".format(
                self.file['paired_1'], self.file['paired_2']))

        single_ends = []
        if self.file['single_1_count']:
            single_ends.append("'{}'".format(self.file['single_1']))
        if self.file['single_2_count']:
            single_ends.append("'{}'".format(self.file['single_2']))
        if self.file['single_any_count']:
            single_ends.append("'{}'".format(self.file['single_any']))
        if single_ends:
            cmd.append("-short {}".format(' '.join(single_ends)))

        if self.file['long_reads'] and not self.args['no_long_reads']:
            cmd.append("-long '{}'".format(self.file['long_reads']))

        return ' '.join(cmd)

    def velvetg(self):
        """Build the velvetg for the second assembly step."""
        cmd = ['velvetg',
               '{}'.format(self.work_path()),
               '-ins_length {}'.format(self.args['ins_length']),
               '-exp_cov {}'.format(self.args['exp_coverage']),
               '-min_contig_lgth {}'.format(self.args['min_contig_length'])]

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output."""
        src = self.iter_file('contigs.fa')
        shutil.move(src, self.file['output'])
