"""Wrapper for the Velvet assembler."""

import shutil

from .base import BaseAssembler


class VelvetAssembler(BaseAssembler):
    """Wrapper for the Velvet assembler."""

    def __init__(self, args, cxn, log):
        """Build the assembler."""
        super().__init__(args, cxn, log)
        self.steps = [self.velveth, self.velvetg]

    @staticmethod
    def parse_contig_id(header):
        """Given a fasta header line return the contig ID."""
        return header

    def velveth(self):  # noqa
        """Build the velveth for the first assembly step."""
        cmd = ['velveth',
               '{}'.format(self.work_path()),
               '{}'.format(self.args['velvet_kmer']),
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

        if self.file['long_reads'] and not self.args['velvet_no_long']:
            cmd.append("-long '{}'".format(self.file['long_reads']))

        return ' '.join(cmd)

    def velvetg(self):
        """Build the velvetg for the second assembly step."""
        cmd = ['velvetg',
               '{}'.format(self.work_path()),
               '-ins_length {}'.format(self.args['velvet_ins_length']),
               '-exp_cov {}'.format(self.args['velvet_exp_cov']),
               '-min_contig_lgth {}'.format(
                   self.args['velvet_min_contig_lgth'])]

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output."""
        src = self.iter_file('contigs.fa')
        shutil.move(src, self.file['output'])

    @staticmethod
    def command_line_args(parser):
        """Add command-line arguments for this assembler."""
        group = parser.add_argument_group('optional Velvet arguments')

        group.add_argument(
            '--velvet-no-long', action='store_true',
            help="""Do not use long reads during assembly. for the assemblers
                is to use long reads. So this argument will stop the following:
                -long '<file>'.""")

        group.add_argument(
            '--velvet-kmer', type=int, default=31,
            help="""k-mer size (default 31).""")

        group.add_argument(
            '--velvet-exp-cov', type=int, default=30,
            help="""The expected coverage of the region. -exp_cov <int>
                (default 30).""")

        group.add_argument(
            '--velvet-ins-length', type=int, default=300,
            help="""The size of the fragments used in the short-read library.
                -ins_length <int> (default 300).""")

        group.add_argument(
            '--velvet-min-contig-lgth', type=int, default=100,
            help="""The minimum contig length used by the assembler itself.
                -min_contig_lgth <int> (default is 100).""")
