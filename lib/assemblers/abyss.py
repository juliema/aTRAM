"""Wrapper for the Abyss assembler."""

from os.path import realpath
from shutil import copyfile

from .base import BaseAssembler


class AbyssAssembler(BaseAssembler):
    """Wrapper for the Abyss assembler."""

    def __init__(self, args, cxn, log):
        """Build the assembler."""
        super().__init__(args, cxn, log)
        self.steps = [self.abyss]

    def abyss(self):
        """Build the command for assembly."""
        cmd = ['abyss-pe',
               "-C '{}'".format(self.work_path()),
               'E=0',
               'k={}'.format(self.args['abyss_kmer'])]

        if self.args.get('abyss_p') is not None:
            cmd.append('p={}'.format(self.args['abyss_p']))

        if self.args.get('abyss_n') is not None:
            cmd.append('n={}'.format(self.args['abyss_n']))

        if self.args.get('abyss_N') is not None:
            cmd.append('N={}'.format(self.args['abyss_N']))

        cmd.append("name='{}'".format(self.file['output']))

        if self.args.get('abyss_np'):
            cmd.append('np={}'.format(self.args['abyss_np']))

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

        if (self.file['long_reads']
                and not self.args.get('abyss_no_long')):
            cmd.append("long='LONGREADS'")
            cmd.append("LONGREADS='{}'".format(self.file['long_reads']))

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output into the temp directory."""
        src = realpath(self.file['output'] + '-unitigs.fa')

        copyfile(src, self.file['output'])

    @staticmethod
    def command_line_args(parser):
        """Add command-line arguments for this assembler."""
        group = parser.add_argument_group('optional Abyss arguments')

        group.add_argument(
            '--abyss-no-long', action='store_true',
            help="""Do not use long reads during assembly. for the assemblers
                is to use long reads. So this argument will stop the following:
                long='LONGREADS' LONGREADS='<file>'.""")

        group.add_argument(
            '--abyss-kmer', type=int, default=64,
            help="""k-mer size. I passes k=<int> (default 64).""")

        group.add_argument(
            '--abyss-np', type=int,
            help="""Abyss must have been compiled to use np.
                It passes np=<integer>.""")

        group.add_argument(
            '--abyss-paired-ends', action='store_true',
            help="""Abyss: Normally, we put all of the input sequences in to
                the -se argument. If this is true then we will put paired end
                sequences into the -in argument with any residual single ends
                into the -se argument.""")

        group.add_argument(
            '--abyss-p', type=int,
            help="""Minimum sequence identity of a bubble. Enter a
                single digit integer [0-9].""")

        group.add_argument(
            '--abyss-n', type=int,
            help="""Minimum number of pairs required for building
                contigs [10].""")

        group.add_argument(
            '--abyss-N', type=int,
            help="""Minimum number of pairs required for building
                scaffolds [n].""")
