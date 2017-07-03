"""Wrappers for the various assember programs."""

import os
import shutil
import lib.log as log
import lib.file_util as file_util


class Assembler:
    """A factory class for building the assembers."""

    @staticmethod
    def factory(args):
        """Return the assembler based upon the configuration options."""

        if args.assembler.lower() == 'abyss':
            return AbyssAssembler(args)
        elif args.assembler.lower() == 'trinity':
            return TrinityAssembler(args)
        elif args.assembler.lower() == 'velvet':
            return VelvetAssembler(args)
        elif args.assembler.lower() == 'spades':
            return SpadesAssembler(args)

    def __init__(self, args):
        self.args = args             # Parsed command line arguments
        self.steps = []              # Assembler steps setup by the assembler
        self.is_paired = False       # Did we find paired end sequences?
        self.output_file = None      # Write to this file
        self.long_reads_file = None  # Long-reads file
        self.end_1_file = None       # Sequneces for end 1 reads
        self.end_2_file = None       # Sequences for end 2 reads
        self.iteration = 0           # Current iteration

    @property
    def iter_dir(self):
        """Get the work directory for the current iteration."""

        return file_util.temp_iter_dir(self.args.temp_dir, self.iteration)

    @property
    def work_path(self):
        """The output directory name may have unique requirements."""

        return self.iter_dir

    def iter_file(self, file_name):
        """Build a temporary file name honoring the current iteration
        directory.
        """

        return os.path.join(self.iter_dir, file_name)

    def assemble(self):
        """Use the assembler to build up the contigs. We take and array of
        subprocess steps and execute them in order. We bracket this with
        pre and post assembly steps.
        """

        for step in self.steps:
            cmd = step()
            log.subcommand(cmd, self.args.temp_dir, self.args.timeout)

        self.post_assembly()

    def post_assembly(self):
        """Assembers have unique post assembly steps."""

    def path(self, file_name):
        """Files will go into the temp dir."""

        blast_db = os.path.basename(self.args.blast_db)
        file_name = '{}.{:02d}.{}'.format(blast_db, self.iteration, file_name)
        rel_path = self.iter_file(file_name)

        return os.path.abspath(rel_path)

    def initialize_iteration(self, iteration):
        """Files used by the assembler. Do this at the start of each
        iteration.
        """

        self.iteration = iteration
        self.output_file = self.path('output.fasta')
        self.end_1_file = self.path('paired_end_1.fasta')
        self.end_2_file = self.path('paired_end_2.fasta')

    @staticmethod
    def parse_contig_id(header):
        """Given a fasta header line return the contig ID."""

        return header.split()[0]


class AbyssAssembler(Assembler):
    """Wrapper for the Abyss assembler."""

    def __init__(self, args):
        super().__init__(args)
        self.steps = [self.abyss]

    def abyss(self):
        """Build the command for assembly."""

        cmd = ['abyss-pe',
               "-C '{}'".format(self.work_path),
               'E=0',
               'k={}'.format(self.args.kmer),
               "name='{}'".format(self.output_file)]

        if self.args.mpi:
            cmd.append('np={}'.format(self.args.cpus))

        if self.is_paired:
            cmd.append("in='{} {}'".format(self.end_1_file, self.end_2_file))
        else:
            cmd.append("se='{}'".format(self.end_1_file))

        if self.long_reads_file and not self.args.no_long_reads:
            cmd.append("long='LONGREADS'")
            cmd.append("LONGREADS='{}'".format(self.long_reads_file))

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output."""

        src = os.path.realpath(self.output_file + '-unitigs.fa')

        shutil.copyfile(src, self.output_file)


class TrinityAssembler(Assembler):
    """Wrapper for the trinity assembler."""

    @property
    def work_path(self):
        """The output directory name has unique requirements."""

        return os.path.join(self.iter_dir, 'trinity')

    def __init__(self, args):
        super().__init__(args)
        self.steps = [self.trinity]

    def trinity(self):
        """Build the command for assembly."""

        cmd = ['Trinity',
               '--seqType fa',
               '--max_memory {}G'.format(self.args.max_memory),
               '--CPU {}'.format(self.args.cpus),
               "--output '{}'".format(self.work_path),
               '--full_cleanup']

        if not self.args.bowtie2:
            cmd.append('--no_bowtie')

        if self.is_paired:
            cmd.append("--left '{}'".format(self.end_1_file))
            cmd.append("--right '{}'".format(self.end_2_file))
        else:
            cmd.append("-single '{}'".format(self.end_1_file))
            cmd.append('--run_as_paired')

        if self.long_reads_file and not self.args.no_long_reads:
            cmd.append("--long_reads '{}'".format(self.long_reads_file))

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output."""

        src = os.path.join(self.iter_dir, 'trinity.Trinity.fasta')
        shutil.move(src, self.output_file)


class VelvetAssembler(Assembler):
    """Wrapper for the Velvet assembler."""

    def __init__(self, args):
        super().__init__(args)
        self.steps = [self.velveth, self.velvetg]

    @staticmethod
    def parse_contig_id(header):
        """Given a fasta header line return the contig ID."""

        return header

    def velveth(self):
        """Build the velveth for the first assembly step."""

        cmd = ['velveth'
               '{}'.format(self.work_path),
               '{}'.format(self.args.kmer),
               '-fasta']

        if self.is_paired:
            cmd.append("-shortPaired '{}' '{}'".format(
                self.end_1_file, self.end_2_file))
        else:
            cmd.append("-shortPaired '{}'".format(self.end_1_file))

        if self.long_reads_file and not self.args.no_long_reads:
            cmd.append("-long '{}'".format(self.long_reads_file))

        return ' '.join(cmd)

    def velvetg(self):
        """Build the velvetg for the second assembly step."""

        cmd = ['velvetg',
               '{}'.format(self.work_path),
               '-ins_length {}'.format(self.args.ins_length),
               '-exp_cov {}'.format(self.args.exp_coverage),
               '-min_contig_lgth {}'.format(self.args.min_contig_length)]

        return ' '.join(cmd)

    def post_assembly(self):
        """Copy the assembler output."""

        src = self.iter_file('contigs.fa')
        shutil.move(src, self.output_file)


class SpadesAssembler(Assembler):
    """Wrapper for the Spades assembler."""

    def __init__(self, args):
        super().__init__(args)
        self.steps = [self.spades]

    def spades(self):
        """Build the command for assembly."""
        # spades.py --only-assembler --threads 1 --cov-cutoff 8
        # -1 lib_I9.01.paired_end_1.fasta -2 lib_I9.01.paired_end_2.fasta
        # -o spades

        cmd = ['spades.py --only-assembler']

        return ' '.join(cmd)
