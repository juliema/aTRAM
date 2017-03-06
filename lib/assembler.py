"""Wrappers for the various assember programs."""

import os
import shutil
import subprocess


class Assembler:
    """A factory class for building the assembers."""

    @staticmethod
    def factory(args):
        """Return the assembler based upon the configuration options."""

        if args.assembler.lower() == 'trinity':
            return TrinityAssembler(args)
        elif args.assembler.lower() == 'velvet':
            return VevetAssembler(args)
        elif args.assembler.lower() == 'abyss':
            return AbyssAssembler(args)

    def __init__(self, args):
        self.args = args
        self.is_paired = False
        self.output_file = None
        self.long_reads_file = None
        self.single_end_file = None
        self.paired_end_1_file = None
        self.paired_end_2_file = None

    @property
    def work_path(self):
        """The output directory name may have unique requirements."""

        return self.args.work_dir

    def command(self):
        """Build the command for assembly."""

        raise NotImplementedError()

    def assemble(self):
        """Use the assembler to build up the contigs."""
        cmd = self.command()
        subprocess.check_call(cmd, shell=True)
        self.post_assembly()

    def post_assembly(self):
        """Assembers have unique post assembly steps."""

    def path(self, temp_dir, file_name, iteration=0):
        """Files will go into the temp dir."""

        file_name = '{}.{:02d}.{}'.format(
            self.args.blast_db, iteration, file_name)

        return os.path.join(self.args.work_dir, temp_dir, file_name)

    def iteration_files(self, temp_dir, iteration):
        """Files used by the assembler. Do this at the start of each
        iteration.
        """

        self.output_file = self.path(temp_dir, 'output.fasta', iteration)
        self.single_end_file = self.path(
            temp_dir, 'single_end.fasta', iteration)
        self.paired_end_1_file = self.path(
            temp_dir, 'paired_end_1.fasta', iteration)
        self.paired_end_2_file = self.path(
            temp_dir, 'paired_end_2.fasta', iteration)


class AbyssAssembler(Assembler):
    """Wrapper for the Abyss assembler."""

    def command(self):
        """Build the command for assembly."""

        cmd = ['abyss-pe']
        cmd.append('v=-v')
        cmd.append('E=0')
        cmd.append('k={}'.format(self.args.kmer))
        # cmd.append('np={}'.format(self.args.cpus))
        cmd.append("name='{}'".format(self.output_file))

        if self.is_paired:
            cmd.append("in='{} {}'".format(self.paired_end_1_file,
                                           self.paired_end_2_file))
        else:
            cmd.append("se='{}'".format(self.paired_end_1_file))

        if self.long_reads_file:
            cmd.append("long='{}'".format(self.long_reads_file))

        return ' '.join(cmd)

    def post_assembly(self):
        """This assember has a unique post assembly step."""

        src = os.path.realpath(self.output_file + '-unitigs.fa')
        dst = self.output_file

        # shutil.move(src['raw_contigs'].name)
        with open(src) as in_file, open(dst, 'w') as out_file:
            for line in in_file:
                out_file.write(line)


class TrinityAssembler(Assembler):
    """Wrapper for the trinity assembler."""

    @property
    def work_path(self):
        """The output directory name has unique requirements."""

        return 'trinity'

    def command(self):
        """Build the command for assembly."""

        cmd = ['Trinity']
        cmd.append('--seqType fa')
        cmd.append('--max_memory {}'.format(self.args['max_memory']))
        cmd.append('--CPU {}'.format(self.args['cpus']))
        cmd.append("--output '{}'".format(self.work_path))
        cmd.append('--full_cleanup')

        if self.is_paired:
            cmd.append("--left '{}'".format(self.paired_end_1_file))
            cmd.append("--right '{}'".format(self.paired_end_2_file))
        else:
            cmd.append("-single '{}'".format(self.paired_end_1_file))
            cmd.append('--run_as_paired')

        if self.long_reads_file:
            cmd.append("--long_reads '{}'".format(self.long_reads_file))

        return ' '.join(cmd)

    def post_assembly(self):
        """This assember has a unique post assembly step."""

        shutil.move('trinity.Trinity.fasta', self.output_file)


class VevetAssembler(Assembler):
    """Wrapper for the Velvet assembler."""

    def command(self):
        """Build the command for assembly."""
