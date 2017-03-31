"""Wrappers for the various assember programs."""

import os
import shutil
import logging
import subprocess


class Assembler:  # pylint: disable=too-many-instance-attributes
    """A factory class for building the assembers."""

    @staticmethod
    def factory(args, temp_dir):
        """Return the assembler based upon the configuration options."""

        if args.assembler.lower() == 'abyss':
            return AbyssAssembler(args, temp_dir)
        elif args.assembler.lower() == 'trinity':
            return TrinityAssembler(args, temp_dir)
        elif args.assembler.lower() == 'velvet':
            return VelvetAssembler(args, temp_dir)

    def __init__(self, args, temp_dir):
        self.args = args            # Parsed command line arguments
        self.steps = []             # Assembler steps setup by the assembler
        self.temp_dir = temp_dir    # Temp directory used for storing files
        self.is_paired = False      # Did we find paired end sequences?
        self.output_file = None     # Write to this file
        self.long_reads_file = None  # Long-reads file
        self.ends_1_file = None     # Sequneces for end 1 or single end reads
        self.ends_2_file = None     # Sequences for end 2 reads
        self.cwd = None             # Some assemblers need a directory change

    @property
    def work_path(self):
        """The output directory name may have unique requirements."""

        return self.args.work_dir

    def assemble(self):
        """Use the assembler to build up the contigs. We take and array of
        subprocess steps and execute them in order. We bracket this with
        pre and post assembly steps.
        """

        try:

            self.pre_assembly()

            for step in self.steps:
                cmd = step()
                logging.info(cmd)
                subprocess.check_call(cmd, shell=True)

            self.post_assembly()

        except Exception as exn:  # pylint: disable=broad-except
            raise exn
        finally:
            if self.cwd:
                os.chdir(self.cwd)

    def pre_assembly(self):
        """Assembers have unique pre assembly steps."""

    def post_assembly(self):
        """Assembers have unique post assembly steps."""

    def path(self, file_name, iteration=0):
        """Files will go into the temp dir."""

        file_name = '{}.{:02d}.{}'.format(
            self.args.blast_db, iteration, file_name)

        return os.path.join(self.temp_dir, file_name)

    def iteration_files(self, iteration):
        """Files used by the assembler. Do this at the start of each
        iteration.
        """

        self.output_file = self.path('output.fasta', iteration)
        self.ends_1_file = self.path('paired_end_1.fasta', iteration)
        self.ends_2_file = self.path('paired_end_2.fasta', iteration)


class AbyssAssembler(Assembler):
    """Wrapper for the Abyss assembler."""

    def __init__(self, args, temp_dir):
        super().__init__(args, temp_dir)
        self.steps = [self.abyss]

    def pre_assembly(self):
        """Abyss needs to work from the current directory which will be the
        temporary directory here.
        """

        self.cwd = os.getcwd()
        os.chdir(self.temp_dir)

    def abyss(self):
        """Build the command for assembly."""

        cmd = ['abyss-pe']
        cmd.append('v=-v')
        cmd.append('E=0')
        cmd.append('k={}'.format(self.args.kmer))
        # cmd.append('np={}'.format(self.args.cpus))
        cmd.append("name={}".format(os.path.basename(self.output_file)))

        if self.is_paired:
            cmd.append("in='{} {}'".format(os.path.basename(self.ends_1_file),
                                           os.path.basename(self.ends_2_file)))
        else:
            cmd.append("se={}".format(os.path.basename(self.ends_1_file)))

        if self.long_reads_file and not self.args.no_long_reads:
            cmd.append("long={}".format(
                os.path.basename(self.long_reads_file)))

        return ' '.join(cmd)

    def post_assembly(self):
        """This assember has unique post assembly steps."""

        os.chdir(self.cwd)

        src = os.path.realpath(self.output_file + '-unitigs.fa')

        shutil.copyfile(src, self.output_file)


class TrinityAssembler(Assembler):
    """Wrapper for the trinity assembler."""

    @property
    def work_path(self):
        """The output directory name has unique requirements."""

        return os.path.join(self.args.work_dir, 'trinity')

    def __init__(self, args, temp_dir):
        super().__init__(args, temp_dir)
        self.steps = [self.trinity]

    def trinity(self):
        """Build the command for assembly."""

        cmd = ['Trinity']
        cmd.append('--seqType fa')
        cmd.append('--max_memory {}'.format(self.args.max_memory))
        cmd.append('--CPU {}'.format(self.args.cpus))
        cmd.append("--output '{}'".format(self.work_path))
        cmd.append('--full_cleanup')

        if not self.args.bowtie2:
            cmd.append('--no_bowtie')

        if self.is_paired:
            cmd.append("--left '{}'".format(self.ends_1_file))
            cmd.append("--right '{}'".format(self.ends_2_file))
        else:
            cmd.append("-single '{}'".format(self.ends_1_file))
            cmd.append('--run_as_paired')

        if self.long_reads_file and not self.args.no_long_reads:
            cmd.append("--long_reads '{}'".format(self.long_reads_file))

        return ' '.join(cmd)

    def post_assembly(self):
        """This assember has a unique post assembly step."""

        src = os.path.join(self.args.work_dir, 'trinity.Trinity.fasta')
        shutil.move(src, self.output_file)


class VelvetAssembler(Assembler):
    """Wrapper for the Velvet assembler."""

    def __init__(self, args, temp_dir):
        super().__init__(args, temp_dir)
        self.steps = [self.velveth, self.velvetg]

    def velveth(self):
        """Build the command for assembly."""

        cmd = ['velveth']
        cmd.append('{}'.format(self.temp_dir))
        cmd.append('{}'.format(self.args.kmer))
        cmd.append('-fasta')

        if self.is_paired:
            cmd.append("-shortPaired '{}' '{}'".format(
                self.ends_1_file, self.ends_2_file))
        else:
            cmd.append("-shortPaired '{}'".format(self.ends_1_file))

        if self.long_reads_file and not self.args.no_long_reads:
            cmd.append("-long '{}'".format(self.long_reads_file))

        return ' '.join(cmd)

    def velvetg(self):
        """Build the command for assembly."""

        cmd = ['velvetg']
        cmd.append('{}'.format(self.temp_dir))
        cmd.append('-ins_length {}'.format(self.args.ins_length))
        cmd.append('-exp_cov {}'.format(self.args.exp_coverage))
        cmd.append('-min_contig_lgth {}'.format(self.args.min_contig_length))

        return ' '.join(cmd)

    def post_assembly(self):
        """This assember has a unique post assembly step."""

        src = os.path.join(self.temp_dir, 'contigs.fa')
        shutil.move(src, self.output_file)
