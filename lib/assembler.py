"""Wrappers for the various assember programs."""

import os
import shutil
import subprocess
from lib.filer import Filer


class Assembler:
    """A factory class for building the assembers."""

    def __init__(self, config, iteration):
        self.config = config
        self.filer = Filer(work_dir=config.work_dir,
                           file_prefix=config.file_prefix,
                           iteration=iteration)

    @property
    def work_path(self):
        """The output directory name may have unique requirements."""
        return self.config.work_dir

    def command(self, iteration, paired):
        """Build the command for assembly."""

    def assemble(self, iteration, paired):
        """Use the assembler to build up the contigs."""
        cmd = self.command(iteration, paired)
        subprocess.check_call(cmd, shell=True)
        self.post_assembly(iteration, paired)

    def post_assembly(self, iteration, paired):
        """Some assembers have unique post assembly steps."""

    @staticmethod
    def factory(config, iteration):
        """Return the assembler based upon the configuration options."""
        if config['assembler'].lower() == 'trinity':
            return TrinityAssembler(config, iteration)
        elif config['assembler'].lower() == 'velvet':
            return VevetAssembler(config, iteration)
        elif config['assembler'].lower() == 'abyss':
            return AbyssAssembler(config, iteration)


class TrinityAssembler(Assembler):
    """Wrapper for the trinity assembler."""

    @property
    def work_path(self):
        """The output directory name has unique requirements."""

        return os.path.join(self.config['work_dir'], 'trinity')

    def command(self, iteration, paired):
        """Build the command for assembly."""

        cmd = ['Trinity']
        cmd.append('--seqType fa')
        cmd.append('--max_memory {}'.format(self.config['max_memory']))
        cmd.append('--CPU {}'.format(self.config['cpus']))
        cmd.append("--output '{}'".format(self.work_path))

        if paired:
            cmd.append("--left '{}'".format(self.filer.paired_end_file('1')))
            cmd.append("--right '{}'".format(self.filer.paired_end_file('2')))
        else:
            cmd.append("-single '{}'".format(self.filer.paired_end_file('1')))
            cmd.append('--run_as_paired')

        return ' '.join(cmd)

    def post_assembly(self, iteration, paired):
        """This assember has a unique post assembly step."""
        old_file = os.path.join(self.work_path, 'Trinity.fasta')
        new_file = self.filer.contig_unfiltered_file()
        shutil.move(old_file, new_file)  # Save the file for further processing
        shutil.rmtree(self.work_path)    # Remove so other iterations will work


class VevetAssembler(Assembler):
    """Wrapper for the Velvet assembler."""

    def command(self, iteration, paired):
        """Build the command for assembly."""


class AbyssAssembler(Assembler):
    """Wrapper for the Abyss assembler."""

    def command(self, iteration, paired):
        """Build the command for assembly."""

        cmd = ['abyss-pe']
        cmd.append('v=-v')
        cmd.append('E=0')
        cmd.append('k={}'.format(self.config['kmer']))
        # cmd.append('np={}'.format(self.config['cpus']))
        cmd.append("name='{}'".format(self.filer.contig_unfiltered_file()))

        if paired:
            cmd.append("in='{} {}'".format(
                self.filer.paired_end_file('1'),
                self.filer.paired_end_file('2')))
        else:
            cmd.append("se='{}'".format(
                self.filer.paired_end_file('1')))

        return ' '.join(cmd)
