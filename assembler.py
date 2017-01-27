"""Wrappers for the various assember programs."""

import os
import shutil
import subprocess
import util


class Assembler:
    """A factory class for building the assembers."""

    def __init__(self, config):
        self.config = config

    @property
    def work_path(self):
        """The output directory name mauy have unique requirements."""
        return self.config['work_dir']

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
    def factory(config):
        """Return the assembler based upon the configuration options."""
        if config['assembler'].lower() == 'trinity':
            return TrinityAssembler(config)
        elif config['assembler'].lower() == 'velvet':
            return VevetAssembler(config)
        elif config['assembler'].lower() == 'abyss':
            return AbyssAssembler(config)


class TrinityAssembler(Assembler):
    """Wrapper for the trinity assembler."""

    @property
    def work_path(self):
        """The output directory name has unique requirements."""
        return os.path.join(self.config['work_dir'], 'trinity')

    def post_assembly(self, iteration, paired):
        """This assember has a unique post assembly step."""
        old_file = os.path.join(self.work_path, 'Trinity.fasta')
        new_file = util.contig_unfiltered_file(self.config, iteration)
        shutil.move(old_file, new_file)  # Save the file for further processing
        shutil.rmtree(self.work_path)    # Remove so other iterations will work

    def command(self, iteration, paired):
        """Build the command for assembly."""

        cmd = ['Trinity']
        cmd.append('--seqType fa')
        cmd.append('--max_memory {}'.format(self.config['max_memory']))
        cmd.append('--CPU {}'.format(self.config['cpus']))
        cmd.append("--output '{}'".format(self.work_path))

        if paired:
            cmd.append("--left '{}'".format(util.paired_end_file(
                self.config, iteration, '1')))
            cmd.append("--right '{}'".format(util.paired_end_file(
                self.config, iteration, '2')))
        else:
            cmd.append("-single '{}'".format(util.paired_end_file(
                self.config, iteration, '1')))
            cmd.append('--run_as_paired')

        return ' '.join(cmd)


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
        cmd.append("name='{}'".format(util.contig_unfiltered_file(
            self.config, iteration)))

        if paired:
            cmd.append("in='{} {}'".format(
                util.paired_end_file(self.config, iteration, '1'),
                util.paired_end_file(self.config, iteration, '2')))
        else:
            cmd.append("se='{}'".format(
                util.paired_end_file(self.config, iteration, '1')))

        return ' '.join(cmd)
