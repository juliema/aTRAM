"""Wrappers for the various assember programs."""

import os
import shutil
import subprocess
from lib.filer import Filer


class Assembler:
    """A factory class for building the assembers."""

    def __init__(self, config):
        self.config = config
        self.filer = Filer(work_dir=config.work_dir,
                           db_prefix=config.db_prefix)

    @property
    def work_path(self):
        """The output directory name may have unique requirements."""

        return self.config.work_dir

    def command(self, files):
        """Build the command for assembly."""

        raise NotImplementedError()

    def assemble(self, files):
        """Use the assembler to build up the contigs."""
        cmd = self.command(files)
        subprocess.check_call(cmd, shell=True)
        self.post_assembly(files)

    def post_assembly(self, files):
        """Assembers have unique post assembly steps."""

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

        return 'trinity'

    def command(self, files):
        """Build the command for assembly."""

        cmd = ['Trinity']
        cmd.append('--seqType fa')
        cmd.append('--max_memory {}'.format(self.config['max_memory']))
        cmd.append('--CPU {}'.format(self.config['cpus']))
        cmd.append("--output '{}'".format(self.work_path))
        cmd.append('--full_cleanup')

        if files['is_paired']:
            cmd.append("--left '{}'".format(files['end_1'].name))
            cmd.append("--right '{}'".format(files['end_2'].name))
        else:
            cmd.append("-single '{}'".format(files['end_1'].name))
            cmd.append('--run_as_paired')

        return ' '.join(cmd)

    def post_assembly(self, files):
        """This assember has a unique post assembly step."""

        shutil.move('trinity.Trinity.fasta', files['raw_contigs'].name)


class VevetAssembler(Assembler):
    """Wrapper for the Velvet assembler."""

    def command(self, files):
        """Build the command for assembly."""


class AbyssAssembler(Assembler):
    """Wrapper for the Abyss assembler."""

    def command(self, files):
        """Build the command for assembly."""

        cmd = ['abyss-pe']
        cmd.append('v=-v')
        cmd.append('E=0')
        cmd.append('k={}'.format(self.config['kmer']))
        # cmd.append('np={}'.format(self.config['cpus']))
        cmd.append("name='{}'".format(files['raw_contigs'].name))

        if files['is_paired']:
            cmd.append("in='{} {}'".format(files['end_1'].name,
                                           files['end_2'].name))
        else:
            cmd.append("se='{}'".format(files['end_1'].name))

        return ' '.join(cmd)

    def post_assembly(self, files):
        """This assember has a unique post assembly step."""

        src = os.path.realpath(files['raw_contigs'].name + '-unitigs.fa')
        dst = files['raw_contigs'].name

        # shutil.move(src, files['raw_contigs'].name)
        with open(src) as in_file, open(dst, 'w') as out_file:
            for line in in_file:
                out_file.write(line)
