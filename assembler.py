"""Wrappers for the various assember programs."""

import util


# pylint: disable=too-few-public-methods
class Assembler:
    """A factory class for building the assembers."""
    @staticmethod
    def factory(config):
        """Return the assembler based upon the configuration options."""
        if config['assembler'].lower() == 'trinity':
            return TrinityAssembler(config)
        elif config['assembler'].lower() == 'velvet':
            return VevetAssembler(config)
        elif config['assembler'].lower() == 'abyss':
            return AbyssAssembler(config)
# pylint: enable=too-few-public-methods


class TrinityAssembler:
    """Wrapper for the trinity assembler."""

    def __init__(self, config):
        self.config = config

    def output_file(self):
        """The output file name has unique requirements."""
        pass

    def command(self, iteration, paired):
        """Build the command for assembly."""

        cmd = ['Trinity']
        cmd.append('--seqType fa')
        cmd.append('--max_memory {}'.format(self.config['max_memory']))
        cmd.append('--CPU {}'.format(self.config['cpu']))

        if not paired:
            cmd.append('--left {}'.format(util.paired_end_file(self.config, iteration, '1')))
            cmd.append('--right {}'.format(util.paired_end_file(self.config, iteration, '2')))
        else:
            cmd.append('--single {}'.format(util.paired_end_file(self.config, iteration, '1')))
            cmd.append('--run_as_paired')

        return ' '.join(cmd)

    def assemble(self, iteration, paired):
        """Use the assembler to build up the contigs."""
        cmd = self.command(iteration, paired)
        print(cmd)


class VevetAssembler:
    """Wrapper for the Velvet assembler."""

    def __init__(self, config):
        self.config = config

    def command(self, iteration, paired):
        """Build the command for assembly."""

    def assemble(self, iteration, paired):
        """Use the assembler to build up the contigs."""


class AbyssAssembler:
    """Wrapper for the Abyss assembler."""

    def __init__(self, config):
        self.config = config

    def command(self, iteration, paired):
        """Build the command for assembly."""

    def assemble(self, iteration, paired):
        """Use the assembler to build up the contigs."""
