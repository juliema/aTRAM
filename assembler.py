"""Wrappers for the various assember programs."""


# pylint: disable=R0903
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
# pylint: enable=R0903


class TrinityAssembler:
    """Wrapper for the trinity assembler."""

    def __init__(self, config):
        self.config = config

    def command(self, iteration):
        """Build the command for assembly."""
        fasta_file = '{}matching_seqs_{}.fasta'.format(self.config['blast_db'], iteration)
        # trinity --seqType fa --single $short_read_file --run_as_paired --JM $jm --output $tempdir
        cmd = 'Trinity --seqType fa --single {} --run_as_paired --output $tempdir'
        cmd = cmd.format(fasta_file)

    def assemble(self, iteration):
        """Use the assembler to build up the contigs."""
        cmd = self.command(iteration)
        print(cmd)


class VevetAssembler:
    """Wrapper for the Velvet assembler."""

    def __init__(self, config):
        self.config = config

    def command(self, iteration):
        """Build the command for assembly."""

    def assemble(self, iteration):
        """Use the assembler to build up the contigs."""


class AbyssAssembler:
    """Wrapper for the Abyss assembler."""

    def __init__(self, config):
        self.config = config

    def command(self, iteration):
        """Build the command for assembly."""

    def assemble(self, iteration):
        """Use the assembler to build up the contigs."""
