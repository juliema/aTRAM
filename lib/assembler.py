"""Base class for the various assembler programs."""
import math
import psutil
import lib.log as log
from lib.assemblers.abyss import AbyssAssembler
from lib.assemblers.spades import SpadesAssembler
from lib.assemblers.trinity import TrinityAssembler
from lib.assemblers.velvet import VelvetAssembler


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


def command_line_args(parser):
    """Add command-line arguments for the assemblers."""

    # TODO: Move these into the appropriate assembler files and assemble here

    group = parser.add_argument_group('optional assembler arguments')

    group.add_argument('--no-long-reads', action='store_true',
                       help='Do not use long reads during assembly. '
                            '(Abyss, Trinity, Velvet)')

    group.add_argument('--kmer', type=int, default=64,
                       help='k-mer size. The default is "64" for Abyss and '
                            '"31" for Velvet. Note: the maximum kmer length '
                            'for Velvet is 31. (Abyss, Velvet)')

    group.add_argument('--mpi', action='store_true',
                       help='Use MPI for this assembler. The assembler '
                            'must have been compiled to use MPI. (Abyss)')

    group.add_argument('--bowtie2', action='store_true',
                       help='Use bowtie2 during assembly. (Trinity)')

    max_mem = max(1.0, math.floor(
        psutil.virtual_memory().available / 1024**3 / 2))
    group.add_argument('--max-memory',
                       default=max_mem, metavar='MEMORY', type=int,
                       help='Maximum amount of memory to use in gigabytes. '
                            'The default is "{}". (Trinity, Spades)'.format(
                                max_mem))

    group.add_argument('--exp-coverage', '--expected-coverage',
                       type=int, default=30,
                       help='The expected coverage of the region. '
                            'The default is "30". (Velvet)')

    group.add_argument('--ins-length', type=int, default=300,
                       help='The size of the fragments used in the short-read '
                            'library. The default is "300". (Velvet)')

    group.add_argument('--min-contig-length', type=int, default=100,
                       help='The minimum contig length used by the assembler '
                            'itself. The default is "100". (Velvet)')

    group.add_argument('--cov-cutoff', default='off',
                       help='Read coverage cutoff value. Must be a positive '
                            'float value, or "auto", or "off". '
                            'The default is "off". (Spades)')


def check_command_line_args(args):
    """Make sure assembler command-line arguments are reasonable."""

    # Check kmer
    if args.assembler == 'velvet' and args.kmer > 31:
        args.kmer = 31

    # Check cov_cutoff
    if args.cov_cutoff not in ['off', 'auto']:
        err = ('Read coverage cutoff value. Must be a positive '
               'float value, or "auto", or "off"')
        try:
            value = float(args.cov_cutoff)
        except ValueError:
            log.fatal(err)
        if value < 0:
            log.fatal(err)
