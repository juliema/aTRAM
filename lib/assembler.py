"""Base class for the various assembler programs."""

import math
import psutil
import lib.log as log
from lib.assemblers.abyss import AbyssAssembler
from lib.assemblers.spades import SpadesAssembler
from lib.assemblers.trinity import TrinityAssembler
from lib.assemblers.velvet import VelvetAssembler
from lib.assemblers.none import NoneAssembler


def factory(args, db_conn):
    """Return the assembler based upon the configuration options."""
    name = args['assembler'].lower()
    if name == 'abyss':
        return AbyssAssembler(args, db_conn)
    elif name == 'trinity':
        return TrinityAssembler(args, db_conn)
    elif name == 'velvet':
        return VelvetAssembler(args, db_conn)
    elif name == 'spades':
        return SpadesAssembler(args, db_conn)
    elif name == 'none':
        return NoneAssembler(args, db_conn)


def command_line_args(parser):
    """Add command-line arguments for the assemblers."""
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


def default_kmer(kmer, assembler):
    """Calculate default kmer argument."""
    if assembler == 'velvet' and kmer > 31:
        kmer = 31

    return kmer


def default_cov_cutoff(cov_cutoff):
    """Calculate default coverage cutoff argument."""
    if cov_cutoff not in ['off', 'auto']:
        err = ('Read coverage cutoff value. Must be a positive '
               'float value, or "auto", or "off"')
        try:
            value = float(cov_cutoff)
        except ValueError:
            log.fatal(err)
        if value < 0:
            log.fatal(err)

    return cov_cutoff
