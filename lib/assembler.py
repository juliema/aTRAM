"""Base class for the various assembler programs."""

import sys
from shutil import which
import textwrap
import psutil
import lib.log as log
from .assemblers.abyss import AbyssAssembler
from .assemblers.spades import SpadesAssembler
from .assemblers.trinity import TrinityAssembler
from .assemblers.velvet import VelvetAssembler
from .assemblers.none import NoneAssembler


ASSEMBLERS = {
    'abyss': AbyssAssembler,
    'trinity': TrinityAssembler,
    'velvet': VelvetAssembler,
    'spades': SpadesAssembler,
    'none': NoneAssembler}


def factory(args, cxn):
    """Return the assembler based upon the configuration options."""
    name = args['assembler'].lower()
    assembler = ASSEMBLERS[name]
    return assembler(args, cxn)


def command_line_args(parser):
    """Add command-line arguments for the assemblers."""
    group = parser.add_argument_group('optional assembler arguments')

    group.add_argument(
        '--no-long-reads', action='store_true',
        help="""Do not use long reads during assembly. for the assemblers is
            to use long reads. So this argument will stop the following:
            Abyss: long='LONGREADS' LONGREADS='<file>';
            Trinity: --long_reads '<file>';
            Velvet: -long '<file>'.""")

    group.add_argument(
        '--kmer', type=int, default=64,
        help="""k-mer size.
            Abyss: k=<int> (default 64);
            Velvet: (default 31).""")

    group.add_argument(
        '--mpi', action='store_true',
        help="""Use MPI for this assembler. The assembler must have been
            compiled to use MPI. Abyss: If this is true then pass --cpus <int>
            into np=<int>.""")

    group.add_argument(
        '--abyss-paired-ends', action='store_true',
        help="""Abyss: Normally, we put all of the input sequences in to the
            -se argument. If this is true then we will put paired end sequences
            into the -in argument with any residual single ends into the -se
            argument.""")

    group.add_argument(
        '--bowtie2', action='store_true',
        help="""Use bowtie2 during assembly. Trinity: This will prevent
            --no_bowtie from being passed to Trinity.""")

    total_mem = psutil.virtual_memory().available >> 30
    max_mem = max(1.0, total_mem >> 1)
    group.add_argument(
        '--max-memory', default=max_mem, metavar='MEMORY', type=int,
        help="""Maximum amount of memory to use in gigabytes. We will use {}
            out of {} GB of free/unused memory.
            Trinity: --max_memory <int>G;
            Spades: --memory <int>.""".format(max_mem, total_mem))

    group.add_argument(
        '--exp-coverage', '--expected-coverage', type=int, default=30,
        help="""The expected coverage of the region. Velvet: -exp_cov <int>
            (default 30).""")

    group.add_argument(
        '--ins-length', type=int, default=300,
        help="""The size of the fragments used in the short-read library.
            Velvet: -ins_length <int> (default 300).""")

    group.add_argument(
        '--min-contig-length', type=int, default=100,
        help="""The minimum contig length used by the assembler itself.
            Velvet: -min_contig_lgth <int> (default is 100).""")

    group.add_argument(
        '--careful', action='store_true',
        help="""Tries to reduce the number of mismatches and short indels.
            Spades: --careful.""")

    group.add_argument(
        '--cov-cutoff', default='off',
        help="""Read coverage cutoff value. Must be a positive float value,
            or "auto", or "off". Spades: --cov-cutoff <keyword or int>.""")

    group.add_argument(
        '--abyss-p', type=int,
        help="""Abyss: Minimum sequence identity of a bubble. Enter a single
            digit integer [0-9].""")


def default_kmer(kmer, assembler):
    """Calculate default kmer argument."""
    if assembler == 'velvet' and kmer > 31:
        kmer = 31

    return kmer


def default_cov_cutoff(cov_cutoff):
    """Calculate default coverage cutoff argument."""
    if cov_cutoff in ['off', 'auto']:
        return cov_cutoff

    err = ('Read coverage cutoff value. Must be a positive '
           'float value, or "auto", or "off"')
    try:
        value = float(cov_cutoff)
    except ValueError:
        log.fatal(err)

    if value < 0:
        log.fatal(err)

    return cov_cutoff


def find_program(assembler_name, program, assembler_arg, option=True):
    """Make sure we can find the programs needed by the assembler."""
    if assembler_arg == assembler_name and option and not which(program):
        err = (textwrap.dedent("""
            We could not find the "{}" program. You either need to
            install it or you need to adjust the PATH environment
            variable with the "--path" option so that aTRAM can
            find it.""")).format(program)
        sys.exit(err)
