"""Base class for the various assembler programs."""

import sys
import textwrap
from shutil import which

from .assemblers.abyss import AbyssAssembler
from .assemblers.none import NoneAssembler
from .assemblers.spades import SpadesAssembler
from .assemblers.trinity import TrinityAssembler
from .assemblers.velvet import VelvetAssembler

ASSEMBLERS = {
    'abyss': AbyssAssembler,
    'trinity': TrinityAssembler,
    'velvet': VelvetAssembler,
    'spades': SpadesAssembler,
    'none': NoneAssembler}


def factory(args, cxn, log):
    """Return the assembler based upon the configuration options."""
    name = args['assembler'].lower()
    assembler = ASSEMBLERS[name]
    return assembler(args, cxn, log)


def command_line_args(parser):
    """Add command-line arguments for the assemblers."""
    AbyssAssembler.command_line_args(parser)
    SpadesAssembler.command_line_args(parser)
    TrinityAssembler.command_line_args(parser)
    VelvetAssembler.command_line_args(parser)


def find_program(assembler_name, program, assembler_arg, option=True):
    """Make sure we can find the programs needed by the assembler."""
    if assembler_arg == assembler_name and option and not which(program):
        err = (textwrap.dedent("""
            We could not find the "{}" program. You either need to
            install it or you need to adjust the PATH environment
            variable with the "--path" option so that aTRAM can
            find it.""")).format(program)
        sys.exit(err)
