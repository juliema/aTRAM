"""Base class for the various assembler programs."""

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
