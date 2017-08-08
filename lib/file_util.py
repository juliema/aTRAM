"""Utility functions."""

import os
from shutil import which
import lib.log as log


def temp_root_dir(temp_dir_arg, temp_dir_default):
    """Make a temporary root directory. Verify that it is clean."""

    temp_dir = temp_dir_default

    if temp_dir_arg:
        temp_dir = os.path.abspath(temp_dir_arg)
        os.makedirs(temp_dir, exist_ok=True)
        if os.listdir(temp_dir):
            log.fatal('"{}" directory is not empty.'.format(temp_dir_arg))

    return temp_dir


def temp_iter_dir(temp_dir, iteration):
    """Make a temp directory for the iteration underneath the root temporary
    directory. Make sure that the iteration directory exists.
    """

    iter_name = 'iteration_{:02d}'.format(iteration)
    iter_dir = os.path.join(temp_dir, iter_name)
    os.makedirs(iter_dir, exist_ok=True)
    return iter_dir


def temp_iter_file(temp_dir, iteration, file_name):
    """Make a temporary file name that is inside of the iteration temporary
    directory.
    """

    iter_dir = temp_iter_dir(temp_dir, iteration)
    return os.path.join(iter_dir, file_name)


def output_file(output_prefix, file_suffix):
    """Build the output file name."""

    return '{}.{}'.format(output_prefix, file_suffix)


def find_programs(args):
    """Make sure we can find the programs needed by the assembler and blast."""

    if not (which('makeblastdb') and which('tblastn') and which('blastn')):
        err = ('We could not find the programs "makeblastdb", "tblastn", or '
               '"blastn". You either need to install them or you need adjust '
               'the PATH environment variable with the "--path" option so '
               'that aTRAM can find it.')
        log.fatal(err)

    if args.assembler == 'abyss' and not which('abyss-pe'):
        err = ('We could not find the "abyss-pe" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)

    if args.assembler == 'abyss' and not args.no_long_reads \
            and not which('bwa'):
        err = ('We could not find the "bwa-mem" program. You either need to '
               'install it, adjust the PATH environment variable '
               'with the "--path" option, or you may use the '
               '"--no-long-reads" option to not use this program.')
        log.fatal(err)

    if args.assembler == 'trinity' and not which('Trinity'):
        err = ('We could not find the "Trinity" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)

    if args.assembler == 'trinity' and args.bowtie2 and not which('bowtie2'):
        err = ('We could not find the "bowtie2" program. You either need to '
               'install it, adjust the PATH environment variable '
               'with the "--path" option, or you may skip using this program '
               'by not using the "--bowtie2" option.')
        log.fatal(err)

    if args.assembler == 'velvet' and \
            not (which('velveth') and which('velvetg')):
        err = ('We could not find either the "velveth" or "velvetg" program. '
               'You either need to install it or you need to adjust the PATH '
               'environment variable with the "--path" option so that aTRAM '
               'can find it.')
        log.fatal(err)

    if args.assembler == 'spades' and not which('spades.py'):
        err = ('We could not find the "Spades" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)
