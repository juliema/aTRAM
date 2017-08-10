"""Utility functions."""

import os
from os.path import abspath, join
from shutil import which
import lib.log as log

VERSION = '2.1'


def temp_root_dir(temp_dir_arg, temp_dir_default):
    """Make a temporary root directory. Verify that it is clean."""
    temp_dir = temp_dir_default

    if temp_dir_arg:
        temp_dir = abspath(temp_dir_arg)
        os.makedirs(temp_dir, exist_ok=True)
        if os.listdir(temp_dir):
            log.fatal('"{}" directory is not empty.'.format(temp_dir_arg))

    return temp_dir


def temp_subdir(temp_dir, subdir):
    """Make a subdirectory inside of the root temporary directory."""
    path = join(temp_dir, subdir)
    os.makedirs(path, exist_ok=True)
    return path


def temp_file(temp_dir, subdir, file_name):
    """Create a temporary file in the temporary subdirectory."""
    return join(temp_dir, subdir, file_name)


def iter_dir_name(temp_dir, iteration):
    """Make a directory name from the iteration."""
    return join(temp_dir, 'iteration_{:02d}'.format(iteration))


def temp_iter_dir(temp_dir, iteration):
    """Make a temp directory for the iteration.

    Make sure that the iteration directory exists.
    """
    path = iter_dir_name(temp_dir, iteration)
    os.makedirs(path, exist_ok=True)
    return path


def temp_iter_file(temp_dir, iteration, file_name):
    """Make a temporary file name that is inside of the iteration directory."""
    return join(iter_dir_name(temp_dir, iteration), file_name)


def output_file(output_prefix, file_suffix):
    """Build the output file name."""
    return '{}.{}'.format(output_prefix, file_suffix)


def find_programs(assembler, no_long_reads, bowtie2):
    """Make sure we can find the programs needed by the assembler and blast."""
    if not (which('makeblastdb') and which('tblastn') and which('blastn')):
        err = ('We could not find the programs "makeblastdb", "tblastn", or '
               '"blastn". You either need to install them or you need adjust '
               'the PATH environment variable with the "--path" option so '
               'that aTRAM can find it.')
        log.fatal(err)

    if assembler == 'abyss' and not which('abyss-pe'):
        err = ('We could not find the "abyss-pe" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)

    if assembler == 'abyss' and not no_long_reads and not which('bwa'):
        err = ('We could not find the "bwa-mem" program. You either need to '
               'install it, adjust the PATH environment variable '
               'with the "--path" option, or you may use the '
               '"--no-long-reads" option to not use this program.')
        log.fatal(err)

    if assembler == 'trinity' and not which('Trinity'):
        err = ('We could not find the "Trinity" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)

    if assembler == 'trinity' and bowtie2 and not which('bowtie2'):
        err = ('We could not find the "bowtie2" program. You either need to '
               'install it, adjust the PATH environment variable '
               'with the "--path" option, or you may skip using this program '
               'by not using the "--bowtie2" option.')
        log.fatal(err)

    if assembler == 'velvet' and not (which('velveth') and which('velvetg')):
        err = ('We could not find either the "velveth" or "velvetg" program. '
               'You either need to install it or you need to adjust the PATH '
               'environment variable with the "--path" option so that aTRAM '
               'can find it.')
        log.fatal(err)

    if assembler == 'spades' and not which('spades.py'):
        err = ('We could not find the "Spades" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        log.fatal(err)
