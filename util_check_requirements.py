#!/usr/bin/env python3
"""Check that we have the minimum requirements for running atram."""

import re
import sys
import subprocess
from functools import reduce
from distutils.version import LooseVersion
from shutil import which


RESULTS = {}


def test_format(name, value):
    """Format test results."""
    RESULTS[name] = value
    value = 'OK' if value else 'FAIL'
    print(name.ljust(40, '.'), value)


def parse_requirements(requirements):
    """Parse a requirement into a module and version parts."""
    reqs = {}
    for req in requirements.split():
        match = re.match(r'^([^>=<]+)([>=<]+)([^>=<]+)$', req)
        module = match.group(1)
        compare = match.group(2)
        version = LooseVersion(match.group(3)).version
        reqs[module] = {'compare': compare, 'version': version}
    return reqs


def check_modules():
    """Get installed python modules."""
    modules = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
    installed_list = parse_requirements(modules.decode('utf-8'))

    with open('requirements.txt') as requirements:
        required_list = parse_requirements(requirements.read())

    for module, required in required_list.items():
        installed = installed_list[module]

        cmp = required['compare']
        i_version = installed['version']
        r_version = required['version']

        if cmp == '==' and i_version != r_version:
            test_format(module, False)
        elif cmp == '>=' and i_version > r_version:
            test_format(module, True)
        elif i_version < r_version:
            test_format(module, False)


def check_programs():
    """Verify that external programs can be found."""
    test_format('makeblastdb', which('makeblastdb'))
    test_format('tblastn', which('tblastn'))
    test_format('blastn', which('blastn'))
    test_format('abyss', which('abyss-pe'))
    test_format('trinity', which('Trinity'))
    test_format('velvet', which('velveth') and which('velvetg'))
    test_format('spades', which('spades.py'))
    test_format('bwa', which('bwa'))
    test_format('bowtie2', which('bowtie2'))
    test_format('exonerate', which('exonerate'))


def requires(module, because, program=None):
    """Show that aTRAM will not work without the given program/module."""
    if not program:
        program = 'atram.py and atram_preprocessor.py'
    if not RESULTS[module]:
        print('  {} will not work because {}'.format(program, because))


def assembler(module, but):
    """Report limited aTRAM functionality."""
    if not RESULTS[module]:
        print('  atram.py will work but {}'.format(but))


def report_results():
    """Show the user what they can and cannot do."""
    print('\nResults:')

    if reduce(lambda a, b: a and b, RESULTS.values()):
        print('  aTRAM is ready and all features are available.\n')
        return

    requires('Python version', 'we need Python version 3.6 or above')
    requires('makeblastdb', "BLAST's makeblastdb is not installed")
    requires('tblastn', "BLAST's tblastn is not installed")
    requires('blastn', "BLAST's blastn is not installed")
    requires('biopython', 'the biopython module is missing (install with pip)')
    requires('psutil', 'the psutil module is missing (install it with pip)')
    requires('numpy', 'the numpy module is missing (install it with pip)')
    requires('exonerate', 'exonerate is not installed', 'atram_stitcher.py')

    assembler('abyss', 'you are missing the abyss assembler')
    assembler('trinity', 'you are missing the trinity assembler')
    assembler('velvet', 'you are missing the velvet assembler')
    assembler('spades', 'you are missing the spades assembler')
    assembler(
        'bowtie2',
        'you are missing bowtie2 and cannot use it with the trinity assembler')
    assembler(
        'bwa',
        ('you are missing bwa and will not be able to use the assemble '
         'long reads with the abyss assembler'))


if __name__ == '__main__':
    test_format(
        'Python version',
        sys.version_info.major == 3 and sys.version_info.minor >= 6)

    check_programs()
    check_modules()
    report_results()
