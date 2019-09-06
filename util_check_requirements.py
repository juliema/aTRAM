"""Check that we have the minimum requirements for running atram."""

import re
import sys
import subprocess
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
        version = match.group(3).split('.')
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

        r_cmp = required['compare']

        for r_part, i_part in zip(required['version'], installed['version']):
            if r_cmp == '==' and i_part != r_part:
                test_format(module, False)
                break
            elif r_cmp == '>=' and i_part > r_part:
                test_format(module, True)
                break
            elif r_part > i_part:
                test_format(module, False)
            break
        else:
            test_format(module, True)


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


if __name__ == '__main__':
    test_format(
        'Python version',
        sys.version_info.major == 3 and sys.version_info.minor >= 5)

    check_programs()
    check_modules()
