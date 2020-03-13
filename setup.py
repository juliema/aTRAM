#!/usr/bin/env python
"""Setup the aTRAM environment."""

# Test pip
# 1) Clean the /dist directory
# 2) python3 setup.py sdist bdist_wheel
# 3) pip install --index-url https://test.pypi.org/simple/
#    --extra-index-url https://pypi.org/simple atram
# 4) twine upload --repository-url https://test.pypi.org/legacy/ dist/*

import re
from setuptools import setup, find_packages


def readme():
    """Get README.md content."""
    with open("README.md", 'r') as f:
        return f.read()


def license_():
    """Get LICENSE.txt content."""
    with open("LICENSE.txt", 'r') as f:
        return f.read()


def find_version():
    """Read version from db.py."""
    regex = r"^ATRAM_VERSION = ['\"]v?([^'\"]*)['\"]"
    with open("./lib/db.py", 'r') as f:
        match = re.search(regex, f.read(), re.M)
        if match:
            return match.group(1)

    raise RuntimeError("Unable to find version string.")


def find_requirements():
    """Read requirements.txt file and returns list of requirements."""
    with open("requirements.txt", 'r') as f:
        return f.read().splitlines()


setup(
    name="atram",
    version=find_version(),
    packages=find_packages(),
    install_requires=find_requirements(),
    description="""atram ("automated target restricted assembly method") is
                   an iterative assembler that performs reference-guided
                   local de novo assemblies using a variety of available
                   methods""",
    long_description=readme(),
    license=license_(),
    url="https://github.com/juliema/aTRAM",
    python_requires='>=3.4',
    scripts=[
        'atram.py',
        'atram_preprocessor.py',
        'atram_stitcher.py',
        'atram_framer.py',
        ])
