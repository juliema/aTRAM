#!/usr/bin/env python

import re
from setuptools import setup, find_packages


def readme():
    """Get README.md content."""
    with open("README.md", 'r') as f:
        return f.read()
    return None


def license():
    """Get LICENSE.txt content."""
    with open("LICENSE.txt", 'r') as f:
        return f.read()
    return None


def find_version():
    """Read version from db.py."""
    rex = r"^ATRAM_VERSION = ['\"]([^'\"]*)['\"]"
    with open("./lib/db.py", 'r') as f:
        match = re.search(rex, f.read(), re.M)
        if match:
            return match.group(1)

    raise RuntimeError("Unable to find version string.")


def find_requirements():
    """Read requirements.txt file and returns list of requirements."""
    with open("requirements.txt", 'r') as f:
        return f.read().splitlines()

    raise RuntimeError("Unable to list requirements.")


setup(
    name="aTRAM",
    version=find_version(),
    packages=find_packages(),
    install_requires=find_requirements(),
    description="""aTRAM ("automated target restricted assembly method") is
                   an iterative assembler that performs reference-guided
                   local de novo assemblies using a variety of available
                   methods""",
    long_description=readme(),
    license=license(),
    url="https://github.com/juliema/aTRAM",
    python_requires='>=3.5',
    scripts=['atram/atram.py', 'atram/atram_preprocessor.py'])
