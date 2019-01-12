"""Setup the test environment."""

import sys
from os.path import dirname, abspath, join


ROOT_DIR = join(dirname(dirname(abspath(__file__))), '.', 'atram')
sys.path.append(ROOT_DIR)
