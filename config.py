"""Handle configuration options."""

import configparser

DEFAULTS = {}


class DictAttrs(dict):
    """Allow dot notation for accessing dict items."""
    __getattr__ = dict.__getitem__


def read_config_file(config, args):
    """Read in the config file and hoist the data into the config dict."""
    if not args.get('config_file', None):
        pass
    parser = configparser.ConfigParser()
    parser.read(args['config_file'])
    return config


def get_configs(args=None):
    if not args:
        args = {}
    config = DEFAULTS

    DictAttrs(config)
# Setup defaults
# Read the config file
# Overwrite command-line args
