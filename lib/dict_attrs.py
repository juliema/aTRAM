"""So we can handle dict keys like object attributes."""


class DictAttrs(dict):
    """So we can handle dict keys like object attributes."""
    __getattr__ = dict.__getitem__
