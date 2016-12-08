"""So we can handle dict keys like object attributes."""


class DictAttrs(dict):
    """Allow dot notation for accessing dict items."""
    __getattr__ = dict.__getitem__
