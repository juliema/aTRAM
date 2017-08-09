"""Null object for the assemblers."""

from lib.assemblers.base import BaseAssembler


class NoneAssembler(BaseAssembler):
    """Null object for the assemblers."""

    def __init__(self, args):
        super().__init__(args)
        self.steps = []
