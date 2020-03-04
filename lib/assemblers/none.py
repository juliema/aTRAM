"""Null object for the assemblers."""

from .. import db_atram
from .. import util
from .base import BaseAssembler


class NoneAssembler(BaseAssembler):
    """Null object for the assemblers."""

    def __init__(self, args, cxn):
        """Build the assembler."""
        super().__init__(args, cxn)
        self.steps = []
        self.blast_only = True  # Used to short-circuit the assembler

    def write_final_output(self, blast_db, query):
        """Output this file if we are not assembling the contigs."""
        prefix = self.final_output_prefix(blast_db, query)

        file_name = '{}.fasta'.format(prefix)

        with open(file_name, 'w') as output_file:
            for row in db_atram.get_sra_blast_hits(self.state['cxn'], 1):
                util.write_fasta_record(
                    output_file, row['seq_name'], row['seq'], row['seq_end'])
