"""Null object for the assemblers."""

import lib.db as db
from lib.assemblers.base import BaseAssembler


class NoneAssembler(BaseAssembler):
    """Null object for the assemblers."""

    def __init__(self, args, db_conn):
        """Build the assembler."""
        super().__init__(args, db_conn)
        self.steps = []
        self.blast_only = True  # Used to short-circuit the assembler

    def write_final_output(self, blast_db, query):
        """Output this file if we are not assembling the contigs."""
        prefix = self.final_output_prefix(blast_db, query)

        file_name = '{}.fasta'.format(prefix)

        with open(file_name, 'w') as output_file:
            for row in db.get_sra_blast_hits(self.state['db_conn'], 1):
                output_file.write('>{}{}\n'.format(
                    row['seq_name'], row['seq_end']))
                output_file.write('{}\n'.format(row['seq']))
