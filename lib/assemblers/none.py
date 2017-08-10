"""Null object for the assemblers."""

import lib.db as db
import lib.file_util as file_util
from lib.assemblers.base import BaseAssembler


class NoneAssembler(BaseAssembler):
    """Null object for the assemblers."""

    def __init__(self, args, db_conn):
        """Build the assembler."""
        super().__init__(args)
        self.steps = []

    @property
    def blast_only(self):
        """Use this to flag if the assembler should exit early."""
        return True

    def write_final_output(self, blast_db, query):
        """Output this file if we are not assembling the contigs."""
        prefix = self.final_output_prefix(blast_db, query)

        file_name = file_util.output_file(prefix, 'blast_only.fasta')

        with open(file_name, 'w') as output_file:
            for row in db.get_sra_blast_hits(self.db_conn, 1):
                output_file.write('>{}{}\n'.format(
                    row['seq_name'], row['seq_end']))
                output_file.write('{}\n'.format(row['seq']))
