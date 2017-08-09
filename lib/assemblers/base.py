"""Base class for the various assembler wrappers."""

import os
import lib.db as db
import lib.log as log
import lib.file_util as file_util


class BaseAssembler:
    """A base class for the assemblers."""

    def __init__(self, args):
        self.args = args    # Parsed command line arguments
        self.steps = []     # Assembler steps setup by the assembler
        self.file = {}      # Files and record counts
        self.iteration = 0  # Current iteration

    @property
    def iter_dir(self):
        """Get the work directory for the current iteration."""

        return file_util.temp_iter_dir(self.args['temp_dir'], self.iteration)

    @property
    def work_path(self):
        """The output directory name may have unique requirements."""

        return self.iter_dir

    def iter_file(self, file_name):
        """Build a temporary file name honoring the current iteration
        directory.
        """

        return os.path.join(self.iter_dir, file_name)

    def assemble(self):
        """Use the assembler to build up the contigs. We take and array of
        subprocess steps and execute them in order. We bracket this with
        pre and post assembly steps.
        """

        for step in self.steps:
            cmd = step()
            log.subcommand(cmd, self.args['temp_dir'], self.args['timeout'])

        self.post_assembly()

    def post_assembly(self):
        """Assemblers have unique post assembly steps."""

    def path(self, file_name):
        """Files will go into the temp dir."""

        blast_db = os.path.basename(self.args['blast_db'])
        file_name = '{}.{:02d}.{}'.format(blast_db, self.iteration, file_name)
        rel_path = self.iter_file(file_name)

        return os.path.abspath(rel_path)

    def initialize_iteration(self, iteration):
        """Files used by the assembler. Do this at the start of each
        iteration.
        """

        self.iteration = iteration

        self.file['output'] = self.path('output.fasta')
        self.file['long_reads'] = ''  # Set up in atram.py for now
        self.file['paired_1'] = self.path('paired_end_1.fasta')
        self.file['paired_2'] = self.path('paired_end_2.fasta')
        self.file['single_1'] = self.path('single_end_1.fasta')
        self.file['single_2'] = self.path('single_end_2.fasta')
        self.file['single_any'] = self.path('single_end_any.fasta')

        self.file['paired_count'] = 0       # paired_1 and paired_2 count
        self.file['single_1_count'] = 0
        self.file['single_2_count'] = 0
        self.file['single_any_count'] = 0

    @staticmethod
    def parse_contig_id(header):
        """Given a fasta header line built by the assembler return the
        contig ID.
        """

        return header.split()[0]

    def write_input_files(self, db_conn):
        """Take the matching blast hits and write the sequence and its matching
        end to the appropriate fasta files.
        """

        log.info(
            'Writing assembler input files: iteration %i' % self.iteration)

        with open(self.file['paired_1'], 'w') as end_1, \
                open(self.file['paired_2'], 'w') as end_2:

            for row in db.get_blast_hits_by_end_count(
                    db_conn, self.iteration, 2):

                self.file['paired_count'] += 1
                out_file = end_1 if row['seq_end'] == '1' else end_2

                out_file.write('>{}/{}\n'.format(
                    row['seq_name'], row['seq_end']))
                out_file.write('{}\n'.format(row['seq']))

        with open(self.file['single_1'], 'w') as end_1, \
                open(self.file['single_2'], 'w') as end_2, \
                open(self.file['single_any'], 'w') as end_any:

            for row in db.get_blast_hits_by_end_count(
                    db_conn, self.iteration, 1):

                if row['seq_end'] == '1':
                    out_file = end_1
                    seq_end = '/1'
                    self.file['single_1_count'] += 1
                elif row['seq_end'] == '2':
                    out_file = end_2
                    seq_end = '/2'
                    self.file['single_2_count'] += 1
                else:
                    out_file = end_any
                    seq_end = ''
                    self.file['single_any_count'] += 1

                out_file.write('>{}{}\n'.format(row['seq_name'], seq_end))
                out_file.write('{}\n'.format(row['seq']))
