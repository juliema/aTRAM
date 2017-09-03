"""Base class for the various assembler wrappers."""

from os.path import abspath, basename, exists, getsize, join, splitext
import datetime
import subprocess
import lib.db as db
import lib.log as log
import lib.bio as bio
import lib.file_util as file_util


class BaseAssembler:
    """A base class for the assemblers."""

    def __init__(self, args, db_conn):
        """Build the assembler."""
        self.args = args      # Parsed command line arguments
        self.steps = []       # Assembler steps setup by the assembler
        self.file = {}        # Files and record counts
        self.iteration = 0    # Current iteration
        self.query_file = ''  # Current query file name
        self.blast_db = ''    # Current blast DB name
        self.db_conn = db_conn  # Save the DB connection
        self.blast_only = False  # Used to short-circuit the assembler

    def iter_dir(self):
        """Get the work directory for the current iteration."""
        return file_util.temp_iter_dir(self.args['temp_dir'],
                                       self.blast_db,
                                       self.query_file,
                                       self.iteration)

    def work_path(self):
        """The output directory name may have unique requirements."""
        return self.iter_dir()

    def iter_file(self, file_name):
        """Build a temporary file in the current iteration directory."""
        return join(self.iter_dir(), file_name)

    def run(self):
        """Try to assemble the input."""
        try:
            log.info('Assembling shards with {}: iteration {}'.format(
                self.args['assembler'], self.iteration))
            self.assemble()
        except TimeoutError:
            msg = 'Time ran out for the assembler after {} (HH:MM:SS)'.format(
                datetime.timedelta(seconds=self.args['timeout']))
            log.fatal(msg)
        except subprocess.CalledProcessError as cpe:
            msg = 'The assembler failed with error: ' + str(cpe)
            log.fatal(msg)

    def no_blast_hits(self):
        """Make sure we have blast hits."""
        if not db.sra_blast_hits_count(self.db_conn, self.iteration):
            log.info('No blast hits in iteration %i' % self.iteration)
            return True

        return False

    def nothing_assembled(self):
        """Make there is assembler output."""
        if not exists(self.file['output']) \
                or not getsize(self.file['output']):
            log.info('No new assemblies in iteration %i' % self.iteration)
            return True

        return False

    def assembled_contigs_count(self, high_score):
        """How many contigs were assembled and are above the thresholds."""
        count = db.assembled_contigs_count(
            self.db_conn,
            self.iteration,
            self.args['bit_score'],
            self.args['contig_length'])

        if not count:
            log.info('No contigs had a bit score greater than {} and are at '
                     'least {} long in iteration {}. The highest score for '
                     'this iteration is {}'.format(
                         self.args['bit_score'],
                         self.args['contig_length'],
                         self.iteration,
                         high_score))
        return count

    def no_new_contigs(self, count):
        """Make the are new contigs in the assembler output."""
        if count == db.iteration_overlap_count(
                self.db_conn,
                self.iteration,
                self.args['bit_score'],
                self.args['contig_length']):
            log.info(
                'No new contigs were found in iteration %i' % self.iteration)
            return True

        return False

    def assemble(self):
        """Use the assembler to build up the contigs.

        We take and array of subprocess steps and execute them in order. We
        bracket this with pre and post assembly steps.
        """
        for step in self.steps:
            cmd = step()
            log.subcommand(cmd, self.args['temp_dir'], self.args['timeout'])

        self.post_assembly()

    def post_assembly(self):
        """The assembler may have unique post assembly steps."""

    def path(self, file_name):
        """Put files into the temp dir."""
        blast_db = basename(self.blast_db)
        file_name = '{}.{}.{:02d}.{}'.format(
            blast_db, self.query_file, self.iteration, file_name)
        rel_path = self.iter_file(file_name)

        return abspath(rel_path)

    def initialize_iteration(self, blast_db, query_file, iteration):
        """Setup file names used by the assembler.

        Do this at the start of each iteration.
        """
        self.blast_db = blast_db
        self.query_file = query_file
        self.iteration = iteration

        self.file['output'] = self.path('output.fasta')
        self.file['long_reads'] = ''  # Set up in atram.py for now
        self.file['paired_1'] = self.path('paired_1.fasta')
        self.file['paired_2'] = self.path('paired_2.fasta')
        self.file['single_1'] = self.path('single_1.fasta')
        self.file['single_2'] = self.path('single_2.fasta')
        self.file['single_any'] = self.path('single_any.fasta')

        self.file['paired_count'] = 0       # paired_1 and paired_2 count
        self.file['single_1_count'] = 0
        self.file['single_2_count'] = 0
        self.file['single_any_count'] = 0

    @staticmethod
    def parse_contig_id(header):
        """Given a fasta header line from the assembler return contig ID."""
        return header.split()[0]

    def write_input_files(self):
        """Write blast hits and matching ends to fasta files."""
        log.info(
            'Writing assembler input files: iteration %i' % self.iteration)

        with open(self.file['paired_1'], 'w') as end_1, \
                open(self.file['paired_2'], 'w') as end_2:

            for row in db.get_blast_hits_by_end_count(
                    self.db_conn, self.iteration, 2):

                self.file['paired_count'] += 1
                out_file = end_1 if row['seq_end'] == '1' else end_2

                out_file.write('>{}/{}\n'.format(
                    row['seq_name'], row['seq_end']))
                out_file.write('{}\n'.format(row['seq']))

        with open(self.file['single_1'], 'w') as end_1, \
                open(self.file['single_2'], 'w') as end_2, \
                open(self.file['single_any'], 'w') as end_any:

            for row in db.get_blast_hits_by_end_count(
                    self.db_conn, self.iteration, 1):

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

    def final_output_prefix(self, blast_db, query):
        """Build the prefix for the name of the final output file."""
        blast_db = basename(blast_db)
        query = splitext(basename(query))[0]
        return '{}.{}_{}'.format(self.args['output_prefix'], blast_db, query)

    def write_final_output(self, blast_db, query):
        """Write the assembler results file.

        In this default case we're writing assembled contigs to fasta files.
        """
        prefix = self.final_output_prefix(blast_db, query)

        self.write_filtered_contigs(prefix)
        self.write_all_contigs(prefix)

    def write_filtered_contigs(self, prefix):
        """Write to the filtered contigs to a final output file."""
        if self.args['no_filter']:
            return

        file_name = file_util.output_file(prefix, 'filtered_contigs.fasta')

        contigs = db.get_all_assembled_contigs(
            self.db_conn, self.args['bit_score'], self.args['contig_length'])

        with open(file_name, 'w') as output_file:
            for contig in contigs:
                self.output_assembled_contig(output_file, contig)

    def write_all_contigs(self, prefix):
        """Write all contigs to a final ouput file."""
        file_name = file_util.output_file(prefix, 'all_contigs.fasta')

        with open(file_name, 'w') as output_file:
            for contig in db.get_all_assembled_contigs(self.db_conn):
                self.output_assembled_contig(output_file, contig)

    @staticmethod
    def output_assembled_contig(output_file, contig):
        """Write one assembled contig to the output fasta file."""
        seq = contig['seq']
        suffix = ''

        if contig['query_strand'] and contig['hit_strand'] and \
                contig['query_strand'] != contig['hit_strand']:
            seq = bio.reverse_complement(seq)
            suffix = '_REV'

        header = '>{}_{}{} iteration={} contig_id={} score={}\n'.format(
            contig['iteration'], contig['contig_id'], suffix,
            contig['iteration'], contig['contig_id'], contig['bit_score'])

        output_file.write(header)
        output_file.write('{}\n'.format(seq))

    def get_single_ends(self):
        """Gather single ends files for the assembly command."""
        single_ends = []
        if self.file['single_1_count']:
            single_ends.append(self.file['single_1'])
        if self.file['single_2_count']:
            single_ends.append(self.file['single_2'])
        if self.file['single_any_count']:
            single_ends.append(self.file['single_any'])
        return single_ends
