"""Base class for the various assembler wrappers."""

from os.path import basename, exists, getsize, join, splitext, abspath
import datetime
from subprocess import CalledProcessError
from .. import db_atram
from .. import log
from .. import bio
from .. import util


class BaseAssembler:  # pylint: disable=too-many-public-methods
    """A base class for the assemblers."""

    def __init__(self, args, cxn):
        """Build the assembler."""
        self.args = args         # Parsed command line arguments
        self.blast_only = False  # Used to short-circuit the assembler
        self.steps = []          # Assembler steps setup by the assembler
        self.file = {}           # Files and record counts

        # We need to pass these variables to child processes.
        # So they cannot be directly attached to an object.
        self.state = {
            'iteration': 0,      # Current iteration
            'query_target': '',  # Original name of the query sequence
            'query_file': '',    # Current query file name
            'blast_db': '',      # Current blast DB name
            'iter_dir': '',      # Name of the temp dir for this iteration
            'cxn': cxn}          # Save the DB connection

    def init_iteration(self, blast_db, query_file, iteration):
        """Make file names used by the assembler."""
        self.state['blast_db'] = blast_db
        self.state['query_file'] = query_file
        self.state['iteration'] = iteration
        if iteration == 1:
            self.state['query_target'] = query_file

    def setup_files(self, iter_dir):
        """Build the file names and counts for the iteration."""
        self.state['iter_dir'] = iter_dir
        self.file['long_reads'] = ''  # Set up in atram.py for now

        names = 'output paired_1 paired_2 single_1 single_2 single_any'.split()
        for name in names:
            self.file[name] = self.iter_file(name + '.fasta')

        # paired = paired_1_count + paired_2_count
        for name in 'paired single_1 single_2 single_any'.split():
            self.file[name + '_count'] = 0

    def file_prefix(self):
        """Build a prefix for the iteration's work directory."""
        return '{}_{}_{:02d}_'.format(
            basename(self.state['blast_db']),
            basename(self.state['query_target']),
            self.state['iteration'])

    def iter_file(self, file_name):
        """Put files into the temp dir."""
        rel_path = join(self.state['iter_dir'], file_name)
        # return rel_path
        return abspath(rel_path)

    def work_path(self):
        """Assembler output directory name may have unique requirements."""
        return self.state['iter_dir']

    def run(self):
        """Try to assemble the input."""
        try:
            log.info('Assembling shards with {}: iteration {}'.format(
                self.args['assembler'], self.state['iteration']))
            self.assemble()
        except TimeoutError:
            msg = 'Time ran out for the assembler after {} (HH:MM:SS)'.format(
                datetime.timedelta(seconds=self.args['timeout']))
            log.error(msg)
            raise TimeoutError(msg)
        except CalledProcessError as cpe:
            msg = 'The assembler failed with error: ' + str(cpe)
            log.error(msg)
            raise RuntimeError(msg)

    def count_blast_hits(self):
        """Make sure we have blast hits."""
        count = db_atram.sra_blast_hits_count(
            self.state['cxn'], self.state['iteration'])
        log.info('{} blast hits in iteration {}'.format(
            count, self.state['iteration']))
        return count

    def nothing_assembled(self):
        """Make there is assembler output."""
        if not exists(self.file['output']) \
                or not getsize(self.file['output']):
            log.info('No new assemblies in iteration {}'.format(
                self.state['iteration']))
            return True
        return False

    def assembled_contigs_count(self, high_score):
        """How many contigs were assembled and are above the thresholds."""
        count = db_atram.assembled_contigs_count(
            self.state['cxn'],
            self.state['iteration'],
            self.args['bit_score'],
            self.args['contig_length'])

        if not count:
            log.info('No contigs had a bit score greater than {} and are at '
                     'least {} bp long in iteration {}. The highest score for '
                     'this iteration is {}'.format(
                         self.args['bit_score'],
                         self.args['contig_length'],
                         self.state['iteration'],
                         high_score))
        return count

    def no_new_contigs(self, count):
        """Make the are new contigs in the assembler output."""
        if count == db_atram.iteration_overlap_count(
                self.state['cxn'],
                self.state['iteration'],
                self.args['bit_score'],
                self.args['contig_length']):
            log.info('No new contigs were found in iteration {}'.format(
                self.state['iteration']))
            return True
        return False

    def assemble(self):
        """
        Use the assembler to build up the contigs.

        We take an array of subprocess steps and execute them in order. We
        bracket this with pre and post assembly steps.
        """
        for step in self.steps:
            cmd = step()
            log.subcommand(cmd, self.args['temp_dir'], self.args['timeout'])
        self.post_assembly()

    def post_assembly(self):
        """Handle unique post assembly steps."""

    @staticmethod
    def parse_contig_id(header):
        """Given a fasta header line from the assembler return contig ID."""
        return header.split()[0]

    def write_input_files(self):
        """Write blast hits and matching ends to fasta files."""
        log.info('Writing assembler input files: iteration {}'.format(
            self.state['iteration']))
        self.write_paired_input_files()
        self.write_single_input_files()

    def write_paired_input_files(self):
        """Write blast hits and matching ends to fasta files."""
        with open(self.file['paired_1'], 'w') as end_1, \
                open(self.file['paired_2'], 'w') as end_2:

            for row in db_atram.get_blast_hits_by_end_count(
                    self.state['cxn'], self.state['iteration'], 2):

                self.file['paired_count'] += 1
                out_file = end_1 if row['seq_end'] == '1' else end_2

                util.write_fasta_record(
                    out_file, row['seq_name'], row['seq'], row['seq_end'])

    def write_single_input_files(self):
        """Write blast hits and matching ends to fasta files."""
        with open(self.file['single_1'], 'w') as end_1, \
                open(self.file['single_2'], 'w') as end_2, \
                open(self.file['single_any'], 'w') as end_any:

            rows = db_atram.get_blast_hits_by_end_count(
                self.state['cxn'], self.state['iteration'], 1)

            for row in rows:
                if row['seq_end'] == '1':
                    out_file = end_1
                    seq_end = '1'
                    self.file['single_1_count'] += 1
                elif row['seq_end'] == '2':
                    out_file = end_2
                    seq_end = '2'
                    self.file['single_2_count'] += 1
                else:
                    out_file = end_any
                    seq_end = ''
                    self.file['single_any_count'] += 1

                util.write_fasta_record(
                    out_file, row['seq_name'], row['seq'], seq_end)

    def final_output_prefix(self, blast_db, query):
        """Build the prefix for the name of the final output file."""
        blast_db = basename(blast_db)
        query = splitext(basename(query))[0]
        return util.prefix_file(
            self.args['output_prefix'], '{}_{}'.format(blast_db, query))

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

        count = db_atram.all_assembled_contigs_count(
            self.state['cxn'],
            self.args['bit_score'],
            self.args['contig_length'])
        if not count:
            return

        file_name = '{}.{}'.format(prefix, 'filtered_contigs.fasta')

        contigs = db_atram.get_all_assembled_contigs(
            self.state['cxn'],
            self.args['bit_score'],
            self.args['contig_length'])

        with open(file_name, 'w') as output_file:
            for contig in contigs:
                self.output_assembled_contig(output_file, contig)

    def write_all_contigs(self, prefix):
        """Write all contigs to a final output file."""
        count = db_atram.all_assembled_contigs_count(self.state['cxn'])
        if not count:
            return

        file_name = '{}.{}'.format(prefix, 'all_contigs.fasta')

        contigs = db_atram.get_all_assembled_contigs(self.state['cxn'])

        with open(file_name, 'w') as output_file:
            for contig in contigs:
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

        header = '{}_{}{} iteration={} contig_id={} score={}'.format(
            contig['iteration'], contig['contig_id'], suffix,
            contig['iteration'], contig['contig_id'], contig['bit_score'])

        util.write_fasta_record(output_file, header, seq)

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

    def simple_state(self):
        """Save the state passed to subprocesses."""
        return {
            'blast_db': self.state['blast_db'],
            'iteration': self.state['iteration'],
            'query_file': self.state['query_file'],
            'query_target': self.state['query_target'],
            'iter_dir': self.state['iter_dir']}
