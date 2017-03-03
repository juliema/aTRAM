"""Handle keeping file name consistent and some file management."""

import os
import sys
import glob
import logging
import tempfile


class Filer:
    """Handle file naming as well as some file creation and deletion."""

    def __init__(self, work_dir='', data_prefix=''):
        self.work_dir = work_dir
        self.data_prefix = data_prefix

    def log_setup(self):
        """Set up the logs for a common format. We need the prefix of the
        output log file name and the command-line arguments for the starting
        message. Both are gotten from the user input.
        """

        file_name = '{}.{}.log'.format(self.data_prefix, sys.argv[0][:-3])

        logging.basicConfig(
            filename=os.path.join(self.work_dir, file_name),
            level=logging.DEBUG,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        logging.info(' '.join(sys.argv))

    def path(self, file_name):
        """Standardize file names with a work directory and a file prefix."""

        file_name = '{}{}'.format(self.data_prefix, file_name)
        return os.path.join(self.work_dir, file_name)

    def db_file_name(self):
        """Create an SQLite3 DB name."""

        return self.path('.sqlite.db')

    def blast_shard_name(self, shard_index):
        """Create the BLAST shard DB names."""

        file_name = '.blast_{}'.format(str(shard_index + 1).zfill(3))
        return self.path(file_name)

    def all_blast_shard_names(self):
        """Get all of the BLAST DB names built by the preprocessor."""

        pattern = self.path('.blast_*.nhr')
        files = glob.glob(pattern)
        if not files:
            print(('No blast shards found. Looking for "{}"\n'
                   'Verify the --work-dir and --file-prefix options.').format(
                       pattern[:-4]))
            sys.exit()

        return sorted([f[:-4] for f in files])

    def contig_score_db(self, iteration):
        """Create the contig blast DB name."""

        file_name = '.tmp_contig_scores_{}'.format(str(iteration).zfill(2))
        return os.path.join(self.work_dir, file_name)

    def target_file(self, iteration):
        """Create a target file name for the next iteration."""

        file_name = '.tmp_target_{}.fasta'.format(str(iteration).zfill(2))
        return os.path.join(self.work_dir, file_name)

    def output_result_name(self, prefix):
        """Build the output file name."""

        file_name = '{}.assembled_contigs.fasta'.format(prefix)
        return os.path.join(self.work_dir, file_name)

    def temp_dir(self):
        """Create a temp dir that will hold files needed by blast and the
        assemblers.
        """

        return tempfile.TemporaryDirectory(dir=self.work_dir)

    def temp_file(self, temp_dir=None):
        """Create temp files for output."""

        if not temp_dir:
            temp_dir = self.work_dir

        return tempfile.NamedTemporaryFile(mode='w', dir=temp_dir)

    def open_assembler_files(self, iteration=None):
        """Open files required by the assembler."""

        return {
            'end_1': self.temp_file('.end_1', iteration=iteration),
            'end_2': self.temp_file('.end_2', iteration=iteration),
            'raw_contigs': self.temp_file('.raw_contigs', iteration=iteration),
            'new_contigs': self.temp_file('.new_contigs', iteration=iteration),
            'is_paired': False}

    @staticmethod
    def close_assembler_files(files):
        """Close files required by the assembler."""

        for key, file in files.items():
            if key not in ['is_paired']:
                file.close()

    @staticmethod
    def remove_all_starting_with(pattern):
        """Remove all files starting with the given pattern. This is
        currently used to remove blast DB files.
        """

        pattern += '.*'
        for file in glob.glob(pattern):
            os.remove(file)

    # def paired_end_file(self, iteration, end):
    #     """Create the file name of the paired end file."""
    #
    #     file_name = 'matching_seqs_{}_{}.fasta'.format(iteration, end)
    #     return self.path(file_name)

    # def contig_blast_file(self):
    #     """Create the file name of the blast DB for the assembled contigs."""
    #
    #     return self.path('blast_contigs_{}')

    # def contig_score_file(self, iteration):
    #     """Create the contig blast result file name."""
    #
    #     return self.contig_score_db(iteration) + '.csv'

    # def contig_unfiltered_file(self, iteration):
    #     """Create the file name for the contigs before they are filtered."""
    #
    #     return self.path('raw_contigs_{}.fasta'.format(
    #         str(iteration).zfill(2)))

    # def contig_filtered_file(self):
    #     """Create the file name for the contigs after they are filtered."""
    #
    #     return self.path('contigs_{}.fasta')

    # @staticmethod
    # def blast_target_against_sra(shard_name, iteration):
    #     """Get the file name of the blast result file."""
    #
    #     return '{}_{}.txt'.format(shard_name, str(iteration).zfill(2))
