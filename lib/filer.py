"""Handle keeping file name consistent and some file management."""

import os
import sys
import glob
import logging
import tempfile


class Filer:
    """Handle file naming as well as some file creation and deletion."""

    def __init__(self, work_dir='', file_prefix='', iteration=None):
        self.work_dir = work_dir
        self.file_prefix = file_prefix
        self.iteration = iteration

    def log_setup(self):
        """
        Set up the logs for a common format. We need the prefix of the output
        log file name and the command-line arguments for the starting message.
        Both are gotten from the user input.
        """

        file_name = '{}{}.log'.format(self.file_prefix, sys.argv[0][:-3])

        logging.basicConfig(
            filename=os.path.join(self.work_dir, file_name),
            level=logging.DEBUG,
            format='%(asctime)s %(levelname)s: %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S')
        logging.info(' '.join(sys.argv))

    def path(self, file_name):
        """Standardize file names with a work directory and a file prefix."""

        if self.iteration:
            file_name = file_name.format(str(self.iteration).zfill(2))

        file_name = '{}{}'.format(self.file_prefix, file_name)
        return os.path.join(self.work_dir, file_name)

    def db_file_name(self):
        """Create an SQLite3 DB name."""

        return self.path('sqlite.db')

    def all_blast_shard_names(self):
        """Get all of the BLAST DB names built by the preprocessor."""

        pattern = self.path('blast_*.nhr')
        files = glob.glob(pattern)
        if not files:
            print(('No blast shards found. Looking for "{}"\n'
                   'Verify the --work-dir and --file-prefix options.').format(
                pattern[:-4]))
            sys.exit()

        return sorted([f[:-4] for f in files])

    def blast_shard_name(self, shard_index):
        """Create the BLAST shard DB names."""

        file_name = 'blast_{}'.format(str(shard_index + 1).zfill(3))
        return self.path(file_name)

    def paired_end_file(self, end):
        """Create the file name of the paired end file."""

        file_name = 'matching_seqs_{}_{}.fasta'.format('{}', end)
        return self.path(file_name)

    def contig_blast_file(self):
        """Create the file name of the blast DB for the assembled contigs."""

        return self.path('blast_contigs_{}')

    def contig_score_db(self):
        """Create the contig blast DB name."""

        return self.path('contig_scores_{}')

    def contig_score_file(self):
        """Create the contig blast result file name."""

        return self.contig_score_db() + '.csv'

    def contig_unfiltered_file(self):
        """Create the file name for the contigs before they are filtered."""

        return self.path('raw_contigs_{}.fasta')

    def contig_filtered_file(self):
        """Create the file name for the contigs after they are filtered."""

        return self.path('contigs_{}.fasta')

    def blast_target_against_sra(self, shard_name):
        """Get the file name of the blast result file."""

        return '{}_{}.txt'.format(shard_name, str(self.iteration).zfill(2))

    @staticmethod
    def temp_file():
        """Create temp files for output. Nest these in a "with" statement."""

        return tempfile.NamedTemporaryFile(mode='w')
