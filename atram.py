"""The aTRAM assembly program."""

import os
import csv
import logging
import subprocess
import multiprocessing
from Bio import SeqIO
import lib.db as db
import lib.bio as bio
from lib.filer import Filer
from lib.configure import Configure
from lib.assembler import Assembler


class Atram:
    """The atram program itself."""

    def __init__(self):
        self.filer = None
        self.config = None
        self.db_conn = None
        self.assembler = None
        self.iteration = 0
        self.shard_list = None
        self.temp_dir = None

    def run(self):
        """Setup  and run atram."""

        self.config = Configure().parse_command_line(
            description=""" """,
            args=('target assembler output_prefix protein work_dir '
                  'data_prefix_atram iterations bit_score evalue '
                  'max_target_seqs genetic_code kmer cpus max_memory'))

        self.filer = Filer(self.config.work_dir, self.config.data_prefix)
        self.filer.log_setup()

        self.shard_list = self.filer.all_blast_shard_names()
        self.assembler = Assembler.factory(self.config)

        self.db_conn = db.connect(self.filer)
        db.create_blast_hits_table(self.db_conn)
        db.create_blast_hits_index(self.db_conn)
        db.create_assembled_contigs_table(self.db_conn)
        db.create_assembled_contigs_index(self.db_conn)

        with self.filer.temp_dir() as temp_dir:
            self.temp_dir = temp_dir.name
            self.main_loop()
            self.output_results()

    def main_loop(self):
        """The main program loop."""

        target = self.config.target
        for self.iteration in range(1, self.config.iterations + 1):

            logging.info('aTRAM iteration %i', self.iteration)

            self.blast_all_targets_against_sra(target)

            assembler_files = self.filer.open_assembler_files(self.iteration)
            self.write_paired_end_files(assembler_files)
            self.get_previous_assemblies(assembler_files)

            cwd = os.getcwd()
            try:
                os.chdir(self.config['work_dir'])  # some assemblers need this
                self.assembler.assemble(assembler_files)
            finally:
                os.chdir(cwd)

            self.filter_contigs(assembler_files)

            target = self.create_targets_from_contigs()

    def get_previous_assemblies(self):
        """Get assembiles from the last iteration. This is easier than
        juggling temp files.
        """

    def output_results(self):
        """Write the assembled contigs to a fasta file."""

        file_name = self.filer.output_result_name(self.config.output_prefix)
        with open(file_name, 'w') as out_file:
            for row in db.get_all_assembled_contigs(self.db_conn):
                header = ('>{}_{} iteration={} contig_id={} '
                          'score={}\n').format(
                    row['iteration'], row['contig_id'],
                    row['iteration'], row['contig_id'], row['bit_score'])
                out_file.write(header)
                out_file.write('{}\n'.format(row['seq']))

    def create_targets_from_contigs(self):
        """Crate a new file with the contigs that will be used as the
        next target.
        """

        target = self.filer.target_file(self.iteration)
        with open(target, 'w') as target_file:
            for row in db.get_assembled_contigs(self.db_conn, self.iteration):
                target_file.write('>{}\n'.format(row[0]))
                target_file.write('{}\n'.format(row[1]))

        return target

    def cleanup_files(self, assembler_files):
        """Cleanup files"""

        # self.filer.close_assembler_files(assembler_files)
        # blast_db = self.filer.contig_score_db(self.iteration)
        # self.filer.remove_all_starting_with(blast_db)

    def blast_all_targets_against_sra(self, target):
        """Blast the targets against the SRA databases. We're using a
        map-reduce strategy here. We map the blasting of the target sequences
        and reduce the output into one fasta file.
        """

        # Pass these arguments to the blast subprocess
        blast_params = {
            'target': target,
            'evalue': self.config.evalue,
            'protein': self.config.protein,
            'work_dir': self.config.work_dir,
            'iteration': self.iteration,
            'data_prefix': self.config.data_prefix,
            'genetic_code': self.config.genetic_code,
            'max_target_seqs': self.config.max_target_seqs}

        with open(blast_params['target']) as in_file:
            lines = in_file.readlines()

        with multiprocessing.Pool(processes=self.config.cpus) as pool:
            results = [pool.apply_async(
                blast_target_against_sra,
                (blast_params, shard_name, self.iteration))
                       for shard_name in self.shard_list]
            _ = [result.get() for result in results]

    def write_paired_end_files(self, assembler_files):
        """Take the matching blast hits and write the sequence and any matching
        end to the appropriate fasta files.
        """

        assembler_files['is_paired'] = 0

        for row in db.get_blast_hits(self.db_conn, self.iteration):
            if row[1] and row[1].endswith(('1', '2')):
                end = row[1][-1]
                key = 'end_' + end
                assembler_files[key].write('>{}/{}\n'.format(row[0], end))
                assembler_files[key].write('{}\n'.format(row[2]))
                assembler_files['is_paired'] += int(end == '2')
            else:
                assembler_files['end_1'].write('>{}\n'.format(row[0]))
                assembler_files['end_1'].write('{}\n'.format(row[2]))

    def filter_contigs(self, assembler_files):
        """Remove junk from the assembled contigs."""

        blast_db = self.filer.contig_score_db(self.iteration)
        raw_contigs = assembler_files['raw_contigs'].name
        self.create_blast_db_from_contigs(blast_db, raw_contigs)
        self.blast_target_against_contigs(blast_db, assembler_files)
        filtered_scores = self.filter_contig_scores(assembler_files)
        self.save_contigs(assembler_files, filtered_scores, raw_contigs)

    def save_contigs(self, assembler_files, filtered_scores, raw_contigs):
        """Save the contigs to the database."""

        batch = []
        with open(raw_contigs) as in_file:
            for contig in SeqIO.parse(in_file, 'fasta'):
                if contig.id in filtered_scores:
                    score = filtered_scores[contig.id]
                    batch.append((
                        self.iteration, contig.id,
                        str(contig.seq), contig.description,
                        score['bit_score'], score['target_start'],
                        score['target_end'], score['contig_start'],
                        score['contig_end'], score['contig_len']))
        db.insert_assembled_contigs_batch(self.db_conn, batch)

    @staticmethod
    def create_blast_db_from_contigs(blast_db, raw_contigs):
        """Create a blast DB from the assembled fragments."""

        cmd = 'makeblastdb -dbtype nucl -in {} -out {}'
        cmd = cmd.format(raw_contigs, blast_db)
        subprocess.check_call(cmd, shell=True)

    def blast_target_against_contigs(self, blast_db, assembler_files):
        """Blast the target sequence against the contings. The blast output
        will have the scores for later processing.
        """

        cmd = []

        if self.config.protein:
            cmd.append('tblastn')
            cmd.append('-db_gencode {}'.format(self.config.genetic_code))
        else:
            cmd.append('blastn')

        cmd.append('-db {}'.format(blast_db))
        cmd.append('-query {}'.format(self.config.target))
        cmd.append('-out {}'.format(assembler_files['new_contigs'].name))
        cmd.append(
            "-outfmt '10 qseqid sseqid bitscore qstart qend sstart send slen'")

        cmd = ' '.join(cmd)

        subprocess.check_call(cmd, shell=True)

    def filter_contig_scores(self, assembler_files):
        """Only save contigs that have scores above the bit score cut-off."""

        # qseqid sseqid bitscore qstart qend sstart send slen
        field_names = ['target_id', 'contig_id', 'bit_score',
                       'target_start', 'target_end',
                       'contig_start', 'contig_end', 'contig_len']

        scores = {}
        with open(assembler_files['new_contigs'].name) as in_file:
            for score in csv.DictReader(in_file, field_names):
                score['bit_score'] = float(score['bit_score'])
                if score['bit_score'] >= self.config.bit_score:
                    for field in field_names[3:]:
                        score[field] = int(score[field])
                    scores[score['contig_id']] = score
        return scores


def blast_target_against_sra(blast_params, shard_name, iteration):
    """Blast the target sequences against an SRA blast DB."""

    filer = Filer(work_dir=blast_params['work_dir'],
                  data_prefix=blast_params['data_prefix'])
    db_conn = db.connect(filer)

    with filer.temp_file() as results_file:
        cmd = create_blast_command(blast_params, shard_name, results_file.name)
        subprocess.check_call(cmd, shell=True)
        insert_blast_hits(db_conn, results_file, iteration, shard_name)


def insert_blast_hits(db_conn, results_file, iteration, shard_name):
    """Put all of the blast hits into the database."""
    batch = []
    with open(results_file.name) as blast_hits:
        for line in blast_hits:
            match = bio.PARSE_BLAST_RESULTS.match(line)
            if match:
                seq_name = match.group(1)
                seq_end = match.group(2)
            else:
                seq_name = line
                seq_end = ''
            batch.append((iteration, seq_end, seq_name, shard_name))
    db.insert_blast_hit_batch(db_conn, batch)


def create_blast_command(blast_params, shard_name, out):
    """Build the blast command to blast the target against the SRA database
    shards.
    """

    cmd = []

    if blast_params['protein'] and blast_params['iteration'] == 1:
        cmd.append('tblastn')
        cmd.append('-db_gencode {}'.format(blast_params['genetic_code']))
    else:
        cmd.append('blastn')
        cmd.append('-evalue {}'.format(blast_params['evalue']))

    cmd.append("-outfmt '10 sseqid'")
    cmd.append('-max_target_seqs {}'.format(blast_params['max_target_seqs']))
    cmd.append('-out {}'.format(out))
    cmd.append('-db {}'.format(shard_name))
    cmd.append('-query {}'.format(blast_params['target']))

    cmd = ' '.join(cmd)

    return cmd


if __name__ == '__main__':

    Atram().run()
