"""The aTRAM assembly program."""

import re
import csv
import logging
import sqlite3
import subprocess
import multiprocessing
from more_itertools import chunked
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline as tblastn
from Bio.Blast.Applications import NcbitblastnCommandline as blastn
from lib.filer import Filer
from lib.configure import Configure
from lib.assembler import Assembler


class Atram:
    """The atram program itself."""

    fragment = re.compile(r'^ ( .* ) ( [\s\/_] [12] )', re.VERBOSE)

    def __init__(self):
        self.filer = None
        self.config = None
        self.db_conn = None
        self.assember = None
        self.iteration = 0
        self.shard_list = None

    def run(self):
        """Setup  and run atram."""
        self.config = Configure().parse_command_line(
            description=""" """,
            args=('target protein iterations cpus evalue max_target_seqs '
                  'assembler max_memory file_prefix work_dir bit_score '
                  'genetic_code kmer'))

        self.filer = Filer(self.config.work_dir, self.config.file_prefix)
        self.filer.log_setup()

        self.assember = Assembler.factory(self.config, self.iteration)
        self.shard_list = self.filer.all_blast_shard_names()

        self.atram()

    def atram(self):
        """The main aTRAM program loop."""

        target = self.config.target
        for self.iteration in range(1, self.config.iterations + 1):
            logging.info('aTRAM iteration %i', self.iteration)
            self.blast_all_targets_against_sra(target)
            sequences = self.get_blast_hits()
            paired = self.write_paired_end_files(sequences)
            # assember.assemble(iteration, paired)
            # filter_contigs(iteration)
            # target = new file
            break

    def blast_all_targets_against_sra(self, target):
        """
        Blast the targets against the SRA databases. We're using a map-reduce
        strategy here. We map the blasting of the target sequences and reduce
        the output into one fasta file.
        """

        # Pass these arguments to the blast subprocess
        blast_params = {
            'target': target,
            'evalue': self.config.evalue,
            'protein': self.config.protein,
            'work_dir': self.config.work_dir,
            'iteration': self.iteration,
            'file_prefix': self.config.file_prefix,
            'genetic_code': self.config.genetic_code,
            'max_target_seqs': self.config.max_target_seqs,
        }

        with multiprocessing.Pool(processes=self.config.cpus) as pool:
            results = [pool.apply_async(blast_target_against_sra,
                                        (blast_params, shard_name))
                       for shard_name in self.shard_list]
            _ = [result.get() for result in results]

    def get_blast_hits(self):
        """Get all of the blast hits"""

        sequences = {}
        for shard_name in self.shard_list:
            hits_file_name = self.filer.blast_target_against_sra(shard_name)
            print(hits_file_name)
            with open(hits_file_name) as blast_hits:
                for line in blast_hits:
                    match = Atram.fragment.match(line)
                    sequences[match.group(1) if match else line] = 1
        return sequences

    def write_paired_end_files(self, sequences):
        """
        Take the matching blast hits and write the sequence and any matching
        end to the appropriate fasta files.
        """

        # db_conn = connect_db(work_dir, file_prefix)
        # fasta_1 = self.filer.paired_end_file(work_dir, file_prefix, iteration, '1')
        # fasta_2 = self.filer.paired_end_file(work_dir, file_prefix, iteration, '2')
        # paired = False
        # with open(fasta_1, 'w') as file_1, open(fasta_2, 'w') as file_2:
        #     for ids in chunked(fragments.keys(), 100):
        #         sql = 'SELECT * FROM frags WHERE frag IN ({})'
        #         sql = sql.format(','.join('?' for _ in ids))
        #         for row in db_conn.execute(sql, ids):
        #             if row[1] and row[1].endswith('2'):
        #                 file_2.write('>{}{}\n'.format(row[0], row[1]))
        #                 file_2.write('{}\n'.format(row[2]))
        #                 paired = True
        #             else:
        #                 file_1.write('>{}{}\n'.format(row[0], row[1]))
        #                 file_1.write('{}\n'.format(row[2]))
        # return paired

    def connect_db(self):
        """Setup the DB for our processing needs."""

        db_path = self.filer.db_file_name()
        db_conn = sqlite3.connect(db_path)
        return db_conn


#
#
#
#
# def create_blast_db(work_dir, file_prefix, iteration):
#     """Create a blast DB from the assembled fragments."""
#
#     blast_db = self.filer.contig_score_db(work_dir, file_prefix, iteration)
#     fasta_file = self.filer.contig_unfiltered_file(
#           work_dir, file_prefix, iteration)
#     cmd = 'makeblastdb -dbtype nucl -in {} -out {}'
#     cmd = cmd.format(fasta_file, blast_db)
#     subprocess.check_call(cmd, shell=True)
#     return blast_db
#
#
# def blast_target_against_contigs(work_dir, file_prefix, iteration,
#                                  target, genetic_code, protein):
#     """
#    Blast the target sequence against the contings. The blast output will have
#     the scores for later processing.
#     """
#
#     blast_db = create_blast_db(work_dir, file_prefix, iteration)
#     scored_contigs = self.filer.contig_score_file(
#                                   work_dir, file_prefix, iteration)
#     blast_args = dict(
#         db=blast_db, query=target, out=scored_contigs,
#         outfmt="'10 qseqid sseqid bitscore qstart qend sstart send slen'")
#     if protein:
#         cmd = str(blastn(
#             cmd='tblastn', db_gencode=genetic_code, **blast_args))
#     else:
#         cmd = str(tblastn(
#             cmd='blastn', task='blastn', **blast_args))
#     print(cmd)
#     subprocess.check_call(cmd, shell=True)
#
#
# def filter_contig_scores(work_dir, file_prefix, iteration, bit_score):
#     """Remove contigs that have scores below the bit score cut-off."""
#
#     # qseqid sseqid bitscore qstart qend sstart send slen
#     field_names = ['target_id', 'contig_id', 'bit_score',
#                    'target_start', 'target_end',
#                    'contig_start', 'contig_end', 'contig_len']
#     contig_score_file = self.filer.contig_score_file(work_dir, file_prefix,
#                                                iteration)
#
#     scores = {}
#     with open(contig_score_file) as in_file:
#         for score in csv.DictReader(in_file, field_names):
#             score['bit_score'] = float(score['bit_score'])
#             if score['bit_score'] >= bit_score:
#                 for field in field_names[3:]:
#                     score[field] = int(score[field])
#                 scores[score['contig_id']] = score
#
#     return scores
#
#
# def filter_contigs(work_dir, file_prefix, iteration):
#     """Remove junk from the assembled contigs."""
#
#     blast_target_against_contigs(work_dir, file_prefix, iteration,
#                                  target, genetic_code, protein)
#     filtered_names = filter_contig_scores(work_dir, file_prefix,
#                                           iteration, bit_score)
#     unfiltered_file = self.filer.contig_unfiltered_file(work_dir, file_prefix,
#                                                   iteration)
#     # filtered_file = self.filer.contig_filtered_file(work_dir, file_prefix,
#     #                                           iteration)
#
#     print(filtered_names)
#     with open(unfiltered_file) as in_file:
#         for contig in SeqIO.parse(in_file, 'fasta'):
#             if contig.id in filtered_names:
#                 print(contig)


def blast_target_against_sra(blast_params, shard_name):
    """
    Blast the target sequences against an SRA blast DB.
    """

    filer = Filer(work_dir=blast_params['work_dir'],
                  iteration=blast_params['iteration'],
                  file_prefix=blast_params['file_prefix'])
    out = filer.blast_target_against_sra(shard_name)

    blast_args = dict(outfmt="'10 sseqid'",
                      max_target_seqs=blast_params['max_target_seqs'],
                      out=out, db=shard_name, query=blast_params['target'])
    if blast_params['protein'] and blast_params['iteration'] == 1:
        cmd = str(blastn(cmd='tblastn',
                         db_gencode=blast_params['genetic_code'],
                         **blast_args))
    else:
        cmd = str(tblastn(
            cmd='blastn', task='blastn', evalue=blast_params['evalue'],
            **blast_args))
    print(cmd)
    subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':

    Atram().run()
