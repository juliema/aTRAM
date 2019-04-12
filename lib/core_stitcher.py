"""
Stitch together exons from targeted assemblies.

It uses amino acid targets and DNA assemblies.
"""

import os
from os.path import abspath, join, basename, splitext
import re
from collections import namedtuple
import subprocess
from glob import glob
from Bio.SeqIO.FastaIO import SimpleFastaParser
import lib.db as db
import lib.log as log
import lib.util as util


ExonerateHeader = namedtuple(
    'ExonerateHeader',
    ['gene', 'taxon', 'query_len', 'query_align_len',
     'query_align_beg', 'query_align_end', 'target_id'])


class Sticher:
    """Container for sticher data."""

    def __init__(self, args):
        """Setup."""
        self.args = args
        self.cxn = None
        self.temp_dir = None

    def stitch(self):
        """Stitch the exons together."""
        log.stitcher_setup(self.args['log_file'])

        with util.make_temp_dir(
                where=self.args['temp_dir'],
                prefix='atram_stitcher_',
                keep=self.args['keep_temp_dir']) as temp_dir:
            self.temp_dir = temp_dir

            self.cxn = db.temp_db(temp_dir, 'atram_stitcher')

            self.prepare_fasta_files()
            self.prepare_reference_genes()
            # self.find_fasta_files_for_genes()
            # self.run_exonerate()
            # self.get_contigs()
            # self.stitch_contigs()
            # self.summary_stats()

    def prepare_fasta_files(self):
        """Prepare fasta files for exonerate."""
        db.create_contigs_table(self.cxn)

        batch = []

        pattern = join(self.args['assemblies_dir'], '*.fasta')
        for fasta_path in sorted(glob(pattern)):
            if os.stat(fasta_path).st_size == 0:
                continue

            log.info('Preparing fasta file: {}'.format(fasta_path))

            tiebreaker = 0

            with open(fasta_path) as fasta_old:
                for seq_name, seq in SimpleFastaParser(fasta_old):
                    tiebreaker += 1
                    seq_name = seq_name.strip()
                    seq_name = re.sub(r'[\s,=@+\[\]:!-]+', '_', seq_name)
                    seq_name = '{}_{}'.format(seq_name, tiebreaker)
                    batch.append([fasta_path, seq_name, seq])

        db.insert_contigs(self.cxn, batch)

    def prepare_reference_genes(self):
        """Prepare reference sequences for exonerate."""
        batch = []
        db.create_reference_table(self.cxn)

        ref_genes = self.args['reference_genes']

        log.info('Preparing reference genes: {}'.format(ref_genes))

        with open(ref_genes) as ref_in:

            for seq_name, seq in SimpleFastaParser(ref_in):
                seq_name = re.sub(r'\s+', '_', seq_name)

                batch.append([seq_name, seq])

        db.insert_references(self.cxn, batch)

    # def find_fasta_files_for_genes(self):
    #     """
    #     Find all assemblies for every gene in the gene list.
    #
    #     Put this list into file for the gene. One file listing all assemblies
    #     for every gene in the master gene list.
    #     """
    #     with open(self.all_ref_genes_path) as all_ref_genes_file:
    #         for ref_gene in all_ref_genes_file:
    #             gene_name = self.gene_name_from_ref_gene(ref_gene)
    #
    #             log.info('Finding fasta files for gene: {}'.format(gene_name))
    #
    #             with open(self.gene_list_path(gene_name), 'w') as gene_list:
    #                 pattern = join(
    #                     self.temp_dir, '*{}*.ed.fasta'.format(gene_name))
    #                 fasta_list = glob(pattern)
    #
    #                 for fasta in fasta_list:
    #                     gene_list.write(fasta)
    #                     gene_list.write('\n')
    #
    # def run_exonerate(self):
    #     """Run exonerate on every reference sequence, taxon combination."""
    #     db.create_exonerate_results_table(self.cxn)
    #
    #     with open(self.all_ref_genes_path) as ref_genes_file:
    #         ref_genes = {x.strip() for x in ref_genes_file}
    #
    #     with open(self.args['taxa_list']) as taxa_list:
    #         taxa = {x.strip() for x in taxa_list}
    #
    #     for ref_gene in sorted(ref_genes):
    #         gene_name = self.gene_name_from_ref_gene(ref_gene)
    #         log.info('Reference file: {}'.format(gene_name))
    #         gene_asm_path = self.gene_list_path(gene_name)
    #         with open(gene_asm_path) as gene_asm_file:
    #             asms = [asm.strip() for asm in gene_asm_file
    #                     if self.taxon_from_asm(gene_name, asm) in taxa]
    #             for asm in sorted(asms):
    #                 taxon = self.taxon_from_asm(gene_name, asm)
    #                 self.exonerate_command(ref_gene, gene_name, taxon, asm)
    #         self.insert_exonerate_results(gene_name)

    def exonerate_command(self, ref_file, gene_name, taxon, asm):
        """Build and run the exonerate program."""
        results_file = self.results_file_path(gene_name)
        cmd = (
            'exonerate --verbose 0 --model protein2genome {ref_file} {asm} '
            '--showvulgar no --showalignment no '
            r'--ryo ">{gene},{taxon},%ql,%qal,%qab,%qae,%ti\n%tcs\n" '
            '>> {results_file};'
        ).format(ref_file=ref_file, asm=asm, gene=gene_name, taxon=taxon,
                 results_file=results_file)

        subprocess.check_call(cmd, shell=True)

    def insert_exonerate_results(self, gene_name):
        """Insert the exonerate results into the database."""
        results_file = self.results_file_path(gene_name)
        batch = []
        with open(results_file) as results_fasta:
            for header, seq in SimpleFastaParser(results_fasta):
                header = header.split(',')
                field = ExonerateHeader(*header)
                batch.append([
                    field.gene,
                    field.taxon,
                    field.query_len,
                    field.query_align_len,
                    field.query_align_beg,
                    field.query_align_end,
                    field.target_id,
                    seq])
        db.insert_exonerate_results(self.cxn, batch)

    def get_contigs(self):
        """
        Sort contigs so they are in the same order they appear in the gene.

        1) Sort by the beginning of the contig.
        2) Look for all sorted results from exonerate.
        3) Create a CSV file for each input with all of the information about
            the contigs from exonerate.
        """

    def ref_gene_file_path(self, seq_name):
        """Build reference gene file name from the sequence name."""
        return join(self.temp_dir, '{}.reference.fasta'.format(seq_name))

    def gene_list_path(self, gene_name):
        """Build file name with list of assemblies for a reference gene."""
        return join(self.temp_dir, '{}.list.txt'.format(gene_name))

    def results_file_path(self, gene_name):
        """Build the file name for the exonerate results."""
        results_file = join(
            self.temp_dir, '{}.results.fasta'.format(gene_name))
        return abspath(results_file)

    @staticmethod
    def gene_name_from_ref_gene(ref_file):
        """Get the gene name from the reference gene file path."""
        ref_file = ref_file.strip()
        gene_name = basename(ref_file)
        return re.sub(r'\.reference\.fasta$', '', gene_name)

    @staticmethod
    def taxon_from_asm(gene_name, asm):
        """Check if the assembly in the taxa set."""
        asm = basename(asm)
        asm = re.sub(r'^{}\.'.format(gene_name), '', asm)
        return re.sub(r'\.ed\.fasta\s*$', '', asm)
