"""
Stitch together exons from targeted assemblies.

It uses amino acid targets and DNA assemblies.
"""

import os
from os.path import abspath, join, split, splitext
import re
import subprocess
from glob import glob

from Bio.SeqIO.FastaIO import SimpleFastaParser

import lib.log as log
import lib.util as util
# import lib.db as db


class Sticher:
    """Container for sticher data."""

    def __init__(self, args):
        """Setup."""
        self.args = args
        self.temp_dir = None

    def stitch(self):
        """Stitch the exons together."""
        log.stitcher_setup(self.args['log_file'])

        with util.make_temp_dir(
                where=self.args['temp_dir'],
                prefix='atram_stitcher_',
                keep=self.args['keep_temp_dir']) as temp_dir:
            self.temp_dir = temp_dir

            self.prepare_fasta_files()
            self.prepare_reference_genes()
            self.find_fasta_files_for_genes()
            self.run_exonerate()
            self.get_gontigs()
            # self.stitch_contigs()
            # self.summary_stats()

    def prepare_fasta_files(self):
        """Prepare fasta files for exonerate."""
        pattern = join(self.args['assemblies_dir'], '*.fasta')

        for fasta_old_path in sorted(glob(pattern)):
            if os.stat(fasta_old_path).st_size == 0:
                continue

            log.info('Preparing fasta file: {}'.format(fasta_old_path))

            fasta_new_path = self.new_fasta_name(fasta_old_path)
            tiebreaker = 0

            with open(fasta_old_path) as fasta_old, \
                    open(fasta_new_path, 'w') as fasta_new:

                for seq_name, seq in SimpleFastaParser(fasta_old):
                    tiebreaker += 1
                    seq_name = seq_name.strip()
                    seq_name = re.sub(r'[\s,=@+\[\]:!-]+', '_', seq_name)
                    seq_name = '{}_{}'.format(seq_name, tiebreaker)

                    util.write_fasta_record(fasta_new, seq_name, seq)

    def prepare_reference_genes(self):
        """
        Format reference sequences for exonerate.

        1) Touch up the sequence name.
        2) Put each sequence into its own file.
        3) Save a list of all of the reference gene files in a master list.
        """
        ref_genes = self.args['reference_genes']

        log.info('Preparing reference genes: {}'.format(ref_genes))

        with open(ref_genes) as ref_in, \
                open(self.all_ref_genes_path, 'w') as all_ref_genes_file:

            for seq_name, seq in SimpleFastaParser(ref_in):
                seq_name = re.sub(r'\s+', '_', seq_name)

                ref_gene_path = self.ref_gene_file_path(seq_name)

                all_ref_genes_file.write(ref_gene_path)
                all_ref_genes_file.write('\n')

                with open(ref_gene_path, 'w') as ref_gene_file:
                    util.write_fasta_record(ref_gene_file, seq_name, seq)

    def find_fasta_files_for_genes(self):
        """
        Find all assemblies for every gene in the gene list.

        Put this list into file for the gene. One file listing all assemblies
        for every gene in the master gene list.
        """
        with open(self.all_ref_genes_path) as all_ref_genes_file:
            for ref_gene in all_ref_genes_file:
                gene_name = self.gene_name_from_ref_gene(ref_gene)

                log.info('Finding fasta files for gene: {}'.format(gene_name))

                with open(self.gene_list_path(gene_name), 'w') as gene_list:
                    pattern = join(
                        self.temp_dir, '*{}*.ed.fasta'.format(gene_name))
                    fasta_list = glob(pattern)

                    for fasta in fasta_list:
                        gene_list.write(fasta)
                        gene_list.write('\n')

    def run_exonerate(self):
        """
        Run exonerate on every reference sequence, taxon combination.

        Get exonerate information for each exon.
        Then sort the results file by taxon and by location of the exon.
        Use exonerate to pull out the exons.
        """
        with open(self.all_ref_genes_path) as ref_genes_file:
            ref_genes = {x.strip() for x in ref_genes_file}

        with open(self.args['taxa_list']) as taxa_list:
            taxa = {x.strip() for x in taxa_list}

        for ref_gene in sorted(ref_genes):
            gene_name = self.gene_name_from_ref_gene(ref_gene)
            log.info('Reference file: {}'.format(gene_name))
            gene_asm_path = self.gene_list_path(gene_name)
            with open(gene_asm_path) as gene_asm_file:
                asms = [asm.strip() for asm in gene_asm_file
                        if self.taxon_from_asm(gene_name, taxa, asm) in taxa]
                for asm in sorted(asms):
                    taxon = self.taxon_from_asm(gene_name, taxa, asm)
                    self.exonerate_cmds(ref_gene, gene_name, taxon, asm)

    def exonerate_cmds(self, ref_file, gene_name, taxon, asm):
        """Build and run the exonerate program."""
        results_file = self.results_file_path(gene_name)
        cmd = (
            'exonerate --verbose 0 --model protein2genome {ref_file} {asm} '
            '--showvulgar no --showalignment no '
            r'--ryo "{gene},{taxon},%ql,%qal,%qab,%qae,%ti\n" '
            '>> {results_file};'
        ).format(ref_file=ref_file, asm=asm, gene=gene_name, taxon=taxon,
                 results_file=results_file)
        subprocess.check_call(cmd, shell=True)

        exons_file = self.exons_file_path(gene_name)
        cmd = (
            'exonerate --verbose 0 --model protein2genome {ref_file} {asm} '
            '--showvulgar no --showalignment no '
            r'--ryo ">{taxon},%ti,%qab\n%tcs\n" >> {exons_file}; '
        ).format(ref_file=ref_file, asm=asm, taxon=taxon,
                 exons_file=exons_file)
        subprocess.check_call(cmd, shell=True)

        sorted_file = self.sorted_file_path(gene_name)
        cmd = (
            'LC_ALL=C sort -t, -k 1,1d -k 2,2d -k 5,5d  '
            '{results_file} > {sorted_file};'
        ).format(results_file=results_file, sorted_file=sorted_file)
        subprocess.check_call(cmd, shell=True)

    @property
    def all_ref_genes_path(self):
        """Return the file that holds the list of all reference gene files."""
        return join(self.temp_dir, 'Gene_list.txt')

    def ref_gene_file_path(self, seq_name):
        """Build reference gene file name from the sequence name."""
        return join(self.temp_dir, '{}.reference.fasta'.format(seq_name))

    def gene_list_path(self, gene_name):
        """Build file name with list of assemblies for a reference gene."""
        return join(self.temp_dir, '{}.list.txt'.format(gene_name))

    def results_file_path(self, gene_name):
        """Build the file name for the exonerate results."""
        results_file = join(self.temp_dir, '{}.results.csv'.format(gene_name))
        return abspath(results_file)

    def exons_file_path(self, gene_name):
        """Build the file name for the exonerate exons."""
        exons_file = join(self.temp_dir, '{}.exons.fasta'.format(gene_name))
        return abspath(exons_file)

    def sorted_file_path(self, gene_name):
        """Build the file name for the exonerate sorted results."""
        sorted_file = join(
            self.temp_dir, '{}.results.sorted.csv'.format(gene_name))
        return abspath(sorted_file)

    @staticmethod
    def gene_name_from_ref_gene(ref_file):
        """Get the gene name from the reference gene file path."""
        ref_file = ref_file.strip()
        _, gene_name = split(ref_file)
        return re.sub(r'\.reference\.fasta$', '', gene_name)

    def new_fasta_name(self, old_fasta_name):
        """Get the new fasta assembly file name from the old one."""
        _, file_part = split(old_fasta_name)
        file_part, _ = splitext(file_part)
        return join(self.temp_dir, '{}.ed.fasta'.format(file_part))

    @staticmethod
    def taxon_from_asm(gene_name, taxa, asm):
        """Check if the assembly in the taxa set."""
        _, asm = split(asm)
        asm = re.sub(r'^{}\.'.format(gene_name), '', asm)
        return re.sub(r'\.ed\.fasta\s*$', '', asm)
