"""
Stitch together exons from targeted assemblies.

It uses amino acid targets and DNA assemblies.
"""

import re
import os
from os.path import join, split, splitext
from glob import glob
# import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
# import lib.db as db
import lib.log as log
import lib.util as util


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
            # stitcher.find_fasta_files_for_genes()
            # stitcher.run_exonerate()

    def prepare_fasta_files(self):
        """Remove characters exonerate cannot handle from assembly names."""
        idx = 0
        pattern = join(self.args['assemblies_dir'], '*.fasta')
        for fasta_in_path in sorted(glob(pattern)):
            if os.stat(fasta_in_path).st_size == 0:
                continue

            log.info('Preparing fasta file: {}'.format(fasta_in_path))

            _, file_part = split(fasta_in_path)
            file_part, _ = splitext(file_part)
            fasta_out_path = join(self.temp_dir, '{}.ed.fasta'.format(file_part))

            with open(fasta_in_path) as fasta_in, \
                    open(fasta_out_path, 'w') as fasta_out:
                for seq_name, seq in SimpleFastaParser(fasta_in):
                    idx += 1

                    seq_name = seq_name.strip()
                    seq_name = re.sub(r'[\s,=@+\[\]:!-]+', '_', seq_name)
                    seq_name = '{}_{}'.format(seq_name, idx)

                    util.write_fasta_record(fasta_out, seq_name, seq)

    def prepare_reference_genes(self):
        """Format reference sequences for exonerate."""
        ref_genes = self.args['reference_genes']

        log.info('Preparing reference genes: {}'.format(ref_genes))

        with open(ref_genes) as ref_in, \
                open(self.gene_list_path, 'w') as ref_list:
            for seq_name, seq in SimpleFastaParser(ref_in):
                seq_name = re.sub(r'\s+', '_', seq_name)

                ref_file = self.ref_file(seq_name)

                ref_list.write(ref_file)
                ref_list.write('\n')

                with open(ref_file, 'w') as ref_out:
                    util.write_fasta_record(ref_out, seq_name, seq)

    def find_fasta_files_for_genes(self):
        """Find genes in the gene list."""
        with open(self.gene_list_path) as gene_list_files:
            for gene_list_file in gene_list_files:
                gene_list_file = gene_list_file.strip()
                _, gene = split(gene_list_file)
                gene = gene.split('.')[0]
                log.info('Finding fasta files for gene {}'.format(gene))
                with open(gene_list_file, 'w') as file_for_gene:
                    print(file_for_gene)
                    pattern = join(self.temp_dir, '*{}*.ed.fasta'.format(gene))
                    fasta_list = glob(pattern)
                    print(fasta_list)
                    for fasta in fasta_list:
                        file_for_gene.write(fasta)
                        file_for_gene.write('\n')

    # def run_exonerate(self):
    #     """
    #     Run exonerate on every reference sequence, taxon combination.
    #
    #     Get exonerate information for each exon.
    #     Then sort the results file by taxon and by location of the exon.
    #     Use exonerate to pull out the exons.
    #     """
    #     taxa_list = self.args['taxa_list']
    #     with open(gene_list_path) as gene_file, open(taxa_list) as taxa_file:
    #         for gene in gene_file:
    #             gene = gene.strip()
    #             taxa_file.seek(0)
    #             for taxon in taxa_file:
    #                 taxon = taxon.strip()
    #                 log.info('{} {}'.format(gene, taxon))
    #                 gene_asm_path = join(self.temp_dir, '{}.list.txt'.format(gene))
    #                 with open(gene_asm_path) as gene_asm_file:
    #                     for asm in gene_asm_file:
    #                         asm = asm.strip()
    #                         self.exonerate_cmds(
    #                             gene_list_path, gene, taxon, asm)
    #
    # def exonerate_cmds(self, ref_file, gene, taxon, asm):
    #     """Build and run the exonerate program."""
    #     results_file = join(self.temp_dir, '{}.results.csv'.format(gene))
    #     cmd = (
    #         'exonerate --verbose 0 --model protein2genome {ref_file} {asm} '
    #         '--showvulgar no --showalignment no '
    #         '--ryo "{gene},{taxon},%ql,%qal,%qab,%qae,%ti\n" >> {results_file};'
    #     ).format(ref_file=ref_file, asm=asm, gene=gene, taxon=taxon,
    #              results_file=results_file)
    #     subprocess.check_call(cmd, shell=True)
    #
    #     exons_file = join(self.temp_dir, '{}.exons.csv'.format(gene))
    #     cmd = (
    #         'exonerate --verbose 0 --model protein2genome {ref_file} {asm} '
    #         '--showvulgar no --showalignment no '
    #         '--ryo ">{taxon},%ti,%qab\n%tcs\n" >> {exons_file}; '
    #     ).format(ref_file=ref_file, asm=asm, taxon=taxon, exons_file=exons_file)
    #     subprocess.check_call(cmd, shell=True)
    #
    #     sorted_file = join(self.temp_dir, '{}.results.sorted.csv'.format(gene))
    #     cmd = (
    #         'LC_ALL=C sort -t, -k 1,1d -k 2,2d -k 5,5d  '
    #         '{results_file} > {sorted_file};'
    #     ).format(results_file=results_file, sorted_file=sorted_file)
    #     subprocess.check_call(cmd, shell=True)

    @property
    def gene_list_path(self):
        """Return gene list file."""
        return join(self.temp_dir, 'Gene_list.txt')

    def ref_file(self, seq_name):
        """Build reference file name."""
        return join(self.temp_dir, '{}.reference.fasta'.format(seq_name))
