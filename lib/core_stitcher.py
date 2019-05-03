"""
Stitch together exons from targeted assemblies.

It uses amino acid targets and DNA assemblies from aTRAM.
"""

import os
from os.path import abspath, join, basename
from collections import namedtuple, OrderedDict
import csv
import subprocess
from glob import glob
from textwrap import shorten
from Bio.SeqIO.FastaIO import SimpleFastaParser
import lib.db as db
import lib.log as log
import lib.util as util


class Sticher:
    """Container for sticher data."""

    def __init__(self, args):
        """Setup."""
        self.args = args
        self.cxn = None
        self.temp_dir = None

    def stitch(self):
        """Stitch the exons together."""
        log.stitcher_setup(self.args.log_file)

        with util.make_temp_dir(
                where=self.args.temp_dir,
                prefix='atram_stitcher_',
                keep=self.args.keep_temp_dir) as temp_dir:
            self.temp_dir = temp_dir

            self.cxn = db.temp_db(temp_dir, 'atram_stitcher')

            self.insert_taxa()
            self.insert_reference_genes()
            self.insert_contigs()
            self.create_reference_files()
            self.create_contig_files()
            self.run_exonerate()
            self.get_contigs()
            # self.stitch_contigs()
            # self.summary_stats()

    def insert_taxa(self):
        """Insert taxa into the database."""
        batch = []

        db.create_taxa_table(self.cxn)

        log.info('Preparing taxa')

        with open(self.args.taxa) as taxa:
            for taxon_name in taxa:
                taxon_name = taxon_name.strip()
                batch.append({'taxon_name': taxon_name})

        db.insert_taxa(self.cxn, batch)

    def insert_reference_genes(self):
        """Prepare reference sequences for exonerate."""
        batch = []

        db.create_reference_genes_table(self.cxn)

        ref_genes = self.args.reference_genes

        log.info('Preparing reference genes: {}'.format(ref_genes))

        with open(ref_genes) as ref_in:

            for ref_name, ref_seq in SimpleFastaParser(ref_in):

                ref_name = util.clean_name(ref_name)

                ref_file = abspath(join(
                    self.temp_dir,
                    '{}.fasta'.format(ref_name)))

                results_file = abspath(join(
                    self.temp_dir,
                    '{}.results.fasta'.format(ref_name)))

                batch.append({
                    'ref_name': ref_name,
                    'ref_seq': ref_seq,
                    'ref_file': ref_file,
                    'results_file': results_file,
                    'input_file': ref_genes})

        db.insert_reference_genes(self.cxn, batch)

    def insert_contigs(self):
        """Prepare fasta files for exonerate."""
        batch = []

        db.create_contigs_table(self.cxn)

        ref_names = set(x['ref_name']
                        for x in db.select_reference_genes(self.cxn))

        taxon_names = set(x['taxon_name'] for x in db.select_taxa(self.cxn))

        pattern = join(self.args.assemblies_dir, '*.fasta')
        for contig_path in sorted(glob(pattern)):

            if os.stat(contig_path).st_size == 0:
                continue

            log.info('Preparing fasta file: {}'.format(contig_path))

            tiebreaker = 0

            with open(contig_path) as contig_old:
                for i, (contig_name, contig_seq) \
                        in enumerate(SimpleFastaParser(contig_old)):

                    contig_file = basename(contig_path)
                    ref_name, taxon_name = contig_file.split('.')[0:2]

                    if ref_name not in ref_names \
                            or taxon_name not in taxon_names:
                        continue

                    tiebreaker += 1
                    contig_name = util.clean_name(contig_name)
                    contig_name = '{}_{}'.format(contig_name, tiebreaker)

                    contig_file = abspath(join(self.temp_dir, contig_file))

                    batch.append({
                        'ref_name': ref_name,
                        'taxon_name': taxon_name,
                        'contig_name': contig_name,
                        'contig_seq': contig_seq,
                        'contig_file': contig_file,
                        'contig_rec': i})

        db.insert_contigs(self.cxn, batch)

    def create_reference_files(self):
        """Create reference gene fasta files for exonerate."""
        log.info('Preparing reference gene files for exonerate')
        for ref in db.select_reference_genes(self.cxn):
            with open(ref['ref_file'], 'w') as ref_file:
                util.write_fasta_record(
                    ref_file,
                    ref['ref_name'],
                    ref['ref_seq'])

    def create_contig_files(self):
        """Create contig fasta files for exonerate."""
        log.info('Preparing contig files for exonerate')

        for ref in db.select_reference_genes(self.cxn):
            ref_name = ref['ref_name']
            for contig_path in \
                    db.select_conting_files(self.cxn, ref_name):
                contig_file = contig_path['contig_file']
                with open(contig_file, 'w') as fasta_file:
                    for contig in db.select_contings_in_file(
                            self.cxn, ref_name, contig_file):
                        util.write_fasta_record(
                            fasta_file,
                            contig['contig_name'],
                            contig['contig_seq'])

    def run_exonerate(self):
        """Run exonerate on every reference sequence, taxon combination."""
        db.create_exonerate_table(self.cxn)

        for ref in db.select_reference_genes(self.cxn):
            log.info('Running exonerate for: {}'.format(ref['ref_name']))

            for contig in db.select_conting_files(self.cxn, ref['ref_name']):
                self.exonerate_command(ref, contig)
                self.insert_exonerate_results(ref)

    @staticmethod
    def exonerate_command(ref, contig):
        """Build and run the exonerate program."""
        cmd = shorten(r"""
            exonerate --verbose 0 --model protein2genome {ref_file}
            {contig_file}
            --showvulgar no --showalignment no
            --ryo ">{ref_name},{taxon_name},%ti,%ql,%qal,%qab,%qae\n%tcs\n"
            >> {results_file};""", 9999).format(
                ref_file=ref['ref_file'],
                contig_file=contig['contig_file'],
                ref_name=ref['ref_name'],
                taxon_name=contig['taxon_name'],
                results_file=ref['results_file'])

        subprocess.check_call(cmd, shell=True)

    def insert_exonerate_results(self, ref):
        """Insert the exonerate results into the database."""
        ExonerateHeader = namedtuple(
            'ExonerateHeader',
            ['ref_name', 'taxon_name', 'contig_name', 'query_len',
             'align_len', 'align_beg', 'align_end'])

        batch = []
        with open(ref['results_file']) as results_fasta:
            for header, target_seq in SimpleFastaParser(results_fasta):
                header = header.split(',')
                field = ExonerateHeader(*header)
                batch.append({
                    'ref_name':    field.ref_name,
                    'taxon_name':  field.taxon_name,
                    'contig_name': field.contig_name,
                    'query_len':   field.query_len,
                    'align_len':   field.align_len,
                    'align_beg':   field.align_beg,
                    'align_end':   field.align_end,
                    'target_seq':  target_seq})
        db.insert_exonerate_results(self.cxn, batch)

    def get_contigs(self):
        """
        Create stats file.

        1) Sort by the beginning of the contig.
        2) Look for all sorted results from exonerate.
        3) Create a CSV file for each input with all of the information about
            the contigs from exonerate.
        """
        for ref in db.select_reference_genes(self.cxn):
            ref_name = ref['ref_name']

            log.info('Creating stats file for: {}'.format(ref_name))

            stats_path = '{}.{}.overlap.{}.contig_list.csv'.format(
                self.args.output_prefix, ref_name, self.args.overlap)

            with open(stats_path, 'w') as stats_file:

                writer = csv.writer(stats_file)

                taxon_count = db.taxon_names_count(self.cxn, ref_name)

                writer.writerow([
                    'Input Gene', ref_name,
                    'In File', ref['input_file'],
                    'Allowing Overlap', self.args.overlap,
                    'Number of Libraries', taxon_count])

                writer.writerow([
                    'Gene', 'Taxon', 'Number of Contigs', 'Gene Length',
                    'Contigs to Keep', 'Total Overlap',
                    'Combined Exon Length', 'Beginning', 'End',
                    'Beginning', 'End', 'Contig Name'])

                for taxon in db.select_exonerate_taxa(self.cxn, ref_name):

                    row = [ref_name,
                           taxon['taxon_name'],
                           taxon['contig_count'],
                           taxon['gene_len']]

                    if taxon['max_cover'] == taxon['gene_len'] \
                            or taxon['contig_count'] == 1:
                        self.use_one_contig(ref, taxon, row)
                    else:
                        self.stitch_multiple_contigs(ref, taxon, row)

                    writer.writerow(row)

    def use_one_contig(self, ref, taxon, row):
        """Use a single contig."""
        top_contig = db.select_top_exonerate_contig(
            self.cxn, ref['ref_name'], taxon['taxon_name'])
        row += [
            1,  # Contigs to keep
            0,  # Total overlap
            top_contig['align_len'],
            top_contig['align_beg'],
            top_contig['align_end'],
            top_contig['contig_name']]

    def stitch_multiple_contigs(self, ref, taxon, row):
        """Stitch multiple contigs together."""
        sum_ = 0
        overlap = 0
        keepers = OrderedDict()

        contigs = db.select_exonerate_stitch(
            self.cxn, ref['ref_name'], taxon['taxon_name'])

        curr = next(contigs)
        for nxt in contigs:
            curr_end = curr['align_end']
            next_beg = nxt['align_beg']
            next_end = nxt['align_end']

            if next_beg > (curr_end - self.args.overlap) \
                    and (next_end > curr_end):

                overlap += (curr['align_end'] - next_beg)
                if not sum_:
                    sum_ = curr['query_len'] + nxt['query_len']
                else:
                    sum_ += nxt['query_len']

                keepers[curr['contig_name']] = curr
                keepers[nxt['contig_name']] = nxt

            elif next_beg <= (curr_end - self.args.overlap):
                pass

            curr = nxt

        row += []

    def stitch_contigs(self):
        """Stitch the contigs together."""
