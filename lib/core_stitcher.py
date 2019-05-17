"""
Stitch together exons from targeted assemblies.

It uses amino acid targets and DNA assemblies from aTRAM.
"""

import os
from os.path import abspath, join, basename
from collections import namedtuple
# import csv
import subprocess
from glob import glob
from Bio.SeqIO.FastaIO import SimpleFastaParser
import lib.stitcher_db as db
import lib.log as log
import lib.util as util


CODON_LEN = 3


class Sticher:
    """Container for sticher data."""

    def __init__(self, args):
        """Setup."""
        self.args = args
        self.cxn = None
        self.temp_dir = None
        self.taxon_names = []

    def stitch(self):
        """Stitch the exons together."""
        log.stitcher_setup(self.args.log_file)

        with util.make_temp_dir(
                where=self.args.temp_dir,
                prefix='atram_stitcher_',
                keep=self.args.keep_temp_dir) as temp_dir:
            self.temp_dir = temp_dir

            self.cxn = db.connect(temp_dir, 'atram_stitcher')

            self.insert_taxa()
            self.insert_reference_genes()
            self.insert_contigs()
            self.create_reference_files()
            self.create_contig_files(run=1)
            self.first_exonerate_run()
            self.first_contig_stitch()
            # self.output_results()

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
        self.taxon_names = sorted([x['taxon_name'] for x in batch])

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

                batch.append({
                    'ref_name': ref_name,
                    'ref_seq': ref_seq,
                    'ref_file': ref_file})

        db.insert_reference_genes(self.cxn, batch)

    def insert_contigs(self):
        """Prepare fasta files for exonerate."""
        log.info('Preparing fasta files: {}'.format(self.args.assemblies_dir))

        batch = []

        db.create_contigs_table(self.cxn)

        ref_names = set(x['ref_name']
                        for x in db.select_reference_genes(self.cxn))

        pattern = join(self.args.assemblies_dir, '*.fasta')
        for contig_path in sorted(glob(pattern)):

            if os.stat(contig_path).st_size == 0:
                continue

            tiebreaker = 0

            with open(contig_path) as contig_old:
                for i, (contig_name, contig_seq) \
                        in enumerate(SimpleFastaParser(contig_old)):

                    contig_file = basename(contig_path)
                    ref_name, taxon_name = contig_file.split('.')[0:2]

                    if ref_name not in ref_names \
                            or taxon_name not in self.taxon_names:
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
                    ref_file, ref['ref_name'], ref['ref_seq'])

    def create_contig_files(self, run=1):
        """Create contig fasta files for exonerate."""
        log.info('Preparing contig files for exonerate')

        for contig_path in db.select_contig_files(self.cxn):
            contig_file = contig_path['contig_file']
            with open(contig_file, 'w') as fasta_file:
                for contig in db.select_contigs_in_file(self.cxn, contig_file):
                    util.write_fasta_record(
                        fasta_file,
                        contig['contig_name'],
                        contig['contig_seq'])

    def first_exonerate_run(self):
        """Run exonerate on every reference sequence, taxon combination."""
        db.create_exonerate_table(self.cxn)

        for ref in db.select_reference_genes(self.cxn):
            log.info('First exonerate run for: {}'.format(ref['ref_name']))

            results_file = abspath(join(
                    self.temp_dir,
                    '{}.results.fasta'.format(ref['ref_name'])))

            for contig_file in db.select_first_exonerate_run(
                    self.cxn, ref['ref_name']):
                self.exonerate_command(ref, contig_file, results_file)

            self.insert_first_exonerate_results(results_file)

    @staticmethod
    def exonerate_command(ref, contig_file, results_file):
        """Build and run the exonerate program."""
        cmd = util.shorten(r"""
            exonerate --verbose 0 --model protein2genome {ref_file}
            {contig_file}
            --showvulgar no --showalignment no
            --ryo ">{ref_name},{taxon_name},%ti,%qab,%qae\n%tcs\n"
            >> {results_file};""").format(
                ref_file=ref['ref_file'],
                contig_file=contig_file['contig_file'],
                ref_name=ref['ref_name'],
                taxon_name=contig_file['taxon_name'],
                results_file=results_file)

        subprocess.check_call(cmd, shell=True)

    def insert_first_exonerate_results(self, results_file):
        """Insert the exonerate results into the database."""
        ExonerateHeader = namedtuple(
                'ExonerateHeader',
                ['ref_name', 'taxon_name', 'contig_name', 'beg', 'end'])

        batch = []
        with open(results_file) as results_fasta:
            for header, seq in SimpleFastaParser(results_fasta):
                header = header.split(',')
                field = ExonerateHeader(*header)
                result = {
                    'ref_name': field.ref_name,
                    'taxon_name': field.taxon_name,
                    'contig_name': field.contig_name,
                    'beg': field.beg,
                    'end': field.end,
                    'seq': seq}
                batch.append(result)

        db.insert_exonerate_results(self.cxn, batch)

    def first_contig_stitch(self):
        """Build one long contig that covers the reference gene.

        Build one long sequence from all of the non-overlapping contigs in the
        exonerate results. We want maximal coverage of the reference gene.
        """
        log.info('First stitch run')

        db.create_stitch_table(self.cxn)

        for stitch in db.select_stitch(self.cxn):
            contigs = []
            position = 0
            prev_contig = {'end': -1}
            curr_contig = None

            while prev_contig:

                if prev_contig['end'] > 0:
                    curr_contig = db.select_overlap(
                        self.cxn,
                        stitch['ref_name'],
                        stitch['taxon_name'],
                        prev_contig['beg'] + 1,
                        prev_contig['end'] - self.args.overlap,
                        prev_contig['end'])

                if not curr_contig:
                    curr_contig = db.select_next(
                        self.cxn,
                        stitch['ref_name'],
                        stitch['taxon_name'],
                        beg=prev_contig['end'])

                if curr_contig:
                    curr_contig = dict(curr_contig)
                    position += 1
                    contigs.append({
                        'ref_name': stitch['ref_name'],
                        'taxon_name': stitch['taxon_name'],
                        'contig_name': curr_contig['contig_name'],
                        'position': position,
                        'seq': curr_contig['seq']})

                prev_contig = curr_contig

            db.insert_stitched_genes(self.cxn, contigs)

    def second_exonerate_run(self):
        """Run exonerate on the stitched genes from the first run."""

    def second_contig_stitch(self):
        """Choose the contigs to cover the reference gene."""
        db.create_stitch_table(self.cxn)

        for ref in db.select_reference_genes(self.cxn):
            ref_name = ref['ref_name']
            ref_len = len(ref['ref_seq']) * CODON_LEN

            log.info('Second stitching run for: {}'.format(ref_name))

            for taxon_name in self.taxon_names:

                contigs = []
                position = 0
                prev_contig = {'end': -1}
                curr_contig = None
                first_time = True

                while prev_contig:
                    seq = ''

                    if not first_time:
                        curr_contig = db.select_overlap(
                            self.cxn,
                            ref_name,
                            taxon_name,
                            prev_contig['beg'] + 1,
                            prev_contig['end'] - self.args.overlap,
                            prev_contig['end'])
                        if curr_contig:
                            curr_contig = dict(curr_contig)
                            prev_end = prev_contig['end']
                            beg = (prev_end - curr_contig['beg']) * CODON_LEN
                            seq = curr_contig['seq'][beg:]

                    if not curr_contig:
                        curr_contig = db.select_next(
                            self.cxn,
                            ref_name,
                            taxon_name,
                            beg=prev_contig['end'])
                        if curr_contig:
                            curr_contig = dict(curr_contig)
                            seq = curr_contig['seq']
                            gap = curr_contig['beg'] - prev_contig['end'] - 1
                            # gap *= CODON_LEN
                            if not first_time and gap > 0:
                                position += 1
                                contigs.append({
                                    'ref_name': ref_name,
                                    'taxon_name': taxon_name,
                                    'contig_name': None,
                                    'position': position,
                                    'seq': 'N' * (gap * CODON_LEN)})

                    if curr_contig:
                        position += 1
                        contigs.append({
                            'ref_name': ref_name,
                            'taxon_name': taxon_name,
                            'contig_name': curr_contig['contig_name'],
                            'position': position,
                            'seq': seq})

                    prev_contig = curr_contig
                    first_time = False

                # Final NNNs
                stitch_len = sum(len(x['seq']) for x in contigs)
                missing = ref_len - stitch_len
                if missing > 0:
                    position += 1
                    contigs.append({
                        'ref_name': ref_name,
                        'taxon_name': taxon_name,
                        'contig_name': None,
                        'position': position,
                        'seq': 'N' * missing})
                db.insert_stitched_genes(self.cxn, contigs)

    def output_results(self):
        """Print results."""

        for ref in db.select_reference_genes(self.cxn):
            ref_name = ref['ref_name']

            out_path = '{}{}.stitched_exons.fasta'.format(
                    self.args.output_prefix, ref_name)

            with open(out_path, 'w') as out_file:

                for taxon_name in self.taxon_names:
                    contig_names = set()
                    first_contig_name = None
                    seqs = []

                    for contig in db.select_stitched_contigs(
                            self.cxn, ref_name, taxon_name):

                        contig_name = contig['contig_name']

                        if contig_name:
                            if not first_contig_name:
                                first_contig_name = contig_name

                            contig_names.add(contig_name)

                        seqs.append(contig['seq'])

                    header = '{}. {}.{}'.format(
                            taxon_name, first_contig_name, len(contig_names))

                    util.write_fasta_record(out_file, header, ''.join(seqs))
