"""
Stitch together exons from targeted assemblies.

It uses amino acid targets and DNA assemblies from aTRAM.
"""

import os
from os.path import abspath, join, basename
from pathlib import Path
from collections import namedtuple
import csv
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
        self.iteration = 0

    def stitch(self):
        """Stitch the exons together."""
        log.stitcher_setup(self.args.log_file)

        with util.make_temp_dir(
                where=self.args.temp_dir,
                prefix='atram_stitcher_',
                keep=self.args.keep_temp_dir) as temp_dir:
            self.temp_dir = temp_dir

            self.cxn = db.connect(temp_dir, 'atram_stitcher')

            self.create_tables()

            self.get_taxa()
            self.insert_reference_genes()
            self.create_reference_files()

            # First iteration gets its data from the input
            self.iteration += 1
            self.get_contigs_from_fasta()
            self.contig_file_write()
            self.run_exonerate()

            # Iterations 2-N get their data from the previous stitch
            for i in range(1, self.args.iterations):
                self.stitch_everything()

                self.iteration += 1
                self.get_contigs_from_previous_stitch()
                self.contig_file_write()
                self.run_exonerate()

            # The final stitch is pickier than the others
            self.stitch_with_gaps()

            log.info('Writing output')
            self.output_stitched_genes()
            self.output_summary_per_gene()
            self.output_summary_per_taxon()

            log.info('Finished.')

    def create_tables(self):
        """Create database tables."""
        db.create_reference_genes_table(self.cxn)
        db.create_exonerate_table(self.cxn)
        db.create_contigs_table(self.cxn)
        db.create_stitch_table(self.cxn)

    def get_taxa(self):
        """Insert taxa into the database."""
        log.info('Preparing taxa')

        with open(self.args.taxa) as taxa:
            for taxon_name in taxa:
                taxon_name = taxon_name.strip()
                self.taxon_names.append(taxon_name)

        self.taxon_names = sorted(self.taxon_names)

    def insert_reference_genes(self):
        """Prepare reference sequences for exonerate."""
        batch = []

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

    def create_reference_files(self):
        """Create reference gene fasta files for exonerate."""
        log.info('Preparing reference gene files for exonerate')
        for ref in db.select_reference_genes(self.cxn):
            with open(ref['ref_file'], 'w') as ref_file:
                util.write_fasta_record(
                    ref_file, ref['ref_name'], ref['ref_seq'])

    def get_contigs_from_fasta(self):
        """Prepare fasta files for exonerate.

        In this iteration we are getting the contigs from the given fasta
        files.
        """
        log.info('{} contig insert: {}'.format(
            util.as_word(self.iteration), self.args.assemblies_dir))

        batch = []

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
                        'contig_rec': i,
                        'iteration': self.iteration})

        db.insert_contigs(self.cxn, batch)

    def get_contigs_from_previous_stitch(self):
        """Prepare fasta files for exonerate.

        In this iteration we are getting all of the contigs from the first
        stitch and combining them into one long contig sequence.
        """
        log.info('{} contig insert'.format(util.as_word(self.iteration)))

        batch = []

        for ref in db.select_reference_genes(self.cxn):
            ref_name = ref['ref_name']

            for taxon_name in self.taxon_names:

                seqs = []
                contig_name = '{}.{}'.format(ref_name, taxon_name)
                contig_file = '{}.stitched.fasta'.format(contig_name)
                contig_file = abspath(join(self.temp_dir, contig_file))

                for contig in db.select_stitched_contigs(
                        self.cxn,
                        ref_name,
                        taxon_name,
                        iteration=self.iteration - 1):

                    seqs.append(contig['seq'])

                batch.append({
                    'ref_name': ref_name,
                    'taxon_name': taxon_name,
                    'contig_name': contig_name,
                    'contig_seq': ''.join(seqs),
                    'contig_file': contig_file,
                    'contig_rec': 1,
                    'iteration': self.iteration})

        db.insert_contigs(self.cxn, batch)

    def contig_file_write(self):
        """Create contig fasta files for exonerate."""
        log.info('{} contig file write'.format(util.as_word(self.iteration)))

        for contig_path in db.select_contig_files(
                self.cxn, iteration=self.iteration):
            contig_file = contig_path['contig_file']
            with open(contig_file, 'w') as fasta_file:
                for contig in db.select_contigs_in_file(
                        self.cxn, contig_file, iteration=self.iteration):
                    util.write_fasta_record(
                        fasta_file,
                        contig['contig_name'],
                        contig['contig_seq'])

    def run_exonerate(self):
        """Run exonerate on every reference sequence, taxon combination."""
        for ref in db.select_reference_genes(self.cxn):
            log.info('{} exonerate run for: {}'.format(
                util.as_word(self.iteration), ref['ref_name']))

            results_file = abspath(join(
                    self.temp_dir,
                    '{}.iteration_{}.results.fasta'.format(
                        ref['ref_name'], self.iteration)))

            Path(results_file).touch()

            for contig_file in db.select_contigs(
                    self.cxn, ref['ref_name'], iteration=self.iteration):

                if util.fasta_file_is_empty(contig_file['contig_file']):
                    continue

                self.exonerate_command(ref, contig_file, results_file)

            self.insert_exonerate_results(results_file)

    def exonerate_command(self, ref, contig_file, results_file):
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

        log.subcommand(cmd, self.temp_dir)

    def insert_exonerate_results(self, results_file):
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
                    'iteration': self.iteration,
                    'seq': seq}
                batch.append(result)

        db.insert_exonerate_results(self.cxn, batch)

    def stitch_everything(self):
        """Build one long contig that covers the reference gene.

        Build one long sequence from all of the non-overlapping contigs in the
        exonerate results. We want maximal coverage of the reference gene.
        """
        log.info('{} stitching run'.format(util.as_word(self.iteration)))

        for stitch in db.select_stitch(self.cxn, iteration=self.iteration):
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
                        prev_contig['end'],
                        iteration=self.iteration)

                if not curr_contig:
                    curr_contig = db.select_next(
                        self.cxn,
                        stitch['ref_name'],
                        stitch['taxon_name'],
                        beg=prev_contig['end'] - self.args.overlap,
                        iteration=self.iteration)

                if curr_contig:
                    curr_contig = dict(curr_contig)
                    position += 1
                    contigs.append({
                        'ref_name': stitch['ref_name'],
                        'taxon_name': stitch['taxon_name'],
                        'contig_name': curr_contig['contig_name'],
                        'position': position,
                        'iteration': self.iteration,
                        'seq': curr_contig['seq']})

                prev_contig = curr_contig

            db.insert_stitched_genes(self.cxn, contigs)

    def stitch_with_gaps(self):
        """Choose the contigs to cover the reference gene.

        In this stitching iteration we are assembling the exon with gaps.
        """
        for ref in db.select_reference_genes(self.cxn):
            ref_name = ref['ref_name']
            ref_len = len(ref['ref_seq']) * CODON_LEN

            log.info('{} stitching run for: {}'.format(
                    util.as_word(self.iteration), ref_name))

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
                            prev_contig['end'],
                            iteration=self.iteration)
                        if curr_contig:
                            curr_contig = dict(curr_contig)
                            beg = prev_contig['end'] - curr_contig['beg'] - 1
                            seq = curr_contig['seq'][beg * CODON_LEN:]

                    if not curr_contig:
                        curr_contig = db.select_next(
                            self.cxn,
                            ref_name,
                            taxon_name,
                            beg=prev_contig['end'] - self.args.overlap,
                            iteration=self.iteration)
                        if curr_contig:
                            curr_contig = dict(curr_contig)
                            seq = curr_contig['seq']
                            gap = curr_contig['beg'] - 1
                            gap -= max(-1, prev_contig['end'])
                            if gap > 0:
                                position += 1
                                contigs.append({
                                    'ref_name': ref_name,
                                    'taxon_name': taxon_name,
                                    'contig_name': None,
                                    'position': position,
                                    'iteration': self.iteration,
                                    'seq': 'N' * (gap * CODON_LEN)})

                    if curr_contig:
                        position += 1
                        contigs.append({
                            'ref_name': ref_name,
                            'taxon_name': taxon_name,
                            'contig_name': curr_contig['contig_name'],
                            'position': position,
                            'iteration': self.iteration,
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
                        'iteration': self.iteration,
                        'seq': 'N' * missing})
                db.insert_stitched_genes(self.cxn, contigs)

    def output_stitched_genes(self):
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

                    hits = db.select_stitched_contig_count(
                        self.cxn,
                        ref_name,
                        taxon_name,
                        iteration=self.iteration)

                    has_hits = 'hits' if hits else 'no hits'
                    log.info('{} {} has {}'.format(
                        ref_name, taxon_name, has_hits))

                    if not hits:
                        continue

                    for contig in db.select_stitched_contigs(
                            self.cxn,
                            ref_name,
                            taxon_name,
                            iteration=self.iteration):

                        contig_name = contig['contig_name']

                        if contig_name:
                            if not first_contig_name:
                                first_contig_name = contig_name

                            contig_names.add(contig_name)

                        seqs.append(contig['seq'])

                    header = '{}. {}.{}'.format(
                            taxon_name, first_contig_name, len(contig_names))

                    util.write_fasta_record(out_file, header, ''.join(seqs))

    def output_summary_per_gene(self):
        """Print per gene summary statistics."""
        out_path = '{}summary_stats_per_gene.csv'.format(
                self.args.output_prefix)
        with open(out_path, 'w') as out_file:
            writer = csv.writer(out_file)
            writer.writerow(
                    ['Locus', 'Taxon', 'Query_Length', 'Target_Length'])

            for stats in db.select_per_gene_stats(self.cxn):
                writer.writerow([
                    stats['ref_name'],
                    stats['taxon_name'],
                    stats['query_len'],
                    stats['target_len']])

    def output_summary_per_taxon(self):
        """Print per taxon summary statistics."""
        out_path = '{}summary_stats_per_taxon.csv'.format(
                self.args.output_prefix)
        with open(out_path, 'w') as out_file:
            writer = csv.writer(out_file)
            writer.writerow([
                'Taxon',
                'Number_Genes',
                'Full_Exons',
                '95%', '90%', '80%', '70%', '50%', '10%'])

            for stats in db.select_per_taxon_stats(self.cxn):
                writer.writerow([
                    stats['taxon_name'],
                    stats['genes'],
                    stats['eq100'],
                    stats['gt95'],
                    stats['ge90'],
                    stats['ge80'],
                    stats['ge70'],
                    stats['ge50'],
                    stats['ge10']])
