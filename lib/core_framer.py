"""Put exons into the correct reading frames."""

import re
import os
from os.path import abspath, join, basename
from collections import defaultdict
from glob import glob
from Bio.SeqIO.FastaIO import SimpleFastaParser
import lib.exonerate as exonerate
import lib.stitcher_db as db
import lib.log as log
import lib.util as util


CODON_LEN = 3


def frame(args):
    """Frame the exons."""
    log.stitcher_setup(args.log_file)
    iteration = 0

    with util.make_temp_dir(
            where=args.temp_dir,
            prefix='atram_framer_',
            keep=args.keep_temp_dir) as temp_dir:
        with db.connect(temp_dir, 'atram_framer') as cxn:
            cxn.row_factory = lambda c, r: {
                col[0]: r[idx] for idx, col in enumerate(c.description)}
            create_tables(cxn)

            taxon_names = get_taxa(args)
            insert_reference_genes(args, temp_dir, cxn)
            check_file_counts(args, cxn, taxon_names)
            create_reference_files(cxn)

            iteration += 1
            get_contigs_from_fasta(args, temp_dir, cxn, taxon_names, iteration)
            contig_file_write(cxn)
            exonerate.run_exonerate(temp_dir, cxn, iteration)
            frame_contigs(cxn)

            log.info('Writing output')

        log.info('Finished.')


def create_tables(cxn):
    """Create database tables."""
    db.create_reference_genes_table(cxn)
    db.create_exonerate_table(cxn)
    db.create_contigs_table(cxn)
    db.create_stitch_table(cxn)


def get_taxa(args):
    """Insert taxa into the database."""
    log.info('Preparing taxa')

    taxon_names = []

    with open(args.taxa) as taxa:
        for taxon_name in taxa:
            taxon_name = taxon_name.strip()
            taxon_names.append(taxon_name)

    return sorted(taxon_names)


def insert_reference_genes(args, temp_dir, cxn):
    """Prepare reference sequences for exonerate."""
    batch = []

    ref_genes = args.reference_genes

    log.info('Preparing reference genes: {}'.format(ref_genes))

    with open(ref_genes) as ref_in:

        for ref_name, ref_seq in SimpleFastaParser(ref_in):

            ref_name = util.clean_name(ref_name)

            ref_file = abspath(join(
                temp_dir,
                '{}.fasta'.format(ref_name)))

            batch.append({
                'ref_name': ref_name,
                'ref_seq': ref_seq,
                'ref_file': ref_file})

    db.insert_reference_genes(cxn, batch)


def create_reference_files(cxn):
    """Create reference gene fasta files for exonerate."""
    log.info('Preparing reference gene files for exonerate')
    for ref in db.select_reference_genes(cxn):
        with open(ref['ref_file'], 'w') as ref_file:
            util.write_fasta_record(ref_file, ref['ref_name'], ref['ref_seq'])


def check_file_counts(args, cxn, taxon_names):
    """Only one contig file may match a reference/taxon pair."""
    ref_names = set(x['ref_name'] for x in db.select_reference_genes(cxn))

    pattern = join(args.assemblies_dir, args.file_filter)
    contig_files = sorted(glob(pattern))

    counts = defaultdict(list)
    for contig_file in contig_files:
        ref_name, taxon_name = parse_contig_file_name(
            ref_names, taxon_names, contig_file)
        if not ref_name or not taxon_name:
            continue
        counts[(ref_name, taxon_name)].append(contig_file)

    counts = {k: v for k, v in counts.items() if len(v) > 1}

    if not counts:
        return

    msg = []
    for key, contig_files in counts.items():
        msg.append(
            'Multiple files match reference {} and taxon {}:'.format(*key))
        msg += contig_files

    log.fatal('\n'.join(msg))


def parse_contig_file_name(ref_names, taxon_names, contig_file):
    """Extract the reference & taxon names from the contig file name."""
    sep = r'[_. ]'

    ref_names = sorted(ref_names, key=lambda x: (len(x), x), reverse=True)
    ref_names = [x + sep for x in ref_names]

    taxon_names = sorted(taxon_names, key=lambda x: (len(x), x), reverse=True)
    taxon_names = [x + sep for x in taxon_names]

    ref_name = [x[:-len(sep)] for x in ref_names if re.search(x, contig_file)]
    taxon_name = [x[:-len(sep)] for x in taxon_names
                  if re.search(x, contig_file)]

    ref_name += [None]
    taxon_name += [None]

    return ref_name[0], taxon_name[0]


def get_contigs_from_fasta(args, temp_dir, cxn, taxon_names, iteration):
    """Prepare fasta files for exonerate.

    In this iteration we are getting the contigs from the given fasta files.
    """
    log.info('{} contig insert: {}'.format(
        util.as_word(iteration), args.assemblies_dir))

    batch = []

    ref_names = set(x['ref_name'] for x in db.select_reference_genes(cxn))

    pattern = join(args.assemblies_dir, args.file_filter)
    for contig_path in sorted(glob(pattern)):

        if os.stat(contig_path).st_size == 0:
            continue

        tiebreaker = 0

        with open(contig_path) as contig_old:
            for i, (contig_name, contig_seq) \
                    in enumerate(SimpleFastaParser(contig_old)):

                contig_file = basename(contig_path)
                ref_name, taxon_name = parse_contig_file_name(
                    ref_names, taxon_names, contig_file)

                if ref_name not in ref_names or taxon_name not in taxon_names:
                    continue

                tiebreaker += 1
                contig_name = '{}@{}_{}_{}'.format(
                    ref_name, taxon_name, contig_name, tiebreaker)
                contig_name = re.sub(r'^[\w@]+', '_', contig_name.strip())

                contig_file = abspath(join(temp_dir, contig_name + '.fasta'))

                batch.append({
                    'ref_name': ref_name,
                    'taxon_name': taxon_name,
                    'contig_name': contig_name,
                    'contig_seq': contig_seq,
                    'contig_file': contig_file,
                    'contig_rec': i,
                    'iteration': iteration})

    db.insert_contigs(cxn, batch)


def contig_file_write(cxn):
    """Create contig fasta files for exonerate."""
    log.info('Write contig files')

    for contig in db.select_all_contigs(cxn):
        with open(contig['contig_file'], 'w') as fasta_file:
            util.write_fasta_record(
                fasta_file,
                contig['contig_name'],
                contig['contig_seq'])


def frame_contigs(cxn):
    """Add NNNs to align the contigs to the reference sequence."""


def output_framed_contigs(args, cxn):
    """Print results."""
    for ref in db.select_reference_genes(cxn):
        ref_name = ref['ref_name']

        out_path = '{}.{}.framed_exons.fasta'.format(
            args.output_prefix, ref_name)

        with open(out_path, 'w') as out_file:
            pass

            # hits = db.select_stitched_contig_count(
            #     cxn,
            #     ref_name,
            #     taxon_name,
            #     iteration=iteration)
            #
            #     # has_hits = 'hits' if hits else 'no hits'
            #     # log.info('{}.{} has {}'.format(
            #     #     ref_name, taxon_name, has_hits))
            #
            #     if not hits:
            #         continue
            #
            #     for contig in db.select_stitched_contigs(
            #             cxn, ref_name, taxon_name, iteration=iteration):
            #
            #         seqs.append(contig['seq'])
            #
            #     header = taxon_name
            #     if args.reference_name:
            #         header = '{}.{}'.format(ref_name, taxon_name)
            #
            #     util.write_fasta_record(out_file, header, ''.join(seqs))
