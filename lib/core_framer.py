"""Put exons into the correct reading frames."""

import csv
from collections import defaultdict
from itertools import product
from . import bio
from . import exonerate
from . import db_stitcher as db
from . import log
from . import util


def frame(args):
    """Frame the exons."""
    log.stitcher_setup(args.log_file, args.log_level)
    iteration = 0

    with util.make_temp_dir(
            where=args.temp_dir,
            prefix='atram_framer_',
            keep=args.keep_temp_dir) as temp_dir:
        with db.connect(temp_dir, 'atram_framer') as cxn:
            cxn.row_factory = lambda c, r: {
                col[0]: r[idx] for idx, col in enumerate(c.description)}
            exonerate.create_tables(cxn)

            taxon_names = exonerate.get_taxa(args)
            exonerate.insert_reference_genes(args, temp_dir, cxn)
            exonerate.check_file_counts(args, cxn, taxon_names)
            exonerate.create_reference_files(cxn)

            iteration += 1
            exonerate.get_contigs_from_fasta(
                args, temp_dir, cxn, taxon_names, iteration)
            exonerate.contig_file_write(cxn)
            exonerate.run_exonerate(temp_dir, cxn, iteration)

            output_contigs(args, cxn)

            log.info('Writing output')
            output_summary_per_gene(args, cxn, taxon_names)
            output_summary_per_taxon(args, cxn, taxon_names)

        log.info('Finished')


def output_contigs(args, cxn):
    """Add NNNs to align the contigs to the reference sequence."""
    log.info('Framing contigs')

    for ref in db.select_reference_genes(cxn):
        ref_name = ref['ref_name']
        ref_len = len(ref['ref_seq']) * bio.CODON_LEN

        names_seen = defaultdict(int)

        out_path = util.prefix_file(
            args.output_prefix, '{}.fasta'.format(ref_name))

        with open(out_path, 'w') as out_file:

            for contig in db.select_exonerate_ref_gene(
                    cxn, ref_name, args.min_length):

                contig_name = exonerate.handle_duplicate_name(
                    contig['contig_name'], names_seen)

                seq = 'N' * (contig['beg'] * bio.CODON_LEN)
                seq += contig['seq']
                seq += 'N' * (ref_len - len(seq))
                util.write_fasta_record(out_file, contig_name, seq)


def output_summary_per_gene(args, cxn, taxon_names):
    """Print per gene summary statistics."""
    longest = max(db.select_longest(cxn), 1)
    lengths = db.select_seq_lengths(cxn)

    counts = {t: {'total': set(), 'long': set()} for t in taxon_names}

    for length in lengths:
        taxon_name = length['taxon_name']
        ref_name = length['ref_name']
        counts[taxon_name]['total'].add(ref_name)
        fraction = length['len'] / longest
        if fraction >= args.long_contig:
            counts[taxon_name]['long'].add(ref_name)

    out_path = util.prefix_file(
        args.output_prefix, 'summary_stats_per_ref_gene.csv')
    with open(out_path, 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(['Taxon',
                         'Total_Genes',
                         'Total_Genes_>={:0.2}'.format(args.long_contig)])
        for taxon, count in counts.items():
            writer.writerow([taxon, len(count['total']), len(count['long'])])


def output_summary_per_taxon(args, cxn, taxon_names):
    """Print per taxon summary statistics."""
    longest = max(db.select_longest(cxn), 1)
    lengths = db.select_seq_lengths(cxn)
    ref_names = [r['ref_name'] for r in db.select_reference_genes(cxn)]

    counts = {c: {'total': 0, 'long': 0}
              for c in product(taxon_names, ref_names)}

    for length in lengths:
        taxon_name = length['taxon_name']
        ref_name = length['ref_name']
        key = (taxon_name, ref_name)
        counts[key]['total'] += 1
        fraction = length['len'] / longest
        if fraction >= args.long_contig:
            counts[key]['long'] += 1

    out_path = util.prefix_file(
        args.output_prefix, 'summary_stats_per_taxon.csv')
    with open(out_path, 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(['Taxon',
                         'Gene',
                         'Total_Contigs',
                         'Total_Contigs_>={:0.2}'.format(args.long_contig)])
        for key, count in counts.items():
            writer.writerow([key[0], key[1], count['total'], count['long']])
