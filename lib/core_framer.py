"""Put exons into the correct reading frames."""

import lib.bio as bio
import lib.exonerate as exonerate
import lib.db_stitcher as db
import lib.log as log
import lib.util as util


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
            output_contigs(args, cxn, taxon_names)

            log.info('Writing output')

        log.info('Finished.')


def output_contigs(args, cxn, taxon_names):
    """Add NNNs to align the contigs to the reference sequence."""
    log.info('Framing contigs')

    for taxon in taxon_names:
        out_path = util.prefix_file(
            args.output_prefix, '{}.fasta'.format(taxon))

        with open(out_path, 'w') as out_file:

            for contig in db.select_exonerate_taxa(
                    cxn, taxon, args.min_length):
                ref_len = len(contig['ref_seq']) * bio.CODON_LEN
                seq = 'N' * (contig['beg'] * bio.CODON_LEN)
                seq += contig["seq"]
                seq += 'N' * (ref_len - len(seq))
                util.write_fasta_record(out_file, contig['contig_name'], seq)
