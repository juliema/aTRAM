"""Put exons into the correct reading frames."""

import lib.exonerate as exonerate
import lib.stitcher_db as db
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
            frame_contigs(cxn)

            log.info('Writing output')

        log.info('Finished.')


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
