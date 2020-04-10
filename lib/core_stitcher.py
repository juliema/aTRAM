"""
Stitch together exons from targeted assemblies.

It uses amino acid targets and DNA assemblies from aTRAM.
"""

from os.path import abspath, join
import csv
from . import bio
from . import db_stitcher as db
from . import log
from . import util
from . import exonerate


def stitch(args):
    """Stitch the exons together."""
    log.stitcher_setup(args.log_file, args.log_level)
    iteration = 0

    with util.make_temp_dir(
            where=args.temp_dir,
            prefix='atram_stitcher_',
            keep=args.keep_temp_dir) as temp_dir:
        with db.connect(temp_dir, 'atram_stitcher') as cxn:
            cxn.row_factory = lambda c, r: {
                col[0]: r[idx] for idx, col in enumerate(c.description)}

            exonerate.create_tables(cxn)

            taxon_names = exonerate.get_taxa(args)
            exonerate.insert_reference_genes(args, temp_dir, cxn)
            exonerate.check_file_counts(args, cxn, taxon_names)
            exonerate.create_reference_files(cxn)

            # First iteration gets its data from the input
            iteration += 1
            exonerate.get_contigs_from_fasta(
                args, temp_dir, cxn, taxon_names, iteration)
            exonerate.contig_file_write(cxn)
            exonerate.run_exonerate(temp_dir, cxn, iteration)
            early_exit_check(cxn)

            # Iterations 2-N get their data from the previous stitch
            for _ in range(1, args.iterations):
                stitch_everything(args, cxn, iteration)

                iteration += 1
                get_contigs_from_previous_stitch(
                    temp_dir, cxn, taxon_names, iteration)
                exonerate.contig_file_write(cxn)
                exonerate.run_exonerate(temp_dir, cxn, iteration)

            # The final stitch is pickier than the others
            stitch_with_gaps(args, cxn, taxon_names, iteration)

            log.info('Writing output')
            output_stitched_genes(args, cxn, taxon_names, iteration)
            output_summary_per_gene(args, cxn, iteration)
            output_summary_per_taxon(args, cxn, iteration)

        log.info('Finished')


def get_contigs_from_previous_stitch(temp_dir, cxn, taxon_names, iteration):
    """Prepare fasta files for exonerate.

    In this iteration we are getting all of the contigs from the first stitch
    and combining them into one long contig sequence.
    """
    log.info('{} contig insert'.format(util.as_word(iteration)))

    batch = []

    for ref in db.select_reference_genes(cxn):
        ref_name = ref['ref_name']

        for taxon_name in taxon_names:

            seqs = []
            contig_name = '{}.{}'.format(ref_name, taxon_name)
            contig_file = '{}.stitched.fasta'.format(contig_name)
            contig_file = abspath(join(temp_dir, contig_file))

            for contig in db.select_stitched_contigs(
                    cxn, ref_name, taxon_name, iteration=iteration - 1):

                seqs.append(contig['seq'])

            batch.append({
                'ref_name': ref_name,
                'taxon_name': taxon_name,
                'contig_name': contig_name,
                'contig_seq': ''.join(seqs),
                'contig_file': contig_file,
                'contig_rec': 1,
                'iteration': iteration})

    db.insert_contigs(cxn, batch)


def early_exit_check(cxn):
    """Check for empty results after the first iteration."""
    if db.select_exonerate_count(cxn):
        return
    log.fatal(util.shorten("""
        There were no hits. Please make sure your contig file names contain
        both a reference gene name and taxon name in them."""))


def stitch_everything(args, cxn, iteration):
    """Build one long contig that covers the reference gene.

    Build one long sequence from all of the non-overlapping contigs in the
    exonerate results. We want maximal coverage of the reference gene.
    """
    log.info('{} stitching run'.format(util.as_word(iteration)))

    for thread in db.select_stitch(cxn, iteration=iteration):
        contigs = []
        position = 0
        prev_contig = {'end': -1}
        curr_contig = None

        while prev_contig:

            if prev_contig['end'] > 0:
                curr_contig = db.select_overlap(
                    cxn,
                    thread['ref_name'],
                    thread['taxon_name'],
                    prev_contig['beg'] + 1,
                    prev_contig['end'] - args.overlap,
                    prev_contig['end'],
                    iteration=iteration)

            if not curr_contig:
                curr_contig = db.select_next(
                    cxn,
                    thread['ref_name'],
                    thread['taxon_name'],
                    beg=prev_contig['end'] - args.overlap,
                    iteration=iteration)

            if curr_contig:
                curr_contig = dict(curr_contig)
                position += 1
                contigs.append({
                    'ref_name': thread['ref_name'],
                    'taxon_name': thread['taxon_name'],
                    'contig_name': curr_contig['contig_name'],
                    'position': position,
                    'iteration': iteration,
                    'seq': curr_contig['seq']})

            prev_contig = curr_contig

        db.insert_stitched_genes(cxn, contigs)


def stitch_with_gaps(args, cxn, taxon_names, iteration):
    """
    Choose the contigs to cover the reference gene.

    In this stitching iteration we are assembling the exon with gaps.
    """
    for ref in db.select_reference_genes(cxn):
        ref_name = ref['ref_name']
        ref_len = len(ref['ref_seq']) * bio.CODON_LEN

        log.info('{} stitching run for: {}'.format(
            util.as_word(iteration), ref_name))

        for taxon_name in taxon_names:

            contigs = []
            position = 0
            prev_contig = {'end': -1}
            curr_contig = None
            first_time = True

            while prev_contig:
                seq = ''

                if not first_time:
                    curr_contig = db.select_overlap(
                        cxn,
                        ref_name,
                        taxon_name,
                        prev_contig['beg'] + 1,
                        prev_contig['end'] - args.overlap,
                        prev_contig['end'],
                        iteration=iteration)
                    if curr_contig:
                        curr_contig = dict(curr_contig)
                        beg = prev_contig['end'] - curr_contig['beg'] - 1
                        seq = curr_contig['seq'][beg * bio.CODON_LEN:]

                if not curr_contig:
                    curr_contig = db.select_next(
                        cxn,
                        ref_name,
                        taxon_name,
                        beg=prev_contig['end'] - args.overlap,
                        iteration=iteration)
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
                                'iteration': iteration,
                                'seq': 'N' * (gap * bio.CODON_LEN)})

                if curr_contig:
                    position += 1
                    contigs.append({
                        'ref_name': ref_name,
                        'taxon_name': taxon_name,
                        'contig_name': curr_contig['contig_name'],
                        'position': position,
                        'iteration': iteration,
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
                    'iteration': iteration,
                    'seq': 'N' * missing})
            db.insert_stitched_genes(cxn, contigs)


def output_stitched_genes(args, cxn, taxon_names, iteration):
    """Print results."""
    for ref in db.select_reference_genes(cxn):
        ref_name = ref['ref_name']

        out_path = util.prefix_file(
            args.output_prefix, '{}.stitched_exons.fasta'.format(ref_name))

        with open(out_path, 'w') as out_file:

            for taxon_name in taxon_names:
                seqs = []

                hits = db.select_stitched_contig_count(
                    cxn,
                    ref_name,
                    taxon_name,
                    iteration=iteration)

                if not hits:
                    continue

                for contig in db.select_stitched_contigs(
                        cxn, ref_name, taxon_name, iteration=iteration):

                    seqs.append(contig['seq'])

                header = taxon_name
                if args.reference_name:
                    header = '{}.{}'.format(ref_name, taxon_name)

                util.write_fasta_record(out_file, header, ''.join(seqs))


def output_summary_per_gene(args, cxn, iteration):
    """Print per gene summary statistics."""
    out_path = util.prefix_file(
        args.output_prefix, 'summary_stats_per_ref_gene.csv')
    with open(out_path, 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(['Locus', 'Taxon', 'Query_Length', 'Target_Length'])

        for stats in db.select_per_gene_stats(cxn, iteration):
            writer.writerow([
                stats['ref_name'],
                stats['taxon_name'],
                stats['query_len'],
                stats['target_len']])


def output_summary_per_taxon(args, cxn, iteration):
    """Print per taxon summary statistics."""
    out_path = util.prefix_file(
        args.output_prefix, 'summary_stats_per_taxon.csv')
    with open(out_path, 'w') as out_file:
        writer = csv.writer(out_file)
        writer.writerow([
            'Taxon',
            'Number_Genes',
            'Full_Exons',
            '95%', '90%', '80%', '70%', '50%', '10%'])

        for stats in db.select_per_taxon_stats(cxn, iteration):
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
