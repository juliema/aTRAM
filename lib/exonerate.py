"""Common exonerate related functions."""

from os.path import abspath, join
from pathlib import Path
from collections import namedtuple
from Bio.SeqIO.FastaIO import SimpleFastaParser
import lib.stitcher_db as db
import lib.log as log
import lib.util as util


def run_exonerate(temp_dir, cxn, iteration):
    """Run exonerate on every reference sequence, taxon combination."""
    for ref in db.select_reference_genes(cxn):
        log.info('{} exonerate run for: {}'.format(
            util.as_word(iteration), ref['ref_name']))

        results_file = abspath(join(
            temp_dir,
            '{}.iteration_{}.results.fasta'.format(
                ref['ref_name'], iteration)))

        Path(results_file).touch()

        for contig_file in db.select_contigs(
                cxn, ref['ref_name'], iteration=iteration):

            if util.fasta_file_is_empty(contig_file['contig_file']):
                continue

            exonerate_command(temp_dir, ref, contig_file, results_file)

        insert_exonerate_results(cxn, iteration, results_file)


def exonerate_command(temp_dir, ref, contig_file, results_file):
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

    log.subcommand(cmd, temp_dir)


def insert_exonerate_results(cxn, iteration, results_file):
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
                'iteration': iteration,
                'seq': seq}
            batch.append(result)

    db.insert_exonerate_results(cxn, batch)
