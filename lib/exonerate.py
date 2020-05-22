"""Common exonerate related functions."""

import re
import os
from os.path import abspath, join, basename
from collections import defaultdict, namedtuple
from glob import glob
from pathlib import Path
from Bio.SeqIO.FastaIO import SimpleFastaParser
from . import db_stitcher as db
from . import log
from . import util


def run_exonerate(temp_dir, cxn, iteration):
    """Run exonerate on every reference sequence, taxon combination."""
    for ref in db.select_reference_genes(cxn):
        log.info('Exonerate run for: {}'.format(ref['ref_name']))

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


def create_tables(cxn):
    """Create database tables."""
    db.create_reference_genes_table(cxn)
    db.create_exonerate_table(cxn)
    db.create_contigs_table(cxn)
    db.create_stitch_table(cxn)


def get_taxa(args):
    """Insert taxa into the database."""
    log.info('Preparing taxa')
    with open(args.taxa) as taxa:
        taxon_names = [n.strip() for n in taxa]
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

    ref_names = [(x, re.sub(sep, sep, x) + sep) for x in ref_names]
    ref_names = sorted(
            ref_names, key=lambda x: (len(x[1]), x), reverse=True)

    taxon_names = [(x, re.sub(sep, sep, x) + sep) for x in taxon_names]
    taxon_names = sorted(
            taxon_names, key=lambda x: (len(x[1]), x), reverse=True)

    ref_name = [x[0] for x in ref_names if re.search(x[1], contig_file)]
    taxon_name = [x[0] for x in taxon_names if re.search(x[1], contig_file)]

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
    names_seen = defaultdict(int)

    ref_names = set(x['ref_name'] for x in db.select_reference_genes(cxn))

    pattern = join(args.assemblies_dir, args.file_filter)
    for contig_path in sorted(glob(pattern)):

        if os.stat(contig_path).st_size == 0:
            continue

        with open(contig_path) as contig_old:
            for i, (header, contig_seq) \
                    in enumerate(SimpleFastaParser(contig_old)):

                contig_file = basename(contig_path)
                ref_name, taxon_name = parse_contig_file_name(
                    ref_names, taxon_names, contig_file)

                if ref_name not in ref_names or taxon_name not in taxon_names:
                    continue

                contig_name = name_contig(
                    taxon_name, ref_name, header, names_seen)
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


DEFAULT = 0
ITERATION = re.compile(r'iteration.(\d+)', re.IGNORECASE)
COVERAGE = re.compile(r'cov.([\d.]+)', re.IGNORECASE)
SCORE = re.compile(r'score.([\d.]+)', re.IGNORECASE)


def name_contig(taxon_name, ref_name, header, names_seen):
    """Shorten contig names."""
    global DEFAULT  # pylint: disable=global-statement
    DEFAULT += 1

    match = ITERATION.search(header)
    iteration = 'i{}'.format(match[1]) if match else ''

    match = COVERAGE.search(header)
    coverage = 'c{}'.format(round(float(match[1]))) if match else ''

    match = SCORE.search(header)
    score = 's{}'.format(round(float(match[1]))) if match else ''

    contig = '{}{}{}'.format(iteration, coverage, score)
    contig = contig if contig else str(DEFAULT)

    name = '{}@{}_{}'.format(taxon_name, ref_name, contig)
    name = re.sub(r'[^\w@]+', '_', name.strip())

    name = handle_duplicate_name(name, names_seen)

    return name


def handle_duplicate_name(contig_name, names_seen):
    """Add a tiebreaker to a duplicate contig name."""
    name = re.sub(r'_v\d+$', '', contig_name, re.IGNORECASE)

    names_seen[name] += 1

    if names_seen[name] > 1:
        name += '_v{}'.format(names_seen[name])

    return name
