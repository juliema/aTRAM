"""
Stitch together exons from targeted assemblies.

It uses amino acid targets and DNA assemblies.
"""

import re
import os
from os.path import join, split
from glob import glob
from Bio.SeqIO.FastaIO import SimpleFastaParser
import lib.db as db
import lib.log as log
import lib.util as util


def stitch(args):
    """Stitch the exons together."""
    log.stitcher_setup(args['log_file'])

    with util.make_temp_dir(
            where=args['temp_dir'],
            prefix='atram_stitcher_',
            keep=args['keep_temp_dir']) as temp_dir:
        util.update_temp_dir(temp_dir, args)

        with db.temp_db(args['temp_dir'], 'atram_stitcher') as cxn:
            db.create_stitcher_tables(cxn)

            asm_dir = join(temp_dir, 'asm')
            os.makedirs(asm_dir, exist_ok=True)
            edit_assemblies(args, asm_dir)

            ref_dir = join(temp_dir, 'ref')
            os.makedirs(ref_dir, exist_ok=True)
            edit_reference_genes(args, ref_dir)


def edit_assemblies(args, asm_dir):
    """Remove characters exonerate cannot handle from assembly names."""
    idx = 0
    pattern = join(args['assemblies_dir'], '*.fasta')
    for fasta_in in glob(pattern):

        _, file_part = split(fasta_in)
        fasta_out = join(asm_dir, file_part)

        with open(fasta_in) as asm_in, open(fasta_out, 'w') as asm_out:
            for seq_name, seq in SimpleFastaParser(asm_in):
                idx += 1
                seq_name = seq_name.strip()
                seq_name = re.sub(r'[\s,=@+\[\]:!-]+', '_', seq_name)
                seq_name = '{}_{}'.format(seq_name, idx)

                util.write_fasta_record(asm_out, seq_name, seq)


def edit_reference_genes(args, ref_dir):
    """Format reference sequences for exonerate."""
    ref_genes = args['reference_genes']
    gene_list = join(ref_dir, 'gene_list.txt')

    with open(ref_genes) as ref_in, open(gene_list, 'w') as ref_list:
        for seq_name, seq in SimpleFastaParser(ref_in):
            seq_name = re.sub(r'\s+', '_', seq_name)

            ref_list.write(seq_name)
            ref_list.write('\n')

            ref_file = join(ref_dir, seq_name + '.reference.fasta')
            with open(ref_file, 'w') as ref_out:
                util.write_fasta_record(ref_out, seq_name, seq)
