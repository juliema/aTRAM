"""
Format the data so that atram can use it later in atram itself.

It takes sequence read archive (SRA) files and converts them into coordinated
blast and sqlite3 databases.
"""

import multiprocessing
import sys
from os.path import basename, join, splitext

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from . import blast, db, db_preprocessor, log, util


def preprocess(args):
    """Build the databases required by atram."""
    log.setup(args['log_file'], args['log_level'], args['blast_db'])

    with util.make_temp_dir(
            where=args['temp_dir'],
            prefix='atram_preprocessor_',
            keep=args['keep_temp_dir']) as temp_dir:
        util.update_temp_dir(temp_dir, args)

        with db.connect(args['blast_db'], clean=True) as cxn:
            db_preprocessor.create_metadata_table(cxn, args)

            db_preprocessor.create_sequences_table(cxn)
            load_seqs(args, cxn)

            log.info('Creating an index for the sequence table')
            db_preprocessor.create_sequences_index(cxn)

            create_all_blast_shards(args, cxn, args['shard_count'])


def load_seqs(args, cxn):
    """Load sequences from a fasta/fastq files into the atram database."""
    # We have to clamp the end suffix depending on the file type.
    for (ends, clamp) in [('mixed_ends', ''), ('end_1', '1'),
                          ('end_2', '2'), ('single_ends', '')]:
        if args.get(ends):
            for file_name in args[ends]:
                load_one_file(args, cxn, file_name, ends, clamp)


def load_one_file(args, cxn, file_name, ends, seq_end_clamp=''):
    """Load sequences from a fasta/fastq file into the atram database."""
    log.info('Loading "{}" into sqlite database'.format(file_name))

    parser = get_parser(args, file_name)

    with util.open_file(args, file_name) as sra_file:
        batch = []

        for rec in parser(sra_file):
            title = rec[0].strip()
            seq = rec[1]
            seq_name, seq_end = blast.parse_fasta_title(
                title, ends, seq_end_clamp)

            batch.append((seq_name, seq_end, seq))

            if len(batch) >= db.BATCH_SIZE:
                db_preprocessor.insert_sequences_batch(cxn, batch)
                batch = []

        db_preprocessor.insert_sequences_batch(cxn, batch)


def get_parser(args, file_name):
    """Get either a fasta or fastq file parser."""
    is_fastq = util.is_fastq_file(args, file_name)
    return FastqGeneralIterator if is_fastq else SimpleFastaParser


def create_all_blast_shards(args, cxn, shard_count):
    """
    Assign processes to make the blast DBs.

    One process for each blast DB shard.
    """
    log.info('Making blast DBs')
    db_preprocessor.aux_db(cxn, args['temp_dir'])
    db_preprocessor.create_seq_names_table(cxn)

    with multiprocessing.Pool(processes=args['cpus']) as pool:
        results = []
        for shard_idx in range(shard_count):
            fasta_path = fill_blast_fasta(args, cxn, shard_count, shard_idx)
            results.append(pool.apply_async(
                create_one_blast_shard,
                (args, fasta_path, shard_idx)))

        all_results = [result.get() for result in results]
    db_preprocessor.aux_detach(cxn)
    log.info('Finished making all {} blast DBs'.format(len(all_results)))


def fill_blast_fasta(args, cxn, shard_count, shard_index):
    """Fill the shard input files with sequences."""
    exe_name, _ = splitext(basename(sys.argv[0]))
    fasta_name = '{}_{:03d}.fasta'.format(exe_name, shard_index + 1)
    fasta_path = join(args['temp_dir'], fasta_name)

    with open(fasta_path, 'w') as fasta_file:
        for row in db_preprocessor.get_sequences_in_shard(
                cxn, shard_count, shard_index):
            util.write_fasta_record(fasta_file, row[0], row[2], row[1])

    return fasta_path


def create_one_blast_shard(args, fasta_path, shard_index):
    """Create a blast DB from the shard."""
    shard = '{}.{:03d}.blast'.format(args['blast_db'], shard_index + 1)
    blast.create_db(args['temp_dir'], fasta_path, shard)
