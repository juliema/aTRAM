"""
Format the data so that atram can use it later in atram itself.

It takes sequence read archive (SRA) files and converts them into coordinated
blast and sqlite3 databases.
"""

import os
import sys
import multiprocessing
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser
import lib.db as db
import lib.log as log
import lib.blast as blast


def preprocess(args):
    """Build the databases required by atram."""
    log.setup(args['log_file'], args['blast_db'])

    with db.connect(args['blast_db'], clean=True) as cxn:
        db.create_metadata_table(cxn)

        db.create_sequences_table(cxn)
        load_seqs(args, cxn)

        log.info('Creating an index for the sequence table')
        db.create_sequences_index(cxn)

        shard_list = assign_seqs_to_shards(cxn, args['shard_count'])

    create_all_blast_shards(args, shard_list)


def load_seqs(args, cxn):
    """Load sequences from a fasta/fastq files into the atram database."""
    # We have to clamp the end suffix depending on the file type.
    for (ends, clamp) in [('mixed_ends', ''), ('end_1', '1'),
                          ('end_2', '2'), ('single_ends', '')]:
        if args.get(ends):
            for file_name in args[ends]:
                load_one_file(cxn, file_name, ends, clamp)


def load_one_file(cxn, file_name, ends, seq_end_clamp=''):
    """Load sequences from a fasta/fastq file into the atram database."""
    log.info('Loading "{}" into sqlite database'.format(file_name))

    is_fastq = file_name.lower().endswith('q')
    parser = FastqGeneralIterator if is_fastq else SimpleFastaParser

    with open(file_name) as sra_file:
        batch = []

        for rec in parser(sra_file):
            title = rec[0].strip()
            seq = rec[1]
            seq_name, seq_end = blast.parse_fasta_title(
                title, ends, seq_end_clamp)

            batch.append((seq_name, seq_end, seq))

            if len(batch) >= db.BATCH_SIZE:
                db.insert_sequences_batch(cxn, batch)
                batch = []

        db.insert_sequences_batch(cxn, batch)


def assign_seqs_to_shards(cxn, shard_count):
    """Assign sequences to blast DB shards."""
    log.info('Assigning sequences to shards')

    total = db.get_sequence_count(cxn)
    offsets = np.linspace(0, total - 1, dtype=int, num=shard_count + 1)
    cuts = [db.get_shard_cut(cxn, offset) for offset in offsets]

    # Make sure the last sequence gets included
    cuts[-1] = cuts[-1] + 'z'

    # Now organize the list into pairs of sequence names
    pairs = [(cuts[i - 1], cuts[i]) for i in range(1, len(cuts))]

    return pairs


def create_all_blast_shards(args, shard_list):
    """
    Assign processes to make the blast DBs.

    One process for each blast DB shard.
    """
    log.info('Making blast DBs')

    with multiprocessing.Pool(processes=args['cpus']) as pool:
        results = []
        for idx, shard_params in enumerate(shard_list, 1):
            results.append(pool.apply_async(
                create_one_blast_shard,
                (args, shard_params, idx)))

        all_results = [result.get() for result in results]
    log.info('Finished making blast all {} DBs'.format(len(all_results)))


def create_one_blast_shard(args, shard_params, shard_index):
    """
    Create a blast DB from the shard.

    We fill a fasta file with the appropriate sequences and hand things off
    to the makeblastdb program.
    """
    shard = '{}.{:03d}.blast'.format(args['blast_db'], shard_index)
    fasta_name = '{}_{:03d}.fasta'.format(os.path.basename(sys.argv[0]),
                                          shard_index)
    fasta_path = os.path.join(args['temp_dir'], fasta_name)

    fill_blast_fasta(args['blast_db'], fasta_path, shard_params)

    blast.create_db(args['temp_dir'], fasta_path, shard)


def fill_blast_fasta(blast_db, fasta_path, shard_params):
    """
    Fill the fasta file used as input into blast.

    Use sequences from the sqlite3 DB. We use the shard partitions passed in to
    determine which sequences to get for this shard.
    """
    with db.connect(blast_db) as cxn:
        limit, offset = shard_params

        with open(fasta_path, 'w') as fasta_file:
            for row in db.get_sequences_in_shard(cxn, limit, offset):
                seq_end = '/{}'.format(row[1]) if row[1] else ''
                fasta_file.write('>{}{}\n'.format(row[0], seq_end))
                fasta_file.write('{}\n'.format(row[2]))
