"""The aTRAM assembly program."""

import re
import logging
import sqlite3
# import subprocess
import multiprocessing
from more_itertools import chunked
from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
# from Bio.Blast.Applications import NcbitblastxCommandline
import configure
import util

FRAGMENT = re.compile(r'^ ( .* ) ( [\s\/_] [12] )', re.VERBOSE)


def blast(config, target, iteration, shard):
    """Blast the target sequences against an SRA blast DB."""
    out = util.blast_result_file(shard, iteration)
    blast_args = dict(outfmt="'6 sseqid'", max_target_seqs=config['max_target_seqs'],
                      out=out, db=shard, query=target)
    if config['protein'] and iteration == 1:
        cmd = str(NcbitblastnCommandline(cmd='tblastn', **blast_args))
    else:
        cmd = str(NcbiblastnCommandline(cmd='blastn', task='blastn',
                                        evalue=config['evalue'], **blast_args))
    print(cmd)
    # subprocess.check_call(cmd, shell=True)


def blast_sra(config, iteration, shards, target):
    """
    Blast the targets against the SRA databases. We're using a map-reduce strategy here.
    We map the blasting of the target sequences and reduce the output into one fasta file.
    """
    with multiprocessing.Pool(processes=config['cpu']) as pool:
        results = [pool.apply_async(blast, (config, target, iteration, shard)) for shard in shards]
        _ = [result.get() for result in results]


def connect_db(config):
    """Setup the DB for our processing needs."""
    db_path = '{}sqlite.db'.format(config['blast_db'])
    db_conn = sqlite3.connect(db_path)
    return db_conn


def get_matching_ends(config, iteration, shards):
    """Take all of the blast hits and append any matching ends that are not already found."""
    db_conn = connect_db(config)
    frags = {}
    for shard in shards:
        file_name = util.blast_result_file(shard, iteration)
        with open(file_name, 'r') as match_file:
            for line in match_file:
                match = FRAGMENT.match(line)
                frags[match.group(1)] = 1

    fasta_file = '{}matching_seqs_{}.fasta'.format(config['blast_db'], iteration)
    batches = chunked(frags.keys(), 100)
    with open(fasta_file, 'w') as out_file:
        for ids in batches:
            sql = 'SELECT * FROM frags WHERE frag IN ({})'.format(','.join('?' for _ in ids))
            for row in db_conn.execute(sql, ids):
                out_file.write('>{}{}\n'.format(row[0], row[1]))
                out_file.write('{}\n'.format(row[2]))


def assemble_hits_trinity(config, iteration):
    """Use an assembler to build up the contigs."""
    fasta_file = '{}matching_seqs_{}.fasta'.format(config['blast_db'], iteration)
    # trinity --seqType fa --single $short_read_file --run_as_paired --JM $jm --output $tempdir
    cmd = 'Trinity --seqType fa --single {} --run_as_paired --output $tempdir'
    cmd = cmd.format(fasta_file)
    print(cmd)


def assemble_hits(config, iteration):
    """Use an assembler to build up the contigs."""
    assemble_hits_trinity(config, iteration)


def filter_contigs():
    """Remove junk from the assembled contigs."""


def atram(config):
    """The main aTRAM program loop."""
    shards = util.get_blast_shards(config)
    target = config['target']
    for iteration in range(1, config['iterations'] + 1):
        logging.info('aTRAM iteration %i', iteration)
        blast_sra(config, iteration, shards, target)
        get_matching_ends(config, iteration, shards)
        assemble_hits(config, iteration)
        # filter_contigs()
        # target = new file
        break


if __name__ == '__main__':
    ARGS = configure.parse_command_line(
        description=""" """,
        args=['blast_db', 'target', 'protein', 'iterations', 'cpu', 'evalue', 'max_target_seqs'])
    util.log_setup(ARGS)
    atram(ARGS)
