"""The aTRAM assembly program."""

import re
import logging
import sqlite3
import subprocess
import multiprocessing
from more_itertools import chunked
from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
import configure
import util
from assembler import Assembler

FRAGMENT = re.compile(r'^ ( .* ) ( [\s\/_] [12] )', re.VERBOSE)


def blast(config, target, iteration, shard):
    """Blast the target sequences against an SRA blast DB."""
    out = util.blast_result_file(shard, iteration)
    blast_args = dict(outfmt="'6 sseqid'", max_target_seqs=config['max_target_seqs'],
                      out=out, db=shard, query=target)
    if config['protein'] and iteration == 1:
        cmd = str(NcbitblastnCommandline(cmd='tblastn', db_gencode=config['genetic_code'],
                                         **blast_args))
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
    db_path = util.db_file(config)
    db_conn = sqlite3.connect(db_path)
    return db_conn


def get_matching_fragments(iteration, shards):
    """Get all of the blast hits."""
    frags = {}
    for shard in shards:
        file_name = util.blast_result_file(shard, iteration)
        with open(file_name, 'r') as match_file:
            for line in match_file:
                match = FRAGMENT.match(line)
                frags[match.group(1)] = 1
    return frags


def write_sequences(config, iteration, fragments):
    """
    Take the matching blast hits and write the sequence and any matching end to the appropriate
    fasta files.
    """
    db_conn = connect_db(config)
    fasta_1 = util.paired_end_file(config, iteration, '1')
    fasta_2 = util.paired_end_file(config, iteration, '2')
    paired = False
    with open(fasta_1, 'w') as file_1, open(fasta_2, 'w') as file_2:
        for ids in chunked(fragments.keys(), 100):
            sql = 'SELECT * FROM frags WHERE frag IN ({})'.format(','.join('?' for _ in ids))
            for row in db_conn.execute(sql, ids):
                if row[1] and row[1].endswith('2'):
                    file_2.write('>{}{}\n'.format(row[0], row[1]))
                    file_2.write('{}\n'.format(row[2]))
                    paired = True
                else:
                    file_1.write('>{}{}\n'.format(row[0], row[1]))
                    file_1.write('{}\n'.format(row[2]))
    return paired


def create_blast_db(config, iteration):
    """Create a blast DB from the assembled fragments."""
    blast_db = util.blast_contig_file(config, iteration)
    fasta_file = util.contig_file(config, iteration)
    cmd = 'makeblastdb -dbtype nucl -in {} -out {}'.format(fasta_file, blast_db)
    subprocess.check_call(cmd, shell=True)
    return blast_db


def filter_contigs(config, iteration):
    """Remove junk from the assembled contigs."""
    blast_db = create_blast_db(config, iteration)
    scored_contigs = util.contig_score_file(config, iteration)
    blast_args = dict(db=blast_db, query=config['target'], out=scored_contigs,
                      outfmt="'6 qseqid sseqid bitscore qstart qend sstart send qlen'")
    if config['protein']:
        cmd = str(NcbitblastnCommandline(cmd='tblastn', db_gencode=config['genetic_code'],
                                         **blast_args))
    else:
        cmd = str(NcbiblastnCommandline(cmd='blastn', task='blastn', **blast_args))
    print(cmd)
    # subprocess.check_call(cmd, shell=True)


def atram(config):
    """The main aTRAM program loop."""
    shards = util.get_blast_shards(config)
    assember = Assembler.factory(config)
    target = config['target']
    for iteration in range(1, config['iterations'] + 1):
        logging.info('aTRAM iteration %i', iteration)
        blast_sra(config, iteration, shards, target)
        fragments = get_matching_fragments(iteration, shards)
        paired = write_sequences(config, iteration, fragments)
        assember.assemble(iteration, paired)
        filter_contigs(config, iteration)
        # target = new file
        break


if __name__ == '__main__':
    ARGS = configure.parse_command_line(
        description=""" """,
        args=['target', 'protein', 'iterations', 'cpu', 'evalue', 'max_target_seqs', 'assembler',
              'max_memory', 'file_prefix', 'work_dir', 'bit_score', 'genetic_code'])
    util.log_setup(ARGS)
    atram(ARGS)
