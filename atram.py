"""The aTRAM assembly program."""

import re
import csv
import logging
import sqlite3
import subprocess
import multiprocessing
from more_itertools import chunked
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
import configure
import util
from assembler import Assembler

FRAGMENT = re.compile(r'^ ( .* ) ( [\s\/_] [12] )', re.VERBOSE)


def blast(config, target, iteration, shard):
    """Blast the target sequences against an SRA blast DB."""

    out = util.blast_result_file(shard, iteration)
    blast_args = dict(outfmt="'10 sseqid'",
                      max_target_seqs=config['max_target_seqs'],
                      out=out, db=shard, query=target)
    if config['protein'] and iteration == 1:
        cmd = str(NcbitblastnCommandline(
            cmd='tblastn', db_gencode=config['genetic_code'], **blast_args))
    else:
        cmd = str(NcbiblastnCommandline(
            cmd='blastn', task='blastn', evalue=config['evalue'],
            **blast_args))
    print(cmd)
    subprocess.check_call(cmd, shell=True)


def blast_sra(config, iteration, shards, target):
    """
    Blast the targets against the SRA databases. We're using a map-reduce
    strategy here. We map the blasting of the target sequences and reduce
    the output into one fasta file.
    """

    with multiprocessing.Pool(processes=config['cpus']) as pool:
        results = [pool.apply_async(blast, (config, target, iteration, shard))
                   for shard in shards]
        _ = [result.get() for result in results]


def connect_db(work_dir, file_prefix):
    """Setup the DB for our processing needs."""

    db_path = util.db_file(work_dir, file_prefix)
    db_conn = sqlite3.connect(db_path)
    return db_conn


def get_matching_fragments(iteration, shards):
    """Get all of the blast hits."""

    frags = {}
    for shard in shards:
        file_name = util.blast_result_file(shard, iteration)
        print(file_name)
        with open(file_name) as match_file:
            for line in match_file:
                match = FRAGMENT.match(line)
                frags[match.group(1) if match else line] = 1
    return frags


def write_sequences(work_dir, file_prefix, iteration, fragments):
    """
    Take the matching blast hits and write the sequence and any matching end to
    the appropriate fasta files.
    """

    db_conn = connect_db(work_dir, file_prefix)
    fasta_1 = util.paired_end_file(work_dir, file_prefix, iteration, '1')
    fasta_2 = util.paired_end_file(work_dir, file_prefix, iteration, '2')
    paired = False
    with open(fasta_1, 'w') as file_1, open(fasta_2, 'w') as file_2:
        for ids in chunked(fragments.keys(), 100):
            sql = 'SELECT * FROM frags WHERE frag IN ({})'
            sql = sql.format(','.join('?' for _ in ids))
            for row in db_conn.execute(sql, ids):
                if row[1] and row[1].endswith('2'):
                    file_2.write('>{}{}\n'.format(row[0], row[1]))
                    file_2.write('{}\n'.format(row[2]))
                    paired = True
                else:
                    file_1.write('>{}{}\n'.format(row[0], row[1]))
                    file_1.write('{}\n'.format(row[2]))
    return paired


def create_blast_db(work_dir, file_prefix, iteration):
    """Create a blast DB from the assembled fragments."""

    blast_db = util.contig_score_db(work_dir, file_prefix, iteration)
    fasta_file = util.contig_unfiltered_file(work_dir, file_prefix, iteration)
    cmd = 'makeblastdb -dbtype nucl -in {} -out {}'
    cmd = cmd.format(fasta_file, blast_db)
    subprocess.check_call(cmd, shell=True)
    return blast_db


def blast_target_against_contigs(work_dir, file_prefix, iteration,
                                 target, genetic_code, protein):
    """
    Blast the target sequence against the contings. The blast output will have
    the scores for later processing.
    """

    blast_db = create_blast_db(work_dir, file_prefix, iteration)
    scored_contigs = util.contig_score_file(work_dir, file_prefix, iteration)
    blast_args = dict(
        db=blast_db, query=target, out=scored_contigs,
        outfmt="'10 qseqid sseqid bitscore qstart qend sstart send slen'")
    if protein:
        cmd = str(NcbitblastnCommandline(
            cmd='tblastn', db_gencode=genetic_code, **blast_args))
    else:
        cmd = str(NcbiblastnCommandline(
            cmd='blastn', task='blastn', **blast_args))
    print(cmd)
    subprocess.check_call(cmd, shell=True)


def filter_contig_scores(work_dir, file_prefix, iteration, bit_score):
    """Remove contigs that have scores below the bit score cut-off."""

    # qseqid sseqid bitscore qstart qend sstart send slen
    field_names = ['target_id', 'contig_id', 'bit_score',
                   'target_start', 'target_end',
                   'contig_start', 'contig_end', 'contig_len']
    contig_score_file = util.contig_score_file(work_dir, file_prefix,
                                               iteration)

    scores = {}
    with open(contig_score_file) as in_file:
        for score in csv.DictReader(in_file, field_names):
            score['bit_score'] = float(score['bit_score'])
            if score['bit_score'] >= bit_score:
                for field in field_names[3:]:
                    score[field] = int(score[field])
                scores[score['contig_id']] = score

    return scores


def filter_contigs(work_dir, file_prefix, iteration):
    """Remove junk from the assembled contigs."""

    blast_target_against_contigs(work_dir, file_prefix, iteration,
                                 target, genetic_code, protein)
    filtered_names = filter_contig_scores(work_dir, file_prefix,
                                          iteration, bit_score)
    unfiltered_file = util.contig_unfiltered_file(work_dir, file_prefix,
                                                  iteration)
    # filtered_file = util.contig_filtered_file(work_dir, file_prefix,
    #                                           iteration)

    print(filtered_names)
    with open(unfiltered_file) as in_file:
        for contig in SeqIO.parse(in_file, 'fasta'):
            if contig.id in filtered_names:
                print(contig)


def atram(config):
    """The main aTRAM program loop."""

    shards = util.shard_db_names(config)
    assember = Assembler.factory(config)
    target = config['target']
    for iteration in range(1, config['iterations'] + 1):
        logging.info('aTRAM iteration %i', iteration)
        blast_sra(config, iteration, shards, target)
        frag = get_matching_fragments(iteration, shards)
        paired = write_sequences(config, iteration, frag)
        assember.assemble(iteration, paired)
        filter_contigs(config, iteration)
        # target = new file
        break


if __name__ == '__main__':
    CONFIG = configure.parse_command_line(
        description=""" """,
        args=('target protein iterations cpus evalue max_target_seqs '
              'assembler max_memory file_prefix work_dir bit_score '
              'genetic_code kmer'))
    util.log_setup(CONFIG.work_dir, CONFIG.file_prefix)
    atram(CONFIG)
