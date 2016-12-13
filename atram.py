"""The aTRAM assembly program."""

import argparse
import logging
import subprocess
import multiprocessing
from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
# from Bio.Blast.Applications import NcbitblastxCommandline
import configure
import util


def blast(config, target, iteration, shard):
    """Blast the target sequences against an SRA blast DB."""
    out = util.blast_result_file(shard, iteration)
    blast_args = dict(outfmt="'6 sseqid'", max_target_seqs=config.max_target_seqs,
                      out=out, db=shard, query=target)
    if config.protein and iteration == 1:
        cmd = str(NcbitblastnCommandline(cmd='tblastn', **blast_args))
    else:
        cmd = str(NcbiblastnCommandline(cmd='blastn', task='blastn',
                                        evalue=config.evalue, **blast_args))
    # print(cmd)
    subprocess.check_call(cmd, shell=True)


def blast_sra(config, iteration, shards, target):
    """
    Blast the targets against the SRA databases. We're using a map-reduce strategy here.
    We map the blasting of the target sequences and reduce the output into one fasta file.
    """
    with multiprocessing.Pool(processes=config.processes) as pool:
        for shard in shards:
            proc = pool.Process(target=blast, args=(config, target, iteration, shard))
            proc.start()
    pool.join()


def retrieve_paired_ends(config, iteration, shard):
    """Spin thru the file and grab matching sequences from the DB"""
    file_name = util.blast_result_file(shard, iteration)
    print(file_name)


def get_matching_ends(config, iteration, shards):
    """Take all of the blast hits and append any matching ends that are not already found."""
    with multiprocessing.Pool(processes=config.processes) as pool:
        for shard in shards:
            proc = pool.Process(target=retrieve_paired_ends, args=(config, iteration, shard))
            proc.start()
    pool.join()


def assemble_hits():
    """Use an assembler to build up the contigs."""


def filter_contigs():
    """Remove junk from the assembled contigs."""


def atram(config):
    """The main aTRAM program loop."""
    shards = util.get_blast_shards(config)
    target = config.target
    for iteration in range(1, config.iterations + 1):
        logging.info('aTRAM iteration %i', iteration)
        blast_sra(config, iteration, shards, target)
        get_matching_ends(config, iteration, shards)
        assemble_hits()
        filter_contigs()
        # target = new file
        break


def parse_args():
    """Parse the input arguments and assign defaults."""
    parser = argparse.ArgumentParser(description=''' ''')
    configure.add_argument(parser, 'out')
    configure.add_argument(parser, 'blast_db')
    configure.add_argument(parser, 'target')
    configure.add_argument(parser, 'protein')
    configure.add_argument(parser, 'iterations')
    configure.add_argument(parser, 'processes')
    configure.add_argument(parser, 'evalue')
    configure.add_argument(parser, 'max_target_seqs')
    config = configure.parse_args(parser)
    return config


if __name__ == '__main__':
    ARGS = parse_args()
    util.log_setup(ARGS)
    atram(ARGS)
