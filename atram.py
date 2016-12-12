"""The aTRAM assembly program."""

import argparse
import logging
# import subprocess
# import multiprocessing
from Bio.Blast.Applications import NcbiblastnCommandline, NcbitblastnCommandline
# from Bio.Blast.Applications import NcbitblastxCommandline
import log
import configure


def blast(config, targets, blast_db, iteration):
    """Blast the target sequences against an SRA blast DB."""
    out = '{}_{}'.format(blast_db, iteration)  # ???
    blast_args = dict(outfmt="'6 sseqid'", max_target_seqs=config.max_target_seqs,
                      out=out, db=blast_db, query=targets)
    if config.protein:
        cmd = str(NcbitblastnCommandline(cmd='tblastn', **blast_args))
    else:
        cmd = str(NcbiblastnCommandline(cmd='blastn', task='blastn',
                                        evalue=config.evalue, **blast_args))
    print(cmd)
    # subprocess.check_call(cmd, shell=True)


def blast_sra(config, iteration):
    """
    Blast the targets against the SRA databases. We're using a map-reduce strategy here.
    We map the blasting of the target sequences and reduce the output into one fasta file.
    """
    blast(config, [], '', iteration)
    # with multiprocessing.Pool(processes=config.processes) as pool:
    #     for shard in config.shards:
    #         targets = shard   # ???
    #         blast_db = shard  # ???
    #         proc = pool.Process(target=blast, args=(config, targets, blast_db, shard, iteration))
    #         proc.start()
    # pool.join()


def get_matching_ends():
    """Take all of the blast hits and append any matching ends that are not already found."""


def assemble_hits():
    """Use an assembler to build up the contigs."""


def filter_contigs():
    """Remove junk from the assembled contigs."""


def atram(config):
    """The main aTRAM program loop."""
    for iteration in range(1, config.iterations + 1):
        logging.info('aTRAM iteration %i', iteration)
        blast_sra(config, iteration)
        get_matching_ends()
        assemble_hits()
        filter_contigs()


def parse_args():
    """Parse the input arguments and assign defaults."""
    parser = argparse.ArgumentParser(description=''' ''')
    configure.add_argument(parser, 'out')
    configure.add_argument(parser, 'blast_db')
    configure.add_argument(parser, 'protein')
    configure.add_argument(parser, 'iterations')
    configure.add_argument(parser, 'evalue')
    configure.add_argument(parser, 'max_target_seqs')
    config = configure.parse_args(parser)
    return config


if __name__ == '__main__':
    ARGS = parse_args()
    log.setup(ARGS)
    atram(ARGS)
