#!/usr/bin/python3

"""The aTRAM assembly program."""

import argparse
import logging
import subprocess
import multiprocessing
import util


def blast(targets, blast_db):
    """Blast the target sequences against an SRA blast DB."""
    subprocess.check_call(("blastn -task blastn -evalue {} -max_target_seqs {} -db {} -query {} "
                           "-outfmt '6 sseqid' -out {}").format(
                               targets, blast_db, targets, targets, targets), shell=True)


def blast_sra(config):
    """
    Get the output from blast into memory. We're using a map-reduce strategy here.
    We map the blasting of the target sequences and reduce the output into one fasta file.
    """
    with multiprocessing.Pool(processes=config.processes) as pool:
        for i, shard in enumerate(config.shards):
            proc = pool.Process(target=blast, args=(config, shard, i))
            proc.start()
    pool.join()


def get_matching_ends():
    """Take all of the blast hits and append any matching ends that are not already found."""


def get_assemblies():
    """Get the ouput from the assembly program into memory."""


def assemble_hits():
    """Use an assembler to build up the contigs."""
    return get_assemblies()


def create_contig_db():
    """Create a blast DB from the assembled contigs."""


def blast_contigs():
    """Blast the target sequences against the contig DB"""


def filter_contigs():
    """Remove junk from the assembled contigs."""
    create_contig_db()
    blast_contigs()


def atram(config):
    """The main aTRAM program loop."""
    for i in range(config.iterations):
        logging.info('aTRAM iteration %i', i)
        blast_sra(config)
        get_matching_ends()
        assemble_hits()
        filter_contigs()


def parse_args():
    """Parse the input arguments and assign defaults."""
    parser = argparse.ArgumentParser(description=''' ''')
    parser.add_argument('-o', '--out',
                        help=('output aTRAM files with this prefix. '
                              'May include a directory in the prefix.'))
    config = parser.parse_args()
    return config


if __name__ == '__main__':
    ARGS = parse_args()
    util.setup_log(ARGS)
    atram(ARGS)
