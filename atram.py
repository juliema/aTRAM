"""The aTRAM assembly program."""

import argparse
import logging
import subprocess
import multiprocessing
from Bio.Blast.Applications import NcbiblastnCommandline
import log
import configure


def blast(config, targets, blast_db, iteration):
    """Blast the target sequences against an SRA blast DB."""
    # push @pids, fork_cmd (get_bin("blastn"), "-task blastn -evalue $evalue
    # -max_target_seqs $max_target_seqs -db $atram_db.$s.db -query $search_fasta
    # -outfmt '6 sseqid' -out $current_shard");
    out = '{}_{}'.format(blast_db, iteration)  # ???
    cmd = NcbiblastnCommandline(
        cmd='blastn', outfmt='6 sseqid', evalue=config.evalue, db=blast_db,
        max_target_seqs=config.max_target_seqs, out=out, query=targets)
    subprocess.check_call(cmd, shell=True)


def blast_sra(config, iteration):
    """
    Blast the targets against the SRA databases. We're using a map-reduce strategy here.
    We map the blasting of the target sequences and reduce the output into one fasta file.
    """
    with multiprocessing.Pool(processes=config.processes) as pool:
        for shard in config.shards:
            targets = shard   # ???
            blast_db = shard  # ???
            proc = pool.Process(target=blast, args=(config, targets, blast_db, shard, iteration))
            proc.start()
    pool.join()


def get_matching_ends():
    """Take all of the blast hits and append any matching ends that are not already found."""


def assemble_hits():
    """Use an assembler to build up the contigs."""


def filter_contigs():
    """Remove junk from the assembled contigs."""


def atram(config):
    """The main aTRAM program loop."""
    for iteration in range(config.iterations):
        logging.info('aTRAM iteration %i', iteration)
        blast_sra(config, iteration)
        get_matching_ends()
        assemble_hits()
        filter_contigs()


def parse_args():
    """Parse the input arguments and assign defaults."""
    parser = argparse.ArgumentParser(description=''' ''')
    configure.add_argument(parser, 'out')
    config = parser.parse_args()
    return config


if __name__ == '__main__':
    ARGS = parse_args()
    log.setup(ARGS)
    atram(ARGS)
