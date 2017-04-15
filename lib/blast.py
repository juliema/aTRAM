"""All blast commands used by aTRAM."""

import os
import re
import sys
import glob
import lib.log as log


# Try to get the sequence name and which end it is from the fasta header
PARSE_HEADER = re.compile(r'^ [>@] \s* ( .* ) ( [\s/._] [12] ) \s* $',
                          re.VERBOSE)

# Parse blast hits file
PARSE_RESULTS = re.compile(r'^ ( .* ) ( [\s\/_] [12] )', re.VERBOSE)


def create_db(temp_dir, fasta_file, blast_db):
    """Create a blast DB."""

    cmd = 'makeblastdb -dbtype nucl -in {} -out {}'
    cmd = cmd.format(fasta_file, blast_db)
    log.subcommand(cmd, temp_dir)


def against_sra(args, blast_db, query, hits_file, iteration):
    """Blast the query sequences against an SRA blast DB."""

    cmd = []

    if args['protein'] and iteration == 1:
        cmd.append('tblastn')
        cmd.append('-db_gencode {}'.format(args['db_gencode']))
    else:
        cmd.append('blastn')
        cmd.append('-evalue {}'.format(args['evalue']))

    cmd.append("-outfmt '10 sseqid'")
    cmd.append('-max_target_seqs {}'.format(args['max_target_seqs']))
    cmd.append('-out {}'.format(hits_file))
    cmd.append('-db {}'.format(blast_db))
    cmd.append('-query {}'.format(query))

    command = ' '.join(cmd)
    log.subcommand(command, args['temp_dir'])


def against_contigs(args, blast_db, query, hits_file):
    """Blast the query sequence against the contings. The blast output
    will have the scores for later processing.
    """

    cmd = []

    if args.protein:
        cmd.append('tblastn')
        cmd.append('-db_gencode {}'.format(args.db_gencode))
    else:
        cmd.append('blastn')

    cmd.append('-db {}'.format(blast_db))
    cmd.append('-query {}'.format(query))
    cmd.append('-out {}'.format(hits_file))
    cmd.append(
        "-outfmt '10 qseqid sseqid bitscore qstart qend sstart send slen'")

    command = ' '.join(cmd)
    log.subcommand(command, args.temp_dir)


def shard_path(blast_db, shard_index):
    """Create the BLAST shard DB names."""

    file_name = '{}.blast_{:03d}'.format(blast_db, shard_index)
    return file_name


def all_shard_paths(blast_db):
    """Get all of the BLAST DB names built by the preprocessor."""

    pattern = '{}.blast_*.nhr'.format(blast_db)

    files = glob.glob(pattern)
    if not files:
        err = ('No blast shards found. Looking for "{}"\n'
               'Verify the --work-dir and --file-prefix options.').format(
                   pattern[:-4])
        sys.exit(err)

    return sorted([f[:-4] for f in files])


def output_file(temp_dir, shrd_path, iteration):
    """Create a file name for blast results."""

    shard_name = os.path.basename(shrd_path)
    file_name = '{}.{}.results'.format(shard_name, iteration)

    return os.path.join(temp_dir, file_name)


def temp_db(temp_dir, blast_db, iteration):
    """Generate a name for the temp DB used to filter the contigs."""

    file_name = '{}.blast.{:02d}'.format(blast_db, iteration)
    return os.path.join(temp_dir, file_name)
