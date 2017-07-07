"""All blast commands used by aTRAM."""

import os
import re
import glob
import json
import lib.log as log
import lib.file_util as file_util


# Try to get the sequence name and which end it is from the fasta header
PARSE_HEADER = re.compile(r'^ [>@] \s* ( [^\s/._]+? ) [\s/._] ( [12] )',
                          re.VERBOSE)

# Parse blast hits file
PARSE_RESULTS = re.compile(r'^ ( [^\s/._]+? ) [\s/._] ( [12] )', re.VERBOSE)


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

    cmd.append('-outfmt 15')
    cmd.append('-max_target_seqs {}'.format(args['max_target_seqs']))
    cmd.append('-out {}'.format(hits_file))
    cmd.append('-db {}'.format(blast_db))
    cmd.append('-query {}'.format(query))

    command = ' '.join(cmd)
    log.subcommand(command, args['temp_dir'])


def against_contigs(args, blast_db, query, hits_file):
    """Blast the query sequence against the contigs. The blast output
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
    cmd.append('-outfmt 15')

    command = ' '.join(cmd)
    log.subcommand(command, args.temp_dir)


def shard_path(blast_db, shard_index):
    """Create the BLAST shard DB names."""

    file_name = '{}.{:03d}.blast'.format(blast_db, shard_index)
    return file_name


def all_shard_paths(blast_db):
    """Get all of the BLAST DB names built by the preprocessor."""

    pattern = '{}.*.blast.nhr'.format(blast_db)

    files = glob.glob(pattern)
    if not files:
        err = ('No blast shards found. Looking for "{}"\n'
               'Verify the --work-dir and --file-prefix options.').format(
                   pattern[:-4])
        log.fatal(err)

    return sorted([f[:-4] for f in files])


def output_file_name(temp_dir, shrd_path, iteration):
    """Create a file name for blast results."""

    shard_name = os.path.basename(shrd_path)
    file_name = '{}.{:02d}.results.json'.format(shard_name, iteration)
    return file_util.temp_iter_file(temp_dir, iteration, file_name)


def temp_db_name(temp_dir, blast_db, iteration):
    """Generate a name for the temp DB used to filter the contigs."""

    file_name = os.path.basename(blast_db)
    file_name = '{}.{:02d}'.format(file_name, iteration)
    return file_util.temp_iter_file(temp_dir, iteration, file_name)


def hits(json_file):
    """Extract the blast hits from the blast json output file."""

    with open(json_file) as blast_file:
        raw = blast_file.read()
        obj = json.loads(raw)

    hits_list = []
    for raw in obj['BlastOutput2'][0]['report']['results']['search']['hits']:
        for i, desc in enumerate(raw['description']):
            hit = dict(desc)
            hit['len'] = raw['len']
            hit.update(raw['hsps'][i])
            hits_list.append(hit)

    return hits_list


def command_line_args(parser):
    """Add optional blast arguments to the command-line parser."""

    group = parser.add_argument_group('optional blast arguments')

    group.add_argument('--db-gencode', type=int, default=1,
                       metavar='CODE',
                       help='The genetic code to use during blast runs. '
                            'The default is "1".')

    group.add_argument('--evalue', type=float, default=1e-10,
                       help='The default evalue is "1e-10".')

    group.add_argument('--max-target-seqs', type=int, default=100000000,
                       metavar='MAX',
                       help='Maximum hit sequences per shard. '
                            'Default is calculated based on the available '
                            'memory and the number of shards.')


def check_command_line_args(args, temp_dir):
    """Make sure optional blast arguments are reasonable."""

    # Calculate the default max_target_seqs per shard
    if not args.max_target_seqs:
        all_shards = all_shard_paths(args.blast_db)
        args.max_target_seqs = int(2 * args.max_memory / len(all_shards)) * 1e6

    file_util.temp_root_dir(args, temp_dir)
