"""The aTRAM assembly program."""

import sys
import re
import os
from os.path import basename, splitext, join
from shutil import which
import argparse
import textwrap
import tempfile
import multiprocessing
from Bio import SeqIO
import lib.db as db
import lib.log as log
import lib.bio as bio
import lib.util as util
import lib.blast as blast
import lib.assembler as assembly


__all__ = ('assemble', )


def assemble(args):
    """Loop thru every blast/query pair and run an assembly for each one."""
    queries = split_queries(args)

    for blast_db in args['blast_db']:

        with db.connect(blast_db, check_version=True) as db_conn:
            for query in queries:
                db.aux_db(db_conn, args['temp_dir'], blast_db, query)
                clean_database(db_conn)

                log.setup(args['log_file'], blast_db, query)

                assembler = assembly.factory(args, db_conn)

                try:
                    assembly_loop(assembler, blast_db, query)
                    assembler.write_final_output(blast_db, query)
                except (TimeoutError, RuntimeError):
                    pass
                except Exception as err:  # pylint: disable=broad-except
                    log.error('Exception: {}'.format(err))

                db.aux_detach(db_conn)


def assembly_loop(assembler, blast_db, query):
    """Iterate the assembly processes."""
    for iteration in range(1, assembler.args['iterations'] + 1):
        log.info('aTRAM blast DB = "{}", query = "{}", iteration {}'.format(
            blast_db, query, iteration))

        assembler.initialize_iteration(blast_db, query, iteration)

        os.makedirs(assembler.iter_dir(), exist_ok=True)

        blast_query_against_all_shards(assembler)

        if assembler.blast_only or assembler.no_blast_hits():
            break

        assembler.write_input_files()

        assembler.run()

        if assembler.nothing_assembled():
            break

        high_score = filter_contigs(assembler)

        count = assembler.assembled_contigs_count(high_score)

        if not count:
            break

        if assembler.no_new_contigs(count):
            break

        query = create_query_from_contigs(assembler)

    else:
        log.info('All iterations completed')


def split_queries(args):
    """
    Create query target for every query and query-split file.

    We put each query record into its own file for blast queries.
    """
    if not args.get('query_split'):
        return args['query'][:]

    queries = []

    path = join(args['temp_dir'], 'queries')
    os.makedirs(path, exist_ok=True)

    for query_path in args['query_split']:

        query_name = splitext(basename(query_path))[0]

        with open(query_path) as query_file:
            for i, rec in enumerate(SeqIO.parse(query_file, 'fasta'), 1):

                query_id = re.sub(r'\W+', '_', rec.id)

                query_file = join(
                    path,
                    '{}_{}_{}.fasta'.format(query_name, query_id, i))

                write_query_seq(query_file, rec.id, str(rec.seq))

                queries.append(query_file)

    return queries


def write_query_seq(file_name, seq_id, seq):
    """Write the sequence to a fasta file."""
    with open(file_name, 'w') as query_file:
        query_file.write('>{}\n'.format(seq_id))
        query_file.write('{}\n'.format(seq))


def clean_database(db_conn):
    """Create database tables for an atram run."""
    db.create_sra_blast_hits_table(db_conn)
    db.create_contig_blast_hits_table(db_conn)
    db.create_assembled_contigs_table(db_conn)


def blast_query_against_all_shards(assembler):
    """
    Blast the query against the SRA databases.

    We're using a map-reduce strategy here. We map the blasting of the query
    sequences and reduce the output into one fasta file.
    """
    log.info('Blasting query against shards: iteration {}'.format(
        assembler.state['iteration']))

    all_shards = shard_fraction(assembler)

    with multiprocessing.Pool(processes=assembler.args['cpus']) as pool:
        results = [pool.apply_async(
            blast_query_against_one_shard,
            (assembler.args, assembler.simple_state(), shard))
                   for shard in all_shards]  # noqa
        _ = [result.get() for result in results]  # noqa


def shard_fraction(assembler):
    """
    Get the shards we are using.

    We may not want the entire DB for highly redundant libraries.
    """
    all_shards = blast.all_shard_paths(assembler.state['blast_db'])
    last_index = int(len(all_shards) * assembler.args['fraction'])
    return all_shards[:last_index]


def blast_query_against_one_shard(args, state, shard):
    """
    Blast the query against one blast DB shard.

    Then write the results to the database.
    """
    temp_dir = util.iter_dir(
        args['temp_dir'],
        state['blast_db'],
        state['query_target'],
        state['iteration'])

    output_file = blast.output_file_name(temp_dir, shard)

    blast.against_sra(args, state, output_file, shard)

    with db.connect(state['blast_db']) as db_conn:
        db.aux_db(
            db_conn,
            args['temp_dir'],
            state['blast_db'],
            state['query_target'])

        shard = basename(shard)

        batch = []

        hits = blast.hits(output_file)
        for hit in hits:
            match = blast.PARSE_RESULTS.match(hit['title'])
            if match:
                seq_name = match.group(1)
                seq_end = match.group(2)
            else:
                seq_name = hit['title']
                seq_end = ''
            batch.append((state['iteration'], seq_end, seq_name, shard))
        db.insert_blast_hit_batch(db_conn, batch)


def filter_contigs(assembler):
    """Remove junk from the assembled contigs."""
    log.info('Saving assembled contigs: iteration {}'.format(
        assembler.state['iteration']))

    blast_db = blast.temp_db_name(
        assembler.iter_dir(), assembler.state['blast_db'])

    hits_file = blast.output_file_name(
        assembler.iter_dir(), assembler.state['blast_db'])

    blast.create_db(
        assembler.iter_dir(), assembler.file['output'], blast_db)

    blast.against_contigs(
        blast_db,
        assembler.state['query_target'],
        hits_file,
        protein=assembler.args['protein'],
        db_gencode=assembler.args['db_gencode'],
        temp_dir=assembler.args['temp_dir'])

    save_blast_against_contigs(assembler, hits_file)

    all_hits = {row['contig_id']: row
                for row
                in db.get_contig_blast_hits(
                    assembler.state['db_conn'],
                    assembler.state['iteration'])}

    return save_contigs(assembler, all_hits)


def save_blast_against_contigs(assembler, hits_file):
    """Save all of the blast hits."""
    batch = []

    for hit in blast.hits(hits_file):
        contig_id = assembler.parse_contig_id(hit['title'])
        batch.append((
            assembler.state['iteration'],
            contig_id,
            hit['title'],
            hit['bit_score'],
            hit['len'],
            hit['query_from'],
            hit['query_to'],
            hit.get('query_strand', ''),
            hit['hit_from'],
            hit['hit_to'],
            hit.get('hit_strand', '')))

    db.insert_contig_hit_batch(assembler.state['db_conn'], batch)


def save_contigs(assembler, all_hits):
    """Save the contigs to the database."""
    batch = []
    high_score = 0
    with open(assembler.file['output']) as in_file:
        for contig in SeqIO.parse(in_file, 'fasta'):
            contig_id = assembler.parse_contig_id(contig.description)
            if contig_id in all_hits:
                hit = all_hits[contig_id]
                batch.append((
                    assembler.state['iteration'],
                    contig.id,
                    str(contig.seq),
                    contig.description,
                    hit['bit_score'],
                    hit['len'],
                    hit['query_from'],
                    hit['query_to'],
                    hit['query_strand'],
                    hit['hit_from'],
                    hit['hit_to'],
                    hit['hit_strand']))
    db.insert_assembled_contigs_batch(assembler.state['db_conn'], batch)

    return high_score


def create_query_from_contigs(assembler):
    """Crate a new file with the contigs used as the next query."""
    log.info('Creating new query files: iteration {}'.format(
        assembler.state['iteration']))

    query = assembler.iter_file('long_reads.fasta')
    assembler.file['long_reads'] = query

    with open(query, 'w') as query_file:
        for row in db.get_assembled_contigs(
                assembler.state['db_conn'],
                assembler.state['iteration'],
                assembler.args['bit_score'],
                assembler.args['contig_length']):
            query_file.write('>{}\n'.format(row[0]))
            query_file.write('{}\n'.format(row[1]))

    return query


def parse_command_line(temp_dir_default):
    """Process command-line arguments."""
    description = """
        This is the aTRAM script. It takes a query sequence and a blast
        database built with the atram_preprocessor.py script and builds an
        assembly.

        If you specify more than one query sequence and/or more than one blast
        database then aTRAM will build one assembly for each query/blast
        DB pair.

        NOTE: You may use a text file to hold the command-line arguments
        like: @/path/to/args.txt. This is particularly useful when specifying
        multiple blast databases or multiple query sequences.
        """
    parser = argparse.ArgumentParser(
        fromfile_prefix_chars='@',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    parser.add_argument('--version', action='version',
                        version='%(prog)s {}'.format(db.ATRAM_VERSION))

    required_command_line_args(parser)
    optional_command_line_args(parser)
    filter_command_line_args(parser)
    blast.command_line_args(parser)
    assembly.command_line_args(parser)

    args = vars(parser.parse_args())

    # Check query arguments
    if not args['query'] and not args['query_split']:
        err = 'You must have at least one --query or --query-split argument.'
        log.fatal(err)

    # Set defaults and adjust arguments based on other arguments
    args['cov_cutoff'] = assembly.default_cov_cutoff(args['cov_cutoff'])
    args['blast_db'] = blast.touchup_blast_db_names(args['blast_db'])
    args['kmer'] = assembly.default_kmer(args['kmer'], args['assembler'])
    args['max_target_seqs'] = blast.default_max_target_seqs(
        args['max_target_seqs'], args['blast_db'], args['max_memory'])

    # Setup temp dir
    if not args['temp_dir']:
        args['temp_dir'] = temp_dir_default
    else:
        os.makedirs(args['temp_dir'], exist_ok=True)

    if args['no_filter']:
        args['bit_score'] = 0
        args['contig_length'] = 0

    if not args['protein'] and args['query']:
        args['protein'] = bio.fasta_file_has_protein(args['query'])

    # Prepend to PATH environment variable if requested
    if args['path']:
        os.environ['PATH'] = '{}:{}'.format(args['path'], os.environ['PATH'])

    find_programs(args['assembler'], args['no_long_reads'], args['bowtie2'])

    return args


def find_programs(assembler, no_long_reads, bowtie2):
    """Make sure we can find the programs needed by the assembler and blast."""
    if not (which('makeblastdb') and which('tblastn') and which('blastn')):
        err = ('We could not find the programs "makeblastdb", "tblastn", or '
               '"blastn". You either need to install them or you need adjust '
               'the PATH environment variable with the "--path" option so '
               'that aTRAM can find them.')
        sys.exit(err)

    if assembler == 'abyss' and not which('abyss-pe'):
        err = ('We could not find the "abyss-pe" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        sys.exit(err)

    if assembler == 'abyss' and not no_long_reads and not which('bwa'):
        err = ('We could not find the "bwa-mem" program. You either need to '
               'install it, adjust the PATH environment variable '
               'with the "--path" option, or you may use the '
               '"--no-long-reads" option to not use this program.')
        sys.exit(err)

    if assembler == 'trinity' and not which('Trinity'):
        err = ('We could not find the "Trinity" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        sys.exit(err)

    if assembler == 'trinity' and bowtie2 and not which('bowtie2'):
        err = ('We could not find the "bowtie2" program. You either need to '
               'install it, adjust the PATH environment variable '
               'with the "--path" option, or you may skip using this program '
               'by not using the "--bowtie2" option.')
        sys.exit(err)

    if assembler == 'velvet' and not (which('velveth') and which('velvetg')):
        err = ('We could not find either the "velveth" or "velvetg" program. '
               'You either need to install it or you need to adjust the PATH '
               'environment variable with the "--path" option so that aTRAM '
               'can find it.')
        sys.exit(err)

    if assembler == 'spades' and not which('spades.py'):
        err = ('We could not find the "Spades" program. You either need to '
               'install it or you need to adjust the PATH environment '
               'variable with the "--path" option so that aTRAM can find it.')
        sys.exit(err)


def required_command_line_args(parser):
    """Add required arguments to the parser."""
    group = parser.add_argument_group('required arguments')

    group.add_argument('-b', '--blast-db', '--sra', '--db', '--database',
                       required=True, metavar='DB', nargs='+',
                       help='''This needs to match the DB prefix you
                            entered for atram_preprocessor.py. You may repeat
                            this argument to run the --query sequence(s)
                            against multiple blast databases.''')

    group.add_argument('-q', '--query', '--target', '--probe',
                       required=False, nargs='+',
                       help='''The path to the fasta file with sequences of
                            interest. You may repeat this argument. If you do
                            then Each --query sequence  file will be run
                            against every --blast-db.''')

    group.add_argument('-Q', '--query-split', '--target-split',
                       required=False, nargs='+',
                       help='''The path to the fasta file with multiple
                            sequences of interest. This will take every
                            sequence in the fasta file and treat it as if it
                            were its own --query argument. So every sequence in
                            --query-split will be run against every --blast-db.
                            ''')

    group.add_argument('-o', '--output-prefix', required=True,
                       help='''This is the prefix of all of the output files.
                            So you can identify different blast output file
                            sets. You may include a directory as part of the
                            prefix. aTRAM will add suffixes to differentiate
                            ouput files.''')


def optional_command_line_args(parser):
    """Add optional atram arguments to the parser."""
    group = parser.add_argument_group('optional aTRAM arguments')

    group.add_argument('-a', '--assembler', default='none',
                       choices=['abyss', 'trinity', 'velvet', 'spades',
                                'none'],
                       help='''Which assembler to use. Choosing "none" (the
                            default) will do a single blast run and stop before
                            any assembly.''')

    group.add_argument('-i', '--iterations', type=int, default=5, metavar='N',
                       help='''The number of pipeline iterations.
                            The default is "5".''')

    group.add_argument('-p', '--protein', action='store_true',
                       help='''Are the query sequences protein?
                            aTRAM will guess if you skip this argument.''')

    group.add_argument('--fraction', type=float, default=1.0,
                       help='''Use only the specified fraction of the aTRAM
                            database. The default is "1.0"''')

    cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
    group.add_argument('--cpus', '--processes', '--max-processes',
                       type=int, default=cpus,
                       help='''Number of CPU threads to use. This will also be
                            used for the assemblers when possible. On this
                            machine the default is ("{}")'''.format(cpus))

    group.add_argument('--log-file',
                       help='''Log file (full path)."''')

    group.add_argument('--path',
                       help='''If the assembler or blast you want to use is not
                            in your $PATH then use this to prepend
                            directories to your path.''')

    group.add_argument('-t', '--temp-dir', metavar='DIR',
                       help='''You may save intermediate files for debugging
                            in this directory. The directory must be empty.''')

    group.add_argument('-T', '--timeout', metavar='SECONDS', default=300,
                       type=int,
                       help='''How many seconds to wait for an assembler before
                            stopping the run. To wait forever set this to 0.
                            The default is "300" (5 minutes).''')


def filter_command_line_args(parser):
    """Add optional values for blast-filtering contigs to the parser."""
    group = parser.add_argument_group(
        'optional values for blast-filtering contigs')

    group.add_argument('--no-filter', action='store_true',
                       help='''Do not filter the assembled contigs. This will:
                            set both the --bit-score and --contig-length
                            to 0''')

    group.add_argument('--bit-score', type=float, default=70.0,
                       metavar='SCORE',
                       help='''Remove contigs that have a value less than this.
                            The default is "70.0". This is turned off by the
                            --no-filter argument.''')

    group.add_argument('--contig-length', '--length', type=int, default=100,
                       help='''Remove blast hits that are shorter than this
                            length. The default is "100". This is turned
                            off by the --no-filter argument.''')


if __name__ == '__main__':

    with tempfile.TemporaryDirectory(prefix='atram_') as TEMP_DIR_DEFAULT:
        ARGS = parse_command_line(TEMP_DIR_DEFAULT)
        assemble(ARGS)
