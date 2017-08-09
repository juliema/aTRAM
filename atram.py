"""The aTRAM assembly program."""

import os
from os.path import basename, exists, getsize
import argparse
import textwrap
import tempfile
import datetime
import subprocess
import multiprocessing
from Bio import SeqIO
import lib.db as db
import lib.bio as bio
import lib.log as log
import lib.blast as blast
import lib.file_util as file_util
import lib.assembler as assembly


def main(args):

    """Loop thru every blast/query pair and run an assembly for each."""

    for blast_db in args['blast_db']:

        db_conn = db.connect(blast_db)

        all_shards = fraction_of_shards(blast_db, args['fraction'])

        if args['split_queries']:
            queries = split_queries(args)
        else:
            queries = unified_queries(args)

        for query in queries:
            print(query)


def split_queries(args):
    """Handle every record as a separate query taget. This means we also need
    to put each query record into its own file for the blast queries.
    """

    for query_file_name in args['query']:
        with open(query_file_name) as query_file:
            for query_rec in SeqIO.parse(query_file, 'fasta'):
                query_file_name = ''
                # Put the sequence into a separate file
                yield query_file_name


def unified_queries(args):
    """Handle each query file as a query target."""

    for query_file_name in args['query']:
        yield query_file_name


def main_old(args):
    """Setup and run atram for each blast-db/query pair. This loops thru all of
    the blast databases and every sequence in the query fasta files and then
    runs the atram main loop for each."""

    # Loop thru every blast DB setup with atram_preprocessor.py
    for blast_db in args['blast_db']:

        db_conn = db.connect(blast_db)

        # We may not want the entire DB for highly redundant libraries
        all_shards = fraction_of_shards(blast_db, args['fraction'])

        # Now loop thru every query fasta file
        for query_file in args['query']:

            # Set up the log file for the run
            log_file = log.file_name(args['log_file'], blast_db, query_file)
            log.setup(log_file)

            with open(query_file) as in_file:

                # Now loop thru every sequence in the query fasta file
                for query_rec in SeqIO.parse(in_file, 'fasta'):

                    # Get a new assembler for the sequence
                    assembler = None
                    if args['assembler']:
                        assembler = assembly.factory(args)

                    # Setup the query sequence file
                    query = initialize_query(
                        args, db_conn, assembler, str(query_rec.seq))

                    atram_loop(db_conn, assembler, query, all_shards, args)

                    output_prefix = get_output_prefix(
                        args, basename(blast_db), query_file, query_rec.id)
                    atram_output(db_conn, assembler, output_prefix, args)

        db_conn.close()


def atram_loop(db_conn, assembler, query, all_shards, args):
    """The main program loop."""

    for iteration in range(args['start_iteration'], args['iterations'] + 1):
        log.info('aTRAM iteration %i' % iteration, line_break='')

        blast_query_against_all_sras(query, all_shards, iteration, args)

        # If we don't have an assembler then we just want the blast hits
        if not assembler:
            break

        # Exit if there are no blast hits
        if not db.sra_blast_hits_count(db_conn, iteration):
            log.info('No blast hits in iteration %i' % iteration)
            break

        assembler.initialize_iteration(iteration)

        assembler.write_input_files(db_conn)

        try:
            log.info('Assembling shards with {}: iteration {}'.format(
                args['assembler'], iteration))
            assembler.assemble()
        except TimeoutError:
            msg = 'Time ran out for the assembler after {} (HH:MM:SS)'.format(
                datetime.timedelta(seconds=args['timeout']))
            log.fatal(msg)
        except subprocess.CalledProcessError as cpe:
            msg = 'The assembler failed with error: ' + str(cpe)
            log.fatal(msg)

        # Exit if nothing was assembled
        if not exists(assembler.file['output']) \
                or not getsize(assembler.file['output']):
            log.info('No new assemblies in iteration %i' % iteration)
            break

        high_score = filter_contigs(db_conn, assembler, iteration, args)

        count = db.assembled_contigs_count(
            db_conn, iteration, args['bit_score'], args['contig_length'])
        if not count:
            log.info(('No contigs had a bit score greater than {} and are at '
                      'least {} long in iteration {}. The highest score for '
                      'this iteration is {}').format(args['bit_score'],
                                                     args['contig_length'],
                                                     iteration,
                                                     high_score))
            break

        if count == db.iteration_overlap_count(
                db_conn, iteration, args['bit_score'], args['contig_length']):
            log.info('No new contigs were found in iteration %i' % iteration)
            break

        query = create_queries_from_contigs(
            db_conn, assembler, iteration, args)
    else:
        log.info('All iterations completed', line_break='')


def get_output_prefix(args, blast_db, query_file, seq_id):
    """Setup output file names."""

    if args['split_queries']:
        return '{}.{}_{}'.format(
            args['output_prefix'], basename(blast_db), seq_id)

    return '{}.{}_{}'.format(
        args['output_prefix'], basename(blast_db), query_file)


def initialize_query(args, db_conn, assembler, query):
    """Get the first set of query sequences. Handle a restart of atram by
    getting the data from the given iteration and setting up an input file.
    """

    if args['start_iteration'] == 1:
        db.create_sra_blast_hits_table(db_conn)
        db.create_contig_blast_hits_table(db_conn)
        db.create_assembled_contigs_table(db_conn)
    else:
        start_iteration = args['start_iteration'] - 1
        query = create_queries_from_contigs(
            db_conn, assembler, start_iteration, args)
    if getsize(query) < 10:
        err = 'There are no sequences '
        if args['start_iteration'] == 1:
            err += 'in {}'.format(args['query'])
        else:
            err += 'for starting iteration {}'.format(
                args['start_iteration'])
        log.fatal(err)

    return query


def fraction_of_shards(blast_db, fraction):
    """Get the fraction of shards to use. This is so that we can blast against
    a portion of the blast DBs we build. You may use this if the blast library
    is highly redundant.
    """

    all_shards = blast.all_shard_paths(blast_db)
    last_index = int(len(all_shards) * fraction)

    return all_shards[:last_index]


def blast_query_against_all_sras(query, all_shards, iteration, args):
    """Blast the queries against the SRA databases. We're using a
    map-reduce strategy here. We map the blasting of the query sequences
    and reduce the output into one fasta file.
    """

    log.info('Blasting query against shards: iteration %i' % iteration)

    with multiprocessing.Pool(processes=args['cpus']) as pool:
        results = [pool.apply_async(
            blast_query_against_sra,
            (shard_path, query, iteration, dict(vars(args))))
                   for shard_path in all_shards]
        _ = [result.get() for result in results]  # noqa


def blast_query_against_sra(shard_path, query, iteration, args):
    """Blast the query against one blast DB shard. Then write the results to
    the database.
    """

    output_file = blast.output_file_name(
        args['temp_dir'], shard_path, iteration)

    blast.against_sra(args, shard_path, query, output_file, iteration)

    db_conn = db.connect(args['blast_db'])

    shard = basename(shard_path)

    batch = []

    for hit in blast.hits(output_file):
        match = blast.PARSE_RESULTS.match(hit['title'])
        if match:
            seq_name = match.group(1)
            seq_end = match.group(2)
        else:
            seq_name = hit['title']
            seq_end = ''
        batch.append((iteration, seq_end, seq_name, shard))
    db.insert_blast_hit_batch(db_conn, batch)
    db_conn.close()


def filter_contigs(db_conn, assembler, iteration, args):
    """Remove junk from the assembled contigs."""

    log.info('Saving assembled contigs: iteration %i' % iteration)

    blast_db = blast.temp_db_name(
        args['temp_dir'], args['blast_db'], iteration)
    hits_file = blast.output_file_name(
        args['temp_dir'], args['blast_db'], iteration)

    blast.create_db(args['temp_dir'], assembler.file['output'], blast_db)

    blast.against_contigs(blast_db,
                          args['query'],
                          hits_file,
                          protein=args['protein'],
                          db_gencode=args['db_gencode'],
                          temp_dir=args['temp_dir'])

    save_blast_against_contigs(db_conn, assembler, hits_file, iteration)

    all_hits = {row['contig_id']: row
                for row
                in db.get_contig_blast_hits(db_conn, iteration)}

    return save_contigs(db_conn, all_hits, assembler, iteration)


def save_blast_against_contigs(db_conn, assembler, hits_file, iteration):
    """Save all of the blast hits."""

    batch = []

    for hit in blast.hits(hits_file):
        contig_id = assembler.parse_contig_id(hit['title'])
        batch.append((iteration, contig_id, hit['title'],
                      hit['bit_score'], hit['len'],
                      hit['query_from'], hit['query_to'],
                      hit.get('query_strand', ''),
                      hit['hit_from'], hit['hit_to'],
                      hit.get('hit_strand', '')))

    db.insert_contig_hit_batch(db_conn, batch)


def save_contigs(db_conn, all_hits, assembler, iteration):
    """Save the contigs to the database."""

    batch = []
    high_score = 0
    with open(assembler.file['output']) as in_file:
        for contig in SeqIO.parse(in_file, 'fasta'):
            contig_id = assembler.parse_contig_id(contig.description)
            if contig_id in all_hits:
                hit = all_hits[contig_id]
                batch.append((
                    iteration, contig.id, str(contig.seq), contig.description,
                    hit['bit_score'], hit['len'],
                    hit['query_from'], hit['query_to'], hit['query_strand'],
                    hit['hit_from'], hit['hit_to'], hit['hit_strand']))
    db.insert_assembled_contigs_batch(db_conn, batch)

    return high_score


def create_queries_from_contigs(db_conn, assembler, iteration, args):
    """Crate a new file with the contigs that will be used as the
    next query query.
    """

    log.info('Creating new query files: iteration %i' % iteration)

    query = assembler.path('long_reads.fasta')
    assembler.file['long_reads'] = query

    with open(query, 'w') as query_file:
        for row in db.get_assembled_contigs(
                db_conn, iteration, args['bit_score'], args['length']):
            query_file.write('>{}\n'.format(row[0]))
            query_file.write('{}\n'.format(row[1]))

    return query


def atram_output(db_conn, assembler, output_prefix, args):
    """Output the results of the assembly or blast."""

    if assembler:
        output_assembly_results(db_conn, output_prefix, args)
    else:
        output_blast_only_results(db_conn, output_prefix)


def output_assembly_results(db_conn, output_prefix, args):
    """Write the assembled contigs to a fasta file."""

    if not args['no_filter']:
        file_name = file_util.output_file(
            output_prefix, 'filtered_contigs.fasta')
        with open(file_name, 'w') as out_file:
            for row in db.get_all_assembled_contigs(
                    db_conn, args['bit_score'], args['contig_length']):
                output_one_assembly(out_file, row)

    file_name = file_util.output_file(output_prefix, 'all_contigs.fasta')
    with open(file_name, 'w') as out_file:
        for row in db.get_all_assembled_contigs(db_conn):
            output_one_assembly(out_file, row)


def output_one_assembly(out_file, row):
    """Write an assembly to the output fasta file."""

    seq = row['seq']
    suffix = ''
    if row['query_strand'] and row['hit_strand'] and \
            row['query_strand'] != row['hit_strand']:
        seq = bio.reverse_complement(seq)
        suffix = '_REV'

    header = ('>{}_{}{} iteration={} contig_id={} '
              'score={}\n').format(
                  row['iteration'], row['contig_id'], suffix,
                  row['iteration'], row['contig_id'], row['bit_score'])
    out_file.write(header)
    out_file.write('{}\n'.format(seq))


def output_blast_only_results(db_conn, output_prefix):
    """Output this file if we are not assembling the contigs."""

    log.info('Output blast only results')

    file_name = file_util.output_file(output_prefix, 'blast_only.fasta')
    with open(file_name, 'w') as out_file:
        for row in db.get_sra_blast_hits(db_conn, 1):
            out_file.write('>{}{}\n'.format(row['seq_name'], row['seq_end']))
            out_file.write('{}\n'.format(row['seq']))


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
                        version='%(prog)s {}'.format(file_util.VERSION))

    required_command_line_args(parser)
    optional_command_line_args(parser)
    filter_command_line_args(parser)
    blast.command_line_args(parser)
    assembly.command_line_args(parser)

    args = vars(parser.parse_args())

    args['cov_cutoff'] = assembly.default_cov_cutoff(args['cov_cutoff'])
    args['blast_db'] = blast.touchup_blast_db_names(args['blast_db'])
    args['kmer'] = assembly.default_kmer(args['kmer'], args['assembler'])
    args['temp_dir'] = file_util.temp_root_dir(
        args['temp_dir'], temp_dir_default)
    args['max_target_seqs'] = blast.default_max_target_seqs(
        args['max_target_seqs'], args['blast_db'], args['max_memory'])

    if args['no_filter']:
        args['bit_score'] = 0
        args['contig_length'] = 0

    # Check query
    if not args['query'] and not args['start_iteration']:
        log.fatal('We need at least one "--query" sequence.')

    # Prepend to PATH environment variable if requested
    if args['path']:
        os.environ['PATH'] = '{}:{}'.format(args['path'], os.environ['PATH'])

    file_util.find_programs(
        args['assembler'], args['no_long_reads'], args['bowtie2'])

    return args


def required_command_line_args(parser):
    """Add required arguments to the parser."""

    group = parser.add_argument_group('required arguments')

    group.add_argument('-b', '--blast-db', '--sra', '--db', '--database',
                       required=True, metavar='DB', action='append',
                       help='''This needs to match the DB prefix you
                            entered for atram_preprocessor.py. You may repeat
                            this argument to run the --query sequence(s)
                            against multiple blast databases.''')

    group.add_argument('-q', '--query', '--target', '--probe',
                       required=True, action='append',
                       help='''The path to the fasta file with sequences of
                            interest. Required unless you specify a
                            "--start-iteration". You may have multiple query
                            sequences in one file or you may repeat this
                            argument. In either case, each --query sequence
                            will be run against every --blast-db.''')

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
                       help='''Number of cpus to use. This will also be used
                            for the assemblers when possible. Defaults to:
                            Total CPUs - 4 = "{}"'''.format(cpus))

    group.add_argument('--log-file',
                       help='''Log file (full path)."''')

    group.add_argument('--path',
                       help='''If the assembler or blast you want to use is not
                            in your $PATH then use this to prepend
                            directories to your path.''')

    group.add_argument('--split-queries', action='store_true',
                       help='''Split every record in the --query file into
                            separate query sequence, as if they were in
                            separate files.''')

    # group.add_argument('--restart', type='store_true',
    #                    help='''Use this to resume from a previous run.''')

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
        main(ARGS)
