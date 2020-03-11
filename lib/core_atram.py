"""Build assemblies using the aTRAM algorithm.."""

import re
import os
from os.path import basename, split, splitext, join
from multiprocessing import Pool
from Bio import SeqIO
from . import db
from . import db_atram
from . import log
from . import bio
from . import util
from . import blast
from . import assembler as assembly


def assemble(args):
    """Loop thru every blast/query pair and run an assembly for each one."""
    with util.make_temp_dir(
            where=args['temp_dir'],
            prefix='atram_',
            keep=args['keep_temp_dir']) as temp_dir:
        util.update_temp_dir(temp_dir, args)

        queries = split_queries(args)

        for blast_db in args['blast_db']:

            with db.connect(blast_db, check_version=True) as cxn:
                for query in queries:
                    db.aux_db(cxn, args['temp_dir'], blast_db, query)
                    clean_database(cxn)

                    log.setup(
                        args['log_file'], args['log_level'], blast_db, query)

                    assembler = assembly.factory(args, cxn)

                    try:
                        assembly_loop(args, assembler, blast_db, query)
                    except (TimeoutError, RuntimeError):
                        pass
                    except Exception as err:  # pylint: disable=broad-except
                        log.error('Exception: {}'.format(err))
                    finally:
                        assembler.write_final_output(blast_db, query)

                    db.aux_detach(cxn)


def assembly_loop(args, assembler, blast_db, query):
    """Iterate over the assembly processes."""
    for iteration in range(1, assembler.args['iterations'] + 1):
        log.info('aTRAM blast DB = "{}", query = "{}", iteration {}'.format(
            blast_db, split(query)[1], iteration))

        assembler.init_iteration(blast_db, query, iteration)

        with util.make_temp_dir(
                where=args['temp_dir'],
                prefix=assembler.file_prefix(),
                keep=args['keep_temp_dir']) as iter_dir:

            assembler.setup_files(iter_dir)

            query = assembly_loop_iteration(args, assembler)

            if not query:
                break

    else:
        log.info('All iterations completed')


def assembly_loop_iteration(args, assembler):
    """One iteration of the assembly loop."""
    blast_query_against_all_shards(assembler)

    count = assembler.count_blast_hits()
    if assembler.blast_only or count == 0:
        return False

    assembler.write_input_files()

    assembler.run()

    if assembler.nothing_assembled():
        return False

    high_score = filter_contigs(assembler)

    count = assembler.assembled_contigs_count(high_score)

    if not count:
        return False

    if assembler.no_new_contigs(count):
        return False

    return create_query_from_contigs(args, assembler)


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

    if not args.get('protein'):
        args['protein'] = bio.fasta_file_has_protein(queries)

    return queries


def write_query_seq(file_name, seq_id, seq):
    """Write the sequence to a fasta file."""
    with open(file_name, 'w') as query_file:
        util.write_fasta_record(query_file, seq_id, seq)


def clean_database(cxn):
    """Create database tables for an atram run."""
    db_atram.create_sra_blast_hits_table(cxn)
    db_atram.create_contig_blast_hits_table(cxn)
    db_atram.create_assembled_contigs_table(cxn)


def blast_query_against_all_shards(assembler):
    """
    Blast the query against the SRA databases.

    We're using a map-reduce strategy here. We map the blasting of the query
    sequences and reduce the output into one fasta file.
    """
    log.info('Blasting query against shards: iteration {}'.format(
        assembler.state['iteration']))

    all_shards = shard_fraction(assembler)

    with Pool(processes=assembler.args['cpus']) as pool:
        results = [pool.apply_async(
            blast_query_against_one_shard,
            (assembler.args, assembler.simple_state(), shard))
                   for shard in all_shards]
        all_results = [result.get() for result in results]

    insert_blast_results(all_shards, assembler.args, assembler.simple_state())
    log.info('All {} blast results completed'.format(len(all_results)))


def insert_blast_results(all_shards, args, state):
    """Add all blast results to the auxiliary  database."""
    with db.connect(state['blast_db']) as cxn:
        db.aux_db(
            cxn,
            args['temp_dir'],
            state['blast_db'],
            state['query_target'])

        for shard in all_shards:
            shard = basename(shard)

            batch = []
            output_file = blast.output_file_name(state['iter_dir'], shard)

            hits = blast.hits(output_file)
            is_single_end = db.is_single_end(cxn)
            for hit in hits:
                seq_name, seq_end = blast.parse_blast_title(
                    hit['title'], is_single_end)
                batch.append((state['iteration'], seq_end, seq_name, shard))
            db_atram.insert_blast_hit_batch(cxn, batch)

        db.aux_detach(cxn)


def blast_query_against_one_shard(args, state, shard):
    """Blast the query against one blast DB shard."""
    output_file = blast.output_file_name(state['iter_dir'], shard)
    blast.against_sra(args, state, output_file, shard)


def shard_fraction(assembler):
    """
    Get the shards we are using.

    We may not want the entire DB for highly redundant libraries.
    """
    all_shards = blast.all_shard_paths(assembler.state['blast_db'])
    last_index = int(len(all_shards) * assembler.args['fraction'])
    return all_shards[:last_index]


def filter_contigs(assembler):
    """Remove junk from the assembled contigs."""
    log.info('Saving assembled contigs: iteration {}'.format(
        assembler.state['iteration']))

    blast_db = blast.temp_db_name(
        assembler.state['iter_dir'], assembler.state['blast_db'])

    hits_file = blast.output_file_name(
        assembler.state['iter_dir'], assembler.state['blast_db'])

    blast.create_db(
        assembler.state['iter_dir'], assembler.file['output'], blast_db)

    blast.against_contigs(
        blast_db,
        assembler.state['query_target'],
        hits_file,
        protein=assembler.args['protein'],
        db_gencode=assembler.args['db_gencode'],
        temp_dir=assembler.args['temp_dir'],
        timeout=assembler.args['timeout'])

    save_blast_against_contigs(assembler, hits_file)

    all_hits = {row['contig_id']: row
                for row
                in db_atram.get_contig_blast_hits(
                    assembler.state['cxn'],
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

    db_atram.insert_contig_hit_batch(assembler.state['cxn'], batch)


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
    db_atram.insert_assembled_contigs_batch(assembler.state['cxn'], batch)

    return high_score


def create_query_from_contigs(args, assembler):
    """Crate a new file with the contigs used as the next query."""
    log.info('Creating new query files: iteration {}'.format(
        assembler.state['iteration']))

    query_dir = join(args['temp_dir'], 'queries')
    os.makedirs(query_dir, exist_ok=True)

    query_file = assembler.file_prefix() + 'long_reads.fasta'
    query = join(query_dir, query_file)
    assembler.file['long_reads'] = query

    with open(query, 'w') as query_file:
        for row in db_atram.get_assembled_contigs(
                assembler.state['cxn'],
                assembler.state['iteration'],
                assembler.args['bit_score'],
                assembler.args['contig_length']):
            util.write_fasta_record(query_file, row[0], row[1])

    return query
