"""Build assemblies using the aTRAM algorithm.."""

import re
import os
from os.path import basename, splitext, join
import multiprocessing
from tempfile import TemporaryDirectory
from Bio import SeqIO
import lib.db as db
import lib.log as log
import lib.util as util
import lib.blast as blast
import lib.assembler as assembly


def assemble(args):
    """Loop thru every blast/query pair and run an assembly for each one."""
    queries = split_queries(args)

    with TemporaryDirectory(prefix='atram_', dir=args['temp_dir']) as temp_dir:
        util.update_temp_dir(temp_dir, args)

        for blast_db in args['blast_db']:

            with db.connect(blast_db, check_version=True) as cxn:
                for query in queries:
                    db.aux_db(cxn, args['temp_dir'], blast_db, query)
                    clean_database(cxn)

                    log.setup(args['log_file'], blast_db, query)

                    assembler = assembly.factory(args, cxn)

                    try:
                        assembly_loop(assembler, blast_db, query)
                    except (TimeoutError, RuntimeError):
                        pass
                    except Exception as err:  # pylint: disable=broad-except
                        log.error('Exception: {}'.format(err))
                    finally:
                        assembler.write_final_output(blast_db, query)

                    db.aux_detach(cxn)


def assembly_loop(assembler, blast_db, query):
    """Iterate over the assembly processes."""
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


def clean_database(cxn):
    """Create database tables for an atram run."""
    db.create_sra_blast_hits_table(cxn)
    db.create_contig_blast_hits_table(cxn)
    db.create_assembled_contigs_table(cxn)


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
                   for shard in all_shards]
        all_results = [result.get() for result in results]
    log.info('All {} blast results completed'.format(len(all_results)))


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

    with db.connect(state['blast_db']) as cxn:
        db.aux_db(
            cxn,
            args['temp_dir'],
            state['blast_db'],
            state['query_target'])

        shard = basename(shard)

        batch = []

        hits = blast.hits(output_file)
        for hit in hits:
            seq_name, seq_end = blast.parse_blast_title(hit['title'])
            batch.append((state['iteration'], seq_end, seq_name, shard))
        db.insert_blast_hit_batch(cxn, batch)


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

    db.insert_contig_hit_batch(assembler.state['cxn'], batch)


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
    db.insert_assembled_contigs_batch(assembler.state['cxn'], batch)

    return high_score


def create_query_from_contigs(assembler):
    """Crate a new file with the contigs used as the next query."""
    log.info('Creating new query files: iteration {}'.format(
        assembler.state['iteration']))

    query = assembler.iter_file('long_reads.fasta')
    assembler.file['long_reads'] = query

    with open(query, 'w') as query_file:
        for row in db.get_assembled_contigs(
                assembler.state['cxn'],
                assembler.state['iteration'],
                assembler.args['bit_score'],
                assembler.args['contig_length']):
            query_file.write('>{}\n'.format(row[0]))
            query_file.write('{}\n'.format(row[1]))

    return query
