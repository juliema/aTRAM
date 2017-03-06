"""The aTRAM assembly program."""

import os
import csv
import logging
import argparse
import textwrap
import tempfile
import multiprocessing
from Bio import SeqIO
import lib.db as db
import lib.blast as blast
import lib.logger as logger
from lib.assembler import Assembler


def run(args):
    """Setup  and run atram."""

    logger.setup(args.work_dir, args.blast_db)

    all_shards = blast.all_shard_paths(args.work_dir, args.blast_db)
    assembler = Assembler.factory(args)

    db_conn = db.connect(args.work_dir, args.blast_db)
    db.create_blast_hits_table(db_conn)
    db.create_assembled_contigs_table(db_conn)
    db_conn.close()

    with tempfile.TemporaryDirectory(dir=args.work_dir) as temporary_dir:
        temp_dir = os.path.join(args.work_dir, 'temp_dir')  # ###############
        os.makedirs(temp_dir, exist_ok=True)                # ###############
        temp_dir = os.path.abspath(temp_dir)
        # temp_dir = os.path.abspath(temporary_dir)
        atram_loop(args, assembler, all_shards, temp_dir)
        # output_results()


def atram_loop(args, assembler, all_shards, temp_dir):
    """The run program loop."""

    query = args.query

    for iteration in range(1, args.iterations + 1):
        logging.info('aTRAM iteration %i', iteration)

        assembler.iteration_files(temp_dir, iteration)

        blast_target_against_all_sras(
            args, temp_dir, query, all_shards, iteration)

        write_assembler_files(args, assembler, iteration)

        # cwd = os.getcwd()
        # try:
        #     os.chdir(args.work_dir)  # some assemblers need this
        assembler.assemble()
        # finally:
        #     os.chdir(cwd)

        filter_contigs(args, assembler, temp_dir, iteration)

        # query = create_targets_from_contigs()
        # assembler.close_files()


def blast_target_against_all_sras(
        args, temp_dir, query, all_shards, iteration):
    """Blast the targets against the SRA databases. We're using a
    map-reduce strategy here. We map the blasting of the query sequences
    and reduce the output into one fasta file.
    """

    with multiprocessing.Pool(processes=args.cpus) as pool:
        results = [pool.apply_async(
            blast_target_against_sra,
            (dict(vars(args)), shard_path, query, temp_dir, iteration))
                   for shard_path in all_shards]
        _ = [result.get() for result in results]


def blast_target_against_sra(args, shard_path, query, temp_dir, iteration):
    """Blast the target against one blast DB shard. Then write the results to
    the database.
    """
    # NOTE: Because this is called in a child process, the address space is not
    # shared with the parent (caller) hence we cannot use object variables.

    output_file = blast.output_file(temp_dir, shard_path, iteration)

    blast.against_sra(args, shard_path, query, output_file, iteration)

    db_conn = db.connect(args['work_dir'], args['blast_db'])

    batch = []
    with open(output_file) as blast_hits:
        for line in blast_hits:
            match = blast.PARSE_RESULTS.match(line)
            if match:
                seq_name = match.group(1)
                seq_end = match.group(2)
            else:
                seq_name = line
                seq_end = ''
            batch.append((iteration, seq_end, seq_name, shard_path))
    db.insert_blast_hit_batch(db_conn, batch)


def write_assembler_files(args, assembler, iteration):
    """Take the matching blast hits and write the sequence and any matching
    end to the appropriate fasta files.
    """

    db_conn = db.connect(args.work_dir, args.blast_db)

    assembler.is_paired = False

    with open(assembler.paired_end_1_file, 'w') as end_1, \
            open(assembler.paired_end_2_file, 'w') as end_2:

        for row in db.get_blast_hits(db_conn, iteration):

            # NOTE: Some assemblers require a slash delimiter for the seq_end
            seq_end = ''
            if row['seq_end']:
                seq_end = '/{}'.format(row['seq_end'][-1])

            if seq_end.endswith('2'):
                assembler.is_paired = True
                out_file = end_2
            else:
                out_file = end_1

            out_file.write('>{}{}\n'.format(row['seq_name'], seq_end))
            out_file.write('{}\n'.format(row['seq']))

    db_conn.close()


def filter_contigs(args, assembler, temp_dir, iteration):
    """Remove junk from the assembled contigs."""

    blast_db = blast.temp_db(temp_dir, args.blast_db, iteration)
    output_file = blast.output_file(temp_dir, args.blast_db,
                                    'blast_contigs', iteration)

    blast.create_db(assembler.output_file, blast_db)

    blast.against_contigs(args, blast_db, args.query, output_file)

    filtered_scores = filter_contig_scores(args, assembler)
    save_contigs(args, assembler, filtered_scores, iteration)


def filter_contig_scores(args, assembler):
    """Only save contigs that have bit scores above the cut-off."""

    # qseqid sseqid bitscore qstart qend sstart send slen
    field_names = ['target_id', 'contig_id', 'bit_score',
                   'target_start', 'target_end',
                   'contig_start', 'contig_end', 'contig_len']

    scores = {}
    with open(assembler.output_file) as in_file:
        for score in csv.DictReader(in_file, field_names):
            score['bit_score'] = float(score['bit_score'])
            if score['bit_score'] >= args.bit_score:
                for field in field_names[3:]:
                    score[field] = int(score[field])
                scores[score['contig_id']] = score
    return scores


def save_contigs(args, assembler, filtered_scores, iteration):
    """Save the contigs to the database."""

    db_conn = db.connect(args.work_dir, args.blast_db)
    batch = []
    with open(assembler.output_file) as in_file:
        for contig in SeqIO.parse(in_file, 'fasta'):
            if contig.id in filtered_scores:
                score = filtered_scores[contig.id]
                batch.append((
                    iteration, contig.id,
                    str(contig.seq), contig.description,
                    score['bit_score'], score['target_start'],
                    score['target_end'], score['contig_start'],
                    score['contig_end'], score['contig_len']))
    db.insert_assembled_contigs_batch(db_conn, batch)
    db_conn.close()


def create_targets_from_contigs():
    """Crate a new file with the contigs that will be used as the
    next query target.
    """

    # query = filer.target_file(iteration)
    # with open(query, 'w') as target_file:
    #     for row in db.get_assembled_contigs(db_conn, iteration):
    #         target_file.write('>{}\n'.format(row[0]))
    #         target_file.write('{}\n'.format(row[1]))
    #
    # return query


def output_results():
    """Write the assembled contigs to a fasta file."""

    # file_name = filer.output_result_name(args.output_prefix)
    # with open(file_name, 'w') as out_file:
    #     for row in db.get_all_assembled_contigs(db_conn):
    #         header = ('>{}_{} iteration={} contig_id={} '
    #                   'score={}\n').format(
    #                       row['iteration'], row['contig_id'],
    #                       row['iteration'], row['contig_id'], row['bit_score'])
    #         out_file.write(header)
    #         out_file.write('{}\n'.format(row['seq']))


def parse_command_line():
    """Process command-line arguments."""

    description = """This is the aTRAM script. It takes a target sequence and
    a set of blast databases built with the atram_preprocessor.py script and
    builds an assembly.
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    group = parser.add_argument_group('required arguments')

    group.add_argument('-q', '--query', required=True,
                       help='The path to the fasta file with sequences of '
                            'interest.')

    group.add_argument('-d', '--work-dir', default='.', metavar='DIR',
                       required=True,
                       help='Where to find the BLAST_DB files. aTRAM will '
                            'also use this as a place to store temporary '
                            'files.')

    group.add_argument('-b', '--blast-db', required=True,
                       help='This needs to match the BLAST_DB you '
                            'entered for atram_preprocessor.py.')

    group.add_argument('-o', '--output', required=True,
                       help='Output the aTRAM results to this file.')

    group.add_argument('-a', '--assembler', required=True,
                       choices=['abyss', 'trinity', 'velvet'],
                       help='Which assembler to use.')

    group = parser.add_argument_group('optional aTRAM arguments')

    group.add_argument('-p', '--protein', action='store_true',
                       help='Are the query sequences protein? '
                            'The default is False.')

    group.add_argument('-i', '--iterations', type=int, default=5, metavar='N',
                       help='The number of pipline iterations. '
                            'The default is 5.')

    cpus = os.cpu_count() - 4 if os.cpu_count() > 4 else 1
    group.add_argument('-c', '--cpus', type=int, default=cpus,
                       help=('Number of cpus to use. The default is {}. '
                             'This will also be used for the assemblers '
                             'when possible.').format(cpus))

    group = parser.add_argument_group('optional blast arguments')

    group.add_argument('-s', '--bit-score', type=float, default=70.0,
                       metavar='SCORE',
                       help='Remove contigs that have a value less than this. '
                            'The default is 70.0')

    group.add_argument('-e', '--evalue', type=float, default=1e-9,
                       help='The default evalue is 1e-9.')

    group.add_argument('-M', '--max-target-seqs', type=int, default=100000000,
                       metavar='MAX',
                       help='Maximum hit sequences per shard. '
                            'Default is 100000000.')

    group.add_argument('-g', '--db-gencode', type=int, default=1,
                       metavar='CODE',
                       help='The genetic code to use during blast runs. '
                            'The default is 1.')

    group = parser.add_argument_group('optional assembler arguments')

    group.add_argument('-k', '--kmer', type=int, default=31,
                       help='k-mer size for assembers that use it. '
                            'The default is 31. (Abyss)')

    group.add_argument('-m', '--max_memory', default='50G', metavar='MEMORY',
                       help='Maximum amount of memory to use. The default is '
                            '50G. (Trinity)')

    args = parser.parse_args()

    # Set to absolute paths because we change directory during assembly
    args.work_dir = os.path.abspath(args.work_dir)
    args.output = os.path.abspath(args.output)
    args.query = os.path.abspath(args.query)

    return args


if __name__ == '__main__':

    ARGS = parse_command_line()
    run(ARGS)
