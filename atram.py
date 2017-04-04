"""The aTRAM assembly program."""

import os
import sys
import csv
import math
import logging
import argparse
import textwrap
import tempfile
import multiprocessing
from shutil import which
import psutil
from Bio import SeqIO
import lib.db as db
import lib.bio as bio
import lib.blast as blast
from lib.assembler import Assembler


def run(args):
    """Setup  and run atram."""

    logging.basicConfig(
        filename=args.log_file,
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(' '.join(sys.argv))

    all_shards = blast.all_shard_paths(args.work_dir, args.blast_db)

    db_conn = db.connect(args.work_dir, args.blast_db)
    db.create_blast_hits_table(db_conn)
    db.create_assembled_contigs_table(db_conn)

    # TODO: Setup starting iteration

    with tempfile.TemporaryDirectory(dir=args.work_dir) as temporary_dir:
        temp_dir = os.path.join(args.work_dir, temporary_dir)
        temp_dir = os.path.join(args.work_dir, 'temp_dir')  # TODO: Delete me
        os.makedirs(temp_dir, exist_ok=True)                # TODO: Delete me
        assembler = None
        if args.assembler:
            assembler = Assembler.factory(args, temp_dir)
        atram_loop(args, db_conn, assembler, all_shards, temp_dir)
        output_results(args, db_conn)

    db_conn.close()


def atram_loop(args, db_conn, assembler, all_shards, temp_dir):
    """The run program loop."""

    query = args.query

    for iteration in range(1, args.iterations + 1):
        logging.info('aTRAM iteration %i', iteration)

        blast_target_against_all_sras(
            args, temp_dir, query, all_shards, iteration)

        # If we don't have an assembler then we just want the blast hits
        if not args.assembler:
            output_blast_only_results(args, db_conn)
            sys.exit()

        # Exit if there are no blast hits
        if not db.blast_hits_count(db_conn, iteration):
            logging.info('No blast hits in iteration %i', iteration)
            break

        assembler.iteration_files(iteration)

        write_assembler_files(db_conn, assembler, iteration)

        try:
            assembler.assemble()
        except Exception as exn:  # pylint: disable=broad-except
            msg = 'The assembler failed: ' + str(exn)
            logging.error(msg)
            sys.exit(msg)

        # Exit if nothing was assembled
        if not os.path.exists(assembler.output_file) \
                or not os.path.getsize(assembler.output_file):
            logging.info('No new assemblies in iteration %i', iteration)
            break

        high_score = filter_contigs(
            args, db_conn, assembler, temp_dir, iteration)

        count = db.assembled_contigs_count(db_conn, iteration)
        if not count:
            logging.info(
                ('No contigs had a bit score greater than %i and are at '
                 'least %i long in iteration %i. The highest score for this '
                 'iteration is %d'),
                args.bit_score, args.length, iteration, high_score)
            break

        if count == db.iteration_overlap_count(db_conn, iteration):
            logging.info('No new contigs were found in iteration %i',
                         iteration)
            break

        # TODO: Exit if the target was covered

        query = create_targets_from_contigs(db_conn, assembler, iteration)


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
    # shared with the parent (caller) hence we cannot share object variables.

    output_file = blast.output_file(temp_dir, shard_path, iteration)

    blast.against_sra(args, shard_path, query, output_file, iteration)

    db_conn = db.connect(args['work_dir'], args['blast_db'])

    shard = os.path.basename(shard_path)

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
            batch.append((iteration, seq_end, seq_name, shard))
    db.insert_blast_hit_batch(db_conn, batch)
    db_conn.close()


def write_assembler_files(db_conn, assembler, iteration):
    """Take the matching blast hits and write the sequence and any matching
    end to the appropriate fasta files.
    """

    assembler.is_paired = False

    with open(assembler.ends_1_file, 'w') as end_1, \
            open(assembler.ends_2_file, 'w') as end_2:

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


def output_blast_only_results(args, db_conn):
    """Output this file if we are not assembling the contigs."""

    with open(args.output, 'w') as out_file:
        for row in db.get_blast_hits(db_conn, 1):
            out_file.write('>{}{}\n'.format(row['seq_name'], row['seq_end']))
            out_file.write('{}\n'.format(row['seq']))


def filter_contigs(args, db_conn, assembler, temp_dir, iteration):
    """Remove junk from the assembled contigs."""

    blast_db = blast.temp_db(temp_dir, args.blast_db, iteration)
    hits_file = blast.output_file(temp_dir, args.blast_db, iteration)

    blast.create_db(assembler.output_file, blast_db)

    blast.against_contigs(args, blast_db, args.query, hits_file)

    filtered_scores = filter_contig_scores(args, hits_file)
    return save_contigs(db_conn, assembler, filtered_scores, iteration)


def filter_contig_scores(args, hits_file):
    """Only save contigs that have bit scores above the cut-off."""

    # qseqid sseqid bitscore qstart qend sstart send slen
    field_names = ['target_id', 'contig_id', 'bit_score',
                   'target_start', 'target_end',
                   'contig_start', 'contig_end', 'contig_len']

    scores = {}
    with open(hits_file) as in_file:
        for score in csv.DictReader(in_file, field_names):
            score['bit_score'] = float(score['bit_score'])
            if score['bit_score'] >= args.bit_score:
                for field in field_names[3:]:
                    score[field] = int(score[field])
                scores[score['contig_id']] = score
    return scores


def save_contigs(db_conn, assembler, filtered_scores, iteration):
    """Save the contigs to the database."""

    batch = []
    high_score = 0
    with open(assembler.output_file) as in_file:
        for contig in SeqIO.parse(in_file, 'fasta'):
            if contig.id in filtered_scores:
                score = filtered_scores[contig.id]
                if score['bit_score'] > high_score:
                    high_score = score['bit_score']
                batch.append((
                    iteration, contig.id,
                    str(contig.seq), contig.description,
                    score['bit_score'], score['target_start'],
                    score['target_end'], score['contig_start'],
                    score['contig_end'], score['contig_len']))
    db.insert_assembled_contigs_batch(db_conn, batch)

    return high_score


def create_targets_from_contigs(db_conn, assembler, iteration):
    """Crate a new file with the contigs that will be used as the
    next query target.
    """

    query = assembler.path('long_reads.fasta', iteration)
    assembler.long_reads_file = query

    with open(query, 'w') as target_file:
        for row in db.get_assembled_contigs(db_conn, iteration):
            target_file.write('>{}\n'.format(row[0]))
            target_file.write('{}\n'.format(row[1]))

    return query


def output_results(args, db_conn):
    """Write the assembled contigs to a fasta file."""

    seen = {}
    with open(args.output, 'w') as out_file:
        for row in db.get_all_assembled_contigs(db_conn):

            if row['contig_id'] in seen:
                continue
            seen[row['contig_id']] = 1

            header = ('>{}_{} iteration={} contig_id={} '
                      'score={}\n').format(
                          row['iteration'], row['contig_id'],
                          row['iteration'], row['contig_id'], row['bit_score'])
            out_file.write(header)
            out_file.write('{}\n'.format(row['seq']))
            header = ('>{}_{}_REV iteration={} contig_id={} '
                      'score={}\n').format(
                          row['iteration'], row['contig_id'],
                          row['iteration'], row['contig_id'], row['bit_score'])
            out_file.write(header)
            out_file.write('{}\n'.format(bio.reverse_complement(row['seq'])))


def parse_command_line():  # pylint: disable=too-many-statements
    """Process command-line arguments."""

    description = """
        This is the aTRAM script. It takes a target sequence and
        a set of blast databases built with the atram_preprocessor.py script
        and builds an assembly.
        """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(description))

    # required arguments
    group = parser.add_argument_group('required arguments')

    group.add_argument('-b', '--blast-db', '--sra', '--db', '--database',
                       required=True, metavar='DB',
                       help='This needs to match the DB you '
                            'entered for atram_preprocessor.py.')

    group.add_argument('-o', '--output', required=True,
                       help='Output the aTRAM results to this file.')

    group.add_argument('-q', '--query', '--target', required=True,
                       help='The path to the fasta file with sequences of '
                            'interest.')

    # optional aTRAM arguments
    group = parser.add_argument_group('optional aTRAM arguments')

    group.add_argument('-a', '--assembler',
                       choices=['abyss', 'trinity', 'velvet'],
                       help='Which assembler to use. If you do not use this '
                            'argument then aTRAM will do a single blast run '
                            'and stop before assembly.')

    group.add_argument('-d', '--work-dir', default='.', metavar='DIR',
                       help=('Which directory has the files created by '
                             'atram_preprocessor.py. This will also be used '
                             'as a place to store temporary files if TEMP_DIR '
                             'is not specified. Defaults to the current '
                             'directory "{}".'.format(os.getcwd())))

    group.add_argument('-i', '--iterations', type=int, default=5, metavar='N',
                       help='The number of pipline iterations. '
                            'The default is "5".')

    group.add_argument('-p', '--protein', action='store_true',
                       help='Are the query sequences protein? '
                            'The aTRAM will guess if you skip this argument.')

    group.add_argument('--fraction', type=float, default=1.0,
                       help='Use only the specified fraction of the aTRAM '
                            'database. The default is "1.0"')

    group.add_argument('--complete', action='store_true',
                       help='Automatically quit when a complete homolog is '
                            'recovered.')

    cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
    group.add_argument('--cpus', '--processes', '--max-processes',
                       type=int, default=cpus,
                       help=('Number of cpus to use. This will also be used '
                             'for the assemblers when possible. Defaults to: '
                             'Total CPUS  - 4 = "{}"').format(cpus))

    group.add_argument('--log-file',
                       help='Log file (full path). The default is to use the '
                            'DIR and DB arguments to come up with a name like '
                            'so "DIR/DB_atram.log"')

    group.add_argument('--path',
                       help='If the assembler or blast you want to use is not '
                            'in your $PATH then use this to prepend '
                            'directories to your path.')

    group.add_argument('--start-iteration', type=int, default=5, metavar='N',
                       help='If resuming from a previous run, which iteration '
                            'number to start from. The default is "1".')

    group.add_argument('--temp-dir',
                       help='Store temporary files in this directory. '
                            'Temporary files will be deleted if you do not '
                            'specify this argument.')

    # optional values for blast-filtering contigs
    group = parser.add_argument_group(
        'optional values for blast-filtering contigs')

    group.add_argument('--bit-score', type=float, default=70.0,
                       metavar='SCORE',
                       help='Remove contigs that have a value less than this. '
                            'The default is "70.0"')

    group.add_argument('--length', '--contig-length', type=int, default=100,
                       help='Remove blast hits that are shorter than this '
                            'length. The default is "100".')

    # optional blast arguments
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
                            'Default is calulated based on the available '
                            'memory and the number of shards. ')

    # optional assembler arguments
    group = parser.add_argument_group('optional assembler arguments')

    group.add_argument('--no-long-reads', action='store_true',
                       help='Do not use long reads during assembly. '
                            '(Abyss, Trinity, Velvet)')

    group.add_argument('--kmer', type=int, default=64,
                       help='k-mer size. The default is "64" for Abyss and '
                            '"31" for Velvet. Note: the maximum kmer length '
                            'for Velvet is 31. (Abyss, Velvet)')

    group.add_argument('--bowtie2', action='store_true',
                       help='Use bowtie2 during assembly. (Trinity)')

    max_mem = max(1, math.floor(
        psutil.virtual_memory().available / 1024**3 / 2))
    group.add_argument('--max-memory', default=max_mem, metavar='MEMORY',
                       type=int,
                       help=('Maximum amount of memory to use in gigabytes. '
                             'The default is "{}". (Trinity)').format(max_mem))

    group.add_argument('--exp-coverage', '--expected_coverage',
                       type=int, default=30,
                       help='The expected coverage of the region. '
                            'The default is "30". (Velvet)')

    group.add_argument('--ins-length', type=int, default=300,
                       help='The size of the fragments used in the short-read '
                            'library. The default is "300". (Velvet)')

    group.add_argument('--min-contig-length', '--min-contig-len',
                       type=int, default=100,
                       help='The minimum contig length. '
                            'The default is "100". (Velvet)')

    args = parser.parse_args()

    # Check kmer
    if args.assembler == 'velvet' and args.kmer > 31:
        args.kmer = 31

    # Set default log file name
    if not args.log_file:
        file_name = '{}.{}.log'.format(args.blast_db, sys.argv[0][:-3])
        args.log_file = os.path.join(args.work_dir, file_name)

    # If not --protein then probe to see if it's a protein seq
    if not args.protein:
        with open(args.query) as in_file:
            for query in SeqIO.parse(in_file, 'fasta'):
                if bio.is_protein(str(query.seq)):
                    args.protein = True

    # Prepend to PATH environment variable if requested
    if args.path:
        os.environ['PATH'] = '{}:{}'.format(args.path, os.environ['PATH'])

    # Calculate the default max_target_seqs per shard
    if not args.max_target_seqs:
        all_shards = blast.all_shard_paths(args.work_dir, args.blast_db)
        args.max_target_seqs = int(2 * args.max_memory / len(all_shards)) * 1e6

    find_programs(args)

    return args


def find_programs(args):
    """Make sure we can find the programs needed by the assembler and blast."""

    if not (which('makeblastdb') and which('tblastn') and which('blastn')):
        print('We could not find the programs "makeblastdb", "tblastn", or '
              '"blastn". You either need to install them or you need adjust '
              'the PATH environment variable with the "--path" option so that '
              'aTRAM can find it.')
        sys.exit()

    if args.assembler == 'abyss' and not which('abyss-pe'):
        print('We could not find the "abyss-pe" program. You either need to '
              'install it or you need to adjust the PATH environment variable '
              'with the "--path" option so that aTRAM can find it.')
        sys.exit()

    if args.assembler == 'abyss' and not args.no_long_reads \
            and not which('bwa'):
        print('We could not find the "bwa-mem" program. You either need to '
              'install it, adjust the PATH environment variable '
              'with the "--path" option, or you may use the "--no-long-reads" '
              'option to not use this program.')
        sys.exit()

    if args.assembler == 'trinity' and not which('Trinity'):
        print('We could not find the "Trinity" program. You either need to '
              'install it or you need to adjust the PATH environment variable '
              'with the "--path" option so that aTRAM can find it.')
        sys.exit()

    if args.assembler == 'trinity' and args.bowtie2 and not which('bowtie2'):
        print('We could not find the "bowtie2" program. You either need to '
              'install it, adjust the PATH environment variable '
              'with the "--path" option, or you may skip using this program '
              'by not using the "--bowtie2" option.')
        sys.exit()

    if args.assembler == 'velvet' and \
            not (which('velveth') and which('velvetg')):
        print('We could not find either the "velveth" or "velvetg" program. '
              'You either need to install it or you need to adjust the PATH '
              'environment variable with the "--path" option so that aTRAM '
              'can find it.')
        sys.exit()


if __name__ == '__main__':

    ARGS = parse_command_line()
    run(ARGS)
