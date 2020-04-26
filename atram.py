#!/usr/bin/env python3
"""
Start atram.

This wrapper module parses the input arguments and passes them to the module
that does the actual processing (core_atram.py).
"""

import os
import argparse
import textwrap
import lib.db as db
import lib.log as log
import lib.bio as bio
import lib.util as util
import lib.blast as blast
import lib.assembler as assembly
from lib.core_atram import assemble


def parse_command_line():
    """Process command-line arguments."""
    description = """
        This  takes a query sequence and a blast database built with the
        atram_preprocessor.py script and builds assemblies.

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

    group = parser.add_argument_group('required arguments')

    group.add_argument(
        '-b', '--blast-db', '--sra', '--db', '--database',
        required=True, metavar='DB', nargs='+',
        help="""This needs to match the DB prefix you entered for
            atram_preprocessor.py. You may repeat this argument to run the
            --query sequence(s) against multiple blast databases.""")

    group.add_argument(
        '-q', '--query', '--target', '--probe', required=False, nargs='+',
        help="""The path to the fasta file with sequences of interest. You may
            repeat this argument. If you do then Each --query sequence file
            will be run against every --blast-db.""")

    group.add_argument(
        '-Q', '--query-split', '--target-split', required=False, nargs='+',
        help="""The path to the fasta file with multiple sequences of interest.
            This will take every sequence in the fasta file and treat it as if
            it were its own --query argument. So every sequence in
            --query-split will be run against every --blast-db.""")

    group.add_argument(
        '-o', '--output-prefix', required=True,
        help="""This is the prefix of all of the output files. So you can
            identify different blast output file sets. You may include a
            directory as part of the prefix. aTRAM will add suffixes to
            differentiate output files.""")

    group.add_argument(
        '-a', '--assembler', default='none',
        choices=['abyss', 'trinity', 'velvet', 'spades', 'none'],
        help="""Which assembler to use. Choosing "none" (the default) will do
            a single blast run and stop before any assembly.""")

    group.add_argument(
        '-i', '--iterations', type=int, default=5, metavar='N',
        help="""The number of pipeline iterations. The default is "5".""")

    group.add_argument(
        '-p', '--protein', action='store_true',
        help="""Are the query sequences protein? aTRAM will guess if you skip
            this argument.""")

    group.add_argument(
        '--fraction', type=float, default=1.0,
        help="""Use only the specified fraction of the aTRAM database. The
            default is 1.0.""")

    cpus = min(10, os.cpu_count() - 4 if os.cpu_count() > 4 else 1)
    group.add_argument(
        '--cpus', '--processes', '--max-processes', type=int, default=cpus,
        help="""Number of CPU processors to use. This will also be used for
            the assemblers when possible. We will use {} out of {} CPUs.
            """.format(cpus, os.cpu_count()))

    group.add_argument('--log-file', help="""Log file (full path)".""")
    group.add_argument(
        '--log-level', choices=['debug', 'info', 'error'], default='info',
        help="""Log messages of the given level (or above). 'debug' shows the
            most messages and 'error' shows the least. The default is
            'info'""")

    group.add_argument(
        '--path',
        help="""If the assembler or blast you want to use is not in your $PATH\
            then use this to prepend directories to your path.""")

    group.add_argument(
        '-t', '--temp-dir', metavar='DIR',
        help="""Place temporary files in this directory. All files will be
            deleted after aTRAM completes. The directory must exist.""")

    group.add_argument(
        '--keep-temp-dir', action='store_true',
        help="""This flag will keep the temporary files in the --temp-dir
            around for debugging.""")

    group.add_argument(
        '-T', '--timeout', metavar='SECONDS', default=600, type=int,
        help="""How many seconds to wait for an assembler or BLAST before
            stopping the run. To wait forever set this to 0. The default
            is "600" (10 minutes).""")

    group = parser.add_argument_group(
        'optional values for blast-filtering contigs')

    group.add_argument(
        '--no-filter', action='store_true',
        help="""Do not filter the assembled contigs. This will: set both the
            --bit-score and --contig-length to 0""")

    group.add_argument(
        '--bit-score', type=float, default=70.0, metavar='SCORE',
        help="""Remove contigs that have a value less than this. The default
            is "70.0". This is turned off by the --no-filter argument.""")

    group.add_argument(
        '--contig-length', '--length', type=int, default=100,
        help="""Remove blast hits that are shorter than this length. The
            default is "100". This is turned off by the --no-filter argument.
            """)

    blast.command_line_args(parser)
    assembly.command_line_args(parser)

    args = vars(parser.parse_args())

    check_query_args(args)
    blast.check_args(args)

    # Set defaults and adjust arguments based on other arguments
    args['cov_cutoff'] = assembly.default_cov_cutoff(args['cov_cutoff'])
    args['blast_db'] = blast.touchup_blast_db_names(args['blast_db'])
    args['kmer'] = assembly.default_kmer(args['kmer'], args['assembler'])
    args['max_target_seqs'] = blast.default_max_target_seqs(
        args['max_target_seqs'], args['blast_db'], args['max_memory'])

    # Timeout: As always, None != 0
    args['timeout'] = max(0, args['timeout'])
    if not(args['timeout']):
        args['timeout'] = None

    setup_blast_args(args)
    set_protein_arg(args)
    setup_path_arg(args)
    find_programs(args)
    util.temp_dir_exists(args['temp_dir'], args.get('debug_dir'))
    util.set_blast_batch_size(args['batch_size'])

    return args


def setup_path_arg(args):
    """Prepend to PATH environment variable if requested."""
    if args['path']:
        os.environ['PATH'] = '{}:{}'.format(args['path'], os.environ['PATH'])


def setup_blast_args(args):
    """Set up the blast args."""
    if args['no_filter']:
        args['bit_score'] = 0
        args['contig_length'] = 0


def check_query_args(args):
    """Validate the query arguments."""
    if not args.get('query') and not args.get('query_split'):
        err = 'You must have at least one --query or --query-split argument.'
        log.fatal(err)


def set_protein_arg(args):
    """Set up the protein argument."""
    if not args['protein'] and args['query']:
        args['protein'] = bio.fasta_file_has_protein(args['query'])


def find_programs(args):
    """Make sure we can find the programs needed by the assembler and blast."""
    blast.find_program('makeblastdb')
    blast.find_program('tblastn')
    blast.find_program('blastn')

    assembly.find_program(
        'abyss', 'bwa', args['assembler'], not args['no_long_reads'])

    assembly.find_program('trinity', 'Trinity', args['assembler'])
    assembly.find_program(
        'trinity', 'Trinity', args['assembler'], args['bowtie2'])

    assembly.find_program('velvet', 'velveth', args['assembler'])
    assembly.find_program('velvet', 'velvetg', args['assembler'])

    assembly.find_program('spades', 'spades.py', args['assembler'])


if __name__ == '__main__':
    ARGS = parse_command_line()
    assemble(ARGS)
