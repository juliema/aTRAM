"""Handle database functions."""

import sqlite3
import sys
import os
from os.path import basename, join, exists


ATRAM_VERSION = 'v2.2.0'

# DB_VERSION != ATRAM_VERSION
# We don't force DB changes until required.
# Therefore DB_VERSION <= ATRAM_VERSION.
DB_VERSION = '2.0'

BATCH_SIZE = 1e6  # How many sequence records to insert at a time


def connect(blast_db, check_version=False, clean=False):
    """Create DB connection."""
    db_name = get_db_name(blast_db)

    if clean and exists(db_name):
        os.remove(db_name)

    if check_version and not exists(db_name):
        err = 'Could not find the database file "{}".'.format(db_name)
        sys.exit(err)

    cxn = db_setup(db_name)

    if check_version:
        check_versions(cxn)

    return cxn


def get_db_name(db_prefix):
    """Build the SQLite DB name from the prefix."""
    return '{}.sqlite.db'.format(db_prefix)


def aux_db(cxn, temp_dir, blast_db, query_name):
    """Create & attach an temporary database to the current DB connection."""
    db_dir = join(temp_dir, 'db')
    os.makedirs(db_dir, exist_ok=True)

    db_name = '{}_{}_temp.sqlite.db'.format(
        basename(blast_db), basename(query_name))
    db_name = join(db_dir, db_name)

    sql = """ATTACH DATABASE '{}' AS aux""".format(db_name)
    cxn.execute(sql)


def aux_detach(cxn):
    """Detach the temporary database."""
    cxn.execute('DETACH DATABASE aux')


def temp_db(temp_dir, db_prefix):
    """Create a temporary database."""
    db_name = join(temp_dir, get_db_name(db_prefix))
    return db_setup(db_name)


def db_setup(db_name):
    """Database setup."""
    cxn = sqlite3.connect(db_name)
    cxn.execute("PRAGMA page_size = {}".format(2**16))
    cxn.execute("PRAGMA busy_timeout = 10000")
    cxn.execute("PRAGMA journal_mode = WAL")
    return cxn


# ########################### misc functions #################################

def check_versions(cxn):
    """Make sure the database version matches what we built it with."""
    version = get_version(cxn)
    if version != DB_VERSION:
        err = ('The database was built with version {} but you are running '
               'version {}. You need to rebuild the atram database by '
               'running atram_preprocessor.py again.').format(
                   version, DB_VERSION)
        sys.exit(err)


# ########################## metadata table ##################################

def create_metadata_table(cxn):
    """
    Create the metadata table.

    A single record used to tell if we are running atram.py against the
    schema version we built with atram_preprocessor.py.
    """
    cxn.executescript("""
        DROP TABLE IF EXISTS metadata;

        CREATE TABLE metadata (
            label TEXT,
            value TEXT);
        """)

    with cxn:
        sql = """INSERT INTO metadata (label, value) VALUES (?, ?)"""
        cxn.execute(sql, ('version', DB_VERSION))
        cxn.commit()


def get_version(cxn):
    """Get the current database version."""
    sql = """SELECT value FROM metadata WHERE label = ?"""
    try:
        result = cxn.execute(sql, ('version', ))
        return result.fetchone()[0]
    except sqlite3.OperationalError:
        return '1.0'


# ########################## sequences table ##################################

def create_sequences_table(cxn):
    """Create a table to hold the raw input sequences."""
    cxn.executescript("""
        DROP TABLE IF EXISTS sequences;

        CREATE TABLE sequences (
            seq_name TEXT,
            seq_end  TEXT,
            seq      TEXT);
        """)


def create_sequences_index(cxn):
    """
    Create the sequences index after we build the table.

    This speeds up the program significantly.
    """
    sql = 'CREATE INDEX sequences_index ON sequences (seq_name, seq_end)'
    cxn.execute(sql)


def insert_sequences_batch(cxn, batch):
    """Insert a batch of sequence records into the database."""
    if batch:
        sql = """INSERT INTO sequences (seq_name, seq_end, seq)
                    VALUES (?, ?, ?)
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def get_sequence_count(cxn):
    """Get the number of sequences in the table."""
    result = cxn.execute('SELECT COUNT(*) FROM sequences')
    return result.fetchone()[0]


def get_shard_cut(cxn, offset):
    """Get the sequence name at the given offset."""
    sql = 'SELECT seq_name FROM sequences ORDER BY seq_name LIMIT 1 OFFSET {}'
    result = cxn.execute(sql.format(offset))
    cut = result.fetchone()[0]
    return cut


def get_sequences_in_shard(cxn, start, end):
    """Get all sequences in a shard."""
    sql = """
        SELECT seq_name, seq_end, seq
          FROM sequences
         WHERE seq_name >= ?
           AND seq_name < ?
        """
    return cxn.execute(sql, (start, end))


def get_sequence_ends(cxn):
    """Get a list of all seq_ends in the database."""
    return cxn.execute('SELECT DISTINCT seq_end FROM sequences')


def get_all_sequences(cxn):
    """Get a list of all seq_ends in the database."""
    return cxn.execute('SELECT * FROM sequences')


# ######################## sra_blast_hits table ###############################

def create_sra_blast_hits_table(cxn):
    """Create a table to hold the blast hits for all iterations."""
    cxn.executescript("""
        DROP TABLE IF EXISTS aux.sra_blast_hits;

        CREATE TABLE aux.sra_blast_hits (
            iteration INTEGER,
            seq_name  TEXT,
            seq_end   TEXT,
            shard     TEXT);

        CREATE INDEX aux.sra_blast_hits_index
            ON sra_blast_hits (iteration, seq_name, seq_end);
        """)


def insert_blast_hit_batch(cxn, batch):
    """Insert a batch of blast hit records into the database."""
    if batch:
        sql = """
            INSERT INTO aux.sra_blast_hits
                        (iteration, seq_end, seq_name, shard)
                        VALUES (?, ?, ?, ?)
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def sra_blast_hits_count(cxn, iteration):
    """Count the blast hits for the iteration."""
    sql = """
        SELECT COUNT(*) AS count
          FROM aux.sra_blast_hits
         WHERE iteration = ?
        """

    result = cxn.execute(sql, (iteration, ))
    return result.fetchone()[0]


def get_sra_blast_hits(cxn, iteration):
    """Get all blast hits for the iteration."""
    sql = """
        SELECT seq_name, seq_end, seq
          FROM sequences
         WHERE seq_name IN (SELECT DISTINCT seq_name
                              FROM aux.sra_blast_hits
                             WHERE iteration = ?)
      ORDER BY seq_name, seq_end
        """

    cxn.row_factory = sqlite3.Row
    return cxn.execute(sql, (iteration, ))


def get_blast_hits_by_end_count(cxn, iteration, end_count):
    """Get all blast hits for the iteration."""
    sql = """
        SELECT seq_name, seq_end, seq
          FROM sequences
         WHERE seq_name IN (SELECT seq_name
                              FROM sequences
                             WHERE seq_name IN (SELECT DISTINCT seq_name
                                                  FROM aux.sra_blast_hits
                                                 WHERE iteration = ?)
                          GROUP BY seq_name
                            HAVING COUNT(*) = ?)
      ORDER BY seq_name, seq_end
        """

    cxn.row_factory = sqlite3.Row
    return cxn.execute(sql, (iteration, end_count))


# ####################### contig_blast_hits table #############################

def create_contig_blast_hits_table(cxn):
    """Create a table to hold blast hits against the contigs."""
    cxn.executescript("""
        DROP TABLE IF EXISTS aux.contig_blast_hits;

        CREATE TABLE aux.contig_blast_hits (
            iteration    INTEGER,
            contig_id    TEXT,
            description  TEXT,
            bit_score    NUMERIC,
            len          INTEGER,
            query_from   INTEGER,
            query_to     INTEGER,
            query_strand TEXT,
            hit_from     INTEGER,
            hit_to       INTEGER,
            hit_strand   TEXT);

        CREATE INDEX aux.contig_blast_hits_index
            ON contig_blast_hits (iteration, bit_score, len);
        """)


def insert_contig_hit_batch(cxn, batch):
    """Insert a batch of blast hit records into the database."""
    if batch:
        sql = """
            INSERT INTO aux.contig_blast_hits
                        (iteration, contig_id, description, bit_score, len,
                         query_from, query_to, query_strand,
                         hit_from, hit_to, hit_strand)
                 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def get_contig_blast_hits(cxn, iteration):
    """Get all blast hits for the iteration."""
    sql = """
        SELECT iteration, contig_id, description, bit_score, len,
               query_from, query_to, query_strand,
               hit_from, hit_to, hit_strand
          FROM aux.contig_blast_hits
         WHERE iteration = ?
        """

    cxn.row_factory = sqlite3.Row
    return cxn.execute(sql, (iteration, ))


# ####################### assembled_contigs table #############################

def create_assembled_contigs_table(cxn):
    """Create a table to hold the assembled contigs."""
    cxn.executescript("""
        DROP TABLE IF EXISTS aux.assembled_contigs;

        CREATE TABLE aux.assembled_contigs (
            iteration    INTEGER,
            contig_id    TEXT,
            seq          TEXT,
            description  TEXT,
            bit_score    NUMERIC,
            len          INTEGER,
            query_from   INTEGER,
            query_to     INTEGER,
            query_strand TEXT,
            hit_from     INTEGER,
            hit_to       INTEGER,
            hit_strand   TEXT);

        CREATE INDEX aux.assembled_contigs_index
                  ON assembled_contigs (iteration, contig_id);
        """)


def assembled_contigs_count(cxn, iteration, bit_score, length):
    """Count the blast hist for the iteration."""
    sql = """
        SELECT COUNT(*) AS count
          FROM aux.assembled_contigs
         WHERE iteration = ?
           AND bit_score >= ?
           AND len >= ?
        """

    result = cxn.execute(sql, (iteration, bit_score, length))
    return result.fetchone()[0]


def iteration_overlap_count(cxn, iteration, bit_score, length):
    """Count how many assembled contigs match what's in the last iteration."""
    sql = """
        SELECT COUNT(*) AS overlap
          FROM aux.assembled_contigs AS curr_iter
          JOIN aux.assembled_contigs AS prev_iter
            ON (     curr_iter.contig_id = prev_iter.contig_id
                 AND curr_iter.iteration = prev_iter.iteration + 1)
         WHERE curr_iter.iteration = ?
           AND curr_iter.seq = prev_iter.seq
           AND curr_iter.bit_score >= ?
           AND prev_iter.bit_score >= ?
           AND curr_iter.len >= ?
        """
    result = cxn.execute(
        sql, (iteration, bit_score, bit_score, length))
    return result.fetchone()[0]


def insert_assembled_contigs_batch(cxn, batch):
    """Insert a batch of blast hit records into the database."""
    if batch:
        sql = """
            INSERT INTO aux.assembled_contigs
                        (iteration, contig_id, seq, description,
                         bit_score, len,
                         query_from, query_to, query_strand,
                         hit_from, hit_to, hit_strand)
                 VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def get_assembled_contigs(cxn, iteration, bit_score, length):
    """
    Get all assembled contigs for the iteration.

    We will use them as the queries in the next atram iteration.
    """
    sql = """
        SELECT contig_id, seq
          FROM aux.assembled_contigs
         WHERE iteration = ?
           AND bit_score >= ?
           AND len >= ?
         """
    return cxn.execute(sql, (iteration, bit_score, length))


def get_all_assembled_contigs(cxn, bit_score=0, length=0):
    """Get all assembled contigs."""
    sql = """
        SELECT iteration, contig_id, seq, description, bit_score, len,
               query_from, query_to, query_strand,
               hit_from, hit_to, hit_strand
          FROM aux.assembled_contigs
         WHERE bit_score >= ?
           AND len >= ?
      ORDER BY bit_score DESC, iteration
        """

    cxn.row_factory = sqlite3.Row
    return cxn.execute(sql, (bit_score, length))


def all_assembled_contigs_count(cxn, bit_score=0, length=0):
    """Count all assembled contigs."""
    sql = """
        SELECT COUNT(*) AS count
          FROM aux.assembled_contigs
         WHERE bit_score >= ?
           AND len >= ?
        """

    result = cxn.execute(sql, (bit_score, length))
    return result.fetchone()[0]


# ############################################################################
# ############################# stitcher tables ##############################
# ############################################################################

# ############################ stitcher contigs ##############################

def create_contigs_table(cxn):
    """Create a table to hold all of the input fasta files."""
    cxn.executescript("""
        DROP TABLE IF EXISTS contigs;

        CREATE TABLE contigs (
            ref_name    TEXT,
            taxon_name  TEXT,
            contig_name TEXT,
            contig_seq  TEXT,
            contig_file TEXT,
            contig_rec  INTEGER);
        """)


def insert_contigs(cxn, batch):
    """Insert a batch of input contig records into the database."""
    if batch:
        sql = """
            INSERT INTO contigs
                (ref_name, taxon_name, contig_name, contig_seq, contig_file,
                contig_rec)
            VALUES (
                :ref_name, :taxon_name, :contig_name, :contig_seq,
                :contig_file, :contig_rec);
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def select_conting_files(cxn, ref_name):
    """Select all contigs files for a reference gene."""
    sql = """
        SELECT DISTINCT contig_file, taxon_name
          FROM contigs
         WHERE ref_name = ?
         ORDER BY contig_file;"""
    return cxn.execute(sql, (ref_name, ))


def select_contings_in_file(cxn, ref_name, contig_file):
    """Select all contigs for a contig file."""
    sql = """
        SELECT *
          FROM contigs
         WHERE ref_name = ?
           AND contig_file = ?
         ORDER BY contig_rec;"""
    return cxn.execute(sql, (ref_name, contig_file))


def select_taxon_names(cxn, ref_name):
    """Count taxon names for a reference gene."""
    sql = """
        SELECT COUNT(DISTINCT taxon_name) FROM contigs WHERE ref_name = ?;
        """
    result = cxn.execute(sql, (ref_name, ))
    return result.fetchone()[0]


# ############################ stitcher references ###########################

def create_reference_genes_table(cxn):
    """Create a table to hold all reference genes."""
    cxn.executescript("""
        DROP TABLE IF EXISTS reference_genes;

        CREATE TABLE reference_genes (
            ref_name     TEXT,
            ref_seq      TEXT,
            ref_file     TEXT,
            results_file TEXT,
            input_file   TEXT);
        """)


def insert_reference_genes(cxn, batch):
    """Insert a batch of reference gene records into the database."""
    if batch:
        sql = """
            INSERT INTO reference_genes
                (ref_name, ref_seq, ref_file, results_file, input_file)
            VALUES (
                :ref_name, :ref_seq, :ref_file, :results_file, :input_file);"""
        cxn.executemany(sql, batch)
        cxn.commit()


def select_reference_genes(cxn):
    """Select all references."""
    cxn.row_factory = sqlite3.Row
    return cxn.execute('SELECT * FROM reference_genes ORDER BY ref_name;')


# ############################## stitcher taxa ###############################

def create_taxa_table(cxn):
    """Create a table to hold the exonerate taxa."""
    cxn.executescript("""
        DROP TABLE IF EXISTS taxa;

        CREATE TABLE taxa (taxon_name TEXT);
        """)


def insert_taxa(cxn, batch):
    """Insert a batch of taxon records into the database."""
    if batch:
        sql = """
            INSERT INTO taxa (taxon_name) VALUES (:taxon_name);"""
        cxn.executemany(sql, batch)
        cxn.commit()


def select_taxa(cxn):
    """Select all references."""
    cxn.row_factory = sqlite3.Row
    return cxn.execute('SELECT * FROM taxa ORDER BY taxon_name;')


# ######################## sticher exonerate results #########################

def create_exonerate_table(cxn):
    """Create a table to hold the exonerate results."""
    cxn.executescript("""
        DROP TABLE IF EXISTS exonerate;

        CREATE TABLE exonerate (
            gene_name       TEXT,
            taxon_name      TEXT,
            query_len       INTEGER,
            query_align_len INTEGER,
            query_align_beg INTEGER,
            query_align_end INTEGER,
            target_id       TEXT,
            target_seq      TEXT);

        CREATE INDEX exonerate_sort_idx
            ON exonerate (gene_name, taxon_name, query_align_beg);
        """)


def insert_exonerate_results(cxn, batch):
    """Insert a batch of exonerate result records into the database."""
    if batch:
        sql = """
            INSERT INTO exonerate
                (gene_name, taxon_name, query_len, query_align_len,
                 query_align_beg, query_align_end, target_id, target_seq)
            VALUES (
                :gene_name, :taxon_name, :query_len, :query_align_len,
                :query_align_beg, :query_align_end, :target_id, :target_seq);
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def select_exonerate_results(cxn):
    """Select exonerate results ordered by where contig starts."""
    sql = """
        SELECT *
          FROM exonerate
      ORDER BY gene_name, taxon_name, query_align_beg;
        """
    return cxn.execute(sql)
