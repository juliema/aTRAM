"""Handle database functions for the stitcher utility."""

from lib.db import temp_db


def connect(temp_dir, db_prefix):
    """Create DB connection."""
    cxn = temp_db(temp_dir, db_prefix)
    cxn.row_factory = lambda c, r: dict(
            [(col[0], r[idx]) for idx, col in enumerate(c.description)])
    return cxn


# ############################## reference genes #############################

def create_reference_genes_table(cxn):
    """Create a table to hold all reference genes."""
    cxn.executescript("""
        DROP TABLE IF EXISTS reference_genes;

        CREATE TABLE reference_genes (
            ref_name     TEXT,
            ref_seq      TEXT,
            ref_file     TEXT);
        """)


def insert_reference_genes(cxn, batch):
    """Insert a batch of reference gene records into the database."""
    if batch:
        sql = """
            INSERT INTO reference_genes (ref_name, ref_seq, ref_file)
            VALUES (:ref_name, :ref_seq, :ref_file);
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def select_reference_genes(cxn):
    """Select all references."""
    return cxn.execute('SELECT * FROM reference_genes ORDER BY ref_name;')


# ################################## taxa ####################################

def create_taxa_table(cxn):
    """Create a table to hold the exonerate taxa."""
    cxn.executescript("""
        DROP TABLE IF EXISTS taxa;

        CREATE TABLE taxa (taxon_name TEXT);
        """)


def insert_taxa(cxn, batch):
    """Insert a batch of taxon records into the database."""
    sql = 'INSERT INTO taxa (taxon_name) VALUES (:taxon_name);'
    if batch:
        cxn.executemany(sql, batch)
        cxn.commit()


# ################################# contigs ##################################

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


def select_contig_files(cxn):
    """Select all contigs files for a reference gene."""
    sql = """
        SELECT DISTINCT contig_file
          FROM contigs
         ORDER BY contig_file;"""
    return cxn.execute(sql)


def select_contigs_in_file(cxn, contig_file):
    """Select all contigs for a contig file."""
    return cxn.execute(
        """SELECT * FROM contigs 
            WHERE contig_file = ?
         ORDER BY contig_rec;""",
        (contig_file, ))


def select_first_exonerate_run(cxn, ref_name):
    """Select all contig files for a reference name, taxon name combination."""
    return cxn.execute(
        """SELECT DISTINCT taxon_name, contig_file
             FROM contigs
            WHERE ref_name = ?
         ORDER BY taxon_name, contig_file
        """,
        (ref_name, ))


# ############################# exonerate results ############################

def create_exonerate_table(cxn):
    """Create a table to hold the exonerate results."""
    cxn.executescript("""
        DROP TABLE IF EXISTS exonerate;

        CREATE TABLE exonerate (
            ref_name    TEXT,
            taxon_name  TEXT,
            contig_name TEXT,
            beg         INTEGER,
            end         INTEGER,
            seq         TEXT);
        """)


def select_stitch(cxn):
    """Select all reference name, taxon name combination."""
    return cxn.execute(
        """SELECT DISTINCT ref_name, taxon_name
             FROM exonerate
         ORDER BY taxon_name, taxon_name
        """)


def insert_exonerate_results(cxn, batch):
    """Insert a batch of exonerate result records into the database."""
    if batch:
        sql = """
            INSERT INTO exonerate (
                ref_name, taxon_name, contig_name, beg, end, seq)
            VALUES (
                :ref_name, :taxon_name, :contig_name, :beg, :end, :seq);
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def select_next(cxn, ref_name, taxon_name, beg=-1):
    """
    Find the next contig for the assembly.

    It's looking for the closest contig to the given beginning. The tiebreaker
    being the longer contig.
    """
    sql = """
        SELECT *
          FROM exonerate
         WHERE ref_name   = ?
           AND taxon_name = ?
           AND beg        > ?
         ORDER BY beg, end DESC
         LIMIT 1;
        """
    result = cxn.execute(sql, (ref_name, taxon_name, beg))
    return result.fetchone()


def select_overlap(cxn, ref_name, taxon_name, beg_lo, beg_hi, end):
    """
    Find the best overlapping contig for the assembly.

    Find an overlapping contig that starts anywhere between beg_lo & beg_hi.
    Is must also end somewhere after the given end marker. We want the contig
    that extends the stitched sequence by the longest amount so we ORDER BY
    end descending & choose the first one.
    """
    sql = """
        SELECT *
          FROM exonerate
         WHERE ref_name   = ?
           AND taxon_name = ?
           AND end        > ?
           AND beg BETWEEN ? AND ?
         ORDER BY end DESC
         LIMIT 1;
        """
    result = cxn.execute(sql, (ref_name, taxon_name, end, beg_lo, beg_hi))
    return result.fetchone()


# ############################# stitched genes ###############################
def create_stitch_table(cxn):
    """Create a table to hold stitched genes & gap fillers.

    These overlaps are trimmed & the position in the assembled gene is noted.
    """
    cxn.executescript("""
        DROP TABLE IF EXISTS stitched;

        CREATE TABLE stitched (
            ref_name    TEXT,
            taxon_name  TEXT,
            contig_name TEXT,
            position    INTEGER,
            seq         TEXT);
        """)


def insert_stitched_genes(cxn, batch):
    """Insert a batch of stitched contig records into the database."""
    if batch:
        sql = """
            INSERT INTO stitched (
                ref_name, taxon_name, contig_name, position, seq)
            VALUES (
                :ref_name, :taxon_name, :contig_name, :position, :seq);
            """
        cxn.executemany(sql, batch)
        cxn.commit()


def select_stitched_contigs(cxn, ref_name, taxon_name):
    """Select stitched contigs for a reference taxon pair."""
    return cxn.execute(
        """SELECT *
             FROM stitched
            WHERE ref_name   = ?
              AND taxon_name = ?
         ORDER BY position
        """,
        (ref_name, taxon_name))
