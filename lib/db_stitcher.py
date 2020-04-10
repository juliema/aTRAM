"""Handle database functions for the stitcher utility."""

from .db import temp_db


def connect(temp_dir, db_prefix):
    """Create DB connection."""
    cxn = temp_db(temp_dir, db_prefix)
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
    sql = """
        INSERT INTO reference_genes (ref_name, ref_seq, ref_file)
        VALUES (:ref_name, :ref_seq, :ref_file);
        """
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


def select_reference_genes(cxn):
    """Select all references."""
    return cxn.execute('SELECT * FROM reference_genes ORDER BY ref_name;')


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
            contig_rec  INTEGER,
            iteration   INTEGER);
        """)


def insert_contigs(cxn, batch):
    """Insert a batch of input contig records into the database."""
    sql = """
            INSERT INTO contigs
                (ref_name, taxon_name, contig_name, contig_seq, contig_file,
                contig_rec, iteration)
            VALUES (
                :ref_name, :taxon_name, :contig_name, :contig_seq,
                :contig_file, :contig_rec, :iteration);
            """
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


def select_all_contigs(cxn):
    """Select all contigs."""
    sql = """SELECT * FROM contigs;"""
    return cxn.execute(sql)


def select_contig_files(cxn, iteration=0):
    """Select all contigs files for a reference gene."""
    sql = """
        SELECT DISTINCT contig_file
          FROM contigs
         WHERE iteration = ?
         ORDER BY contig_file;"""
    return cxn.execute(sql, (iteration, ))


def select_contigs_in_file(cxn, contig_file, iteration=0):
    """Select all contigs for a contig file."""
    return cxn.execute(
        """SELECT * FROM contigs
            WHERE contig_file = ?
              AND iteration   = ?
         ORDER BY contig_rec;""",
        (contig_file, iteration))


def select_contigs(cxn, ref_name, iteration=0):
    """Select all contig files for a reference name, taxon name combination."""
    return cxn.execute(
        """SELECT DISTINCT taxon_name, contig_file
             FROM contigs
            WHERE ref_name  = ?
              AND iteration = ?
         ORDER BY taxon_name, contig_file
        """,
        (ref_name, iteration))


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
            iteration   INTEGER,
            seq         TEXT);
        """)


def select_exonerate_ref_gene(cxn, ref_name, min_len):
    """Get all exonerate results for a reference gene."""
    return cxn.execute(
        """SELECT *
             FROM exonerate 
            WHERE ref_name = ?
              AND LENGTH(seq) >= ?
         GROUP BY seq;""",
        (ref_name, min_len))


def select_exonerate_count(cxn):
    """Select all reference name, taxon name combination."""
    result = cxn.execute("""SELECT COUNT(*) AS n FROM exonerate;""")
    return result.fetchone()['n']


def select_stitch(cxn, iteration=0):
    """Select all reference name, taxon name combination."""
    return cxn.execute(
        """SELECT DISTINCT ref_name, taxon_name
             FROM exonerate
            WHERE iteration = ?
         ORDER BY taxon_name, taxon_name
        """,
        (iteration, ))


def insert_exonerate_results(cxn, batch):
    """Insert a batch of exonerate result records into the database."""
    sql = """
            INSERT INTO exonerate (
                ref_name, taxon_name, contig_name, beg, end, iteration, seq)
            VALUES (
                :ref_name, :taxon_name, :contig_name, :beg, :end,
                :iteration, :seq);
            """
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


def select_next(cxn, ref_name, taxon_name, beg=-1, iteration=0):
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
           AND iteration  = ?
         ORDER BY beg, end DESC, contig_name
         LIMIT 1;
        """
    result = cxn.execute(sql, (ref_name, taxon_name, beg, iteration))
    return result.fetchone()


def select_longest(cxn):
    """Get the longest contig for the framer reports."""
    sql = """SELECT MAX(LENGTH(seq)) AS max_len FROM exonerate;"""
    result = cxn.execute(sql)
    return result.fetchone()['max_len']


def select_seq_lengths(cxn):
    """Get the sequence lengths for the framer reports."""
    sql = """
        SELECT taxon_name, ref_name, length(seq) AS len
          FROM exonerate;"""
    return cxn.execute(sql)


def select_overlap(
        cxn, ref_name, taxon_name, beg_lo, beg_hi, end, iteration=0):
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
           AND iteration  = ?
           AND end        > ?
           AND beg BETWEEN ? AND ?
         ORDER BY end DESC, contig_name
         LIMIT 1;
        """
    result = cxn.execute(
        sql, (ref_name, taxon_name, iteration, end, beg_lo, beg_hi))
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
            iteration   INTEGER,
            seq         TEXT);
        """)


def insert_stitched_genes(cxn, batch):
    """Insert a batch of stitched contig records into the database."""
    sql = """
        INSERT INTO stitched (
            ref_name, taxon_name, contig_name, position, iteration, seq)
        VALUES (
            :ref_name, :taxon_name, :contig_name, :position, :iteration,
            :seq);
        """
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


def select_stitched_contigs(cxn, ref_name, taxon_name, iteration=0):
    """Select stitched contigs for a reference taxon pair."""
    return cxn.execute(
        """SELECT *
             FROM stitched
            WHERE ref_name   = ?
              AND taxon_name = ?
              AND iteration  = ?
         ORDER BY position
        """,
        (ref_name, taxon_name, iteration))


def select_stitched_contig_count(cxn, ref_name, taxon_name, iteration=0):
    """Select stitched contigs for a reference taxon pair."""
    result = cxn.execute(
        """SELECT COUNT(*) AS hits
             FROM stitched
            WHERE ref_name   = ?
              AND taxon_name = ?
              AND iteration  = ?
              AND contig_name IS NOT NULL;
        """,
        (ref_name, taxon_name, iteration))
    return result.fetchone()['hits']


def select_per_gene_stats(cxn, iteration):
    """Get data for the per gene summary report."""
    return cxn.execute(
        """SELECT ref_name, taxon_name,
                LENGTH(ref_seq)      AS query_len,
                SUM(LENGTH(seq)) / 3 AS target_len
             FROM stitched
             JOIN reference_genes USING (ref_name)
            WHERE contig_name IS NOT NULL
              AND iteration = ?
         GROUP BY ref_name, taxon_name;
        """,
        (iteration, ))


def select_per_taxon_stats(cxn, iteration):
    """Get data for the per taxon summary report."""
    return cxn.execute(
        """
        WITH
        properties AS (
            SELECT taxon_name, ref_name,
                   CAST(SUM(LENGTH(seq)) / 3 AS REAL)
                   / CAST(LENGTH(ref_seq) AS REAL) AS property
             FROM stitched
             JOIN reference_genes USING (ref_name)
            WHERE contig_name IS NOT NULL
              AND iteration = ?
         GROUP BY taxon_name, ref_name),
        thresholds AS (
            SELECT taxon_name, ref_name,
                   property,
                   property  = 1.00  AS eq100,
                   property >= 0.95  AS gt95,
                   property >= 0.90  AS ge90,
                   property >= 0.80  AS ge80,
                   property >= 0.70  AS ge70,
                   property >= 0.50  AS ge50,
                   property >= 0.10  AS ge10,
                   property >  0.10  AS lt10
              FROM properties)
        SELECT taxon_name,
               COUNT(DISTINCT ref_name) AS genes,
               SUM(eq100) AS eq100,
               SUM(gt95)  AS gt95,
               SUM(ge90)  AS ge90,
               SUM(ge80)  AS ge80,
               SUM(ge70)  AS ge70,
               SUM(ge50)  AS ge50,
               SUM(ge10)  AS ge10,
               SUM(lt10)  AS lt10
          FROM thresholds
      GROUP BY taxon_name;
        """,
        (iteration, ))
