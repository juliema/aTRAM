"""Handle database functions."""

import sqlite3


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
    sql = """
        INSERT INTO aux.sra_blast_hits
                    (iteration, seq_end, seq_name, shard)
                    VALUES (?, ?, ?, ?)
        """
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


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


def get_blast_hits(cxn, iteration):
    """Get all blast hits for the iteration."""
    sql = """
        SELECT s.seq_name, s.seq_end, seq
          FROM sequences AS s
          JOIN aux.sra_blast_hits AS h 
               ON (s.seq_name = h.seq_name AND s.seq_end = h.seq_end)
           AND iteration = ?
        """

    cxn.row_factory = sqlite3.Row
    return cxn.execute(sql, (iteration, ))


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
    sql = """
        INSERT INTO aux.contig_blast_hits
                    (iteration, contig_id, description, bit_score, len,
                     query_from, query_to, query_strand,
                     hit_from, hit_to, hit_strand)
             VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


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
    sql = """
         INSERT INTO aux.assembled_contigs
                     (iteration, contig_id, seq, description,
                      bit_score, len,
                      query_from, query_to, query_strand,
                      hit_from, hit_to, hit_strand)
              VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
         """
    if batch:
        with cxn:
            cxn.executemany(sql, batch)


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
      ORDER BY bit_score DESC, len DESC, iteration
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
