# aTRAM

This  takes a query sequence and a blast database built with the 
atram_preprocessor.py script and builds assemblies.

If you specify more than one query sequence and/or more than one blast
database then aTRAM will build one assembly for each query/blast
DB pair.

NOTE: You may use a text file to hold the command-line arguments
like: @/path/to/args.txt. This is particularly useful when specifying
multiple blast databases or multiple query sequences.

 ## Arguments
 
`-h, --help`

Show this help message and exit.

`--version`

Show program's version number and exit.

`-b DB [DB ...], --blast-db DB [DB ...], --sra DB [DB ...], --db DB [DB ...],
--database DB [DB ...]`
                    
This needs to match the DB prefix you entered for
atram_preprocessor.py. You may repeat this argument to
run the --query sequence(s) against multiple blast
databases.
                    
`-q QUERY [QUERY ...], --query QUERY [QUERY ...], --target QUERY [QUERY ...],
--probe QUERY [QUERY ...]`

The path to the fasta file with sequences of interest.
You may repeat this argument. If you do then Each
--query sequence file will be run against every
--blast-db.

`-Q QUERY_SPLIT [QUERY_SPLIT ...], --query-split QUERY_SPLIT [QUERY_SPLIT ...],
--target-split QUERY_SPLIT [QUERY_SPLIT ...]
`
The path to the fasta file with multiple sequences of
interest. This will take every sequence in the fasta
file and treat it as if it were its own --query
argument. So every sequence in --query-split will be
run against every --blast-db.

`-o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX`

This is the prefix of all of the output files. So you
can identify different blast output file sets. You may
include a directory as part of the prefix. aTRAM will
add suffixes to differentiate output files.

`-a {abyss,trinity,velvet,spades,none}, --assembler
{abyss,trinity,velvet,spades,none}`

Which assembler to use. Choosing "none" (the default)
will do a single blast run and stop before any
assembly.

`-i N, --iterations N`

The number of pipeline iterations. The default is "5".

`-p, --protein`

Are the query sequences protein? aTRAM will guess if you skip this argument.

`--fraction FRACTION`

Use only the specified fraction of the aTRAM database. The default is 1.0.

`--cpus CPUS, --processes CPUS, --max-processes CPUS`

Number of CPU processors to use. This will also be
used for the assemblers when possible. We will use 8
out of 12 CPUs.

`--log-file LOG_FILE`

Log file (full path)".

`--path PATH`

If the assembler or blast you want to use is not in
your $PATH then use this to prepend directories to
your path.

`-t DIR, --temp-dir DIR`

Place temporary files in this directory. All files
will be deleted after aTRAM completes. The directory
must exist.
                    
`--keep-temp-dir`

This flag will keep the temporary files in the --temp-dir around for debugging.

`-T SECONDS, --timeout SECONDS`

How many seconds to wait for an assembler before
stopping the run. To wait forever set this to 0. The
default is "300" (5 minutes).

`--no-filter`

Do not filter the assembled contigs. This will: set
both the --bit-score and --contig-length to 0
                    
`--bit-score SCORE`

Remove contigs that have a value less than this. The
default is "70.0". This is turned off by the --no-
filter argument.

`--contig-length CONTIG_LENGTH, --length CONTIG_LENGTH`
                    
Remove blast hits that are shorter than this length.
The default is "100". This is turned off by the --no-
filter argument.

`--db-gencode CODE
`
The genetic code to use during blast runs. The default is "1".

`--evalue EVALUE`

The default evalue is "1e-10".

`--word-size WORD_SIZE
`
Word size for wordfinder algorithm. 'Must be >= 2.

`--max-target-seqs MAX`

Maximum hit sequences per shard. Default is calculated
based on the available memory and the number of
shards.
                    
`--batch-size BATCH_SIZE`
                    
Use this option to control blast memory usage and the
concatenation of queries. Setting this value too low
can degrade performance.

`--no-long-reads`

Do not use long reads during assembly. (Abyss, Trinity, Velvet)
                    
`--kmer KMER`

k-mer size. The default is 64 for Abyss and 31 for
Velvet. Note: the maximum kmer length for Velvet is
31. (Abyss, Velvet)
                    
`--mpi`

Use MPI for this assembler. The assembler 'must have
been compiled to use MPI. (Abyss)

`--bowtie2
`
Use bowtie2 during assembly. (Trinity)

`--max-memory MEMORY`

Maximum amount of memory to use in gigabytes. We will
use 12 out of 24 GB of free/unused memory. (Trinity,
Spades)

`--exp-coverage EXP_COVERAGE, --expected-coverage EXP_COVERAGE`
 
The expected coverage of the region. The default is
"30". (Velvet)
                    
`--ins-length INS_LENGTH`
                    
The size of the fragments used in the short-read
library. The default is "300". (Velvet)
                    
`--min-contig-length MIN_CONTIG_LENGTH`

The minimum contig length used by the assembler
itself. The default is "100". (Velvet)
                    
`--cov-cutoff COV_CUTOFF`
                    
Read coverage cutoff value. Must be a positive float
value, or "auto", or "off". The default is "off".
(Spades)
