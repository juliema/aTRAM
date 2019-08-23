

Fill in the capitalized portions with your options, and fill in the paths.

There are many more options than this, so for reference list them like so:

atram.py -h`

# aTRAM main script arguments

Several arguments have synonyms, given below.

## General:

- `-h, --help`

  - List arguments

- `--version`

  - Give version number

## Required arguments:

- `-b DB, --blast-db DB, --sra DB, --db DB, --database DB`

  - The aTRAM library. match the name you gave for the library prefix in `atram_preprocessor.py`.

- `-o OUTPUT, --output OUTPUT`

  - Give a prefix for output files. You may include a directory path as part of the prefix.

- `-q QUERY, --query QUERY, --target QUERY` or
  `-Q QUERY_SPLIT, --query QUERY_SPLIT, --target QUERY_SPLIT`

  - This specifies the query sequence or sequences.  If one sequence use the "-q" option.  For many query sequences, "-Q".

## Optional aTRAM arguments:

- `a {abyss,trinity,velvet,spades}, --assembler {abyss,trinity,velvet,spades}`

  - Choose which assembler to use from the list. If you do not use this argument then aTRAM will do a single blast run and stop before assembly.

- `-i N, --iterations N`

  - The number of pipeline iterations. The default is "5".

- `-p, --protein`

  - Are the query sequences protein? aTRAM will guess if you skip this argument.

- `--fraction FRACTION`

  - Use only the specified fraction of the aTRAM database. The default is to use all data (=1.0). This option is useful for very large datasets and for high-copy targets such as mitochondria.

- `--cpus CPUS, --processes CPUS, --max-processes CPUS`

  - We default to a number that will not use up all of your cores. You may set this number to use more preocessors, but be aware that aTRAM uses a lot of temporary disk space (usually in /tmp), so you should balance the increased parallelization with the increase in temporary disk space. We should also note that you can use the `--temp-dir` option (below) to use a higher capacity disk.

- `--log-file LOG_FILE`

  - Specifies the full path of the log file (full path). The default is to use the DIR and DB arguments to come up with a name like so: `DIR/DB_atram.log`

- `-T SECONDS, --timeout SECONDS`

  - How many seconds to wait for an assembler before stopping the run. To wait forever set this to 0\. The default is "300" (5 minutes). This option was added to account for assembler module errors we observed in Abyss when long reads are used; under certain conditions Abyss can rarely run indefinitely and must be killed.

## Optional values for blast-filtering contigs:

- `--bit-score SCORE`

  - Remove contigs that have a value less than this. The default is 70.0\. This is turned off by the --no-filter argument. Increasing this arguement is useful if you are getting non-target contigs; if small or divergent targets are missed, reducing it can also be tried.

- `--contig-length CONTIG_LENGTH, --length CONTIG_LENGTH`

  - Remove blast hits that are shorter than this length. The default is 100\. This is turned off by the --no-filter argument.

- `--no-filter`

  - Do not filter the assembled contigs. This will set both the --bit-score and --contig-length to 0

## Optional blast arguments:

- `--db-gencode CODE`

  - The genetic code to use during blast runs. The default is "1", the standard eukaryotic code. For chloroplast, mitochondrial, and bacterial protein targets, refer to [NCBI codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1) for the correct option.

- `--evalue EVALUE`

  - The default evalue for BLAST is 1e-10.

- `--max-target-seqs MAX`

  - Maximum hit sequences per shard. The default is calculated automatically based on the available memory and the number of shards. this could be useful for very high-copy targets.

optional assembler arguments:

- `--no-long-reads`

  - Do not use long reads during assembly (Abyss, Trinity, Velvet). This controls behavior of a new option, where previously recovered contigs are used in assemblies as "long read data." If jobs do not finish due to assembly module problems, refer above to time limits on aTRAM runs.

- `--kmer KMER`

  - k-mer size. The default is 64 for Abyss and 31 for Velvet. Note: the maximum kmer length for Velvet is 31 (Abyss, Velvet).

- `--mpi`

  - Use MPI for this assembler. The assembler must have been compiled to use MPI (Abyss) and mpirun must be available in the path.

- `--bowtie2`

  - Use bowtie2 during assembly (Trinity).

- `--max-memory MEMORY`

  - Maximum amount of memory to use in gigabytes. The default is 30 (Trinity, Spades).

- `--exp-coverage EXP_COVERAGE, --expected-coverage EXP_COVERAGE`

  - The expected coverage of the region. The default is 30 (Velvet). Refer to the Velvet manual for this option.

- `--ins-length INS_LENGTH`

  - The mean size of the fragments used in the short-read library. The default is 300 (Velvet). This can be calculated from short read data given a reference, or refer to your library construction solution.

- `--min-contig-length MIN_CONTIG_LENGTH`

  - The minimum contig length used by the assembler itself. The default is 100 (Velvet).

- `--cov-cutoff COV_CUTOFF`

  - Read coverage cutoff value (Spades). Must be a positive float value, or "auto", or "off". The default value is "off".

# Assembling multiple genes against a library

 aTRAM2.0 can assemble a set of genes against a single library.  Create a single file with multiple fasta-formatted sequences and then simply use `-Q QUERY_SPLIT` where QUERY_SPLIT is the name of the file you created above.


# Example of running a shell loop

In many cases it is convenient to run aTRAM 2 as a loop, assembling a set of genes for a set of taxa. These can be set up in two parts, as shown below.  Note that aTRAM2 has built in functions supporting assembly of many genes against a library, as described just above.

```
# Make aTRAM libraries
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do
  atram_preprocessor.py -c 4 -b path_to_atram_library/lib_${a} path_to_input/${a}_P*.fq
done
```

The part `${a}_P*.fq` will have to be modified to match the name pattern of your input fastq files.

Then, supposing we have a set of genes stored in a single file and wish to use Abyss:

```
# Assemble genes
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do
  atram.py -b path_to_atram_library/lib_${a} -Q file_name -i 5 --cpus 4  --kmer 64 -o path_to_output/lib_${a}.atram2.fasta --log-file path_to_output/lib_${a}.log -a abyss
done
```
