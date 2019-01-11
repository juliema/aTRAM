# Manual for using aTRAM 2.0: automated Target Restricted Assembly Method

# Background

aTRAM 2.0 is a major overhaul of the aTRAM approach to assembling loci from NGS data. The new code has been reimplemented in python, and the approach to short read library construction is completely revamped, resulting in major performance and assembly improvements.

aTRAM ("automated target restricted assembly method") is an iterative assembler that performs reference-guided local de novo assemblies using a variety of available methods. It is well-suited to various tasks where NGS data needs to be queried for gene sequences, such as phylogenomics. The design philosophy is modular and expandable, with support for four de-novo assemblers to date: Velvet, Abyss, Trinity, and Spades.

Please consult the reference below for more information about aTRAM1.0:
`Allen, JM, DI Huang, QC Cronk, KP Johnson. 2015\. aTRAM automated target restricted assembly method a fast method for assembling loci across divergent taxa from next-generation sequencing data. BMC Bioinformatics 16:98 DOI 10.1186/s12859-015-0515-2`

A paper on aTRAM 2.0 is now in press:
`Allen J.M., R. LaFrance, R. A. Folk, K. P. Johnson, and R. P. Guralnick.  In Press.  aTRAM 2.0: An improved, flexible locus assembler for NGS data.  Evolutionary Informatics`

# Installation

You will need to have Python3 installed, as well as pip, a package manager for python. Beyond these, it is easiest to handle aTRAM 2 dependencies by setting up a virtual environment, which is a contained workspace with internally installed python libraries. Run the following code in what you intend to be your working directory:

```
git clone https://github.com/juliema/aTRAM.git
cd path/to/cloned/atram
virtualenv venv -p python3
source venv/bin/activate
pip install -r requirements.txt
```

You should see something like `(venv)` at the beginning of your command prompt after running the second line, indicating the environment is active. Once you have verified that the requirements installed with no errors, only the second line needs to be run before each aTRAM 2 session.

If you choose not to use virtual environments, you will likely have to specify python3.

You will need to install BLAST externally and have it in the path. You will also need one of the supported assembly modules. URLs are given below:

- [BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - version 2.7.1
- [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/)
- [Trinity](http://trinityrnaseq.github.io/) - version 2.5.1
- [Abyss](http://www.bcgsc.ca/platform/bioinfo/software/abyss) - version 2.0.2
- [SPAdes](http://cab.spbu.ru/software/spades/) - version 3.11.1

# Library Preparation

Use `atram_preprocessor.py` for this. Define your new library name with -b (something simple). Then give it your fastq files. You can either list the forward and reverse read files, or put them into one file and use the --mixed-ends option. Under the hood, aTRAM 2 is building a SQLite3 database for rapid read retrieval. Note that aTRAM 2 is not backwards compatible with aTRAM 1 libraries; it is also best to rebuild any libraries after major updates.

```
python path_to_aTRAM/atram_preprocessor.py -c NUMBER_OF_THREADS -b path_to_atram_library/LIBRARY_PREFIX --end-1 path_to_read_1/read_1.fastq --end-2 path_to_read_2/read_2.fastq
```

## Preprocessor arguments:

- `-h, --help`

  - List arguments

- `--version`

  - Give version number

- `-b DB, --blast-db DB, --output DB, --db DB`

  - This is the prefix of all of the blast database files so you can identify different blast database sets and so they can be stored together without resorting to subdirectories. You may include a directory as part of the prefix. The default is `atram_<today's date>`.

- `--cpus CPUS, --processes CPUS, --max-processes CPUS`

  - We default to a number that will not use up all of your cores. You may set this number to use more preocessors, but be aware that aTRAM uses a lot of temporary disk space (usually in /tmp), so you should balance the increased parallelization with the increase in temporary disk space. We should also note that you can use the `--temp-dir` option (below) to use a higher capacity disk.

- `-t DIR, --temp-dir DIR`

  - Place temporary files in this directory. All files will be deleted after aTRAM completes. The directory must exist.

- `-l LOG_FILE, --log-file LOG_FILE`

  - Log file (full path). The default is to use the DB and program name to come up with a name like so: `<DB>_atram_preprocessor.log`

- `-s SHARDS, --shards SHARDS, --number SHARDS`

  - Number of blast DB shards to create. The default is to have each shard contain roughly 250MB of sequence data. Since this parameter affects parallelization and I/O load, it may be worth experimenting with.

# Assembling Loci

```
python path_to_aTRAM/atram.py -b path_to_atram_library/LIBRARY_PREFIX -q path_to_reference_loci/Locus.fasta -i NUMBER_OF_ITERATIONS --cpus NUMBER_OF_THREADS  --kmer KMER_NUMBER -o path_to_output/LIBRARY_PREFIX.Locus.atram2.fasta --log-file path_to_output/LIBRARY_PREFIX.Locus.log -a ASSEMBLER_CHOICE
```

Fill in the capitalized portions with your options, and fill in the paths.

There are many more options than this, so for reference list them like so:

`python path_to_atram/atram.py -h`

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

In many cases it is convenient to run aTRAM 2.0 as a loop, assembling a set of genes for a set of taxa. These can be set up in two parts, as shown below.  Note that aTRAM2.0 has built in functions supporting assembly of many genes against a library, as described just above.  

```
# Make aTRAM libraries
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do
python path_to_aTRAM/atram_preprocessor.py -c 4 -b path_to_atram_library/lib_${a} path_to_input/${a}_P*.fq
done
```

The part `${a}_P*.fq` will have to be modified to match the name pattern of your input fastq files.

Then, supposing we have a set of genes stored in a single file and wish to use Abyss:

```
# Assemble genes
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do
\python ./aTRAM/atram.py -b path_to_atram_library/lib_${a} -Q file_name -i 5 --cpus 4  --kmer 64 -o path_to_output/lib_${a}.atram2.fasta --log-file path_to_output/lib_${a}.log -a abyss
done
```

# Backwards compatibility

For any tools that depend on the output format of aTRAM 1.0, this script will perform the conversion of fasta headers:

```
for i in $(find . -name "*.fasta"); do
sed 's/.* iteration=/>/g' ${i} | sed 's/ contig_id/.0_contigid/g' | sed 's/contigid.*length_//g' | sed 's/_cov.* score=/_/g' | sed 's/\.[0-9]*$//g' > ${i}.aTRAM1.fasta
done
```

For the [exon stitching pipeline](https://github.com/juliema/exon_stitching), output files from aTRAM 2.0 must additionally be named like so: `libraryname_locusname.best.fasta`. If both this and the fasta header conversion are performed, any previously used tools should work.

# Testing suite

For a new aTRAM install, it may be desirable to make sure the install resulted in a fully functional aTRAM. End users can activate our testing suite by running this command in the aTRAM directory:

```pytest tests```
