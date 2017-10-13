## Manual for using aTRAM 2.0: automated Target Restricted Assembly Method

## Background
    1. What does it do?
    2. How it works

aTRAM 2.0 is a major overhaul of the aTRAM approach to assembling loci from NGS data. The new code has been reimplemented in python, and the approach to short read library construction is completely revamped, resulting in major performance and assembly improvements.

aTRAM ("automated target restricted assembly method") is an iterative assembler that performs reference-guided local de novo assemblies using a variety of available methods. It is well-suited to various tasks where NGS data needs to be queried for gene sequences, such as phylogenomics. The design philosophy is modular and expandable, with support for four de-novo assemblers to date: Velvet, Abyss, Trinity, and Spades.

Cite like so:
```Allen, JM, DI Huang, QC Cronk, KP Johnson. 2015. aTRAM automated target restricted assembly method a fast method for assembling loci across divergent taxa from next-generation sequencing data. BMC Bioinformatics 16:98 DOI 10.1186/s12859-015-0515-2```

## Installation
     1. Python 3.0 or greater and a number of dependencies
     2. Dependencies given in requirements.txt

You will need to have Python3 installed, as well as pip, a package manager for python. Beyond these, it is easiest to handle aTRAM 2 dependencies by setting up a virtual environment, which is a contained workspace with internally installed python libraries. Run the following code in what you intend to be your working directory:

```
virtualenv venv -p python3
source venv/bin/activate
pip install -r path_to_atram/requirements.txt 
```

You should see something like `(venv)` at the beginning of your command prompt after running the second line, indicating the environment is active. Once you have verified that the requirements installed with no errors, only the second line needs to be run before each aTRAM 2 session.

If you choose not to use virtual environments, you will likely have to specify python3.

URLs for software:

* BLAST: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
* Velvet: https://www.ebi.ac.uk/~zerbino/velvet/
* Trinity: http://trinityrnaseq.sourceforge.net
* MAFFT: http://mafft.cbrc.jp/alignment/software/
* MUSCLE: http://www.drive5.com/muscle/
* Spades: http://bioinf.spbau.ru/spades

### Library Preparation
Use `atram_preprocessor.py` for this. You can either list the forward and reverse read files, or glob them with wildcards as below. Under the hood, aTRAM 2 is building a SQLite3 database for rapid read retrieval. Note that aTRAM 2 is not backwards compatible with aTRAM 1 libraries; it is also best to rebuild any libraries after major updates. 

``` python path_to_aTRAM/atram_preprocessor.py -c NUMBER_OF_THREADS -b path_to_atram_library/LIBRARY_PREFIX READ_NAME*.fastq ```
  
### Assembling Loci
``` 
python path_to_aTRAM/atram.py -b path_to_atram_library/LIBRARY_PREFIX -q path_to_reference_loci/Locus.fasta -i NUMBER_OF_ITERATIONS --cpus NUMBER_OF_THREADS  --kmer KMER_NUMBER -o path_to_output/LIBRARY_PREFIX.Locus.atram2.fasta --log-file path_to_output/LIBRARY_PREFIX.Locus.log -a ASSEMBLER_CHOICE
```

Fill in the capitalized portions with your options, and fill in the paths. 

There are many more options than this, so for reference list them like so:

```python path_to_atram/atram.py -h```

## Arguments
Several arguments have synonyms, given below.

* `-h, --help`
    * List arguments
    
* `--version`        
    * Give version number

# Required arguments:
* `-b DB, --blast-db DB, --sra DB, --db DB, --database DB`
    * The aTRAM library. match the name you gave for the library prefix in `atram_preprocessor.py`.

* `-o OUTPUT, --output OUTPUT`
    * Give a prefix for output files. You may include a directory path as part of the prefix.
  
* `-q QUERY, --query QUERY, --target QUERY`
    * Path to the fasta file with sequences of interest. Not required if you specify `--start-iteration` to restart a run.

# Optional aTRAM arguments:
* `a {abyss,trinity,velvet,spades}, --assembler {abyss,trinity,velvet,spades}`
    * Choose which assembler to use from the list. If you do not use this argument then aTRAM will do a single blast run and stop before assembly.

* `-i N, --iterations N`  
    * The number of pipeline iterations. The default is "5".

* `-p, --protein`
    * Are the query sequences protein? aTRAM will guess if you skip this argument.

* `--fraction FRACTION`
    * Use only the specified fraction of the aTRAM database. The default is to use all data (=1.0). This option is useful for very large datasets and for high-copy targets such as mitochondria.

* `--cpus CPUS, --processes CPUS, --max-processes CPUS`
    * Number of CPU threads to use. This will also be used for the assemblers when possible. Defaults to: Total system CPUs - 4

* `--log-file LOG_FILE`   
    * Specifies the full path of the log file (full path). The default is to use the DIR and DB arguments to come up with a name like so: `DIR/DB_atram.log`

* `--path PATH`
    * If the assembler or BLAST dependencies you want to use are not in your $PATH then use this to prepend directories to your path.

* `--start-iteration N, --restart N`
    * If resuming from a previous run, which iteration number to start from. The default is 1.

* `-t DIR, --temp-dir DIR`
    * You may save intermediate files for debugging in this directory. The directory must be empty. This option is useful for understanding exactly what aTRAM is coming up with at each step -- BLAST results, etc. These files are always made but by default not kept and handled as OS temporary files.

* `-T SECONDS, --timeout SECONDS`
    * How many seconds to wait for an assembler before stopping the run. To wait forever set this to 0. The default is "300" (5 minutes). This option was added to account for assembler module errors we observed in Abyss when long reads are used; under certain conditions Abyss can rarely run indefinitely and must be killed.

# Optional values for blast-filtering contigs:
* `--bit-score SCORE`
    * Remove contigs that have a value less than this. The default is 70.0. This is turned off by the --no-filter argument. Increasing this arguement is useful if you are getting non-target contigs; if small or divergent targets are missed, reducing it can also be tried.
    
* `--contig-length CONTIG_LENGTH, --length CONTIG_LENGTH`
    * Remove blast hits that are shorter than this length. The default is 100. This is turned off by the --no-filter argument.
    
* `--no-filter`
    * Do not filter the assembled contigs. This will set both the --bit-score and --contig-length to 0

# Optional blast arguments:
* `--db-gencode CODE`     
    * The genetic code to use during blast runs. The default is "1", the standard eukaryotic code. For chloroplast, mitochondrial, and bacterial protein targets, refer to [NCBI codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1) for the correct option.
    
* `--evalue EVALUE`
    * The default evalue for BLAST is 1e-10.

* `--max-target-seqs MAX`
    * Maximum hit sequences per shard. The default is calculated automatically based on the available memory and the number of shards. this could be useful for very high-copy targets.

optional assembler arguments:
* `--no-long-reads`
    * Do not use long reads during assembly (Abyss, Trinity, Velvet). This controls behavior of a new option, where previously recovered contigs are used in assemblies as "long read data." If jobs do not finish due to assembly module problems, refer above to time limits on aTRAM runs.

* `--kmer KMER`
    * k-mer size. The default is 64 for Abyss and 31 for Velvet. Note: the maximum kmer length for Velvet is 31 (Abyss, Velvet).

* `--mpi`
    * Use MPI for this assembler. The assembler must have been compiled to use MPI (Abyss) and mpirun must be available in the path.
    
* `--bowtie2`
    * Use bowtie2 during assembly (Trinity).
    
* `--max-memory MEMORY`
    * Maximum amount of memory to use in gigabytes. The default is 30 (Trinity, Spades).
    
* `--exp-coverage EXP_COVERAGE, --expected-coverage EXP_COVERAGE`
    * The expected coverage of the region. The default is 30 (Velvet). Refer to the Velvet manual for this option.

* `--ins-length INS_LENGTH`
    * The mean size of the fragments used in the short-read library. The default is 300 (Velvet). This can be calculated from short read data given a reference, or refer to your library construction solution.
    
* `--min-contig-length MIN_CONTIG_LENGTH`
    * The minimum contig length used by the assembler itself. The default is 100 (Velvet).
    
* `--cov-cutoff COV_CUTOFF`
    * Read coverage cutoff value (Spades). Must be a positive float value, or "auto", or "off". The default value is "off".


## Example of running a shell loop
In many cases it is convenient to run aTRAM 2 as a loop, assembling a set of genes for a set of taxa. These can be set up in two parts like so:

```
# Make aTRAM libraries
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do 
python path_to_aTRAM/atram_preprocessor.py -c 4 -b path_to_atram_library/lib_${a} path_to_input/${a}_P*.fq
done
```

The part `${a}_P*.fq` will have to be modified to match the name pattern of your input fastq files. Then, supposing we have 300 genes labeled consecutively and wish to use Abyss:

```
# Assemble genes
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do
for (( i=1 ; i<=300; i++ )); # Iterate through locus numbers
do 
python ./aTRAM/atram.py -b path_to_atram_library/lib_${a} -q path_to_reference_loci/Locus_${i}.fasta -i 5 --cpus 4  --kmer 64 -o path_to_output/lib_${a}.Locus_${i}.atram2.fasta --log-file path_to_output/lib_${a}.Locus_${i}.log -a abyss
done
done
```

## Backwards compatibility

For any tools that depend on the output format of aTRAM 1.0, this script will perform the conversion of fasta headers:

```
for i in $(find . -name "*.fasta"); do
sed 's/.* iteration=/>/g' ${i} | sed 's/ contig_id/.0_contigid/g' | sed 's/contigid.*length_//g' | sed 's/_cov.* score=/_/g' | sed 's/\.[0-9]*$//g' > ${i}.aTRAM1.fasta
done
```

For the [exon stitching pipeline](https://github.com/juliema/exon_stitching), output files from aTRAM 2.0 must additionally be named like so: `libraryname_locusname.best.fasta`. If both this and the fasta header conversion are performed, any previously used tools should work.
