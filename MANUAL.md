#     Manual for aTRAM

Here are the required commands for running aTRAM as well as a few things to keep in mind for maximum efficiency.


## Setup:
To determine if aTRAM can run on your computer:

### configure.pl

perl configure.pl 
  
  This script will tell you if you need to download any new programs. If you do make sure they are in your $PATH. You can either add them directly to your /usr/bin directory or add the path to the programs to your $PATH. 
Furthermore, configure.pl will check for the necessary de novo alignment programs including velvet, trinity and SOAPdenovo, you need to have at least one of these available. It will aslo check that you have muscle and blast, these programs are required.

## Running aTRAM

### Input File
The two paired end reads should be in one file either fasta or fasq. Concatenate the two files together if necessary.

### format_SRA.pl

perl format_SRA.pl inputfile.fastq

  This script will create aTRAM blast formatted databases of your Illumina run. 

### aTRAM.pl

perl  aTRAM.pl -reads DatabaseName -target TargetSequence.fasta -fraction Fraction 


aTRAM.pl -reads shortreadfile -target target.fasta [-ins_length int] [-exp_coverage int] [-iterations int] [-start_iteration int] [-log_file filename] [-use_ends] [-output filename]

* -reads:           short read archive (already run through format_sra.pl).
* -target:          fasta file with sequences of interest.
* -output:          optional: the prefix for the pipeline's output files (default name is the same as -reads).
* -ins_length:      optional: the size of the fragments used in the short-read library (default 300).
* -exp_coverage:    optional: the expected coverage of the region for velvetg (default 30).
* -iterations:      optional: the number of pipeline iterations (default 5).
* -start_iteration: optional: if resuming from previous run, which iteration number to start from (default 0).
* -log_file:        optional: a file to store output of the pipeline.
* -use_ends:        optional: if this flag is present, use the first and last $ins_length of long contigs in the search.



## Running on a cluster:

max_processes should equal number of threads

format_sra.pl max processes is set to 4



