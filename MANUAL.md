#     Manual for aTRAM

Here are the required commands for running aTRAM as well as a few things to keep in mind for maximum efficiency. Required commands are indicated and optional commands are in brackets. 


## Setup:
To determine if aTRAM can run on your computer:

### configure.pl

######perl configure.pl 
  
  This script will determine if BLAST, muscle, mafft and other de novo assembly software are available. It will create a text file called configure.txt in the lib folder with the paths to the programs. They can be available in your $PATH by adding them directly to the appropriate folder (e.g. /usr/bin) or add the path to the programs to your $PATH. Alternatively, the configure.txt file can be edited directly with the path to the programs after configure.pl is run. BLAST is required and at least one of the de novo assemblers are required Trinity or Velvet.

## Running aTRAM

### format_SRA.pl

######perl format_SRA.pl -input inputfile.fastq [-output Database Name | -number Number of Shards ]

This script will create an aTRAM blast formatted database of your Illumina run. By default it will split the short read archive into shards for every 250 MB. The final aTRAM database will have roughly (Number of GB of a fasta file)/4 blast formatted databases. You can specify a database name with -output and the number of shards with -number.  We calculated 250 MB based on the memory allocated to BLAST and time necessary to search through each database. 
  
The -input file can be either fastq or fasta, the two paired end reads should be in one file. Concatenate the two files together if necessary.


### aTRAM.pl

######perl  aTRAM.pl -sra DatabaseName -target TargetSequence.fasta  [-fraction Fraction | -complete | -protein | -max_processes | -max_memory | -shards number | -log_file filename | -help | -output_file filename | -temp_name | -save_temp ]

Required Parameters:
  * - sra DatabaseName, the name of the formatted short read archive from format_sra.pl
  * - target  TargetSequence.fasta, the target gene in either DNA or Amino Acid format, if AA turn on -protein

Important Optional Parameters:
  * - complete  [on/off], will autocomplete if the atram output blasts back to both the beginning and end of the target sequence, indicates the full target has been recovered.
  * - fraction number bewteen 0-1,  fraction of the full Illumina database to search, a good starting point would be to aim for 20X coverage
  * - protein  [on/off], turns on if the target sequence is in Amino Acids
  
Other Optional Parameters:
  * - help will give program parameters
  * - log_file  filename,    will create one file for all the log output
  * - output_file filename,  will create the output file
  * - shards number,  number of shards to search
  * - temp_name filename, names a temporary file to store XXXX
  * - save_temp filename, saves the temp file

Running on a Cluster:
  * - max_processes  maximum  number of processes to run for each search
  * - max_memory maximum amount of memory to allocate to aTRAM


Parameters with modifiable default values:    
  * -ins_length number
  * -iterations number, number of iterations before the program stops, default 5
  * -start_iter number, iteration to start on if continuing running aTRAM from a previous run
  * -exp_cov
  * -evalue  number, can allow for different e-values in blast search
  * 

my $ins_length = 300;
my $iterations = 5;
my $start_iter = 1;
my $exp_cov = 30;
my $evalue = 10e-10;
my $max_target_seqs = 100000000;
my $assembler = "Velvet";
my $kmer = 31;



[-ins_length int] [-exp_coverage int] [-iterations int] [-start_iteration int] [-log_file filename] [-use_ends] [-output filename]

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



