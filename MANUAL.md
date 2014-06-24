#     Manual for aTRAM

Here are the required commands for running aTRAM as well as a few things to keep in mind for maximum efficiency. Required commands are indicated and optional commands are in brackets. 


## Setup:
To determine if aTRAM can run on your computer:

### configure.pl

######perl configure.pl 
  
  This script will determine if BLAST, muscle, mafft and other de novo assembly software are available. It will create a text file called configure.txt in the lib folder with the paths to the programs. They can be available in your $PATH by adding them directly to the appropriate folder (e.g. /usr/bin) or add the path to the programs to your $PATH. Alternatively, the configure.txt file can be edited directly with the path to the programs after configure.pl is run. BLAST is required and at least one of the de novo assemblers are required Trinity or Velvet.
  
If you have problems running configure.pl for the first time, delete lib/config.txt.

URLs for software:
* BLAST: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
* Velvet: https://www.ebi.ac.uk/~zerbino/velvet/
* Trinity: http://trinityrnaseq.sourceforge.net
* MAFFT: http://mafft.cbrc.jp/alignment/software/
* muscle: http://www.drive5.com/muscle/

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
  * -exp_coverage number, the expected coverage of the region for velvetg (default 30).
  * -evalue  number, can allow for different e-values in blast search
  * -max_target_seqs, maxumum number of hits BLAST will identify
  * -assembler,  by default this is Velvet but it can be set to run Trinity -assembler = "Trinity"
  * -kmer, kmer number for assemblers


## Pipelines

#### BasicPipeline.pl
Runs aTRAM on a list of genes on a number of aTRAM formatted databases.

######perl BasicPipeline.pl -samples SampleFile.txt -targets -TargetFile.txt [atram parameters -complete -protein]
##### TargetFile.txt  tab delimited file 
genename  path/gene.fasta

##### SampleFile.txt  tab delimited file, will create a directory for each sample  
SampleName  path/atram_db


#### AlignmentPipeline.pl 
aTRAM on a list of genes on a number of aTRAM formatted databases then aligns the aTRAM contigs back to the target

######perl AlignmentPipeline.pl -samples SampleFile.txt -targets -TargetFile.txt [atram parameters -complete -protein]
##### TargetFile.txt  tab delimited file 
genename  path/gene.fasta

##### SampleFile.txt  tab delimited file, will create a directory for each sample  
SampleName  path/atram_db


