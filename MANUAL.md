#Using aTRAM

Here are the required commands for running aTRAM as well as a few things to keep in mind for maximum efficiency. Required commands are indicated and optional commands are in brackets. 


## Setup
To determine if aTRAM can run on your computer:

### configure.pl

```perl configure.pl```
 
This script will determine if BLAST, muscle, mafft and other de novo assembly software are available. It will create a text file called configure.txt in the lib folder with the paths to the programs. They can be available in your $PATH by adding them directly to the appropriate folder (e.g. /usr/bin) or add the path to the programs to your $PATH. Alternatively, the configure.txt file can be edited directly with the path to the programs after configure.pl is run. BLAST is required and at least one of the de novo assemblers,  Trinity or Velvet, are required.
 
If you have problems running configure.pl for the first time, delete lib/config.txt.

URLs for software:

* BLAST: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
* Velvet: https://www.ebi.ac.uk/~zerbino/velvet/
* Trinity: http://trinityrnaseq.sourceforge.net
* MAFFT: http://mafft.cbrc.jp/alignment/software/
* MUSCLE: http://www.drive5.com/muscle/

## Running aTRAM

### format_SRA.pl

```perl format_SRA.pl -input inputfile.fastq [-output Database Name | -number Number of Shards ]```

This script will create an aTRAM-formatted database of your Illumina run. By default, it will split the short read archive into ~250 MB shards. The final aTRAM database will have roughly (number of GB of a fasta file)/4 shards. We chose 250 MB based on the memory allocated to BLAST and time necessary to search through each database. You can specify the number of shards created with `-number`.

`-output` specifies the name of the aTRAM database. 
 
`-input` specifies the fastq or fasta file of short reads. Concatenate the two files together if necessary.

### aTRAM.pl

```perl aTRAM.pl -sra DatabaseName -target TargetSequence.fasta [-fraction Fraction | -complete | -protein | -max_processes | -max_memory | -shards number | -log_file filename | -help | -output_file filename | -temp_name | -save_temp ]```

Required Parameters:
 
 * `-sra DatabaseName`: the name of the formatted short read archive from format_sra.pl
 * `-target TargetSequence.fasta`: the target gene in either DNA or Amino Acid format. (if AA, use the `-protein` flag)

Important Optional Parameters:

 * `-complete` [on/off]: if the atram output blasts back to both the beginning and end of the target sequence, indicates the full target has been recovered.
 * `-fraction [number between 0-1]`: fraction of the full Illumina database to search, a good starting point would be to aim for 20X coverage
 * -protein [on/off], turns on if the target sequence is in Amino Acids
 
Other Optional Parameters:

 * `-help` will give program parameters
 * `-log_file filename` will create one file for all the log output
 * `-output_file filename` specifies the output file prefix
 * `-shards number` specifies the number of shards to search
 * `-temp_name filename` names a temporary file to store XXXX
 * `-save_temp filename` saves the temp file

Running on a Cluster:

 * `-max_processes number`: maximum number of processes to run for each search
 * `-max_memory number`: maximum amount of memory to allocate to aTRAM

Parameters with modifiable default values: 

 * `-ins_length number`: the average insert length of the Illumina library
 * `-iterations number`: the number of iterations before the program stops ( default 5)
 * `-start_iter number`: the iteration to start on if restarting aTRAM from a previous run
 * `-exp_coverage number`: expected coverage of the region for velvetg (default 30).
 * `-evalue number`: e-value for blast search
 * `-max_target_seqs`: maximum number of hits BLAST will identify
 * `-assembler assemblername`: specifies one of the configured de novo assemblers in lib/Assemblers (default Velvet)
 * `-kmer number`: kmer number to pass to the assembler


## Pipelines

### BasicPipeline.pl
Runs aTRAM for a number of target sequences on a number of aTRAM databases.

```perl BasicPipeline.pl -samples SampleFile.txt -targets -TargetFile.txt [atram parameters -complete -protein]```

* `TargetFile.txt`: tab-delimited file; each line contains the target name and the path to its fasta file
* `SampleFile.txt`: tab-delimited file; each line contains the sample name and the path to its aTRAM-formatted database

#### Results:

`pipeline.log` contains the entire console output of the pipeline.

For each sample, there is a directory of results. Each of those directories contains the following:

* `.all.fasta` has all of the contigs pulled out in each iteration
* `.best.fasta` has the best contigs for each iteration
* `.results.txt` is the aTRAM results file for this sample.

### AlignmentPipeline.pl 
Runs aTRAM for a number of target sequences on a number of aTRAM databases, then aligns the resulting contigs back to the target sequences. Target sequences must be DNA, not amino acid.

```perl AlignmentPipeline.pl -samples SampleFile.txt -targets -TargetFile.txt [atram parameters -complete -protein]```

* `TargetFile.txt`: tab-delimited file; each line contains the target name and the path to its fasta file
* `SampleFile.txt`: tab-delimited file; each line contains the sample name and the path to its aTRAM-formatted database

#### Results:

For each target sequence, there is a directory of results. Each of those directories contains the following:

`results.txt` contains the overall summary of the pipeline results in tabular form: each sample's best contig is listed, with summary statistics for its bitscore and percent coverage of the target sequence.

`pipeline.log` contains the entire console output of the pipeline.

The `aTRAM` directory contains the per-sample results. For each sample, there are six files: 

* `.all.fasta` has all of the contigs pulled out in each iteration
* `.best.fasta` has the best contigs for each iteration
* `.complete.fasta` has the contigs that cover both ends of the target sequence
* `.align.fasta` uses MUSCLE to align the contigs from .best.fasta to the target sequence
* `.trimmed.fasta` takes the MUSCLE alignment and trims out any gaps relative to the target sequence.
* `.results.txt` is the aTRAM results file for this sample.

The `alignments` folder contains summary results over all samples. There are two fasta files:

* `.full.fasta` contains the entire sequence of the best-scoring contig for every sample, including introns and assemblies that go beyond the ends of the target sequence.
* `.trimmed.fasta` contains those same contigs aligned to the target sequence and trimmed of gaps relative to the target sequence. If the target sequence was a CDS, this file will not contain introns.


