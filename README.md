##	aTRAM: automated Target Restricted Assembly Method
doi:10.5281/zenodo.10431

aTRAM performs targeted de novo assembly of loci from paired-end Illumina runs. 

When using this software please cite:

```Allen, JM, DI Huang, QC Cronk, KP Johnson. 2015. aTRAM automated target restricted assembly method a fast method for assembling loci across divergent taxa from next-generation sequencing data. BMC Bioinformatics 16:98 DOI 10.1186/s12859-015-0515-2```

###Installation###
aTRAM can be run directly from the downloaded Github repo. To make sure you have all the software required by aTRAM, we provide ```configure.pl``` to check for dependencies.

On OSX, aTRAM is also available through [Homebrew](https://github.com/Homebrew/homebrew-science). Run ```brew tap homebrew/science``` and then ```brew install atram``` to install version 1.04.

if you get an error with gcc try ```brew install homebrew/versions/gcc48``` then ```brew install atram```


###Preparing an aTRAM database###
Given an Illumina paired-end short-read archive in fastq or fasta form, create a master aTRAM database from the short reads:

```format_sra.pl -input my_pe_library.fa|fq -out my_atram_db```

###Finding a homolog with aTRAM###
Given an aTRAM database and a target sequence in fasta form, run the main aTRAM script:

```aTRAM.pl -reads my_atram_db -target target.fasta [-ins_length int] [-exp_coverage int] [-iterations int] [-output filename]```

aTRAM has many command-line options as well:

```
  pipeline parameters:
  -sra|database|db: aTRAM database name (already run through format_sra.pl).
  -target:          fasta file with sequences of interest.
  -output:	        optional: the prefix for the pipeline's output files (default name is the same as -sra).
  -log_file:        optional: a file to store output of the pipeline.
  -tempfiles:       optional: use this name to save the intermediate files from the run.
  -iterations:      optional: the number of pipeline iterations (default 5).
  -start_iteration: optional: if resuming from previous run, which iteration number to start from (default 0).

  optional parameters:
  -protein:         if the target sequence is a protein fasta file (not mandatory, aTRAM will guess).
  -complete:        if specified, automatically quits when a complete homolog is recovered.
  -fraction:        if specified, use only specified fraction of the aTRAM database.
  -processes:       if specified, aTRAM will use no more than this number of processes for multiprocessing.

  optional assembly parameters:
  -assembler:       software to be used for targeted assembly (default is Velvet).
  -ins_length:	    the size of the fragments used in the short-read library (default 300).
  -exp_coverage:    the expected coverage of the region for velvetg (default 30).
  -kmer:            kmer size for assemblers that use it (default 31).

  optional values for searching short reads:
  -evalue:          default value is 10e-10.

  optional values for blast-filtering contigs:
  -bitscore:        default value is 70.
  -length:          default value is 100.
```
