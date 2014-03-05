##	aTRAM: automated Target Restricted Assembly Method

aTRAM does targeted de novo assembly of loci from paired-end Illumina runs and consists of three scripts:

	1. configure.pl:  this script determines that the computer has the proper dependencies for aTRAM.
	2. fomrat_sra.pl: this script creates an aTRAM database from an Illumina paired-end short-read archive.
	3. aTRAM.pl: this script runs aTRAM with a target sequence and the formatted short-read archive.


	aTRAM.pl -sra shortreadfile -target target.fasta [-ins_length int] [-exp_coverage int] [-iterations int] [-start_iteration int] [-log_file filename] [-use_ends] [-output filename]

		* -reads:     	    short read archive (already run through format_sra.pl).
		* -target:          fasta file with sequences of interest.
		* -output:	    optional: the prefix for the pipeline's output files (default name is the same as -reads).
		* -ins_length:	    optional: the size of the fragments used in the short-read library (default 300).
		* -exp_coverage:    optional: the expected coverage of the region for velvetg (default 30).
		* -iterations:      optional: the number of pipeline iterations (default 5).
		* -start_iteration: optional: if resuming from previous run, which iteration number to start from (default 0).
		* -log_file:        optional: a file to store output of the pipeline.
		* -use_ends:        optional: if this flag is present, use the first and last $ins_length of long contigs in the search.
