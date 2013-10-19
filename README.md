##	sTRAM: Simple Target Restricted Assembly Method

sTRAM consists of two scripts:

1.	makelibrary.pl: this script takes a short-read archive and prepares it for sTRAMming.
2.	sTRAM.pl: this script runs the pipeline on a prepared short-read archive.

		sTRAM.pl -reads shortreadfile -target target.fasta [-ins_length int] [-exp_coverage int] [-iterations int] [-start_iteration int] [-log_file filename] [-use_ends] [-output filename]

		* -reads:     		short read archive (already run through makelibrary.pl).
		* -target:          fasta file with sequences of interest.
		* -output:	        optional: the prefix for the pipeline's output files (default name is the same as -reads).
		* -ins_length:	    optional: the size of the fragments used in the short-read library (default 300).
		* -exp_coverage:    optional: the expected coverage of the region for velvetg (default 30).
		* -iterations:      optional: the number of pipeline iterations (default 5).
		* -start_iteration: optional: if resuming from previous run, which iteration number to start from (default 0).
		* -log_file:        optional: a file to store output of the pipeline.
		* -use_ends:        optional: if this flag is present, use the first and last $ins_length of long contigs in the search.
