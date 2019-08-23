# aTRAM Stitcher

This program will find and stitch together exons from targeted
assemblies using amino acid targets and DNA assemblies.

 ## Arguments

`-h, --help`

Show this help message and exit.

`--version`

Show program's version number and exit.

`-T TAXA, --taxa TAXA`

A text file of all of your taxon names.

`-r FASTA, --reference-genes FASTA, --refs FASTA`

Reference amino acid sequences in a FASTA file.

`-a PATH, --assemblies-dir PATH`

The path to the DNA contigs.

`-O OVERLAP, --overlap OVERLAP`

Contigs must overlap by this many codons before it is
considered a real overlap.

`-t DIR, --temp-dir DIR`

Place temporary files in this directory. All files
will be deleted after aTRAM completes. The directory
must exist.

`--keep-temp-dir`

This flag will keep the temporary files in the --temp-
dir around for debugging.

`-l LOG_FILE, --log-file LOG_FILE`

Log file (full path). The default is
"atram_stitcher_<date>.log".

`-i N, --iterations N `

The number of times to run the main stitcher loop.
This must be either 1 or 2, the default is 2.

`-o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX`

This is the prefix of all of the output files. So you
can identify different stitcher output file sets. You
may include a directory as part of the prefix. The
stitcher will add suffixes to differentiate output
files.

`-f FILE_FILTER, --file-filter FILE_FILTER`

Use this to filter files in the assemblies directory.
For example '*filtered*.fasta' will select all fasta
files in the assemblies directory with the word
filtered in them. The default is to select all fasta
files in the assemblies directory '*.fasta'.

`--reference-name`

Prepend the reference name to the final assembled gene
name? if false the gene name in the reference file
with just be the <taxon-name> if you select this then
the assembled gene name will be <reference-
name>.<taxon-name>.
