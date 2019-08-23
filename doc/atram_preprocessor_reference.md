# aTRAM Preprocessor

This script prepares data for use by the atram.py
script. It takes fasta or fastq files of paired-end (or
single-end) sequence reads and creates a set of atram
databases.

You need to prepare the sequence read archive files so that the
header lines contain only a sequence ID with the optional
paired-end suffix at the end of the header line. The separator
for the optional trailing paired-end suffix may be a space,
a slash "/", a dot ".", or an underscore "_".

For example:

    >DBRHHJN1:427:H9YYAADXX:1:1101:10001:77019/1
    GATTAA...
    >DBRHHJN1:427:H9YYAADXX:1:1101:10001:77019/2
    ATAGCC...
    >DBRHHJN1:427:H9YYAADXX:1:1101:10006:63769/2
    CGAAAA...

 ## Arguments

`-h, --help`            

Show this help message and exit.

`--version`

Show program's version number and exit.

`--end-1 FASTA/Q [FASTA/Q ...], -1 FASTA/Q [FASTA/Q ...]`

Sequence read archive files that have only end 1
sequences. The sequence names do not need an end
suffix, we will assume the suffix is always 1. The
files are in fasta or fastq format. You may enter more
than one file or you may use wildcards.

If there are end 2 sequences in the file then you will want to use the
`--mixed-ends` argument instead.

`--end-2 FASTA/Q [FASTA/Q ...], -2 FASTA/Q [FASTA/Q ...]`

Sequence read archive files that have only end 2
sequences. The sequence names do not need an end
suffix, we will assume the suffix is always 2. The
files are in fasta or fastq format. You may enter more
than one file or you may use wildcards.

If there are end 1 sequences in the file then you will want to use the
`--mixed-ends` argument instead.

`--mixed-ends FASTA/Q [FASTA/Q ...], -m FASTA/Q [FASTA/Q ...]`

Sequence read archive files that have a mix of both
end 1 and end 2 sequences (or single ends). The files
are in fasta or fastq format. You may enter more than
one file or you may use wildcards.

The sequence names must have a sequence end suffix like "/1" or "_2".

`--single-ends FASTA/Q [FASTA/Q ...], -0 FASTA/Q [FASTA/Q ...]`

Sequence read archive files that have only unpaired
sequences. Any sequence suffix will be ignored. The
files are in fasta or fastq format. You may enter more
than one file or you may use wildcards.

This option will ignore any sequence ends in the in the sequence name.

`-b DB, --blast-db DB, --db DB`

This is the prefix of all of the blast database files.
So you can identify different blast database sets. You
may include a directory as part of the prefix. The
default is "./atram_<today's date>".

For example, if you want to keep you data base in a directory called:
`/home/my_dir/atram_db/`. And you want to identify that these sequences are for
_Canis lupus_. Then you might prefix this database with
`/home/my_dir/atram_db/canis_lupus`. aTRAM will create a bunch a file in that
directory with names starting with 'canis_lupus'. Like 
'canis_lupus.001.blast.nhr' or 'canis_lupus.sqlite.db', etc. This allows you to
keep many aTRAM databases in the same directory.


`--cpus CPUS, --processes CPUS, --max-processes CPUS`

Number of CPU threads to use. On this machine the
default is to use the number of processes on your computer minus 4. More
threads can improve the speed but if you use too many you will slow down the
computer for any other use.

`-t DIR, --temp-dir DIR`

Place temporary files in this directory. All files
will be deleted after aTRAM completes. The directory
must exist.

Sometimes using the system temporary directory is not the right option because
aTRAM is filling up that disk or it is too slow or you want to debug issues.
Note: If you want to keep the temporary files around for debugging then you
should also use the `--keep-temp-dir` option.

`--keep-temp-dir`

This flag will keep the temporary files in the --temp-dir around for debugging.

`-l LOG_FILE, --log-file LOG_FILE`

Log file (full path). The default is to use the DB and
program name to come up with a name like
"<DB>_atram_preprocessor.log".

`-s SHARDS, --shards SHARDS, --number SHARDS`

Number of blast DB shards to create. The default is to
have each shard contain roughly 250MB of sequence
data.

`--path PATH`

If makeblastdb is not in your $PATH then use this to prepend directories to
your path. For instance, if you installed makeblastdb in 
`/home/my_dir/bin/makeblastdb`, you could use this option to find it.

`--fasta`

Are these fasta files? If you do not specify either --fasta or --fastq then
aTRAM will guess the file type by looking at the last character of the file
name. This option is most useful when you are dealing with compressed files
which will defeat aTRAM's FASTA/Q format guessing algorithm.

`--fastq`

Are these fastq files? If you do not specify either --fasta or --fastq then
aTRAM will guess the file type by looking at the last character of the file
name. This option is most useful when you are dealing with compressed files
which will defeat aTRAM's FASTA/Q format guessing algorithm.

`--gzip`

Are these gzip files? aTRAM does not try to guess if the file is compressed.

`--bzip`

Are these bzip files? aTRAM does not try to guess if the file is compressed.
