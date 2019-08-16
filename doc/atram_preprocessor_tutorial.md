# A tutorial for atram_preprocessor.py

This program reads a series of related FASTA or FASTQ files and builds an aTRAM database. This aTRAM database is actually two or more databases (typically several).

- An SQLite3 database that holds the contents of the FASTA/FASTQ files in a format that can be easily and quickly queried. It takes three pieces of information from the original files are: the sequence name, sequence end (1, 2, or none), and the sequence itself.


- A set of BLAST databases. atram.py uses multiple BLAST databases. This dataset division enables parallelized read queries and greatly improves performance even for serial queries.


- Note: That during the assembly process, in atram.py, other **temporary** SQLite3 and BLAST databases will also be built.

![atram_preprocessor.py](images/atram_preprocessor.png "aTRAM pre-processor")

## The input FASTA files

We have provided some input files that we can use to practice using aTRAM in the `doc/data` directory. For this tutorial, we will be using the `tutorial_end_1.fasta.gz` and `tutorial_end_2.fasta.gz` files.

"Where do we want to store the built aTRAM databases?" I prefer to keep my aTRAM built databases in their own directory on a large disk. So, I create a `atram_db` directory and point aTRAM to that. In Linux Bash this looks like `mkdir -p /path/to/my/dir/atram_db'.


```python

```
