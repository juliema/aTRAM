# Tips

## Creating an argument file

TODO: Add an atram.py example.

There are a lot of arguments to aTRAM and even I don't remember them all. To help with this a lot of people create Bash scripts once they have tuned the arguments for their needs. I prefer to use a slightly different method, an argument file. This is a text file that lists the arguments to a program, one argument in long form per line.

For the atram_preprocessor tutorial I would create a file, let's call it `atram_preprocessor.args`, like so:
```
--blast-db=/path/to/atram_db/tutorial
--end-1=/path/to/doc/data/tutorial_end_1.fasta.gz
--end-2=/path/to/doc/data/tutorial_end_2.fasta.gz
--gzip
```
And then you would use it like this:
```bash
../atram_preprocessor.py @atram_preprocessor.args
```
