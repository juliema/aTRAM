# Tips

- [Assembling multiple genes against a library](#Assembling-multiple-genes-against-a-library)
- [Example of running a shell loop](#Example-of-running-a-shell-loop)
- [Creating an argument file](#Creating-an-argument-file)
- [Backwards compatibility](#Backwards-compatibility)


## Assembling multiple genes against a library

 aTRAM2.0 can assemble a set of genes against a single library. Create a single
 file with multiple fasta-formatted sequences and then simply use `-Q QUERY_SPLIT`
 where QUERY_SPLIT is the name of the file you created above.

## Example of running a shell loop

In many cases it is convenient to run aTRAM 2 as a loop, assembling a set of
genes for a set of taxa. These can be set up in two parts, as shown below.
Note that aTRAM2 has built in functions supporting assembly of many genes
against a library, as described just above.

```bash
# Make aTRAM libraries
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do
  atram_preprocessor.py -c 4 -b path_to_atram_library/lib_${a} path_to_input/${a}_P*.fq
done
```

The part `${a}_P*.fq` will have to be modified to match the name pattern of
your input fastq files.

Then, supposing we have a set of genes stored in a single file and wish to use
Abyss:

```bash
# Assemble genes
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do
  atram.py -b path_to_atram_library/lib_${a} -Q file_name -i 5 --cpus 4  --kmer 64 -o path_to_output/lib_${a}.atram2.fasta --log-file path_to_output/lib_${a}.log -a abyss
done
```

## Creating an argument file

There are a lot of arguments to aTRAM and even I don't remember them all. To
help with this a lot of people create Bash scripts once they have tuned the
arguments for their needs. I prefer to use a slightly different method, an
argument file. This is a text file that lists the arguments to a program, one
argument, in long form, per line.

For the atram_preprocessor tutorial I would create a file, let's call it
`atram_preprocessor.args`, like so:

```
--blast-db=/path/to/atram_db/tutorial
--end-1=/path/to/doc/data/tutorial_end_1.fasta.gz
--end-2=/path/to/doc/data/tutorial_end_2.fasta.gz
--gzip
```

And then you would use it like this:

```bash
atram_preprocessor.py @atram_preprocessor.args
```
You can still add command-line arguments. Like so:

```bash
atram_preprocessor.py @atram_preprocessor.args --cpus=8
```

## Backwards compatibility

For any tools that depend on the output format of aTRAM 1.0, this script will
perform the conversion of fasta headers:

```
for i in $(find . -name "*.fasta"); do
  sed 's/.* iteration=/>/g' ${i} | sed 's/ contig_id/.0_contigid/g' | sed 's/contigid.*length_//g' | sed 's/_cov.* score=/_/g' | sed 's/\.[0-9]*$//g' > ${i}.aTRAM1.fasta
done
```
