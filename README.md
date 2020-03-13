# aTRAM: automated Target Restricted Assembly Method [![Build Status](https://travis-ci.org/juliema/aTRAM.svg?branch=master)](https://travis-ci.org/juliema/aTRAM)

- [Background](#Background)
- [Installation](#Installation)
- [Quick start](#Quick-start)
- [Documentation](#Documentation)

## Background

aTRAM ("automated target restricted assembly method") is an iterative assembler
that performs reference-guided local de novo assemblies using a variety of
available methods. It is well-suited to various tasks where Next-Generation
Sequence (NGS) data needs to be queried for gene sequences, such as
phylogenomics. The design philosophy is modular and expandable, with support
for four de-novo assemblers to date: Velvet, Abyss, Trinity, and Spades.

aTRAM 2 is a major overhaul of the aTRAM approach to assembling loci from (NGS)
data. The new code has been reimplemented in Python, and the approach to short
read library construction is completely revamped, resulting in major
performance and assembly improvements.

Please consult the reference below for more information about aTRAM1.0:
`Allen, JM, DI Huang, QC Cronk, KP Johnson. 2015. aTRAM automated target
restricted assembly method a fast method for assembling loci across divergent
taxa from next-generation sequencing data. BMC Bioinformatics 16:98
DOI 10.1186/s12859-015-0515-2`

The reference for aTRAM 2.0:
`Allen J.M., R. LaFrance, R. A. Folk, K. P. Johnson, and R. P. Guralnick
In Press.  aTRAM 2.0: An improved, flexible locus assembler for NGS data.
Evolutionary Informatics`

## Installation

You will need to have Python3 installed, as well as pip, a package manager for
Python.

```bash
git clone https://github.com/juliema/aTRAM.git
pip install --user --requirement atram/requirements.txt
```

### aTRAM uses these programs so you need to install them.

You will need to use a locally installed BLAST:

- [BLAST](
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download),
version 2.7.1

You will also need at least one of the supported assembly modules:

- [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/)
- [Trinity](http://trinityrnaseq.github.io/), version 2.5.1
- [Abyss](http://www.bcgsc.ca/platform/bioinfo/software/abyss), version 2.0.2
- [SPAdes](http://cab.spbu.ru/software/spades/), version 3.11.1

If you want to use the atram stitcher you will need to install exonerate:

- [exonerate](
https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-user-guide)

### Installation using [`conda`](https://www.anaconda.com/distribution/):

Alternatively, you can install both dependencies and `aTRAM` by using `conda`. Inside the `aTRAM` directory, run the following:

```bash
conda env create -f environment.yml
conda activate aTRAM
```

## Quick start

Note: aTRAM 2 is not backwards compatible with aTRAM 1. It is also best to
rebuild any libraries after major updates.

### Library Preparation

Use `atram_preprocessor.py` for this.

- Define your new library name with the --blast-db option. Which consists of a
path and the library prefix itself. This program will add suffixes to
differentiate different database files.

- Then give it your fastq files. You can either list the forward and reverse
read files, or put them into one file and use the --mixed-ends option.

 Under the hood, aTRAM is building BLAST databases and an SQLite3 database for
 rapid read retrieval.

```bash
atram_preprocessor.py \
  --blast-db=path_to_atram_library/LIBRARY_PREFIX \
  --end-1=path_to_reads/read_1.fastq \
  --end-2=path_to_reads/read_2.fastq
```

### Assembling Loci

`atram.py` uses the databases built by `atram_preprocessor.py` to assemble
loci.

- You need to give it the same --blast-db option from the preprocessor.
- You also need to give it a query sequence. The query sequence is a FASTA
file.
- An assembler choice. The assembler choice is one of the assemblers mentioned
above (velvet, trinity, abyss, or spades).
- And an output prefix. The `--output-prefix` works just like the
`--blast-db-prefix` with the directory part and the library prefix itself.

```bash
atram.py \
  --blast-db=path_to_atram_library/LIBRARY_PREFIX \
  --query=path_to_reference_loci/Locus.fasta \
  --assembler=ASSEMBLER_CHOICE \
  --output-prefix=path_to_output/OUTPUT_PREFIX
```

### Stitching genes from assembled loci

`atram_stitcher.py` Takes the output assemblies from `atram.py` and reference
amino acid targets and then stitches them together using an iterative process.

- Give it a directory containing the assemblies.
- A set of reference amino acid sequences in a FASTA file.
- A list of taxon names. One taxon per line.

```bash
atram_stitcher.py \
  --assemblies-dir=path_to_assemblies \
  --reference-genes=path_to_genes/ref_genes.fasta \
  --taxa=path_to/taxon_list.txt
```

## Documentation

- [aTRAM Overview](doc/introduction.md)
- [Tips for using aTRAM](doc/tips.md)
- [Troubleshooting](doc/troubleshooting.md)

### Tutorials
- [Pre-processing tutorial](doc/atram_preprocessor_tutorial.md)
- [aTRAM tutorial](doc/atram_tutorial.md)
- [Stitcher tutorial](doc/atram_stitcher_tutorial.md)

### Program Reference
- [atram_preprocessor.py](doc/atram_preprocessor_reference.md)
- [atram.py](doc/atram_reference.md)
- [atram_stitcher.py](doc/atram_stitcher_reference.md)
