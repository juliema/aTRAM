## Manual for using aTRAM 2.0: automated Target Restricted Assembly Method

## Background
    1. What does it do?
    2. How it works

aTRAM 2.0 is a major overhaul of the aTRAM approach to assembling loci from phylogenomic data. The new code has been reimplemented in python, and the approach to short read library construction is completely revamped, resulting in major performance improvements.

aTRAM ("automated target restricted assembly method") is an iterative assembler that performs reference-guided local de novo assemblies using a variety of available methods. It is well-suited to various tasks where NGS data needs to be queried for gene sequences, such as phylogenomics. The design philosophy is modular and expandable, with support to date for four de-novo assemblers to date: Velvet, Abyss, SOAP, and Spades.

## Installation
     1. Python 3.0 or greater and a number of dependencies
     2. Dependencies given in requirements.txt

You will need to have Python3 installed, as well as pip, a package manager for python. Beyond these, it is easiest to handle aTRAM 2 dependencies by setting up a virtual environment, which is a contained workspace with internally installed python libraries. Run the following code in what you intend to be your working directory:

```
virtualenv venv -p python3
source venv/bin/activate
pip install -r requirements.txt 
```

You should see something like `(venv)` at the beginning of your command prompt after running the second line, indicating the environment is active. Once you have verified that the requirements installed with no errors, only the second line needs to be run before each aTRAM 2 session.

If you choose not to use virtual environments, you will likely have to specify python3.

URLs for software:

* BLAST: http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
* Velvet: https://www.ebi.ac.uk/~zerbino/velvet/
* Trinity: http://trinityrnaseq.sourceforge.net
* MAFFT: http://mafft.cbrc.jp/alignment/software/
* MUSCLE: http://www.drive5.com/muscle/
* Spades: http://bioinf.spbau.ru/spades

### Library Preparation
Use atram_preprocessor.py for this. You can either list the forward and reverse read files, or glob them with wildcards as below. Note that aTRAM 2 is not backwards compatible with aTRAM 1 libraries; it is also best to rebuild any libraries after major updates. 

``` python path_to_aTRAM/atram_preprocessor.py -c NUMBER_OF_THREADS -b path_to_atram_library/LIBRARY_PREFIX READ_NAME*.fastq ```
  
### Assembling Loci

``` 
python path_to_aTRAM/atram.py -b path_to_atram_library/LIBRARY_PREFIX -q path_to_reference_loci/Locus.fasta -i NUMBER_OF_ITERATIONS --cpus NUMBER_OF_THREADS  --kmer KMER_NUMBER -o path_to_output/LIBRARY_PREFIX.Locus.atram2.fasta --log-file path_to_output/LIBRARY_PREFIX.Locus.log -a ASSEMBLER_CHOICE
```

Fill in the capitalized portions with your options, and fill in the paths. 

There are many more options than this, so for reference list them like so:

```python path_to_atram/atram.py -h```

## Example of running a shell loop

In many cases it is convenient to run aTRAM 2 as a loop, assembling a set of genes for a set of taxa. These can be set up in two parts like so:

```
# Make aTRAM libraries
array=(sample1 sample2 sample3)

for a in "${array[@]}"; # Iterate through samples
do 
python path_to_aTRAM/atram_preprocessor.py -c 4 -b path_to_atram_library/lib_${a} path_to_input/${a}_P*.fq
done
```


