#     Manual for aTRAM

Here are the required commands for running aTRAM as well as a few things to keep in mind for maximum efficiency.


## Setup:
To determine if aTRAM can run on your computer:

### configure.pl

perl configure.pl 
  
  This script will tell you if you need to download any new programs. If you do make sure they are in your $PATH. You can either add them directly to your /usr/bin directory or add the path to the programs to your $PATH. 
Furthermore, configure.pl will check for the necessary de novo alignment programs including velvet, trinity and SOAPdenovo, you need to have at least one of these available. It will aslo check that you have muscle and blast, these programs are required.

## Running aTRAM

### Input File
The two paired end reads should be in one file either fasta or fasq. Concatenate the two files together if necessary.

### format_SRA.pl

perl format_SRA.pl inputfile.fastq

  This script will create aTRAM blast formatted databases of your Illumina run. 

### aTRAM.pl

perl  aTRAM.pl 




## Running on a cluster:

max_processes should equal number of threads

format_sra.pl max processes is set to 4



