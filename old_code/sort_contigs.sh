# sorts an inputted fasta file of Illumina multiplexed short reads, such as the one from 3-getpairedend.pl.
# execute with the shell of your choice, i.e. "bash 3.5-sort_fasta.sh file.fasta" gives you file.fasta.sorted.fasta

infile=$1
outfile=$1.sorted.tab
gawk '{if (NF == 0) next; s = ""; for (i=2;i<=NF;i++) s = s$i; print length(s)","$1","s}' RS=">" $infile | sort -n -r | gawk '{ print $1 "\t" $2 "\t" $3}' FS="," > $outfile
