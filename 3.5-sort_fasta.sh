# sorts an inputted fasta file of Illumina multiplexed short reads, such as the one from 3-getpairedend.pl.
# execute with the shell of your choice, i.e. "bash 3.5-sort_fasta.sh file.fasta" gives you file.fasta.sorted.fasta

infile=$1
outfile=$1.sorted.fasta
gawk '{if (NF==0) next; sub(/lcl\|/,""); s = ""; for (i=2;i<=NF;i++) s = s$i; print $1","s}' RS=">" $infile | sort | gawk '{print ">" $1 "\n" $2}' FS="," > $outfile
