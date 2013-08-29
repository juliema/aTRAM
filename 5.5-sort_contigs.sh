# sorts an inputted fasta file of Illumina multiplexed short reads, such as the one from 3-getpairedend.pl.
# execute with the shell of your choice, i.e. "bash 3.5-sort_fasta.sh file.fasta" gives you file.fasta.sorted.fasta

infile=$1
outfile=$1.sorted.tab
gawk '{if ($0 - /^$/) next; s = ""; for (i=2;i<=NF;i++) s = s$i; print $1","s}' RS=">" $infile | gawk '{s = ""; for (i=1;i<=NF;i++) s = s"_"$i;print $4 "," s}' FS="_" | sort -n -r | gawk '{sub ("_","",$2); print $1 "\t" $2 "\t" $3}' FS="," > $outfile

# gawk '{sub(/[[:digit:]]+_/,"",$1); print ">" $1 "\n" $2}' FS="," > $outfile
