# sorts an inputted fasta file of Illumina multiplexed short reads, such as the one from 3-getpairedend.pl.
# execute with the shell of your choice, i.e. "bash 3.5-sort_fasta.sh file.fasta" gives you file.fasta.sorted.fasta

infile=$1
outfile=$1.sorted.fasta
tempfile=$1.temp
tempfile2=$1.temp2

tempdir=$(dirname $1)
gawk '{if (NF==0) next; sub(/lcl\|/,""); s = ""; for (i=2;i<=NF;i++) s = s$i; print $1","s}' RS=">" $infile > $tempfile
echo "wrote $tempfile"
date
sort -t',' -k 1 -T $tempdir $tempfile > $tempfile2
echo "wrote $tempfile2"
date
rm $tempfile
gawk '{print ">" $1 "\n" $2}' FS="," $tempfile2 > $outfile
echo "wrote $outfile"
date
rm $tempfile2
