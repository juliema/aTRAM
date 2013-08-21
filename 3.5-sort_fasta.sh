# sorts an inputted fasta file of Illumina multiplexed short reads, such as the one from 3-getpairedend.pl.
# execute with the shell of your choice, i.e. "bash file.fasta > sorted_file.fasta"

gawk -e '/#/ {sub(/lcl\|/,""); print $1 "," $2 $3 }' RS=">" $1 | sort | gawk -e '/#/ {print ">" $1 "\n" $2}' FS=","
