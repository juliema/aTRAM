use strict;

my $short_read_archive = shift;

#initially, run the pipeline on the inputted gene:
my $search_fasta = shift;

# for (my $i=0; $i<5; $i++) {
my $i=0;
# 	print ("interation $i starting...\n");
# 	print ("\t blastn...\n");
# 	system("blastn -db $short_read_archive.db -query $search_fasta -outfmt 6 -num_threads 8 -out $short_read_archive.blast.$i");
# 	print ("\t retrieving reads...\n");
# 	system ("perl ~/TRAM/2.5-sequenceretrieval.pl $short_read_archive.1.fasta $short_read_archive.2.fasta $short_read_archive.blast.$i");
# 	print ("\t velveth...\n");
	system ("velveth $short_read_archive.velvet 31 -fasta -shortPaired $short_read_archive.blast.$i.sorted.fasta");
	print ("\t velvetg...\n");
	system ("velvetg $short_read_archive.velvet -ins_length 300 -exp_cov 30 -min_contig_lgth 200");
	$search_fasta = "$short_read_archive.$1.contigs.fa";
	system ("mv $short_read_archive.velvet/contigs.fa $search_fasta");
# }
