use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

my $short_read_archive = 0;
my $search_fasta = 0;
my $ins_length = 300;
my $iterations = 5;
my $start_iter = 0;
my $help = 0;
my $log_file = 0;

GetOptions ('reads=s' => \$short_read_archive,
            'target=s' => \$search_fasta,
            'insert_length|ins_length=i' => \$ins_length,
            'iterations=i' => \$iterations,
            'start_iteration=i' => \$start_iter,
            'log_file=s' => \$log_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive and $search_fasta) {
    pod2usage(-msg => "Must specify a short read archive (that has already been prepared with 0-prepare_files.pl) and a target gene in fasta form.");
}

print $runline;

my $executing_path = dirname(__FILE__);

my $log_fh;
if ($log_file) {
	open $log_fh, ">", $log_file or die "couldn't open $log_file\n";
} else {
	$log_fh = *STDOUT;
}

my $cmd;
if ($start_iter > 0) {
	$search_fasta = "$short_read_archive.$start_iter.contigs.fa";
}

open OUT_FH, ">>", "$short_read_archive.all.fasta";

for (my $i=$start_iter; $i<$iterations; $i++) {
	print ("interation $i starting...\n");
	$cmd = "blastn -db $short_read_archive.db -query $search_fasta -outfmt 6 -num_threads 8 -out $short_read_archive.blast.$i";
	print $log_fh ("\t$cmd\n");
	system($cmd);

	$cmd = "perl ~/TRAM/2.5-sequenceretrieval.pl $short_read_archive.1.fasta $short_read_archive.2.fasta $short_read_archive.blast.$i";
	print $log_fh ("\t$cmd\n");
	system($cmd);

	$cmd = "velveth $short_read_archive.velvet 31 -fasta -shortPaired $short_read_archive.blast.$i.sorted.fasta";
	print $log_fh ("\t$cmd\n");
	system($cmd);

	$cmd = "velvetg $short_read_archive.velvet -ins_length $ins_length -exp_cov 30 -min_contig_lgth 200";
	print $log_fh ("\t$cmd\n");
	system($cmd);
	$search_fasta = "$short_read_archive.$i.contigs.fa";
	system ("mv $short_read_archive.velvet/contigs.fa $search_fasta");

	print OUT_FH `cat $search_fasta | gawk '{sub(/>/,">$1"); print $0}'`;

	$cmd = "bash $executing_path/5.5-sort_contigs.sh $search_fasta";
	system($cmd);

	open FH, "<", "$search_fasta.sorted.tab";
	my @contigs = <FH>;
	close FH;
	my @new_contigs = ();

	for (my $j=0; $j<@contigs; $j++) {
		my ($len,$name,$seq) = split (/\t/,$contigs[$j]);
		if ($len > ($ins_length * 5)) {
			my $name_start = $name. "_start";
			my $name_end = $name . "_end";
			$seq =~ /^(.{$ins_length}).*(.{$ins_length})$/;
			my $seq_start = $1;
			my $seq_end = $2;
			push @new_contigs, "$ins_length\t$name_start\t$seq_start";
			push @new_contigs, "$ins_length\t$name_end\t$seq_end";
		} else {
			push @new_contigs, $contigs[$j];
		}
	}

	open FH, ">", "$search_fasta";
	foreach my $line (@new_contigs) {
		chomp $line;
		$line =~ /(.*?)\t(.*?)\t(.*)/;
		print FH ">$2\n$3";
	}
	close FH;
# 	system ("mv $search_fasta.sorted.tab $search_fasta");

}

close $log_fh;


__END__

=head1 NAME

pipeline.pl

=head1 SYNOPSIS

pipeline.pl -reads -target [-ins_length]

=head1 OPTIONS

  -reads:     		short read archive (already run through 0-prepare_files.pl).
  -target:          fasta file with sequences of interest.
  -ins_length:	    optional: if specified, the size of the fragments used in the short-read library (default 300).

=head1 DESCRIPTION

Takes a fasta file and finds aligned regions in each sequence in the fasta file that
match the reference sequence(es). Returns a fasta file of aligned regions of similarity.
Uses BLASTN to find regions of similarity.

=cut

