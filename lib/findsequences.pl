# given two fasta files corresponding to paired ends of short reads, find the specified sequences and then output them into a single output file.
use strict;
use File::Temp qw/ tempfile tempdir /;

if (@ARGV < 3) {
	die "Usage: 6.5-findsequences.pl fastafile sequencelist outfile\n";
}

my $fastafile = shift;
my $sequencelist = shift;
my $outfile = shift;

unless (-e $sequencelist) {
	die "File $sequencelist does not exist.\n";
}

unless (-e $fastafile) {
	die "File $fastafile does not exist.\n";
}

my (undef, $seq_names) = tempfile(UNLINK => 1);
my (undef, $fasta_sort) = tempfile(UNLINK => 1);

system ("cat $sequencelist | sort | uniq -u > $seq_names");
system ("gawk '{if (NF == 0) next; s = \"\"; for (i=2;i<=NF;i++) s = s\$i; print \$1\",\"s}' RS=\">\" $fastafile | sort | gawk '{ print \">\"\$1\"\\n\"\$2}' FS=\",\" > $fasta_sort");

open LIST_FH, "<", "$seq_names";
open FA_FH, "<", "$fasta_sort";
open OUT_FH, ">", "$outfile";

my $curr_name = readline LIST_FH;

while ($curr_name) {
	chomp $curr_name;
	my $fa_seq = (readline FA_FH) . (readline FA_FH);

	if ($fa_seq eq "") { last; }

	if ($fa_seq =~ /$curr_name/) {
		print OUT_FH "$fa_seq";
		$curr_name = readline LIST_FH;
	}
}

close LIST_FH;
close FA_FH;
close OUT_FH;
