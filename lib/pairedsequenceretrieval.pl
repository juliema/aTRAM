# given two fasta files corresponding to paired ends of short reads, find the specified sequences and then output them into a single output file.
#!/usr/bin/perl
use strict;
use File::Temp qw/ tempfile tempdir /;

if (@ARGV < 3) {
	die "Usage: 2.5-pairedsequenceretrieval.pl shortreads.#.fasta sequencelist outfile\n";
}

my $fastafile = shift;
my $sequencelist = shift;
my $outfile = shift;

if ($fastafile !~ /#/) {
	die "fasta file must have '#' in name, to be replaced by 1 or 2 for the paired end files.\n";
}

unless (-e $sequencelist) {
	die "File $sequencelist does not exist.\n";
}

my $fastafile_1 = "$fastafile";
$fastafile_1 =~ s/#/1/;
my $fastafile_2 = "$fastafile";
$fastafile_2 =~ s/#/2/;

unless (-e $fastafile_2) {
	die "Files $fastafile_1 and $fastafile_2 do not exist.\n";
}


my (undef, $seq_names) = tempfile(UNLINK => 1);

system ("gawk '{sub(\"/1\",\"\");print \$2;}' $sequencelist | sort | uniq -u > $seq_names");

open LIST_FH, "<", "$seq_names";
open FA1_FH, "<", "$fastafile_1";
open FA2_FH, "<", "$fastafile_2";
open OUT_FH, ">", "$outfile";

my $curr_name = readline LIST_FH;

while ($curr_name) {
	chomp $curr_name;
	my $fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
	my $fa_seq2 = (readline FA2_FH) . (readline FA2_FH);

	if ($fa_seq1 eq "") { last; }
	if ($fa_seq2 eq "") { last; }

	if ($fa_seq1 =~ /$curr_name/) {
		print OUT_FH "$fa_seq1$fa_seq2";
		$curr_name = readline LIST_FH;
	}
}

close LIST_FH;
close FA1_FH;
close FA2_FH;
close OUT_FH;
