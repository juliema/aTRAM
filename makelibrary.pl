#!/usr/bin/perl
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;


if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $short_read_archive = "";
my $output_file = "";
my $numlibraries = 8;
my $help = 0;

GetOptions ('input=s' => \$short_read_archive,
            'output=s' => \$output_file,
            'number=i' => \$numlibraries,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive) {
    pod2usage(-msg => "Must specify a short read archive in fasta or fastq form.");
}

unless ($output_file) {
    $output_file = $short_read_archive;
}

my $tempdir = dirname ("$short_read_archive");

# sort fasta short-read file
print "" . timestamp() . ": separating fasta/fastq file.\n";

open SEARCH_FH, "<", $short_read_archive;

my @out1_fhs = ();
my @out2_fhs = ();

for (my $i=0; $i<$numlibraries; $i++) {
	open my $out1_fh, ">", "$short_read_archive.temp.$i.1" or exit_with_msg ("couldn't create result file");
	open my $out2_fh, ">", "$short_read_archive.temp.$i.2" or exit_with_msg ("couldn't create result file");
	push @out1_fhs, $out1_fh;
	push @out2_fhs, $out2_fh;
}

my $name = "";
my $seq = "";
my $seqlen = 0;
while (my $line = readline SEARCH_FH) {
	chomp $line;
	if ($line =~ /[@>](.*?)([\s\/])([12])/) {
		if ($name ne "") {
			if ($name =~ /\/1/) {
				print {$out1_fhs[(hash_to_bucket($name))]} ">$name,$seq\n";
			} elsif ($name =~ /\/2/) {
				print {$out2_fhs[(hash_to_bucket($name))]} ">$name,$seq\n";
			}
		}
		$name = "$1\/$3";
		$seq = "";
	} elsif ($line =~ /^\+/){
		# is this a fastq quality line? eat chars to the length of the full sequence.
		while ($seqlen > 0) {
			$line = readline SEARCH_FH;
			chomp $line;
			$seqlen = $seqlen - length($line);
		}
	} else {
		$seq .= $line;
		$seqlen = length ($seq);
	}
}

if ($name =~ /\/1/) {
	print {$out1_fhs[(hash_to_bucket($name))]} ">$name,$seq\n";
} elsif ($name =~ /\/2/) {
	print {$out2_fhs[(hash_to_bucket($name))]} ">$name,$seq\n";
}
close SEARCH_FH;

for (my $i=0; $i<$numlibraries; $i++) {
	close $out1_fhs[$i];
	close $out2_fhs[$i];
}

for (my $i=0; $i<$numlibraries; $i++) {
	my $tempsort1 = "$short_read_archive.sort.$i.1";
	my $tempsort2 = "$short_read_archive.sort.$i.2";
	print "" . timestamp() . ": running sort on $short_read_archive.sort.$i.1.\n";
	system ("sort -t',' -k 1 -T $tempdir $short_read_archive.temp.$i.1 > $short_read_archive.sort.$i.1");
	print "" . timestamp() . ": running sort on $short_read_archive.sort.$i.2.\n";
	system ("sort -t',' -k 1 -T $tempdir $short_read_archive.temp.$i.2 > $short_read_archive.sort.$i.2");
	print "" . timestamp() . ": sorted.\n";
	system ("rm $short_read_archive.temp.$i.1");
	system ("rm $short_read_archive.temp.$i.2");
	if (-z "$short_read_archive.sort.$i.1") {
		exit_with_msg ("Sort failed. Are you sure $short_read_archive exists?");
	}
}

for (my $i=0; $i<$numlibraries; $i++) {
	print "" . timestamp() . ": Making $short_read_archive.$i.1.fasta.\n";
	open my $in_fh, "<", "$short_read_archive.sort.$i.1" or exit_with_msg ("couldn't create result file");
	open my $out_fh, ">", "$short_read_archive.$i.1.fasta" or exit_with_msg ("couldn't create result file");
	while (my $line = readline $in_fh) {
		chomp $line;
		my ($name, $seq) = split(/,/,$line);
		print $out_fh "$name\n$seq\n";
	}
	close $in_fh;
	close $out_fh;
	system ("rm $short_read_archive.sort.$i.1");

	print "" . timestamp() . ": Making $short_read_archive.$i.2.fasta.\n";
	open my $in_fh, "<", "$short_read_archive.sort.$i.2" or exit_with_msg ("couldn't create result file");
	open my $out_fh, ">", "$short_read_archive.$i.2.fasta" or exit_with_msg ("couldn't create result file");
	while (my $line = readline $in_fh) {
		chomp $line;
		my ($name, $seq) = split(/,/,$line);
		print $out_fh "$name\n$seq\n";
	}
	close $in_fh;
	close $out_fh;
	system ("rm $short_read_archive.sort.$i.2");

	# make the blast db from the first of the paired end files
	print "" . timestamp() . ": Making blastdb from $short_read_archive.$i.1.fasta.\n";
	system ("makeblastdb -in $short_read_archive.$i.1.fasta -dbtype nucl -out $short_read_archive.$i.db");
	system ("rm $short_read_archive.$i.1.fasta");
}

print "" . timestamp() . ": Finished\n";

sub exit_with_msg {
	my $msg = shift;
	print STDERR "$msg\n";
	exit 1;
}

sub timestamp {
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $mon++;
    $mon = sprintf("%02d", $mon);
    $min = sprintf("%02d", $min);
    $sec = sprintf("%02d", $sec);
    $hour = sprintf("%02d", $hour);
    $mday = sprintf("%02d", $mday);

    $year -= 100;
    my $time = "$hour:$min:$sec";
    my $date = "$year$mon$mday";
    return "$date $time";
}

sub hash_to_bucket {
	my $key = shift;
	my $bucket;

	if ($key =~ /(.+)\/\d/) {
		$key = $1;
	}
	$key =~ s/\D//g;
	$key =~ s/.*(.{3})$/$1/;

	$bucket = $key % $numlibraries;
	return $bucket;
}

__END__

=head1 NAME

sTRAM.pl

=head1 SYNOPSIS

makelibrary.pl -input short_read_archive [-output library_name] [-number int]

Takes a fasta or fastq file of paired-end short reads and prepares it for sTRAM.

=head1 OPTIONS

 -input:   short read archive.
 -output:  optional: prefix of output library (default is the same as -input).
 -number:  optional: the number of partitioned libraries to make.

=cut

