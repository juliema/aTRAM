#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw(tempfile cleanup);
use File::Copy qw(copy);
use File::Path qw(make_path);
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/../lib";
use System;
use Configuration;
use Mapreduce;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $atram_db = "";
my $output_file = "";
my $help = 0;
my $debug = 0;
my $rename = 0;
my $numshards = 0;
my $max_processes = 4;

GetOptions ('input=s' => \$atram_db,
            'output=s' => \$output_file,
            'number=i' => \$numshards,
            'debug' => \$debug,
            'rename' => \$rename,
			'max_processes|processes=i' => \$max_processes,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($atram_db) {
    pod2usage(-msg => "Must specify a short read archive in fasta or fastq form.");
}

unless ($output_file) {
    $output_file = $atram_db;
}

# find the path specified in the .atram file, if provided.
if ($atram_db =~ /\.atram$/) {
	open ATRAM_FH, "<:crlf", $atram_db;
	$atram_db = readline ATRAM_FH;
	chomp $atram_db;
	close ATRAM_FH;
}

unless($atram_db) {
    pod2usage(-msg => "Must specify a short read archive (that has already been prepared with format_sra.pl).");
}

# check to make sure that the specified short read archive exists:
unless ((-e "$atram_db.0.1.fasta") && (-e "$atram_db.0.2.fasta")) {
	pod2usage(-msg => "Short read archive does not seem to be in the format made by format_sra.pl. Did you specify the name correctly?");
}

my $original_numshards = get_total_shards("$atram_db");
if ($rename != 0) {
	# we just want to rename the atram_db to the new output name
	$numshards = $original_numshards;
}

if ($numshards == 0) {
    pod2usage(-msg => "No number of shards specified.");
}

# make the output file an absolute path, just to be safe.
$output_file = File::Spec->rel2abs($output_file);
my $output_path = dirname ($output_file);

# check to see if the path for the output_file exists; if not, create it.
unless (-d $output_path) {
	make_path ($output_path);
}

# Look in the config.txt file to find the correct paths to binaries.
Configuration::initialize();

my $log_file = "$output_file.log";
set_log ($log_file);

# making a redirect file to make it easier for users to have something to specify.
my $db_file = "$output_file.atram";
open DB_FH, ">", $db_file;
print DB_FH "$output_file\n";
close DB_FH;

my @pids = ();
my @tempfiles = ();
my @keys = ();

if ($original_numshards != $numshards) {
	my $srasize = (-s $atram_db);
	my $srasizeMB = $srasize / 1e6;
	$srasizeMB =~ s/(\d*)\.(\d{2}).*/$1.$2/;

	set_multiplier ($srasize);

	printlog ("making $numshards shards from $atram_db.");

	# declare how many shards we'll be making.
	set_total_shards ($numshards);
	my $max_shard = get_max_shard();

	for (my $i=0;$i<$numshards; $i++) {
		$keys[$i] = 0;
	}


	my @out1_fhs = ();
	my @out2_fhs = ();
	my @out1_bucketfiles = ();
	my @out2_bucketfiles = ();
	my @out1_sortedfiles = ();
	my @out2_sortedfiles = ();

	# setting up tempfiles for sorting:
	for (my $i=0; $i<$numshards; $i++) {
		my ($fh, $filename) = tempfile("$output_file.bucket.$i.1.XXXX", UNLINK => 1);
		push @out1_bucketfiles, $filename;
		push @out1_fhs, $fh;

		($fh, $filename) = tempfile("$output_file.bucket.$i.2.XXXX", UNLINK => 1);
		push @out2_bucketfiles, $filename;
		push @out2_fhs, $fh;

		($fh, $filename) = tempfile("$output_file.sorted.$i.1.XXXX", UNLINK => 1);
		push @out1_sortedfiles, $filename;

		($fh, $filename) = tempfile("$output_file.sorted.$i.2.XXXX", UNLINK => 1);
		push @out2_sortedfiles, $filename;
	}

	# Divide fasta/fastq short reads into buckets for sorting.
	printlog ("Dividing fasta/fastq file into buckets for sorting.");
	my $name = "";
	my $seq = "";
	my $seqlen = 0;
	my $shard = 0;

	for (my $i=0; $i < $original_numshards; $i++) {
		my $short_read_archive = "$atram_db.$i.1.fasta";
		open SEARCH_FH, "<:crlf", $short_read_archive;
		while (my $line = readline SEARCH_FH) {
			chomp $line;
			if ($line =~ /^[@>](.*?)([\s\/])([12])/) {
				if ($name ne "") {
					if ($name =~ /\/1/) {
						print {$out1_fhs[$shard]} ">$name,$seq\n";
					}
				}
				$name = "$1\/$3";
				$shard = map_to_shard($name);
				$keys[$shard]++;
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
			print {$out1_fhs[$shard]} ">$name,$seq\n";
		}
		close SEARCH_FH;

		$short_read_archive = "$atram_db.$i.2.fasta";
		open SEARCH_FH, "<:crlf", $short_read_archive;
		while (my $line = readline SEARCH_FH) {
			chomp $line;
			if ($line =~ /^[@>](.*?)([\s\/])([12])/) {
				if ($name ne "") {
					if ($name =~ /\/2/) {
						print {$out2_fhs[$shard]} ">$name,$seq\n";
					}
				}
				$name = "$1\/$3";
				$shard = map_to_shard($name);
				$keys[$shard]++;
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

		if ($name =~ /\/2/) {
			print {$out2_fhs[$shard]} ">$name,$seq\n";
		}
		close SEARCH_FH;
	}
	printlog ("starting sort.");
	#
	for (my $i=0; $i<$numshards; $i++) {
		close $out1_fhs[$i];
		close $out2_fhs[$i];
		push @pids, fork_cmd ("sort", "-t',' -k 1 -T $output_path $out1_bucketfiles[$i] > $out1_sortedfiles[$i]");
		if (@pids >= ($max_processes - 1)) {
			# don't spawn off too many threads at once.
			wait_for_forks(\@pids);
		}
	}
	wait_for_forks(\@pids);
	#
	for (my $i=0; $i<$numshards; $i++) {
		push @pids, fork_cmd ("sort", "-t',' -k 1 -T $output_path $out2_bucketfiles[$i] > $out2_sortedfiles[$i]");
		if (@pids >= ($max_processes - 1)) {
			# don't spawn off too many threads at once.
			wait_for_forks(\@pids);
		}
	}
	wait_for_forks(\@pids);
	printlog ("sorted.");
	#
	for (my $i=0; $i<$numshards; $i++) {
		printlog ("Making $output_file.$i.1.fasta.");
		open my $in_fh, "<", "$out1_sortedfiles[$i]" or exit_with_msg ("couldn't read $out1_sortedfiles[$i]");
		open my $out_fh, ">", "$output_file.$i.1.fasta" or exit_with_msg ("couldn't create $output_file.$i.1.fasta");
		while (my $line = readline $in_fh) {
			chomp $line;
			my ($name, $seq) = split(/,/,$line);
			print $out_fh "$name\n$seq\n";
		}
		close $in_fh;
		close $out_fh;
	#
		printlog ("Making $output_file.$i.2.fasta.");
		open $in_fh, "<", "$out2_sortedfiles[$i]" or exit_with_msg ("couldn't read $out2_sortedfiles[$i]");
		open $out_fh, ">", "$output_file.$i.2.fasta" or exit_with_msg ("couldn't create $output_file.$i.2.fasta");
		while (my $line = readline $in_fh) {
			chomp $line;
			my ($name, $seq) = split(/,/,$line);
			print $out_fh "$name\n$seq\n";
		}
		close $in_fh;
		close $out_fh;
	}

} else {
	printlog ("Number of shards requested is the same as original number of shards in $atram_db.");
	for (my $i=0; $i<$numshards; $i++) {
		copy ("$atram_db.$i.1.fasta", "$output_file.$i.1.fasta");
		copy ("$atram_db.$i.2.fasta", "$output_file.$i.2.fasta");
	}
}

printlog ("Making blastdbs.");
for (my $i=0; $i<$numshards; $i++) {
	# make the blast db from the first of the paired end files
	push @pids, fork_cmd (get_bin("makeblastdb"), "-in $output_file.$i.1.fasta -dbtype nucl -out $output_file.$i.db");
	if (@pids >= ($max_processes - 1)) {
		# don't spawn off too many threads at once.
		wait_for_forks(\@pids);
	}
}

wait_for_forks(\@pids);

cleanup();
foreach my $tempfile (@tempfiles) {
	printlog ("removing $tempfile");
	system ("rm $tempfile");
}

open DB_FH, ">>", $db_file;
for (my $i=0;$i<$numshards; $i++) {
	print DB_FH "shard $i has\t$keys[$i] keys\n";
}
close DB_FH;

set_total_shards(0); # reset total shards to force a file count.
printlog ("Finished: aTRAM database $output_file has " . (get_total_shards($output_file)) . " shards.");


__END__

=head1 NAME

format_sra.pl

=head1 SYNOPSIS

format_sra.pl -input atram_db [-output aTRAM_db_name] [-number int]

Takes a fasta or fastq file of paired-end short reads and creates an aTRAM database for use by aTRAM.pl.

=head1 OPTIONS

 -input:   short read archive.
 -output:  optional: prefix of aTRAM database (default is the same as -input).
 -number:  optional: number of shards to create (default is however many are required for each to be ~500MB).

=cut

