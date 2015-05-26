#!/usr/bin/env perl
use strict;
use File::Basename;
use File::Temp qw(tempfile cleanup);
use File::Path qw (make_path);
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use System;
use Configuration;
use Mapreduce;

if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "Running " . basename($0) . " " . join (" ", @ARGV) . ", " . get_version() . "\n";

my $short_read_archive = "";
my $short_read_1 = "";
my $short_read_2 = "";
my $output_file = "";
my $log_file = "";
my $help = 0;
my $debug = 0;
my $numshards = 0;
my $max_processes = 4;

GetOptions ('input=s' => \$short_read_archive,
			'1input=s' => \$short_read_1,
			'2input=s' => \$short_read_2,
            'output=s' => \$output_file,
            'number=i' => \$numshards,
            'debug' => \$debug,
			'max_processes|processes=i' => \$max_processes,
            'log_file=s' => \$log_file,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

if ($log_file eq "") {
	$log_file = "$output_file.log";
}

set_log($log_file);
set_debug($debug);
printlog ("$runline");

my $is_fastq = 0;
my $srasize = 0;
unless(-e $short_read_archive) {
	# check to see if 1input and 2input were specified:
	unless ((-e $short_read_1) && (-e $short_read_2)) {
    	pod2usage(-msg => "Must specify a short read archive in fasta or fastq form.");
	} else {
		printlog "Two input files: $short_read_1 and $short_read_2\n";
		if ($short_read_1 =~ /\.f.*q/) {
			$is_fastq = 1;
		}
		$srasize = (-s $short_read_1);
	}
} else {
	printlog "One input file: $short_read_archive\n";
	if ($short_read_archive =~ /\.f.*q/) {
		$is_fastq = 1;
	}
	$srasize = (-s $short_read_archive);
}

# if it's a fastq file, the file is twice the size that it would be if it were a fasta.
if ($is_fastq) {
	$srasize = $srasize/2;
}

unless ($output_file) {
    $output_file = $short_read_archive;
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

# making a redirect file to make it easier for users to have something to specify.
my $db_file = "$output_file.atram";
open DB_FH, ">", $db_file;
print DB_FH "$output_file\n";
close DB_FH;

my $srasizeMB = $srasize / 1e6;
$srasizeMB =~ s/(\d*)\.(\d{2}).*/$1.$2/;

set_multiplier ($srasize);

# if the user didn't specify how many shards to make, we should make as many as we need so that they average 250MB each.
if ($numshards == 0) {
	$numshards = int($srasize / 2.5e8);

	# but if it's smaller than the minimum, make at least one shard.
	if ($numshards == 0) {
		$numshards = 1;
	}
	printlog ("Input is $srasizeMB MB; we will make $numshards shards.");
} else {
	printlog ("making $numshards shards.");
}

# declare how many shards we'll be making.
set_total_shards ($numshards);
my $max_shard = get_max_shard();

my @keys = ();
for (my $i=0;$i<$numshards; $i++) {
	$keys[$i] = 0;
}

my @tempfiles = ();

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

if ($short_read_archive ne "") {
	open SEARCH_FH, "<:crlf", $short_read_archive;
	while (my $line = readline SEARCH_FH) {
		chomp $line;
		if ($line =~ /^[@>](.*?)([\s\/])([12])/) {
			if ($name ne "") {
				if ($name =~ /\/1/) {
					print {$out1_fhs[$shard]} ">$name,$seq\n";
				} elsif ($name =~ /\/2/) {
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

	if ($name =~ /\/1/) {
		print {$out1_fhs[$shard]} ">$name,$seq\n";
	} elsif ($name =~ /\/2/) {
		print {$out2_fhs[$shard]} ">$name,$seq\n";
	}
	close SEARCH_FH;
} else {
	# there are two files: process the first and then the second.
	my @sra_files = ($short_read_1, $short_read_2);
	my @out_fhs = (\@out1_fhs, \@out2_fhs);
	
	for (my $i=0; $i<2; $i++) {
		$name = "";
		$seq = "";
		open SEARCH_FH, "<:crlf", $sra_files[$i];
		while (my $line = readline SEARCH_FH) {
			chomp $line;
			if ($line =~ /^[@>](.+?)(\/[12])*$/) {
				if ($name ne "") {
					print {@{$out_fhs[$i]}[$shard]} ">$name,$seq\n";
				}
				$name = "$1\/" . ($i+1);
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
		print {@{$out_fhs[$i]}[$shard]} ">$name,$seq\n";
		close SEARCH_FH;
	}	
}

my @pids = ();
printlog ("starting sort.");

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

for (my $i=0; $i<$numshards; $i++) {
    push @pids, fork_cmd ("sort", "-t',' -k 1 -T $output_path $out2_bucketfiles[$i] > $out2_sortedfiles[$i]");
	if (@pids >= ($max_processes - 1)) {
		# don't spawn off too many threads at once.
		wait_for_forks(\@pids);
	}
}
wait_for_forks(\@pids);
printlog ("sorted.");

for (my $i=0; $i<$numshards; $i++) {
	printlog ("Making $output_file.$i.1.fasta.");
	open my $in_fh, "<", "$out1_sortedfiles[$i]" or exit_with_msg ("couldn't read $out1_sortedfiles[$i]");
	open my $out_fh, ">", "$output_file.$i.1.fasta" or exit_with_msg ("couldn't create $output_file.$i.1.fasta");
	my $prevname = "";
	while (my $line = readline $in_fh) {
		chomp $line;
		my ($name, $seq) = split(/,/,$line);

		# if the new sequence is the same name as the last sequence, something is wrong.
		if ($prevname eq $name) {
			printlog ("Error detected in $short_read_archive: Sequence $name occurred at least twice.");
			exit 2;
		}
		$prevname = $name;
		print $out_fh "$name\n$seq\n";
	}
	close $in_fh;
	close $out_fh;

	printlog ("Making $output_file.$i.2.fasta.");
	open $in_fh, "<", "$out2_sortedfiles[$i]" or exit_with_msg ("couldn't read $out2_sortedfiles[$i]");
	open $out_fh, ">", "$output_file.$i.2.fasta" or exit_with_msg ("couldn't create $output_file.$i.2.fasta");
	$prevname = "";
	while (my $line = readline $in_fh) {
		chomp $line;
		my ($name, $seq) = split(/,/,$line);
		# if the new sequence is the same name as the last sequence, something is wrong.
		if ($prevname eq $name) {
			printlog ("Error detected in $short_read_archive: Sequence $name occurred at least twice.");
			exit 2;
		}
		$prevname = $name;
		print $out_fh "$name\n$seq\n";
	}
	close $in_fh;
	close $out_fh;
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

format_sra.pl [-input short_read_archive | -1input sra1 -2input sra2] [-output aTRAM_db_name] [-number int]

Takes a fasta or fastq file of paired-end short reads and creates an aTRAM database for use by aTRAM.pl.

=head1 OPTIONS

 -input:   short read archive in fasta or fastq format. Use this option if both ends are in a single file.
 -1input, -2input: short read archives for 1st and 2nd paired-end reads, in fasta or fastq format.
 -output:  optional: prefix of aTRAM database (default is the same as -input).
 -number:  optional: number of shards to create (default is however many are required for each to be ~250MB).

=cut

