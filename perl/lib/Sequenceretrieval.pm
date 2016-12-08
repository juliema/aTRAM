#!/usr/bin/env perl
package Sequenceretrieval;
use strict;
use File::Temp qw/ tempfile /;
use System;
use Parsing;

BEGIN {
	require Exporter;
	# set the version for version checking
	our $VERSION     = 1.00;
	# Inherit from Exporter to export functions and variables
	our @ISA         = qw(Exporter);
	# Functions and variables which are exported by default
	our @EXPORT      = qw(fork_pair_retrieval pairedsequenceretrieval findsequences);
	# Functions and variables which can be optionally exported
	our @EXPORT_OK   = qw();
}

sub fork_pair_retrieval {
	my $fastafile = shift;
	my $sequencelist = shift;
	my $outfile = shift;

    my $child_pid = fork();
    unless ($child_pid) { #child process
		pairedsequenceretrieval ($fastafile, $sequencelist, $outfile);
		exit;
    } else { #parent process
        return $child_pid;
    }
}

sub pairedsequenceretrieval {
	my $fastafile = shift;
	my $sequencelist = shift;
	my $outfile = shift;

	if ($fastafile !~ /#/) {
		print "fasta file must have '#' in name, to be replaced by 1 or 2 for the paired end files.\n";
		return 0;
	}

	unless (-e $sequencelist) {
		print "File $sequencelist does not exist.\n";
		return 0;
	}

	my $fastafile_1 = "$fastafile";
	$fastafile_1 =~ s/#/1/;
	my $fastafile_2 = "$fastafile";
	$fastafile_2 =~ s/#/2/;

	unless (-e $fastafile_2) {
		print "Files $fastafile_1 and $fastafile_2 do not exist.\n";
		return 0;
	}


	my (undef, $seq_names) = tempfile(UNLINK => 1);

	system ("sort $sequencelist | uniq -u > $seq_names");

	open LIST_FH, "<:crlf", "$seq_names";
	open FA1_FH, "<:crlf", "$fastafile_1";
	open FA2_FH, "<:crlf", "$fastafile_2";
	open OUT_FH, ">", "$outfile";

	my $line = readline LIST_FH;
	my $fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
	my $fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
	while (1) {
		# break if any of these files are done.
		if (!(defined $fa_seq1)) { last; }
		if (!(defined $fa_seq2)) { last; }
		if (!(defined $line)) { last; }
		$line =~ /(.*?)\/1/;
		my $curr_name = $1;
		$fa_seq1 =~ />(.*?)\/1/;
		my $name1 = $1;
		$fa_seq2 =~ />(.*?)\/2/;
		my $name2 = $1;
		if ($name1 =~ /$curr_name/) {
			# we've gotten to the corresponding entry in fa1.
			# check the entry in fa2:
			if (($name2 cmp $curr_name) < 0) {
				# is it smaller? (fa2 has more entries than fa1)
				# then we need to move to the next fa2 entry and loop.
				$fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
				next;
			} elsif (($name2 cmp $curr_name) == 0) {
				# if they're equal, print them all out to OUT_FH, advance iterators, and loop.
				print OUT_FH "$fa_seq1$fa_seq2";
				$line = readline LIST_FH;
				$fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
				$fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
				next;
			} else {
				# if $name2 is bigger and we didn't find $curr_name, that means that $name2 doesn't exist.
				# advance the first two and loop.
				$line = readline LIST_FH;
				$fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
				next;
			}
		}
		$fa_seq1 = (readline FA1_FH) . (readline FA1_FH);

	}

	close LIST_FH;
	close FA1_FH;
	close FA2_FH;
	close OUT_FH;
	return 1;
}

sub findsequences {
	my $fastafile = shift;
	my $names = shift;

	unless (-e $fastafile) {
		printlog ("File $fastafile does not exist.");
		die;
	}

	my $hashed_seqs = {};
	my ($taxa, $taxanames) = parsefasta ($fastafile);

	foreach my $name (@$names) {
		if (exists $taxa->{$name}) {
			$hashed_seqs->{$name} = $taxa->{$name};
		}
	}
	return $hashed_seqs;
}

return 1;
