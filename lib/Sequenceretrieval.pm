#!/usr/bin/env perl
package Sequenceretrieval;
use strict;
use File::Temp qw/ tempfile /;
use System;
use Parsing;

### github

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
	print "Line 29 $child_pid child pid fork pair retrieval my files are  fastafile $fastafile , squencelist $sequencelist , outfile $outfile \n\n";
		exit;
    } else { #parent process
        return $child_pid;
    }
}

sub pairedsequenceretrieval {
    my $read1=0;
    my $read2=0;
	  
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
	print "Line 65 sequence names are $seq_names\n\n";

	open LIST_FH, "<:crlf", "$seq_names";
	open FA1_FH, "<:crlf", "$fastafile_1";
	open FA2_FH, "<:crlf", "$fastafile_2";
	open OUT_FH, ">", "$outfile";
	print "Line 71 fastafile1 $fastafile_1 fastafle2 $fastafile_2 outfile $outfile\n";

	my $line = readline LIST_FH;
	print "line 73  line $line  this is reading the line of LIST_FH which is the $seq_names\n";
	my $fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
	my $fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
	print "Line 75 $fa_seq1  seq1  $fa_seq2 seq2 \n";
	while (1) {
#	    $read1++;
		# break if any of these files are done.
		if (!(defined $fa_seq1)) { last; }
		if (!(defined $fa_seq2)) { last; }
		if (!(defined $line)) { last; }
		$line =~ /(.*?)\/1/;
		my $curr_name = $1;
#		print "Line 84 current name $curr_name\n";
		$fa_seq1 =~ />(.*?)\/1/;
		my $name1 = $1;
		$fa_seq2 =~ />(.*?)\/2/;
		my $name2 = $1;
#	    print "$name1 $name2 \n";  
		if ($name1 =~ /$curr_name/) {
#		    print " we've gotten to the corresponding entry in fa1.\n";
#		    print " check the entry in fa2.\n";
			if ($name2 ne $curr_name) {
#			    print "$name2 cmp $curr_name\n";
#			    print "# is it smaller? (fa2 has more entries than fa1)\n";
#			    print "# then we need to move to the next fa2 entry and loop.\n";
				$fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
				next;
			} 
			elsif ($name2 eq  $curr_name) {
                            $read2++;
			    print " if they're equal, print them all out to OUT_FH, advance iterators, and loop.\n";
			    print OUT_FH "$fa_seq1$fa_seq2";
			    print "what is printing to the output file\n $fa_seq1$fa_seq2\n";
			    print "$name2    $curr_name\n";
			    $line = readline LIST_FH;
			    $fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
			    $fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
			    next;
			}
		}
			    
#			} else {
#				$fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
#				next;



#############################################################################################################################
#			if (($name2 cmp $curr_name) < 0) {
#			    print "$name2 cmp $curr_name\n";
#			    print "# is it smaller? (fa2 has more entries than fa1)\n";
#			    print "# then we need to move to the next fa2 entry and loop.\n";
#				$fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
#				next;
#			} elsif (($name2 cmp $curr_name) == 0) {
#			    print " if they're equal, print them all out to OUT_FH, advance iterators, and loop.\n";
#				print OUT_FH "$fa_seq1$fa_seq2";
#				print "what is printing to the output file\n $fa_seq1$fa_seq2\n";
#				$line = readline LIST_FH;
#				$fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
#				$fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
#				next;
##			} else {
##				$fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
##				next;
##			}
#			} else {
#				# if $name2 is bigger and we didn't find $curr_name, that means that $name2 doesn't exist.
#				# advance the first two and loop.
#				$line = readline LIST_FH;
#				$fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
#				next;
#			}
#		}
		$fa_seq1 = (readline FA1_FH) . (readline FA1_FH);

	}
    print "read1 total $read1  read2 total $read2\n";
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
