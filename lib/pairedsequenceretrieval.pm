# given two fasta files corresponding to paired ends of short reads, find the specified sequences and then output them into a single output file.
#!/usr/bin/perl
use strict;
use File::Temp qw/ tempfile /;

sub fork_pair_retrieval {
	my $fastafile = shift;
	my $sequencelist = shift;
	my $outfile = shift;

    my $child_pid = fork();
    unless ($child_pid) { #child process
		return pairedsequenceretrieval ($fastafile, $sequencelist, $outfile);
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
	return 1;
}

return 1;
