# given two fasta files corresponding to paired ends of short reads, find the specified sequences and then output them into a single output file.
#!/usr/bin/perl
use strict;
use File::Temp qw/ tempfile /;
require Subfunctions;

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

	open LIST_FH, "<", "$seq_names";
	open FA1_FH, "<", "$fastafile_1";
	open FA2_FH, "<", "$fastafile_2";
	open OUT_FH, ">", "$outfile";

	my ($fa_seq1, $fa_seq2, $line) = "";
	$line = readline LIST_FH;
	$fa_seq1 = (readline FA1_FH) . (readline FA1_FH);
	$fa_seq2 = (readline FA2_FH) . (readline FA2_FH);
	while (1) {
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
		# break if any of these files are done.
		if ($fa_seq1 eq "") { last; }
		if ($fa_seq2 eq "") { last; }
		if ($line eq "") { last; }

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
		die "File $fastafile does not exist.\n";
	}

	my $hashed_seqs = {};
	my @sorted_names = sort @$names;

	my (undef, $fasta_sort) = tempfile(UNLINK => 1);
	sortfasta ($fastafile, $fasta_sort, "#");

	open FA_FH, "<", "$fasta_sort";
	my $curr_name = shift @sorted_names;
	my $fa_seq = readline FA_FH;
	while (1) {
		$fa_seq =~ />(.*?)#(.*)/;
		my $name = $1;
		my $seq = $2;
		if ($name =~ /$curr_name/) {
			# fa_seq is the seq we're looking for.
			$hashed_seqs->{$curr_name} = $seq;
			$curr_name = shift @sorted_names;
		} else {
			# we need to go on to the next seq.
			$fa_seq = readline FA_FH;
		}
		# break if any of these files are done.
		if ($fa_seq eq "") { last; }
		if ($curr_name eq undef) { last; }
	}

	close FA_FH;

	return $hashed_seqs;
}

return 1;
