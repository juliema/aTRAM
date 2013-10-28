#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;

my $debug = 0;
if (@ARGV == 0) {
    pod2usage(-verbose => 1);
}

my $runline = "running " . basename($0) . " " . join (" ", @ARGV) . "\n";

#required: input file names
my $short_read_archive = 0;
my $search_fasta = 0;

#optional parameters:
my $help = 0;
my $log_file = 0;
my $use_ends = 0;
my $protein = 0;
my $blast_name = 0;
my $complete = 0;
my $velvet = 0;
my $bitscore = 70;

#parameters with modifiable default values
my $output_file = 0;
my $ins_length = 300;
my $iterations = 5;
my $start_iter = 1;
my $exp_cov = 30;
my $evalue = 10e-10;
my $max_target_seqs = 100000000;

GetOptions ('reads=s' => \$short_read_archive,
            'target=s' => \$search_fasta,
            'insert_length|ins_length=i' => \$ins_length,
            'exp_coverage|expected_coverage=i' => \$exp_cov,
            'iterations=i' => \$iterations,
            'start_iteration=i' => \$start_iter,
            'log_file=s' => \$log_file,
            'use_ends' => \$use_ends,
            'output=s' => \$output_file,
            'protein' => \$protein,
            'blast=s' => \$blast_name,
            'debug' => \$debug,
            'velvet' => \$velvet,
            'bitscore=i' => \$bitscore,
            'max_target_seqs=i' => \$max_target_seqs,
            'evalue=f' => \$evalue,
            'complete' => \$complete,
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive and $search_fasta) {
    pod2usage(-msg => "Must specify a short read archive (that has already been prepared with makelibrary.pl) and a target gene in fasta form.");
}

# check to make sure that the specified short read archive exists:
unless ((-e "$short_read_archive.1.fasta") && (-e "$short_read_archive.2.fasta")) {
    pod2usage(-msg => "Short read archive does not seem to be in the format made by makelibrary.pl. Did you specify the name correctly?");
}

print $runline;

my $executing_path = dirname(__FILE__);
my $cmd;

my $log_fh;
if ($log_file) {
	open $log_fh, ">", $log_file or die "couldn't open $log_file\n";
} else {
	open $log_fh, ">", "$output_file.log" or die "couldn't open $output_file.log\n";
}

print $log_fh $runline;

unless ($output_file) {
    $output_file = $short_read_archive;
}

open CONTIGS_FH, ">", "$output_file.all.fasta";
truncate CONTIGS_FH, 0;
close CONTIGS_FH;

my $start_seq = "";
my $end_seq = "";
my @hit_matrices = ();

# process the target sequence file to look for the start seq and end seq.
my @target_seqs = ();

open SEARCH_FH, "<", $search_fasta;
my $name = "";
my $seq = "";
while (my $line=readline SEARCH_FH) {
	if ($line =~ />(.*)/) {
		if ($name ne "") {
			push @target_seqs, "$name,$seq";
		}
		$name = $1;
		chomp $name;
		$seq = "";
	} else {
		$seq .= $line;
		chomp $seq;
	}
}
push @target_seqs, "$name,$seq";

close SEARCH_FH;

my $len = 100;
if ($protein==1) {
	$len = 30;
}

# check first target seq, see if we need to make just the first bit separated for checking completeness.
my $firstseq = shift @target_seqs;
if ($firstseq =~ /(.*),(.{$len})(.*)/) {
	# split into two smaller seqs:
	$start_seq = "sTRAM_target_start";
	unshift @target_seqs, "$1,$3";
	unshift @target_seqs, "$start_seq,$2";
} else {
	# it's short enough to just go as-is
	$firstseq =~ /(.*),/;
	$start_seq = $1;
	unshift @target_seqs, $firstseq;
}

# check last target seq, see if we need to make just the first bit separated for checking completeness.
my $lastseq = pop @target_seqs;
if ($lastseq =~ /(.*),(.*)(.{$len})/) {
	# split into two smaller seqs:
	$end_seq = "sTRAM_target_end";
	push @target_seqs, "$1,$2";
	push @target_seqs, "$end_seq,$3";
} else {
	# it's short enough to just go as-is
	$lastseq =~ /(.*),/;
	$end_seq = $1;
	push @target_seqs, $lastseq;
}

# okay, let's re-assemble the file for the target fasta.
my ($TARGET_FH, $target_fasta) = tempfile(UNLINK => 1);
my @targets = ();
foreach my $line (@target_seqs) {
	$line =~ /(.*),(.*)/;
	print $TARGET_FH ">$1\n$2\n";
	push @targets, $1;
}

# make a database from the target so that we can compare contigs to the target.
my (undef, $targetdb) = tempfile(UNLINK => 1);
if ($protein == 1) {
	$cmd = "makeblastdb -in $target_fasta -dbtype prot -out $targetdb.db -input_type fasta";
} else {
	$cmd = "makeblastdb -in $target_fasta -dbtype nucl -out $targetdb.db -input_type fasta";
}
system_call ($cmd);

if ($start_iter > 1) {
	my $x = $start_iter-1;
	$search_fasta = "$output_file.$x.contigs.fa";
}

for (my $i=$start_iter; $i<=$iterations; $i++) {
	print ("interation $i starting...\n");
	print $log_fh ("interation $i starting...\n");

	my %hit_matrix = ();
	push @hit_matrices, \%hit_matrix;

	if ($velvet==0) {
	# 1. blast to find any short reads that match the target.
	print "\tblasting short reads...\n";
		if (($protein == 1) && ($i == 1)) {
			system_call ("tblastn -max_target_seqs $max_target_seqs -db $short_read_archive.db -query $search_fasta -outfmt 6 -num_threads 8 -out $output_file.blast.$i");
		} else {
			system_call ("blastn -task blastn -evalue $evalue -max_target_seqs $max_target_seqs -db $short_read_archive.db -query $search_fasta -outfmt 6 -num_threads 8 -out $output_file.blast.$i");
		}

		# did we not find any reads? Go ahead and quit.
		if ((-s "$output_file.blast.$i") == 0) {
			die "No similar reads were found.\n";
		}

		# 2 and 3. find the paired end of all of the blast hits.
		print "\tgetting paired ends...\n";
		system_call("perl $executing_path/lib/pairedsequenceretrieval.pl $short_read_archive.#.fasta $output_file.blast.$i $output_file.blast.$i.fasta");
	}
	# 4. assemble the reads using velvet
	print "\tassembling with Velvet...\n";
	system_call("velveth $output_file.velvet 31 -fasta -shortPaired $output_file.blast.$i.fasta");
	system_call("velvetg $output_file.velvet -ins_length $ins_length -exp_cov $exp_cov -min_contig_lgth 200");

	# 5. now we filter out the contigs to look for just the best ones.
	print "\tfiltering contigs...\n";

	my $blast_file = "";
	unless ($blast_name) {
		(undef, $blast_file) = tempfile(UNLINK => 1);
	} else {
		$blast_file = "$blast_name.$i";
	}

	if ($protein == 1) {
		system_call ("blastx -db $targetdb.db -query $output_file.velvet/contigs.fa -out $blast_file -outfmt '6 qseqid sseqid bitscore qstart qend sstart send'");
	} else {
		system_call ("tblastx -db $targetdb.db -query $output_file.velvet/contigs.fa -out $blast_file -outfmt '6 qseqid sseqid bitscore qstart qend sstart send'");
	}
	# 6. we want to keep the contigs that have a bitscore higher than 100.
	open BLAST_FH, "<", $blast_file;
	while (my $line = readline BLAST_FH) {
		$line =~ /(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*?)\s+(.*)\s/;
		my $contig = $1;
		my $baitseq = $2;
		my $score = $3;
		my $qstart = $4;
		my $qend = $5;
		my $sstart = $6;
		my $send = $7;
		my $strand = (($qend-$qstart) / ($send-$sstart));
		my $currscore = $hit_matrix{$contig}->{$baitseq};
		if ($currscore == undef) {
			$hit_matrix{$contig}->{$baitseq} = $strand * $score;
		} else {
			if ($currscore < $score) {
			$hit_matrix{$contig}->{$baitseq} = $strand * $score;
			}
		}
	}
	close BLAST_FH;

	foreach my $contig (keys %hit_matrix) {
		my $contig_high_score = 0;
		$hit_matrix{$contig}->{"strand"} = 1;
		foreach my $baitseq (@targets) {
			my $score = abs($hit_matrix{$contig}->{$baitseq});
			if ($score > $contig_high_score) {
				$contig_high_score = $score;
				$hit_matrix{$contig}->{"strand"} = ($hit_matrix{$contig}->{$baitseq})/$score;
				$hit_matrix{$contig}->{$baitseq} = $score;
			}
			if ($score < 70) {
				delete $hit_matrix{$contig}->{$baitseq};
			}
		}
		if ($contig_high_score < $bitscore) {
			delete $hit_matrix{$contig};
		}
	}

	my ($SORT_FH, $sort_results) = tempfile(UNLINK => 1);
	foreach my $contig (keys %hit_matrix ) {
		print $SORT_FH $contig . "\n";
	}
	close $SORT_FH;

	# we'll use the resulting contigs as the query for the next iteration.
	$search_fasta = "$output_file.$i.contigs.fa";

	system_call("perl $executing_path/lib/findsequences.pl $output_file.velvet/contigs.fa $sort_results $search_fasta");

	# SHUTDOWN CHECK:
	if (-z $search_fasta) {
		die ("No contigs had a bitscore greater than $bitscore; quitting at iteration $i.");
	}

	# revcomping contigs with negative strand directions:
	my @contigs = ();

	open SEARCH_FH, "<", $search_fasta;
	my $name = "";
	my $seq = "";
	while (my $line=readline SEARCH_FH) {
		if ($line =~ />(.*)/) {
			if ($name ne "") {
				if ($hit_matrix{$name}->{"strand"} < 0) {
					$seq = reverse_complement($seq);
				}
				push @contigs, ">$i"."_$name\n$seq";
			}
			$name = $1;
			chomp $name;
			$seq = "";
		} else {
			$seq .= $line;
			chomp $seq;
		}
	}
	if ($hit_matrix{$name}->{"strand"} < 0) {
		$seq = reverse_complement($seq);
	}
	push @contigs, ">$i"."_$name\n$seq";
	close SEARCH_FH;

	# save off these resulting contigs to the ongoing contigs file.
	open CONTIGS_FH, ">>", "$output_file.all.fasta";
	print CONTIGS_FH join("\n",@contigs) . "\n";
	close CONTIGS_FH;

	# SHUTDOWN CHECK:
	if ($complete == 1) {
		foreach my $contig (keys %hit_matrix) {
			my $contigname = "$i"."_$contig";
			if (($hit_matrix{$contig}->{$start_seq} > 0) && ($hit_matrix{$contig}->{$end_seq} > 0)) {
				die ("Contig $contigname contains both the start and end of the target sequence; quitting at iteration $i.");
			}
		}
	}

	# if we flagged to use just the ends of the contigs, clean that up.
	if ($use_ends != 0) {
		system_call("bash $executing_path/lib/sort_contigs.sh $search_fasta");

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
			print FH ">$2\n$3\n";
		}
		close FH;
	}
}

print "Finished!\n\nContigs by target coverage:\n";

open CONTIGS_FH, ">", "$output_file.results.txt";
print CONTIGS_FH "contig\t" . join ("\t",@targets) . "\n";

my @complete_contigs = ();

for (my $i=0; $i<@hit_matrices; $i++) {
	my $hitmatrix = @hit_matrices[$i];
	foreach my $contig (keys $hitmatrix) {
		my $contigname = "".($i+1)."_$contig";
		print CONTIGS_FH "$contigname\t";
		foreach my $target (@targets) {
			if ($hitmatrix->{$contig}->{$target} == undef) {
				print CONTIGS_FH "-\t";
			} else {
				print CONTIGS_FH "".$hitmatrix->{$contig}->{$target}."\t";
			}
		}
		print CONTIGS_FH "\n";
		if (($hitmatrix->{$contig}->{$start_seq} > 0) && ($hitmatrix->{$contig}->{$end_seq} > 0)) {
			push @complete_contigs, $contigname;
		}
	}
}
close CONTIGS_FH;

system ("cat $output_file.results.txt");
print "\nContigs containing the entire target sequence:\n\t" . join("\n\t", @complete_contigs) . "\n";
print $log_fh "Contig coverage in $output_file.results.txt\n";
print $log_fh "\nContigs containing the entire target sequence:\n\t" . join("\n\t", @complete_contigs) . "\n";

close $log_fh;

sub exit_with_msg {
	my $msg = shift;
	print STDERR "$msg\n";
	exit 1;
}

sub system_call {
	my $cmd = shift;
	print $log_fh ("\t$cmd\n");
	my ($saveout, $saveerr);
	if ($debug == 0) {
		open $saveout, ">&STDOUT";
		open $saveerr, ">&STDERR";
		open STDOUT, '>', File::Spec->devnull();
		open STDERR, '>', File::Spec->devnull();
	}
	my $exit_val = eval {
		system ($cmd);
	};

	if ($debug == 0) {
		open STDOUT, ">&", $saveout;
		open STDERR, ">&", $saveerr;
	}

	if ($exit_val != 0) {
		print "System call \"$cmd\" exited with $exit_val\n";
		exit;
	}
	return $exit_val;
}

sub reverse_complement {
	my $charstr = shift;

	# reverse the DNA sequence
	my $revcomp = reverse($charstr);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
	return $revcomp;
}

__END__

=head1 NAME

sTRAM.pl

=head1 SYNOPSIS

sTRAM.pl -reads shortreadfile -target target.fasta [-ins_length int] [-exp_coverage int] [-iterations int] [-start_iteration int] [-log_file filename] [-use_ends] [-output filename]

=head1 OPTIONS

  -reads:     		short read archive (already run through makelibrary.pl).
  -target:          fasta file with sequences of interest.
  -output:	        optional: the prefix for the pipeline's output files (default name is the same as -reads).
  -ins_length:	    optional: the size of the fragments used in the short-read library (default 300).
  -exp_coverage:    optional: the expected coverage of the region for velvetg (default 30).
  -iterations:      optional: the number of pipeline iterations (default 5).
  -start_iteration: optional: if resuming from previous run, which iteration number to start from (default 0).
  -log_file:        optional: a file to store output of the pipeline.
  -use_ends:        optional: if this flag is present, use the first and last $ins_length of long contigs in the search.
  -protein:         optional: if the target sequence is a protein fasta file.


=head1 DESCRIPTION

Takes a fasta file and finds aligned regions in each sequence in the fasta file that
match the reference sequence(es). Returns a fasta file of aligned regions of similarity.
Uses BLASTN to find regions of similarity.

=cut

