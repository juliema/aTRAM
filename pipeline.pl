use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use IPC::System::Simple qw(run system capture EXIT_ANY);
# use autodie qw(:all);
use File::Temp qw/ tempfile tempdir /;

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

#parameters with modifiable default values
my $output_file = 0;
my $ins_length = 300;
my $iterations = 5;
my $start_iter = 0;
my $exp_cov = 30;

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
            'help|?' => \$help) or pod2usage(-msg => "GetOptions failed.", -exitval => 2);

if ($help) {
    pod2usage(-verbose => 1);
}

unless($short_read_archive and $search_fasta) {
    pod2usage(-msg => "Must specify a short read archive (that has already been prepared with 0-prepare_files.pl) and a target gene in fasta form.");
}

# check to make sure that the specified short read archive exists:
unless ((-e "$short_read_archive.1.fasta") && (-e "$short_read_archive.2.fasta")) {
    pod2usage(-msg => "Short read archive does not seem to be in the format made by 0-prepare_files.pl. Did you specify the name correctly?");
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

my (undef, $targetdb) = tempfile(UNLINK => 1);
my (undef, $blast_file) = tempfile(UNLINK => 1);
my (undef, $sort_file) = tempfile(UNLINK => 1);

# make a database from the target so that we can compare contigs to the target.
if ($protein == 1) {
	$cmd = "makeblastdb -in $search_fasta -dbtype prot -out $targetdb.db";
} else {
	$cmd = "makeblastdb -in $search_fasta -dbtype nucl -out $targetdb.db";
}
print $log_fh ("$cmd\n");
capture(EXIT_ANY, $cmd);


if ($start_iter > 0) {
	$search_fasta = "$output_file.$start_iter.contigs.fa";
}

open CONTIGS_FH, ">>", "$output_file.all.fasta";

for (my $i=$start_iter; $i<$iterations; $i++) {
	print ("interation $i starting...\n");
	print $log_fh("interation $i starting...\n");

	# 1. blast to find any short reads that match the target.
	if (($protein == 1) && ($i == 0)) {
		$cmd = "tblastn -max_target_seqs 100000000 -db $short_read_archive.db -query $search_fasta -outfmt 6 -num_threads 8 -out $output_file.blast.$i";
	} else {
		$cmd = "blastn -task blastn -evalue 10e-10 -max_target_seqs 100000000 -db $short_read_archive.db -query $search_fasta -outfmt 6 -num_threads 8 -out $output_file.blast.$i";
	}
	print $log_fh ("\t$cmd\n");
	capture(EXIT_ANY, $cmd);

	# 2 and 3. find the paired end of all of the blast hits.
	$cmd = "perl $executing_path/2.5-pairedsequenceretrieval.pl $short_read_archive.#.fasta $output_file.blast.$i $output_file.blast.$i.fasta";
	print $log_fh ("\t$cmd\n");
	capture (EXIT_ANY, $cmd);

	# 4. assemble the reads using velvet
	$cmd = "velveth $output_file.velvet 31 -fasta -shortPaired $output_file.blast.$i.fasta";
	print $log_fh ("\t$cmd\n");
	capture (EXIT_ANY, $cmd);

	$cmd = "velvetg $output_file.velvet -ins_length $ins_length -exp_cov $exp_cov -min_contig_lgth 200";
	print $log_fh ("\t$cmd\n");
	capture (EXIT_ANY, $cmd);

	# we'll use the resulting contigs as the query for the next iteration.
	$search_fasta = "$output_file.$i.contigs.fa";
	capture ("mv $output_file.velvet/contigs.fa $search_fasta");

	# 5. now we filter out the contigs to look for just the best ones.
	if ($protein == 1) {
		$cmd = "blastx -db $targetdb.db -query $search_fasta -out $blast_file -outfmt '6 qseqid bitscore'";
	} else {
		$cmd = "tblastx -db $targetdb.db -query $search_fasta -out $blast_file -outfmt '6 qseqid bitscore'";
	}
	print $log_fh ("\t$cmd\n");
	capture (EXIT_ANY, $cmd);

	# 6. we want to keep the contigs that have a bitscore higher than 100.
	$cmd = "gawk '{print \$2\"\\t\"\$1;}' $blast_file | sort -n -r | gawk '{if (\$1 > 100) print \$2;}' | sort | uniq > $sort_file";
	print $log_fh ("\t$cmd\n");
	capture (EXIT_ANY, $cmd);

	$cmd = "perl $executing_path/6.5-findsequences.pl $search_fasta $sort_file $search_fasta";
	print $log_fh ("\t$cmd\n");
	capture (EXIT_ANY, $cmd);

	# save off these resulting contigs to the ongoing contigs file.
	print CONTIGS_FH `cat $search_fasta | gawk '{sub(/>/,">$i\_"); print \$0}'`;

	# if we flagged to use just the ends of the contigs, clean that up.
	if ($use_ends != 0) {

		$cmd = "bash $executing_path/5.5-sort_contigs.sh $search_fasta";
		capture (EXIT_ANY, $cmd);

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

close $log_fh;


__END__

=head1 NAME

pipeline.pl

=head1 SYNOPSIS

pipeline.pl -reads shortreadfile -target target.fasta [-ins_length int] [-exp_coverage int] [-iterations int] [-start_iteration int] [-log_file filename] [-use_ends] [-output filename]

=head1 OPTIONS

  -reads:     		short read archive (already run through 0-prepare_files.pl).
  -target:          fasta file with sequences of interest.
  -output:	        optional: the prefix for the pipeline's output files (default name is the same as -reads).
  -ins_length:	    optional: the size of the fragments used in the short-read library (default 300).
  -exp_coverage:    optional: the expected coverage of the region for velvetg (default 30).
  -iterations:      optional: the number of pipeline iterations (default 5).
  -start_iteration: optional: if resuming from previous run, which iteration number to start from (default 0).
  -log_file:        optional: a file to store output of the pipeline.
  -use_ends:        optional: if this flag is present, use the first and last $ins_length of long contigs in the search.


=head1 DESCRIPTION

Takes a fasta file and finds aligned regions in each sequence in the fasta file that
match the reference sequence(es). Returns a fasta file of aligned regions of similarity.
Uses BLASTN to find regions of similarity.

=cut

