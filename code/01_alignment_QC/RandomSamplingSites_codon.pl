#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:s:n:', \%opts);
my $input_file = $opts{'i'} or die "use: perl $0 -i input_seq -o output_seq -s sample_length -n sample_number\n";
my $output_file = $opts{'o'} or die "use: perl $0 -i input_seq -o output_seq -s sample_length -n sample_number\n";
my $sample_length = $opts{'s'} or die "use: perl $0 -i input_seq -o output_seq -s sample_length -n sample_number\n";
my $simple_number = $opts{'n'} or die "use: perl $0 -i input_seq -o output_seq -s sample_length -n sample_number\n";

my @sites;
my %seq = &read_fasta("$input_file");
my @seq_name = sort keys %seq;
my $length = length($seq{$seq_name[0]})/3;
my $output_folder = "${output_file}_sample";
if (!(-e $output_folder)){
	mkdir "$output_folder" or die "Cannot create $output_folder: $!\n";
} else {
	print "$output_folder exists!\n";
}

for my $n (1..$simple_number){
	my @number;
	for my $i (1..$sample_length){
		my $random = int(rand ($length));
		push @number, $random;
	}

	@number = &uniq(@number);
	print "sample $n sampled $#number sites\n";

	open SAM, ">./$output_folder/${output_file}_$n.fas" or die "Cannot create ${output_file}_$n.fas: $!\n";
	my %rand_seq;
	foreach my $i (@seq_name){
		my @cds_seq = split("", $seq{$i});
		my @sample_sites;
		foreach my $j (@number){
			my $codon = $cds_seq[3*$j].$cds_seq[3*$j+1].$cds_seq[3*$j+2];
			push @sample_sites, $codon;
		}
		#print "Site Length is $#sample_sites\n";
		print SAM ">$i\n";
		print SAM @sample_sites;
		print SAM "\n";
	}
	close SAM;
}


sub read_fasta{
	open READ_SEQ, "@_" or die "Can not read Sequence file: $!\n";
	open SEQ_INFO, ">>../Sequence_Info.txt" or die "Can not read Sequence file: $!\n";
	my %sequences;
	my $id;
	my $count = 0;
	my $seq_length;
	while (<READ_SEQ>){
	chomp;
	if ($_ =~ /^>(?<seq_name>.+$)/){
		$id  = $+{seq_name};
		$count += 1;
	}else{
		$sequences{$id} .= $_;
		$seq_length = length ($sequences{$id});
		}
	}
	print "Successfuly read @_, in total $count sequences";
	print SEQ_INFO "Successfuly read @_, in total $count sequences with $seq_length bp(aa).\n";
	close SEQ_INFO;
	print " with $seq_length bp.\n";
	%sequences;
}


sub uniq{
	my %seen;
	return grep {!$seen{$_}++} @_;
}

