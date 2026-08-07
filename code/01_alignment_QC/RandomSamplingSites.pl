#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:s:', \%opts);
my $input_file = $opts{'i'} or die "use: perl $0 -i input_seq -o output_seq -s sample_length\n";
my $output_file = $opts{'o'} or die "use: perl $0 -i input_seq -o output_seq -s sample_length\n";
my $sample_length = $opts{'s'} or die "use: perl $0 -i input_seq -o output_seq -s sample_length\n";

my %table = (
	"AAA" => "K", 
	"AAC" => "N", 
	"AAG" => "K", 
	"AAT" => "N", 
	"ACA" => "T", 
	"ACC" => "T", 
	"ACG" => "T", 
	"ACT" => "T", 
	"AGA" => "R", 
	"AGC" => "S", 
	"AGG" => "R", 
	"AGT" => "S", 
	"ATA" => "I", 
	"ATC" => "I", 
	"ATG" => "M", 
	"ATT" => "I", 
	"CAA" => "Q", 
	"CAC" => "H", 
	"CAG" => "Q", 
	"CAT" => "H", 
	"CCA" => "P", 
	"CCC" => "P", 
	"CCG" => "P", 
	"CCT" => "P", 
	"CGA" => "R", 
	"CGC" => "R", 
	"CGG" => "R", 
	"CGT" => "R", 
	"CTA" => "L", 
	"CTC" => "L", 
	"CTG" => "L", 
	"CTT" => "L", 
	"GAA" => "E", 
	"GAC" => "D", 
	"GAG" => "E", 
	"GAT" => "D", 
	"GCA" => "A", 
	"GCC" => "A", 
	"GCG" => "A", 
	"GCT" => "A", 
	"GGA" => "G", 
	"GGC" => "G", 
	"GGG" => "G", 
	"GGT" => "G", 
	"GTA" => "V", 
	"GTC" => "V", 
	"GTG" => "V", 
	"GTT" => "V", 
	"TAA" => "*", 
	"TAC" => "Y", 
	"TAG" => "*", 
	"TAT" => "Y", 
	"TCA" => "S", 
	"TCC" => "S", 
	"TCG" => "S", 
	"TCT" => "S", 
	"TGA" => "*", 
	"TGC" => "C", 
	"TGG" => "W", 
	"TGT" => "C", 
	"TTA" => "L", 
	"TTC" => "F", 
	"TTG" => "L", 
	"TTT" => "F"
	);


my @sites;
my %seq = &read_fasta("$input_file");
my @seq_name = sort keys %seq;
my $length = length($seq{$seq_name[0]})/3;
my @number;
for my $i (1..$sample_length){
	my $random = int(rand ($length));
	push @number, $random;
}

@number = &uniq(@number);
print "sampled $#number sites\n";

open SAM, ">./$output_file" or die "Cannot create Sampled_seq.fas: $!\n";
my %rand_seq;
foreach my $i (@seq_name){
	my @sites = split("", $seq{$i});
	my @sample_sites;
	foreach my $j (@number){
		my $codon = $sites[3*$j].$sites[3*$j+1].$sites[3*$j+2];
		push @sample_sites, $codon;
	}
	#print "Site Length is $#sample_sites\n";
	print SAM ">$i\n";
	print SAM @sample_sites;
	print SAM "\n";
}

close SAM;


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
