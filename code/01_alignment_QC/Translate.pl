#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder -o output_folder\n";
my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_folder\n";

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


chdir "./$input_folder" or die "Cannot enter cds folder: $!\n";
my @cds = glob ("*.fas");
chdir "../";
#print "@cds\n";
mkdir "./$output_folder" or die "Cannot make folder $output_folder: $!\n";

foreach my $seq (@cds){
	my %cds = &read_fasta("./$input_folder/$seq");
	my $seq0 = $seq;
	$seq0 =~ s/(\w+)\..*/$1/;
	#$seq =~ s/\.cds\.align\.fas\.fas//g;
	open PEP, ">./$output_folder/$seq0.pep.fas" or die "Can not make file $seq.pep.fas\n";
	my @cds_name = sort keys %cds;
	foreach my $i (@cds_name){
		my @cds_seq = split(//,$cds{$i});
		#my $codon_pos = $cds_seq[3].$cds_seq[4].$cds_seq[5];
		#print "$codon_pos\n";
		my $pep_length = int (($#cds_seq + 1)/3);
		my $pep_seq;
		foreach my $j (0..($pep_length-1)){
			my $codon_pos = $cds_seq[3*$j].$cds_seq[3*$j+1].$cds_seq[3*$j+2];
			$codon_pos = uc $codon_pos;
			if (exists $table{$codon_pos}){
				$pep_seq .= $table{$codon_pos};
			} elsif ($codon_pos eq "---") {
				#print "$codon_pos\n";
				$pep_seq .= "-";
			}	else {
				$pep_seq .= "-"
			}
		}
		my @pep_site = split (//, $pep_seq);
		if ($pep_site[-1] eq "*"){
			pop (@pep_site);
		}
		print PEP ">$i\n";
		print PEP @pep_site;
		print PEP "\n";
	}
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
	s/\r//;
	if ($_ =~ /^>(?<seq_name>.+$)/){
		$id  = $+{seq_name};
		$count += 1;
	}else{
		$sequences{$id} .= $_;
		$seq_length = length ($sequences{$id});
		}
	}
	print "Successfuly read @_, in total $count sequences with $seq_length bp.\n";
	print SEQ_INFO "Successfuly read @_, in total $count sequences with $seq_length bp.\n";
	close SEQ_INFO;
	#print "\n";
	%sequences;
}


