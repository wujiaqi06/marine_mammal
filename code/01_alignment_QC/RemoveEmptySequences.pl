#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
# Usage example: perl RemoveGap.pl -i 43mammal_pro_align -o DelAmbSites_43mammal_pro_align -d 0.3 -t phy
# Usage example2: perl RemoveGap.pl -i 43mammal_pro_align -o DelAmbSites_43mammal_pro_align_phy -d 0.3 -t fas
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:m:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder -o output_folder -m maximum_missing_value\n";
my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_folder -m maximum_missing_value\n";
my $maximum_missing_value = $opts{'m'} or die "use: perl $0 -i Genome_folder -o output_folder -m maximum_missing_value\n";

mkdir "$output_folder" or die "Can not create folder $output_folder:$!\n";
open RMI, ">./$output_folder.rmInfo.txt" or die "Cannot create $output_folder.rmInfo.txt: $!\n";
open RMC, ">./$output_folder.rmCount.txt" or die "Cannot create $output_folder.rmCount.txt: $!\n";
open RMS, ">./$output_folder.keptNumber.txt" or die "Cannot create $output_folder.keptNumber.txt $!\n";
print RMI "gene\tremoved_seq\n";
print RMC "gene\tNumber_of_rm_seq\n";
print RMS "gene\tSeq_Number\n";

chdir "./$input_folder";
my @all_seq_file = glob '*.fas';
chdir "../";

foreach my $seq (@all_seq_file){
	my %read_sequences = &read_fasta("./$input_folder/$seq");
	my @sequences = values %read_sequences;
	my @sequence_name = sort keys %read_sequences;
	my $sequence_number = $#sequence_name +1;
	my $sequence_length = length($sequences[0]);
	my $rm_ount = 0;
	my $kepted_number = $sequence_number;
	foreach my $i (@sequence_name){
		my $seq0 = $read_sequences{$i};
		#print "$seq0\n";
		my $miss = &missing_value_gene($seq0);
		if ($miss >= $maximum_missing_value){
			delete $read_sequences{$i};
			print "Removed One Sequence: $i, $miss\n";
			print RMI "$seq\t$i\n";
			$rm_ount ++;
			$kepted_number --;
		}
		if ($miss == 1){
			print "Find an empyt seq: $i, $miss\n";
		}
		
	}
	print RMC "$seq\t$rm_ount\n";
	print RMS "$seq\t$kepted_number\n";
	
	open OUT, ">./$output_folder/$seq";
	foreach my $i (sort keys %read_sequences){
		print OUT ">$i\n$read_sequences{$i}\n";
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

sub missing_value_gene{
	my $seq = $_[0];
	#print "$seq\n";
	my @sites = split ("", $seq);
	my $gene_length = length ($seq);
	my $missing = 0;
	foreach my $i (@sites){
		if ($i eq "-"){
			$missing ++;
		}
	}
	my $missing_rate = $missing/$gene_length;
	return $missing_rate;
}
