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
getopts('i:o:d:t:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder -o output_folder -d deleting_rate -t fasta\/phy\n";
my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_folder -d deleting_rate -t fasta\/phy\n";
my $del_missing_value = $opts{'d'} or die "use: perl $0 -i Genome_folder -o output_folder -d deleting_rate -t fasta\/phy\n";
my $output_type = $opts{'t'} or die "use: perl $0 -i Genome_folder -o output_folder -d deleting_rate -t fasta\/phy\n";

# my $input_folder = "test3";
# my $output_folder = "test4";
# my $del_missing_value = 0.3;
# my $output_type = "fas";

mkdir "$output_folder" or die "Can not create folder $output_folder:$!\n";
open DELIND, ">./$input_folder.Delete_info.txt";
chdir "./$input_folder";
my @all_seq_file = glob '*.fas';
chdir "../";

## get %table, the codon table.
open TABLE, "Standard.txt" or die "Cannot open Standard.txt: $!\n";
my %table;
while (<TABLE>){
	chomp;
	my @lines = split;
	$table{$lines[0]} = $lines[1];
}


foreach my $seq (@all_seq_file){
	my %read_sequences = &read_fasta("./$input_folder/$seq");
	my @sequences = values %read_sequences;
	my @sequence_name = sort keys %read_sequences;
	my $sequence_number = $#sequence_name +1;
	my $sequence_length = length($sequences[0]);
	my $missing_value = 0 x int ($sequence_length/3);
	my @missing = split(//, $missing_value);
	## Get Pep Missing Info
	foreach my $i (@sequences) {
		my $temp_missing = "";
		my @seq_site = &translate($i);
		#print "Seq Length ".@seq_site."\t";
		#say @seq_site;
		foreach my $j (@seq_site){
			if (($j eq "-") | ($j eq "*")){
				$temp_missing .= 0;
			} else {
				$temp_missing .= 1;
				}
		}
		my @temp_array = split(//, $temp_missing);
		foreach my $j (0..$#missing){
			$missing[$j] = $missing[$j] + $temp_array[$j];
		}
		#print "temp_missing is @temp_array\n";
		#my $test1 = scalar (@temp_array);
		#my $test2 = scalar (@missing);
		#print "$seq test1 $test1, test2 $test2\n";
	}
	#print "missing_value is @missing\n";


	my $cut_number = int (($#sequence_name+1) * $del_missing_value);
	#print "$cut_number\n";
	my @del_index;
	foreach my $i (0..$#missing){
		if ($missing[$i] < $cut_number){
			push @del_index, $i;
		}
	}
	
	my $delete_length = ($#del_index + 1) * 3;
	my $remaining = $sequence_length - $delete_length;
	print "Delete $delete_length sites, remaining $remaining sites\n\n";
	print DELIND "$seq: original length: $sequence_length bp, deleted $delete_length bp, remaining $remaining bp.\n";

	## Remove Gap
	@del_index = reverse (@del_index);
	if ($output_type =~ /fas(ta)?/){
		open OUTPUT, ">./$output_folder/$seq\n";
		foreach my $seq_name (@sequence_name){
			my @site = split(//, $read_sequences{$seq_name});
			#my @site_del = splice @site, @del_index;
			foreach my $del (@del_index){
				splice(@site, $del*3, 3);
			}
			#print "length of site is $#site\n";
			#my $new_seq = join(/:/,@site);
			print OUTPUT ">$seq_name\n";
			print OUTPUT @site;
			print OUTPUT "\n";
		}
		close OUTPUT;
	} elsif($output_type eq "phy"){
		open OUTPUT, ">./$output_folder/$seq.phy\n";
		print OUTPUT "$sequence_number $remaining\n";
		foreach my $seq_name (@sequence_name){
			my @site = split(//, $read_sequences{$seq_name});
			#my @site_del = splice @site, @del_index;
			foreach my $del (@del_index){
				splice(@site, $del*3, 3);
			}
			#print "length of site is $#site\n";
			#my $new_seq = join(/:/,@site);
			print OUTPUT "$seq_name\n";
			print OUTPUT @site;
			print OUTPUT "\n";
		}
			close OUTPUT;
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

sub translate{
	my $input = $_[0];
	#print "$input\n";
	my @cds_seq = split (//,$input);
	#print "cds length". @cds_seq . "\n";
	my $pep_seq;
	my $pep_length = int (($#cds_seq+1)/3);
	foreach my $j (0..($pep_length-1)){
		my $codon_pos = $cds_seq[3*$j].$cds_seq[3*$j+1].$cds_seq[3*$j+2];
		if (exists $table{$codon_pos}){
			$pep_seq .= $table{$codon_pos};
		} else {
			$pep_seq .= "-";
			#print "Codon is $codon_pos\n";
		}
	}
	my @pep_site = split (//, $pep_seq);
	# if ($pep_site[-1] eq "*"){
	# 	pop (@pep_site);
	# }
	@pep_site;
}


