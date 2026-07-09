#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
# 2018.4.26
# Usage example: perl RemoveGap.pl -i 43mammal_pro_align -o DelAmbSites_43mammal_pro_align -d 0.3 -t phy
# Usage example2: perl RemoveGap.pl -i 43mammal_pro_align -o DelAmbSites_43mammal_pro_align_phy -d 0.3 -t fas
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:o:f:g:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder -o output_folder -f fragment_size -g gap_length\n";
my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_folder -f fragment_size -g gap_length\n";
my $fragment_size = $opts{'f'} or die "use: perl $0 -i Genome_folder -o output_folder -f fragment_size -g gap_length\n";
my $gap_length = $opts{'g'} or die "use: perl $0 -i Genome_folder -o output_folder -f fragment_size -g gap_length\n";
#my $output_type = $opts{'t'} or die "use: perl $0 -i Genome_folder -o output_folder -f fragment_size -g gap_length -t fasta\/phy\n";

## make outpup folder
open DELIND, ">./Delete_info.txt" or die "Cannot create file Delete_info.txt :$!\n";
open GAP_INFO, ">./gap_info.txt" or die "Cannot create file gap_info.txt: $!\n";
open EMPTY, ">./empty_sequence_info.txt" or die "Cannot create file ${input_folder}_empty_sequence_info.txt: $!\n";
mkdir "$output_folder" or die "Can not create folder $output_folder:$!\n";
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
	my @seq_name = sort keys %read_sequences;
	my @sequences = values %read_sequences;
	foreach my $i (@seq_name){
		####remove empty sequences
		if ($read_sequences{$i} =~ /^\-+$/){
			print EMPTY "$seq\t$i\n";
			print "$seq\t$i\tEmapy\n";
			delete $read_sequences{$i};
			next;
		}
		####remove N in sequences
		$read_sequences{$i} =~ s/N/\-/ig;
		my @protein = &translate($read_sequences{$i});
		#print "Seq is @protein\n";
		#my $count_seq = 0;
		my %fre_start_end;
		my %fre_length;
		my $fre_start = 1;
		my $gap_count = 0;
		my %fre_gap_before;
		my %fre_gap_after;
		my $fre_end;
		my $count_fre = 0;
		foreach my $j (0..($#protein-1)){
			if (($protein[$j] eq "-") & ($protein[$j+1] ne "-")){
				$fre_gap_after{$fre_start} = $gap_count + 1;
				$fre_start = $j+2;
				$fre_end = $j+1;
				$count_fre ++;
				$fre_gap_before {$fre_start} = $gap_count + 1;
			}  elsif (($protein[$j] ne "-") & ($protein[$j+1] ne "-") & ($j != ($#protein-1))){
				$count_fre ++;
				$fre_end++;
				$gap_count = 0;
			}  elsif (($protein[$j] ne "-") & ($protein[$j+1] eq "-")){
				$fre_end ++;
				$fre_start_end{$fre_start} = $fre_end;
				$fre_length{$fre_start} = $count_fre ++;
				$count_fre = 0;
			}  elsif (($protein[$j] eq "-") & ($protein[$j+1] eq "-")){
				$count_fre = 0;
				$gap_count++;
			} 

			if(($protein[$j] ne "-") & ($protein[$j+1] ne "-") & ($j == ($#protein-1))){
				$count_fre ++;
				$fre_end ++;
				$fre_start_end{$fre_start} = $fre_end + 1;
				$fre_length{$fre_start} = $count_fre ++;
				$count_fre = 0;
				$fre_gap_after{$fre_start} = 0;
			}
			
			if ($fre_start == 1){
				$fre_gap_before{$fre_start} = 0;
			}

			if (($protein[$j] eq "-") & ($protein[$j+1] eq "-") & ($j == ($#protein-1))){
				$fre_gap_after{$fre_start} = $#protein - $fre_end +1;
				#print "fre_end is $fre_end\n";
			}

			if (($protein[$j] ne "-") & ($protein[$j+1] eq "-") & ($j == ($#protein-1))){
				$fre_gap_after{$fre_start} = 1;
			}

		}
		my $replace_sum = 0;
		foreach my $t (sort keys %fre_start_end){
				#print "$seq, $i, start: $t, end $fre_start_end{$t}, length $fre_length{$t} ,gap_before $fre_gap_before{$t}, gep_after $fre_gap_after{$t}\n";
				#print "t is $t, fre_start_end is $fre_start_end{$t}, fre_gap_before{t} is $fre_gap_before{$t}, fre_gap_after{t} is $fre_gap_after{$t}\n";
				if ($t==1){
					$fre_length{$t} = $fre_length{$t} + 1;
				}
				if (($fre_gap_before{$t} >= $gap_length) | ($fre_gap_after{$t} >= $gap_length)){
					if ($fre_length{$t} <= $fragment_size){
						print GAP_INFO "$seq, $i, start: $t, end $fre_start_end{$t}, length $fre_length{$t} ,gap_before $fre_gap_before{$t}, gep_after $fre_gap_after{$t}\n";
						my $replace = "---" x $fre_length{$t};
						my $replace_start = $t*3-3;
						my $replace_length = 3 * $fre_length{$t};
						$replace_sum += $fre_length{$t};
						#print"$replace\n";
						substr ($read_sequences{$i}, $replace_start, $replace_length, $replace);
					}
				}
			}
		if ($replace_sum != 0){
				print DELIND "$seq, sequence $i, removed $replace_sum animo acid sites\n";
			}
	}
	my $name = $seq;
	$name =~ s/(\w+)\..*/$1/;
	#print "value is $_\n";
	open OUT, ">./$output_folder/$name.fas" or die "Cannot create file $output_folder/$name.fas: $!\n";
	foreach my $i (sort keys %read_sequences){
		if (!($read_sequences{$i} =~ /^\-+$/)){
			print OUT ">$i\n$read_sequences{$i}\n";
		} else {
			print "Extra empty sequence: $seq\t$i\n";
			print EMPTY "$seq\t$i\n";
		}
	}
	close OUT;
}
close DELIND;

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

