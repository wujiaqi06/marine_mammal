#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;
getopts('i:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i Genome_folder\n";
#my $output_file = $opts{'o'} or die "use: perl $0 -i Genome_folder -o output_file\n";
chdir "./$input_folder";
my @all_seq_file = glob '*.fas';
chdir "../";
open GAP_INFO, ">./$input_folder.gap_info.txt" or die "Cannot make file gap_info.txt: $!\n";
print GAP_INFO "File Name\tSeq Name\tFragment Start\tFragment End\tFragment Length\tGap length before Fragment\tGap length after Fragment\n";
my %gap_info_sum;

print "Step1, summarizing gap and sequence informaton...\n";
foreach my $seq (@all_seq_file){
	my %read_sequences = &read_fasta("./$input_folder/$seq");
	my @seq_name = sort keys %read_sequences;
	my @sequences = values %read_sequences;
	foreach my $i (@seq_name){
		my @seq_site = split(//,$read_sequences{$i});
		#my @seq_site = &translate($read_sequences{$i});
		#print "Seq is @seq_site\n";
		#my $count_seq = 0;
		my %fre_start_end;
		my %fre_length;
		my $fre_start = 1;
		my $gap_count = 0;
		my %fre_gap_before;
		my %fre_gap_after;
		my $fre_end;
		my $count_fre = 0;
		foreach my $j (0..($#seq_site-1)){
			if (($seq_site[$j] eq "-") & ($seq_site[$j+1] ne "-")){
				$fre_gap_after{$fre_start} = $gap_count + 1;
				$fre_start = $j+2;
				$fre_end = $j+1;
				$count_fre ++;
				$fre_gap_before {$fre_start} = $gap_count + 1;
			}  elsif (($seq_site[$j] ne "-") & ($seq_site[$j+1] ne "-") & ($j != ($#seq_site-1))){
				$count_fre ++;
				$fre_end++;
				$gap_count = 0;
			}  elsif (($seq_site[$j] ne "-") & ($seq_site[$j+1] eq "-")){
				$fre_end ++;
				$fre_start_end{$fre_start} = $fre_end;
				$fre_length{$fre_start} = $count_fre ++;
				$count_fre = 0;
			}  elsif (($seq_site[$j] eq "-") & ($seq_site[$j+1] eq "-")){
				$count_fre = 0;
				$gap_count++;
			} 

			if(($seq_site[$j] ne "-") & ($seq_site[$j+1] ne "-") & ($j == ($#seq_site-1))){
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

			if (($seq_site[$j] eq "-") & ($seq_site[$j+1] eq "-") & ($j == ($#seq_site-1))){
				$fre_gap_after{$fre_start} = $#seq_site - $fre_end +1;
			}

			if (($seq_site[$j] ne "-") & ($seq_site[$j+1] eq "-") & ($j == ($#seq_site-1))){
				$fre_gap_after{$fre_start} = 1;
			}

		}
		my $replace_sum = 0;
		foreach my $t (sort keys %fre_start_end){
				print GAP_INFO "$seq\t$i\t$t\t$fre_start_end{$t}\t$fre_length{$t}\t$fre_gap_before{$t}\t$fre_gap_after{$t}\n";
				$gap_info_sum{$seq} .= "$seq\t$i\t$t\t$fre_start_end{$t}\t$fre_length{$t}\t$fre_gap_before{$t}\t$fre_gap_after{$t}\n";
			}
	}
	#print "Gap Info Sum: \n$gap_info_sum{$seq}\n";
}
print "Gap and sequence information are summarized in $input_folder.gap_info.txt\n\n";
close GAP_INFO;




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

sub max{
	if ($_[0] >= $_[1]){
		$_[0];
	} else {
		$_[1];
	}
}

sub min{
	if ($_[0] <= $_[1]){
		$_[0];
	} else {
		$_[1];
	}
}

sub overlap{
	my $start1 = $_[0];
	my $end1 = $_[1];
	my $start2 = $_[2];
	my $end2 = $_[3];
	my $start = &max($start1, $start2);
	my $end = &min($end1, $end2);
	if ($start < $end){
		return ($start,$end);
	} else {
		return (-1,-1);
	}
}

sub p_distance{
	my @seq1 = split ("",$_[0]);
	my @seq2 = split ("",$_[1]);
	my $start = $_[2]-1;
	my $end = $_[3]-1;
	my $identical = 0;
	my $difference = 0;
	foreach my $i ($start..$end){
		if ($seq1[$i] eq $seq2[$i]){
			$identical ++;
		} else {
			$difference ++;
		}
	}
	my $p_distance = $difference /($identical + $difference);
	return ($identical, $difference, $p_distance);
}










