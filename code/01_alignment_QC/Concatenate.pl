#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;

my %opts;

getopts('i:o:f:', \%opts);
my $input_folder = $opts{'i'} or die "use: $0 -i input_folder -o output_file -f Output_format(fas\/phy)\n";
my $output_file = $opts{'o'} or die "use: $0 -i input_folder -o output_file -f Output_format(fas\/phy)\n";
my $Output_format = $opts{'f'} or die "use: $0 -i input_folder -o output_file -f Output_format(fas\/phy)\n";


print "\nThis software is by Jiaqi Wu in The University of Tokyo(2015.12.6).\nRunning...\n\n";
chdir "$input_folder";
my @all_seq_file = glob '*.fas*';
print "@all_seq_file\n";


my %concatenate_sequences = &read_fasta($all_seq_file[0]);
my $first_file = $all_seq_file[0];
shift (@all_seq_file);
my $length;


open ERROR, ">../ERROR.txt" or die "Could not make new :$!\n";
open ANNO_INFO, ">../Annotation_Info.txt" or die "Could not make new :$!\n";

## Write the annotation of first file in Annotation_Info.txt
my @seq_fisrt_keys = sort keys %concatenate_sequences;
my $random_first_seq = $seq_fisrt_keys[rand @seq_fisrt_keys];
my $random_first_seq_length = length ($concatenate_sequences{$random_first_seq});
print ANNO_INFO "$first_file = 1-$random_first_seq_length;\n";

foreach my $seq_file (@all_seq_file){
	my %read_seq = &read_fasta($seq_file);
	my @seq_keys = sort keys %read_seq;
	my $random_seq = $seq_keys[rand @seq_keys];
	my $random_seq_length = length ($read_seq{$random_seq});
#	print "Randomly Picked up:$random_seq with length $random_seq_length.\n";

#	get all sequences;
	my @all_keys = &uniq(keys %concatenate_sequences, @seq_keys);
#	print "Total Sequences are ($#all_keys+1)\n";
	
	foreach my $seq (@all_keys){
		##### If sequences is OK in both concatenated file and new file, concatenate!
		##### Check alignment is correct or not;
		if((exists $read_seq{$seq})&&(exists $concatenate_sequences{$seq})){
			my $add_length = length ($read_seq{$seq});
		#my $length0 = length ($concatenate_sequences{$seq});
			if (($add_length == $random_seq_length)&&(exists $read_seq{$seq})&&(exists $concatenate_sequences{$seq})){
				$concatenate_sequences{$seq} .= $read_seq{$seq};
				$length = length ($concatenate_sequences{$seq});
			}
			#### Alignemnt check. If sequence length is different among sequences in new file, give wornning and exit.
			if ($add_length != $random_seq_length){
				print "Alignment Error in reading $seq_file, $seq: Sequence lengths are not equal! Please check your alignment!\n";
				print ERROR "Alignment Error in reading $seq_file, $seq: Sequence lengths are not equal! Please check your alignment!\n";
				print "Programm Quit due to fatal error!\n";
				print ERROR "Programm Quit due to fatal error!\n";
			exit();
			}
		}
		
		#### If the sequence is missing in new file, give wornning and replace this sequence with "-";
		if ((!(exists $read_seq{$seq}))){
			print "Error! $seq_file missing sequence: $seq!\n";
			print "Please check sequences! Continue...\n";
			print ERROR "File $seq_file missing sequence: $seq!\n";
			$concatenate_sequences{$seq} .= "-"x$random_seq_length;
			#exit();
		}
		#### If sequence are not find in concatenated file, give warnning and exit.
		if (!(exists $concatenate_sequences{$seq})){
			print "Unknown Sequences! $seq_file has one sequences:$seq missing in other files!\nPlease check your sequence file!\n";
			print ERROR "Unknown Sequences! $seq_file has one sequences:$seq missing in other files!\nPlease check your sequence file!\n";			
			print "Programm Quit due to fatal error!\n";
			print ERROR "Programm Quit due to fatal error!\n";
			exit();
		}
	}
			print "Current Concatenated Sequence length is $length bp.\n\n\n";
			my $start = $length-$random_seq_length+1;
			print ANNO_INFO "$seq_file = $start-$length;\n";
}

##Write Concatenated sequences

my $number_seq = keys %concatenate_sequences;
	
if ($Output_format eq "phy"){
	open OUTFILE, ">../$output_file\.phy" or die "Can not write file:$output_file*$!\n";
	print OUTFILE "$number_seq $length\n";
	foreach my $seq (sort keys %concatenate_sequences){
		print OUTFILE "$seq     $concatenate_sequences{$seq}\n";
	}
}

		if($Output_format eq "fas"){
			open OUTFILE, ">../$output_file\.fas" or die "Can not write file:$output_file*$!\n";
			foreach my $seq (sort keys %concatenate_sequences){
				print OUTFILE ">$seq\n$concatenate_sequences{$seq}\n";
			}
		}
print "Finished! Including $number_seq sequences with length $length bp.\n\n";

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
	print SEQ_INFO "Successfuly read @_, in total $count sequences with $seq_length bp.\n";
	close SEQ_INFO;
	print " with $seq_length bp.\n";
	%sequences;
}

sub uniq{
	my %seen;
	return grep {!$seen{$_}++} @_;
}

## Second version of unif subroutine.
#sub uniq{
#	my @unique;
#	my %seen;
#	foreach my $value (@_){
#		if (!$seen{$value}){
#			push @unique, $value;
#			$seen{$value} = 1;
#		}
#	}
#	@unique;
#}



