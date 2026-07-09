#!/usr/bin/env perl 
# (Copyright) Jiaqi Wu
#use diagnostics;
use 5.010;
#use strict;
use Cwd;
use Getopt::Std;


my %opts;
getopts('c:p:o:', \%opts);
my $input_cds = $opts{'c'} or die "use: perl $0 -c input_cds_folder -p input_pep_align_folder -o output_fodler\n";
my $input_pep_align = $opts{'p'} or die "use: perl $0 -c input_cds_folder -p input_pep_align_folder -o output_fodler\n";
my $output_folder = $opts{'o'} or die "use: perl $0 -c input_cds_folder -p input_pep_align_folder -o output_fodler\n";



# my $input_cds = "test";
# my $input_pep_align = "test2";
# my $output_folder = "test3";

open TABLE, "Standard.txt" or die "Cannot open Standard.txt: $!\n";
my %table;
while (<TABLE>){
	chomp;
	my @lines = split;
	$table{$lines[0]} = $lines[1];
}
close TABLE;


chdir "$input_cds" or die "Cannot open folder $input_cds: $!\n";
my @folder1 = glob "*.fas";
chdir "../$input_pep_align" or die "Cannot open folder2 $input_pep_align: $!\n";
my @floder2 = glob "*.fas";
chdir "../";
mkdir "./$output_folder" or die "Cannot make folder $output_folder: $!\n";

foreach my $i (@folder1){
	open LENGTH, ">./CDS_align_seq_length.txt" or die "Cannot make file CDS_align_seq_length.txt: $!\n";
	my $seq;
	if ($i =~ /(.*)\.fas/){
		$seq = $1;
	}
	open OUTFILE, ">./$output_folder/$seq.cds.align.fas" or die "Cannot make file $seq.cds.align.fas: $!\n";
	my %seq_ori = &read_fasta("./$input_cds/$seq.fas");
	my %seq_align = &read_fasta("./$input_pep_align/$seq.pep.align.fas");
	my @seq_ori_name = sort keys %seq_ori;
	my @seq_length;
	foreach my $j (@seq_ori_name){
		my @translate_seq = &translate($seq_ori{$j});
		#print ">$j\n";
		#print @translate_seq;
		#print "\n";
		my @align_seq = split (//,$seq_align{$j});
		my @cds_seq = &cut_cds($seq_ori{$j});
		for (my $t = 0; $t <= $#align_seq; $t++){
			#my $codon = $cds_seq[($t+1-$count)*3].$cds_seq[($t+1-$count)*3+1].$cds_seq[($t+1-$count)*3+2]; 
			if ($align_seq[$t] ne $translate_seq[$t]){
				#$count ++;
				splice @translate_seq, $t, 0, "-";
				#my $codon = $cds_seq[($t+1-$count)*3].$cds_seq[($t+1-$count)*3+1].$cds_seq[($t+1-$count)*3+2];
				#print "$codon \n";
				splice @cds_seq, $t*3, 0, "-";
				splice @cds_seq, $t*3+1, 0, "-";
				splice @cds_seq, $t*3+2, 0, "-";
			}
		}
		push @seq_length, scalar (@cds_seq);
		#print "translate_seq, @translate_seq\n";
		#print "align_seq, @align_seq\n";
		print OUTFILE ">$j\n";
		print OUTFILE @cds_seq;
		print OUTFILE "\n";
	}
	print LENGTH "$i\t@seq_length\n";
}
close LENGTH;
close OUTFILE;





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
		$sequences{$id} .= uc $_;
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
	my $pep_seq;
	my $pep_length = int (($#cds_seq+1)/3);
	foreach my $j (0..($pep_length-1)){
		my $codon_pos = $cds_seq[3*$j].$cds_seq[3*$j+1].$cds_seq[3*$j+2];
		$codon_pos = uc $codon_pos;
		if ((exists $table{$codon_pos}) & ($table{$codon_pos} ne "*")){
			$pep_seq .= $table{$codon_pos};
		} 
	}
	my @pep_site = split (//, $pep_seq);
	## Remove stop codon
	if ($pep_site[-1] eq "*"){
		pop (@pep_site);
	}
	@pep_site;
}

sub cut_cds{
	my $input = $_[0];
	#print "$input\n";
	my @cds_seq = split (//,$input);
	my $cds_cut_seq = "";
	my $pep_length = int (($#cds_seq+1)/3);
	foreach my $j (0..($pep_length-1)){
		my $codon_pos = $cds_seq[3*$j].$cds_seq[3*$j+1].$cds_seq[3*$j+2];
		if ((exists $table{$codon_pos}) & ($table{$codon_pos} ne "*")){
			$cds_cut_seq .= $codon_pos;
		} 
	}
	my @cds_site = split (//, $cds_cut_seq);
	@cds_site;
}

