#!/usr/bin/env perl -w
# (Copyright) Jiaqi Wu
use diagnostics;
use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use File::Copy qw(move);

my %opts;
getopts('i:', \%opts);
my $input_folder = $opts{'i'} or die "use: perl $0 -i input_sequence_folder\n";
#my $type = $opts{'t'} or die "use: $0 -i Genome_folder -t cdna\/pep -o output_folder\n";
#my $output_folder = $opts{'o'} or die "use: perl $0 -i Genome_folder -t cdna\/pep -o output_folder\n";

#my $input_folder = "seq";
my $output_folder = "BaseMLResult";
mkdir "./$output_folder" or die "Can not make folder $output_folder:$!\n";
mkdir "./rst_sum" or die "Can not make folder rst:$!\n";
my @sequence_file = glob "./$input_folder/*.phy";

foreach my $seq (@sequence_file){
	my %sequ = &read_pamlnuc($seq);
	open SPEC, ">./species.sub.txt" or die $!;
	my $species_names = join("\n", sort keys %sequ);
	print SPEC "$species_names\n";
	close SPEC;

	system "Rscript SubTree.R";

	open CTL, ">>./baseml.ctl" or die "Can not create baseml.ctl file: $!\n";
	## Get the sequence name;
	#my $sequence_root;
	my $sequence_name;
	if ($seq =~ /\/(\w+).phy/){
		#$sequence_root = $1;
		$sequence_name = $1;
		#print "$2\n"
	}
	##generate the baseml.ctl file
	open CTL0, "_baseml.ctl" or die "Can not open _baseml.ctl file: $!\n";
	while (<CTL0>){
		#my $test = $_;
		s/sequence_file/$seq/;
		s/output_file/\.\/$output_folder\/$sequence_name.out/;
		my $control_file = $_;
		print CTL "$control_file";
	}
	#print "\n\n\n";
	#close CTL;

	##run BaseML
	system "baseml";
	rename "rst", "$sequence_name.rst";
	move "$sequence_name.rst","./rst_sum/$sequence_name.rst";

	system "rm baseml.ctl";
}





sub read_pamlnuc{
	open READ_SEQ, "@_" or die "Can not read Sequence file: $!\n";
	#open SEQ_INFO, ">>Sequence_Info.txt" or die "Can not read Sequence file: $!\n";
	my %sequences;
	my $id;
	my $count = 0;
	my $seq_length;
	my $line = 0;
	while (<READ_SEQ>){
		$line ++;
		chomp;
		s/\r//;
		if ($line == 1){
		}elsif (($line % 2) == 0){
			$id  = $_;
			$count += 1;
		}elsif (($line % 2) == 1){
			$sequences{$id} .= $_;
			$seq_length = length ($sequences{$id});
		}
	}
	print "Successfuly read @_, in total $count sequences.\n";
	#print SEQ_INFO "Successfuly read @_, in total $count sequences.\n";
	#close SEQ_INFO;
	print "\n";
	%sequences;
}
