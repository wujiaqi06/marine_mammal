perl RemoveEmptySequences.pl -i hg38_mammals_species_name -o hg38_noEmpty0.9 -m 0.9

perl RemoveGap_cds.pl -i hg38_noEmpty0.9 -o hg38_noEmpty0.9_cut0.7 -d 0.7 -t fas

perl RemoveFragmentSeq_cds.pl -i hg38_noEmpty0.9_cut0.7 -o hg38_noEmpty0.9_cut0.7f10g5 -f 10 -g 5

perl RemoveEmptySequences.pl -i hg38_noEmpty0.9_cut0.7f10g5 -o hg38_noEmpty0.9_cut0.7f10g5_noEmpty0.9 -m 0.9


"hg38_species_name_cleaned" is the same with "hg38_noEmpty0.9_cut0.7f10g5_noEmpty0.9"