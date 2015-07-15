#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email:  jiang_li@affymetrix.com
# Date:   Mon Nov 17 13:25:14 2014
###################################
use strict;
use warnings;

my $gene = $ARGV[0];

print <<HERE;
## 1: make bed file
Bed_maker.py --gene_name ${gene} --output_file ${gene}.Exons.bed

## 2. make possible homs
hom_designer.py --input_bed ${gene}.Exons.bed --output_file ${gene}.homs.tsv --output_fasta ${gene}.homs.fa --output_bed ${gene}.homs.bed

## 3. Run xhyb to search hits of homs genome wide
## This should be run on cluster1 (ssh cluster1)
## The go to the ${gene}_pbs folder and loop the .sh files
# for i in `ls *.sh`; do echo \$i; qsub \$i; done 
CreatePBS_xhyb.pl $gene

## 4. Add homs
hom_hit.py --input_file ${gene}.homs.tsv --output_file ${gene}.design.homs.tsv --hom_hit_dir ${gene}_results/ --output_pickle_h1 ${gene}.xhyb.h1.locations --output_pickle_h2 ${gene}.xhyb.h2.locations 

## 5. Select MIPs
select_MIPs.py --gene_name $gene --input_bed ${gene}.Exons.bed --input_design ${gene}.design.homs.tsv --potential_pair_pickle ${gene}.pair_file --h1_genome_hits ${gene}.xhyb.h1.locations --h2_genome_hits ${gene}.xhyb.h2.locations --hom_disqualification_log ${gene}.hom.disualifications.log --mip_disqualification_log ${gene}.mip.disualifications.log --output_bed ${gene}.mip.bed --output_design ${gene}.MIP.design.tsv --H1_H2_fasta_database ${gene}.blast_db.fa --H1_H2_fasta_query ${gene}.blast_query.fa --allow_snps_in_homs 1 --allow_snps_in_homs_distance 5 --allow_snps_in_homs_maxNum 100 --force 2>\&1

HERE
