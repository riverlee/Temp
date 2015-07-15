#!/bin/env python

"""
Take a bed file of all the exons from a given gene and then find all the possible homlogous 
which cover the exons regions

"""

import argparse

## Customized modules/functions
from local_snp_finder import local_snp_id
from nib_frag import nibFragger
from mip_tm_clip import optimize_seq
from dust import score_dust
from mylog import info, error,warn, debug


GENE_NAME = "MYC"
parser = argparse.ArgumentParser(description='Identifies all potential hom sequences within a given region.')

parser.add_argument("--input_bed", dest="input_bed", help="The regions to design for - bed file format", default="{}.Exons.bed".format(GENE_NAME))
parser.add_argument("--output_file", dest="output_file", help="File to store the results in", default="{}.homs.tsv".format(GENE_NAME))
parser.add_argument("--output_fasta", dest="output_fasta", help="File containing fasta sequences of the homs", default="{}.homs.fa".format(GENE_NAME))
parser.add_argument("--output_bed", dest="output_bed", help="File containing regions of the homs in bed format, used to check whether all the possible homs cover all the exons", default="{}.homs.bed".format(GENE_NAME))

#Default options
parser.add_argument('--tm_min', dest="tm_min", help="MIPs with Tms lower than this are ignored", type=float, default=55)
parser.add_argument('--tm_max', dest="tm_max", help="MIPs with Tms higher than this are ignored", type=float, default=62)
parser.add_argument('--min_hom_length', dest="min_hom_length", help="Min hom length (default 17)", type=int, default=17)
parser.add_argument('--max_hom_length', dest="max_hom_length", help="Max hom length (default 28)", type=int, default=28)
parser.add_argument('--gc_threshold_min', dest="gc_threshold_min", help="Low GC content threshold (default .15)", type=float, default=.15)
parser.add_argument('--gc_threshold_max', dest="gc_threshold_max", help="High GC content threshold (default .85)", type=float, default=.85)
parser.add_argument("--mip_offset", dest="mip_offset", help="how far back from the region start to begin MIP design (default 65)", type=int, default=65)

args = parser.parse_args()

# Assign to variables
input_bed, output_file, output_fasta = args.input_bed, args.output_file, args.output_fasta
mip_offset, tm_min, tm_max = args.mip_offset, args.tm_min, args.tm_max
gc_threshold_min, gc_threshold_max = args.gc_threshold_min, args.gc_threshold_max
min_hom_length, max_hom_length = args.min_hom_length, args.max_hom_length
output_bed = args.output_bed


info("Input: {} | Outputs: {}, {}, {}".format(input_bed,output_file,output_fasta,output_bed))


#Parse the regions
target_regions = []
for file_line in open(input_bed):
    if file_line.startswith("track") or file_line.startswith("#"): continue
    file_parts = file_line.rstrip().split("\t")
    chrom, start, stop = file_parts[0], int(file_parts[1]), int(file_parts[2])
    target_regions.append((chrom, start, stop, ))

snp_finder = local_snp_id(chrom)

output = open(output_file, "w")
output.write("hom_name\tregion_index\tchrom\tregion_start\tregion_stop\thom_start\thom_stop\tseq\tseq_tm\tgc_count\tgc_pct\tdust_scoreH1\tdust_scoreH2\tdust_pct_H1\tdust_pct_H2\thp_run\tSNPs\tSMs\n")

out_fasta = open(output_fasta, "w")
out_bed = open(output_bed, "w")

info("There are total of {} exons".format(str(len(target_regions))))

for region_index, (chrom, start, stop) in enumerate(target_regions):
    #For every position in the index, design a hom - extend mip to the right
    i = region_index+1
    bp = str(abs(stop-start+1))
    info("Finding homs on exon {} ({}-{}| {}bp)".format(str(i),str(start),str(stop),bp))

    for hom_position in range(start - mip_offset - 20, stop + mip_offset):
        hom_seq = nibFragger(chrom.replace("chr",""), hom_position, 35)
        opt_seq, opt_tm = optimize_seq(hom_seq)
        gc_count = opt_seq.count("C") + opt_seq.count("G")
        gc_content = gc_count / float(len(opt_seq))
        (hp_run, dust_score_H1, dust_score_H2, dust_pct_H1, dust_pct_H2) = score_dust(opt_seq)

        #Disqualify homs based on thresholds
        if len(opt_seq) > max_hom_length or len(opt_seq) < min_hom_length: continue
        if opt_tm < tm_min or opt_tm > tm_max: continue
        if gc_content > gc_threshold_max or gc_content < gc_threshold_min: continue
        SNPs = snp_finder.local_snps(hom_position, hom_position + len(opt_seq) - 1)
        SMs = snp_finder.local_sms(hom_position, hom_position + len(opt_seq) - 1)
        end_position = hom_position + len(opt_seq) - 1

        output.write("{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}\t{}\t{:.2f}\t{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\t{}\n".format(region_index, hom_position, region_index, chrom, start, stop, hom_position, end_position, opt_seq, opt_tm, \
                gc_count, gc_content, dust_score_H1, dust_score_H2, dust_pct_H1, dust_pct_H2, hp_run, SNPs, SMs))
        out_fasta.write(">{}-{}\n{}\n".format(region_index, hom_position, opt_seq))
        out_bed.write("{}\t{}\t{}\t{}-{}\n".format("chr"+chrom,hom_position,end_position,region_index,hom_position))

output.close()
out_fasta.close()
out_bed.close()
info("Finished")

