#!/bin/env python
import argparse
from collections import defaultdict
from os import listdir
from os.path import isfile, join
import re

try:
    import cPickle as pickle
except:
    import pickle

## Customized modules/functions
from mylog import info, error, warn, debug

GENE_NAME = "EGFR"

parser = argparse.ArgumentParser(description='Adds columns to the design file for H1 genome hits and H2 genome hits.')
parser.add_argument("--input_file", dest="input_file", help="File containing the design", default="../designs/{}.design.tsv".format(GENE_NAME))
parser.add_argument("--output_file", dest="output_file", help="Where to write the output file", default="../designs/{}.design.hom.tsv".format(GENE_NAME))
parser.add_argument("--hom_hit_dir", dest="hom_hit_dir", help="The output location of the xhyb files", default="{}_results/".format(GENE_NAME))
#parser.add_argument("--xhyb_chromosomes", dest="xhyb_chromosomes", help="The location of the datagram used to create xhyb results", default="script/xhyb_check.datagrams.lst")
parser.add_argument("--hom_hit_mask_threshold", dest="hom_hit_mask_threshold", help="Min mask value to count as a hit (mask: [15,15,10,10,5,5,5,1,1,1,1,1,1,1,1,1,1]) used for hom hits against genome", type=int, default=65)
parser.add_argument("--hit_loc_mask_threshold", dest="hit_loc_mask_threshold", help="Min mask value to count as a hit (mask: [15,15,10,10,5,5,5,1,1,1,1,1,1,1,1,1,1]) used for gap fill checks", type=int, default=65)
parser.add_argument("--output_pickle_h1", dest="output_pickle_h1", help="Where to write the output pickle", default="../designs/{}_xhyb.h1.locations".format(GENE_NAME))
parser.add_argument("--output_pickle_h2", dest="output_pickle_h2", help="Where to write the output pickle", default="../designs/{}_xhyb.h2.locations".format(GENE_NAME))

args = parser.parse_args()

input_file, output_file, hom_hit_dir = args.input_file, args.output_file, args.hom_hit_dir
#xhyb_chromosomes, hom_hit_mask_thresh, hit_loc_mask_thresh = args.xhyb_chromosomes, args.hom_hit_mask_threshold, args.hit_loc_mask_threshold
hom_hit_mask_thresh, hit_loc_mask_thresh = args.hom_hit_mask_threshold, args.hit_loc_mask_threshold
output_pickle_h1, output_pickle_h2 = args.output_pickle_h1, args.output_pickle_h2

info("Inputs: {},{}| Outputs: {}".format(input_file,hom_hit_dir,output_file))

def align_score(alignment):
    H1_alignment_score, H2_alignment_score = 0, 0
    mask = [15,15,10,10,5,5,5,1,1,1,1,1,1,1,1,1,1]
    H1_sequence, H2_sequence = alignment[:17], alignment[::-1][:17]
    for mask_score, align_character in zip(mask, H1_sequence):
        if align_character != ".": H1_alignment_score += mask_score
    for mask_score, align_character in zip(mask, H2_sequence):
        if align_character != ".": H2_alignment_score += mask_score
    return H1_alignment_score, H2_alignment_score

#chromosome_list = [x.rstrip() for x in open(xhyb_chromosomes)]
chromosome_list = [ f.split(".")[1] for f in listdir(hom_hit_dir) if isfile(join(hom_hit_dir,f)) and re.match("fasta_all.*tab",f) ]

mip_information_H1, mip_information_H2 = defaultdict(int), defaultdict(int)
h1_hom_location_data, h2_hom_location_data = defaultdict(dict), defaultdict(dict)

for chromosome in chromosome_list:
    info("Doing {}".format(str(chromosome)))

    for xhyb_line in open("{}/fasta_all.{}.matches.tab".format(hom_hit_dir, chromosome)):
        (chrom, strand, start, stop, length, hom_name, hom_seq, mismatches, length2, alignment) = xhyb_line.rstrip('\n').split("\t")
        H1_alignment_score, H2_alignment_score = align_score(alignment)

        #Hom hits for H1 and H2
        if (H1_alignment_score > hom_hit_mask_thresh): mip_information_H1[hom_name] += 1
        if (H2_alignment_score > hom_hit_mask_thresh): mip_information_H2[hom_name] += 1

        #Store positions for gap fill analysis
        if H1_alignment_score >= hit_loc_mask_thresh:
            if not chromosome in h1_hom_location_data[hom_name]:
                h1_hom_location_data[hom_name][chromosome] = {}
            if not strand in h1_hom_location_data[hom_name][chromosome]:
                h1_hom_location_data[hom_name][chromosome][strand] = []
            h1_hom_location_data[hom_name][chromosome][strand].append(int(start))

        if H2_alignment_score >= hit_loc_mask_thresh:
            if not chromosome in h2_hom_location_data[hom_name]:
                h2_hom_location_data[hom_name][chromosome] = {}

            if not strand in h2_hom_location_data[hom_name][chromosome]:
                h2_hom_location_data[hom_name][chromosome][strand] = []
            h2_hom_location_data[hom_name][chromosome][strand].append(int(stop))

output = open(output_file, "w")
for design_line in open(input_file):
    if design_line.startswith("hom_name"):
        output_line = "{}\tH1_genome_hits\tH2_genome_hits\n".format(design_line.rstrip('\n'))
    else:
        file_parts = design_line.rstrip('\n').split("\t")
        hom_name = file_parts[0]
        output_line = "{}\t{}\t{}\n".format(design_line.rstrip('\n'), mip_information_H1[hom_name], mip_information_H2[hom_name])
    output.write(output_line)
output.close()

pickle.dump(h1_hom_location_data, open(output_pickle_h1, "wb"), -1)
pickle.dump(h1_hom_location_data, open(output_pickle_h2, "wb"), -1)
info("Finished")
