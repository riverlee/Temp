#!/bin/env python
import sys
import math
import argparse
import os.path
import pandas as pd
import numpy as np
from collections import defaultdict
try:
    import cPickle as pickle
except:
    import pickle

## Customized modules/functions
from mylog import info, error, warn, debug
from local_repeatMasker_finder import local_repeatMasker_id

parser = argparse.ArgumentParser(description='Create MIPs from a list of homs.')

GENE_NAME = "FTH1"

#Input files
parser.add_argument("--gene_name", dest="gene_name", help="The name of the gene we are designing", default="{}".format(GENE_NAME))
parser.add_argument("--input_bed", dest="input_bed", help="The regions to design for - bed file format", default="{}.Exons.bed".format(GENE_NAME))
parser.add_argument("--input_design", dest="input_design", help="input design file of homs", default="{}.design.hom.tsv".format(GENE_NAME))
#parser.add_argument("--repeat_masked_file", dest="repeat_masked_file", help="Repeat masked file of the gene", default="rm_files/{}.RM.seq".format(GENE_NAME))
parser.add_argument("--homs_to_omit", dest="homs_to_omit", help="Hom probes to avoid picking (results of self-blast)", default="{}.blast.bad_homs".format(GENE_NAME))
parser.add_argument("--potential_pair_pickle", dest="potential_pair_pickle", help="Serialized potential pairs file if you give a new file name it creates the data", default="{}.pair_file".format(GENE_NAME))
parser.add_argument("--h1_genome_hits", dest="h1_genome_hits", help="pickle file containing hom genome locations for h1", default="{}_xhyb.h1.locations".format(GENE_NAME))
parser.add_argument("--h2_genome_hits", dest="h2_genome_hits", help="pickle file containing hom genome locations for h2", default="{}_xhyb.h2.locations".format(GENE_NAME))

#output file
parser.add_argument("--hom_disqualification_log", dest="hom_disqualification_log", help="location to log hom disualifications", default="{}.hom.disqualification.log".format(GENE_NAME))
parser.add_argument("--mip_disqualification_log", dest="mip_disqualification_log", help="location to log mip disualifications", default="{}.snp.disqualification.log".format(GENE_NAME))
parser.add_argument("--output_bed", dest="output_bed", help="Create a bed file of the design", default="{}.SNPs.bed".format(GENE_NAME))
parser.add_argument("--output_design", dest="output_design", help="The design file", default="{}.MIP.design".format(GENE_NAME))
parser.add_argument("--H1_H2_fasta_database", dest="H1_H2_fasta_database", help="creates fasta database for self-blast", default="{}.blast_db.fa".format(GENE_NAME))
parser.add_argument("--H1_H2_fasta_query", dest="H1_H2_fasta_query", help="creates fasta query for self-blast", default="{}.blast_query.fa".format(GENE_NAME))

#Default options
parser.add_argument("--ideal_gap_length", dest="ideal_gap_length", help="the ideal size of a gap to fill (default 100)", type=int, default=100)
parser.add_argument("--gap_length_min", dest="gap_length_min", help="the minimum acceptable gap fill length (default 90)", type=int, default=90)
parser.add_argument("--gap_length_max", dest="gap_length_max", help="the maximum acceptable gap fill length (default 110)", type=int, default=110)
parser.add_argument("--overlap_step", dest="overlap_step", help="The min size to offset mips from each other (default 50)", type=int, default=50)
parser.add_argument("--cover_exon_upstream", dest="cover_exon_upstream", help="How many bases upstream of an exon to cover (default 65)", type=int, default=65)
#GC
parser.add_argument("--gc_min", dest="gc_min", help="the minimum acceptable gc percent (default 25)", type=float, default=.25)
parser.add_argument("--gc_max", dest="gc_max", help="the maximum acceptable gc percent (default 75)", type=float, default=.75)
#General
parser.add_argument("--hp_run_max", dest="hp_run_max", help="The maximum amount of homopolymer runs to allow (default 0)", type=int, default=0)
parser.add_argument("--hom_hit_max", dest="hom_hit_max", help="The max times a hom can hit the genome (default 1500 - based on mask score 65)", type=int, default=1500)
#DUST / RM thresholds
parser.add_argument("--dust_max", dest="dust_max", help="the max dust percent allowed (default .5)", type=float, default=.5)
parser.add_argument("--rm_max", dest="rm_max", help="the max repeat percent allowed (default .5)", type=float, default=.5)
#Off-target hits
parser.add_argument("--off_target_gap_threshold", dest="off_target_gap_threshold", help="off target gep threshold length max", type=int, default=250)
#Whether allows SNPs in homs (0, no snps in homs, 1 or others numbers allow snps in homs)
parser.add_argument("--allow_snps_in_homs", dest="allow_snps_in_homs", help="Whether allows SNPs in homs", type=int, default=0)
#What's the mininal distance to allow SNPs to the end of upstream and start of downstream
parser.add_argument("--allow_snps_in_homs_distance", dest="allow_snps_in_homs_distance", help="What's the mininal distance to allow snps in homs to the end of upsteam and start of downstream",type=int,default=10)
parser.add_argument("--allow_snps_in_homs_maxNum", dest="allow_snps_in_homs_maxNum", help="Maxinum allowed number of SNPs in homs", type=int, default=5)

args = parser.parse_args()

input_bed, input_design, potential_pair_pickle = args.input_bed, args.input_design, args.potential_pair_pickle
output_bed, H1_H2_fasta_database, H1_H2_fasta_query, output_design = args.output_bed, args.H1_H2_fasta_database, args.H1_H2_fasta_query, args.output_design
hom_disqualification_log, mip_disqualification_log = args.hom_disqualification_log, args.mip_disqualification_log
ideal_gap_length, gl_min, gl_max, overlap_step = args.ideal_gap_length, args.gap_length_min, args.gap_length_max, args.overlap_step
cover_exon_upstream = args.cover_exon_upstream
gc_min, gc_max, dust_max, rm_max = args.gc_min, args.gc_max, args.dust_max, args.rm_max
hp_run_max, hom_hit_max = args.hp_run_max, args.hom_hit_max
h1_genome_hit_file, h2_genome_hit_file, off_target_gap_threshold = args.h1_genome_hits, args.h2_genome_hits, args.off_target_gap_threshold
homs_to_omit = args.homs_to_omit

GENE_NAME=args.gene_name
allow_snps_in_homs = args.allow_snps_in_homs
allow_snps_in_homs_distance = args.allow_snps_in_homs_distance
allow_snps_in_homs_maxNum = args.allow_snps_in_homs_maxNum

info("Inputs:{},{},{},{} | Outputs: {},{},{},{},{},{} ".format(input_bed,input_design,h1_genome_hit_file,h2_genome_hit_file,hom_disqualification_log,mip_disqualification_log,output_bed,output_design,H1_H2_fasta_database,H1_H2_fasta_query))


#QC CHECKS
#TODO: check to make sure genome hits file exist
h1_genome_hits, h2_genome_hits = pickle.load(open(h1_genome_hit_file, 'rb')), pickle.load(open(h2_genome_hit_file, 'rb'))

#FUNCTIONS
def positions_in_range(position, threshold, position_list, strand="+"):
    match_count = 0
    for test_position in position_list:
        if strand == "+":
            if test_position >= position and test_position <= position + threshold:
                match_count += 1
        else:
            if test_position <= position and test_position >= position - threshold:
                match_count += 1
    return match_count

def gap_fills(hom_data_h1, hom_data_h2, max_fill):
    gap_fill_matches = 0
    for chromosome in (set(hom_data_h1) & set(hom_data_h2)):
        if '+' in hom_data_h1[chromosome] and '+' in hom_data_h2[chromosome]:
            h1_positions_pos = hom_data_h1[chromosome]["+"]
            h2_positions_pos = hom_data_h2[chromosome]["+"]
            for pos_position in h2_positions_pos:
                gap_fill_matches += positions_in_range(pos_position, max_fill, h1_positions_pos, "+")

        if '-' in hom_data_h1[chromosome] and '-' in hom_data_h2[chromosome]:
            h1_positions_neg = hom_data_h1[chromosome]["-"]
            h2_positions_neg = hom_data_h2[chromosome]["-"]
            for neg_position in h2_positions_neg:
                gap_fill_matches += positions_in_range(neg_position, max_fill, h1_positions_neg, "-")
    return gap_fill_matches

def fetch_positions(region_index, hom_name):
    hom_index, hom_position = hom_name.split("-")
    if int(hom_index) == region_index:
        hom_start = int(design_source.hom_start[(design_source.region_index == region_index) & (design_source.hom_name == hom_name)].values[0])
        hom_stop = int(design_source.hom_stop[(design_source.region_index == region_index) & (design_source.hom_name == hom_name)].values[0])
        return hom_start, hom_stop
    else:
        return -1, -1
def fetch_snps(region_index,hom_name,start,stop,flag):
    hom_index, hom_position = hom_name.split("-")
    if int(hom_index) == region_index:
        snps = design_source.SNPs[ (design_source.region_index == region_index) & (design_source.hom_name == hom_name)].values[0]
        minD=sys.maxsize
        if pd.isnull(snps):
            return minD,0
        else:
            snpsPos = [ int(s.split(":")[1]) for s in snps.split("|") ]
            # get minimal distance
            for pos in snpsPos:
                d=1000
                # flag 1 means downstream
                # flag 0 means upstream
                if flag:
                    d= pos - start
                else:
                    d= stop - pos
                
                if d<minD: minD=d
            
            return minD, len(snpsPos)
    else:
        return minD,0

def fetch_sequence(hom_name):
    sequence = design_source.seq[design_source.hom_name == hom_name].values[0]
    return sequence

def revcomp(dna, reverse=True):
    bases = 'ATGCTACG'
    complement_dict = {bases[i]:bases[i+4] for i in range(4)}
    if reverse:
        dna = reversed(dna)
    result = [complement_dict[base] for base in dna]
    return ''.join(result)

if os.path.isfile(potential_pair_pickle):
    make_pairs = False
else:
    make_pairs = True

#load the design
info("Loading design file - {}".format(input_design))
design_source = pd.read_table(input_design)
chromosome = int(design_source.chrom.head(1).values)


## Old version, replace by Jiang at line 159
"""
# Use repeatMasker program
if os.path.isfile(rm_file):
    for seq_line in open(rm_file):
        if seq_line.startswith(">"):
            (chrom, start, stop) = seq_line[1:].rstrip().split("-")
            current_index = int(start)
        else:
            for position, base in enumerate(seq_line.rstrip()):
                if base == "N": rm_data.append(int(position) + current_index)
            current_index += position
print (rm_data)
"""

#Homs to censor
if os.path.isfile(homs_to_omit):
    info("Loading hom omit file - {}".format(homs_to_omit))
    homs_to_omit = [x.rstrip() for x in open(homs_to_omit)]


#Get the region set
info("Loading exons region - {}".format(input_bed))
region_start = sys.maxsize
region_end = -sys.maxsize
region_list = []
for region_line in open(input_bed):
    if region_line.startswith("track"): continue
    region_parts = region_line.rstrip().split("\t")
    start, stop = int(region_parts[1]), int(region_parts[2])
    region_list.append((start, stop,))
    if start<region_start:
        region_start=start
    
    if stop>region_end:
        region_end = stop

#Repeat Masked Positions
info("Loading repeatMaskers")
rm_finder = local_repeatMasker_id("chr{}".format(str(chromosome)))
rm_data = rm_finder.local_repeatMaskers(region_start-cover_exon_upstream,region_end+cover_exon_upstream)

#Storage of the design
Final_MIP_Set = defaultdict(list)

#Generate all pairs valid for a given region
if make_pairs:
    hom_fail_log = open(hom_disqualification_log, "w")
    hom_fail_log.write("""hom_name\tH1_H2_both\tseq\tstart\tstop\tgc_pct\trm_pct\thomopolymer_run\tdust_pct_h1\tdust_pct_h2\tgenome_hits_h1\tgenome_hits_h2\tinvalidGC\tinvalidRM\tinvalidDustH1\tinvalidDustH2\tinvalidHomHitH1\tinvalidHomHitH2\thasSNPs\thasSMs\tSelf_Blast_Fail\n""")

    info("Make possible hom pairs and do hom filter | do filter")
    #TODO: print out header information about the settings
    #Make a list of all the acceptable homs - this reduces the amount of searching for MIP evaluation
    valid_H1_indexes, valid_H2_indexes, hom_count = [], [], 0
    for probe in design_source.iterrows():
        hom_count += 1
        hom_index, hom_name, hom_seq, start, stop = probe[0], probe[1]["hom_name"], probe[1]["seq"], int(probe[1]["hom_start"]), int(probe[1]["hom_stop"])
        hom_gc, hom_dust_H1, hom_dust_H2, hom_hit_H1, hom_hit_H2 = probe[1]["gc_pct"], probe[1]["dust_pct_H1"], probe[1]["dust_pct_H2"], probe[1]["H1_genome_hits"], probe[1]["H2_genome_hits"],
        hom_SNPS, hom_SMS, hom_hp = probe[1]["SNPs"], probe[1]["SMs"], probe[1]["hp_run"]

        gcFail, H1dustFail, H2dustFail, H1hitFail, H2hitFail, hpFail, rmFail, SNPFail, SMFail, blastFail = False, False, False, False, False, False, False, False, False, False

        if hom_gc < gc_min or hom_gc > gc_max: gcFail = True
        if hom_hp > hp_run_max: hpFail = True
        if hom_dust_H1 > dust_max: H1dustFail = True
        if hom_dust_H2 > dust_max: H2dustFail = True
        if hom_hit_H1 > hom_hit_max: H1hitFail = True
        if hom_hit_H2 > hom_hit_max: H2hitFail = True
        if not pd.isnull(hom_SNPS): SNPFail = True
        if not pd.isnull(hom_SMS): SMFail = True

        #Repeat mask
        rm_count = 0
        #for masked_base in rm_data:
        #    if masked_base > start and masked_base < stop:
        #        print (masked_base, start, stop)
        #        rm_count += 1
        
        for i in range(start,stop):
            if i in rm_data:
                rm_count+=1

        rm_pct = rm_count / float(len(hom_seq))
        if rm_pct >= rm_max:
            rmFail = True

        #If bad self blast
        if hom_name in homs_to_omit:
            blastFail = True

        validH1, validH2 = True, True
        if allow_snps_in_homs:
            if rmFail or blastFail or gcFail or hpFail or H1dustFail or H1hitFail or SMFail: validH1 = False
            if rmFail or blastFail or gcFail or hpFail or H2dustFail or H2hitFail or SMFail: validH2 = False
        else:
            if SNPFail or rmFail or blastFail or gcFail or hpFail or H1dustFail or H1hitFail or SMFail: validH1 = False
            if SNPFail or rmFail or blastFail or gcFail or hpFail or H2dustFail or H2hitFail or SMFail: validH2 = False

        if validH1: valid_H1_indexes.append(hom_index)
        if validH2: valid_H2_indexes.append(hom_index)

        failType = ""
        if not validH1 and not validH2: failType = "Both"
        elif not validH1: failType = "H1"
        elif not validH2: failType = "H2"
        if failType != "":
            hom_fail_log.write("""{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}\t{20}\n""".format(hom_name,failType,hom_seq,start,stop,hom_gc,rm_pct,hom_hp,
                hom_dust_H1,hom_dust_H2,hom_hit_H1,hom_hit_H2,gcFail,rmFail,H1dustFail,H2dustFail,H1hitFail,H2hitFail,SNPFail,SMFail,blastFail))
    hom_fail_log.close()
    print ("Filtered homs!\nStarting hom count: {}\nH1_probes_remaining: {}\tH2 probes remaining: {}\n".format(hom_count, len(valid_H1_indexes), len(valid_H2_indexes)))

    info("Make possible hom pairs and do hom filter | make pairs")
    #Now make the pairs
    mip_fail_log = open(mip_disqualification_log, "w")
    #mip_fail_log.write("h2_hom_name\th1_hom_name\toff_target_hits\n")
    mip_fail_log.write("h2_hom_name\th1_hom_name\toff_target_hits\n")
    for r_index, (start, stop) in enumerate(region_list):
        Final_MIP_Set[str(r_index)] = []
        start, stop, = start - cover_exon_upstream, stop + cover_exon_upstream
        for H2_hom in design_source[(design_source.region_index == r_index) & (design_source.index.isin(valid_H2_indexes))].iterrows():
            H2_hom_index, H2_hom_name, H2_start, H2_stop = H2_hom[0], H2_hom[1]["hom_name"], H2_hom[1]["hom_start"], H2_hom[1]["hom_stop"]

            for H1_hom in design_source[(design_source.hom_start >= H2_stop + gl_min) & (design_source.hom_start <= H2_stop + gl_max) & (design_source.region_index == r_index) & (design_source.index.isin(valid_H1_indexes))].iterrows():
                H1_hom_index, H1_hom_name, H1_start, H1_stop = H1_hom[0], H1_hom[1]["hom_name"], H1_hom[1]["hom_start"], H1_hom[1]["hom_stop"]

                gap_length = H1_start - H2_stop
                gap_fill_count = gap_fills(h1_genome_hits[H1_hom_name], h2_genome_hits[H2_hom_name], off_target_gap_threshold)

                if gap_fill_count > 1:
                    mip_fail_log.write("{}\t{}\t{}\n".format(H2_hom_name, H1_hom_name, gap_fill_count))
                else:
                    Final_MIP_Set[str(r_index)].append((H2_hom_name, H1_hom_name,))

    pickle.dump(Final_MIP_Set, open(potential_pair_pickle, "wb"), -1)
    mip_fail_log.close()
else:
    info("Loading homs pairs")
    Final_MIP_Set = pickle.load(open(potential_pair_pickle, "rb"))

print("Created pairs!\nRegion\tPosition\tSize\tAvailable Pairs\n")
for r_index, (start, stop) in enumerate(region_list):
    print ("{}\t{}-{}\t{}\t{}".format(r_index, start, stop, stop-start, len(Final_MIP_Set[str(r_index)])))
print ("\n")

#Select the MIPs
info("Make final homs pairs to make MIP")
exon_pairs = defaultdict(list)
for region_index, (start, stop) in enumerate(region_list):
    valid_mips = Final_MIP_Set[str(region_index)]
    last_H2_end, last_H1_start = start - cover_exon_upstream, stop + cover_exon_upstream
    last_H1_end = last_H1_start
    for H2_mip, H1_mip in valid_mips:
        H2_start, H2_stop = fetch_positions(region_index, H2_mip)
        H1_start, H1_stop = fetch_positions(region_index, H1_mip)
        
        upstream_snp_minD,upstream_snp_num = fetch_snps(region_index,H2_mip,H2_start,H2_stop,0)
        downstream_snp_minD,downstream_snp_num = fetch_snps(region_index,H1_mip,H1_start,H2_stop,1)

        #print (H2_start, H2_stop, H1_start, H1_stop)
        if (H2_start >= last_H2_end + overlap_step) and (H2_stop < last_H1_start):
            if allow_snps_in_homs and upstream_snp_minD > allow_snps_in_homs_distance and downstream_snp_minD > allow_snps_in_homs_distance and upstream_snp_num < allow_snps_in_homs_maxNum and downstream_snp_num <allow_snps_in_homs_maxNum:
                exon_pairs[str(region_index)].append((H2_mip, H1_mip,))
                last_H2_end, last_H1_start, last_H1_end = H2_stop, H1_start, H1_stop
        elif (H2_start > last_H1_end):
            #In case there is a hole in the coverage
            if allow_snps_in_homs and upstream_snp_minD > allow_snps_in_homs_distance and downstream_snp_minD > allow_snps_in_homs_distance and upstream_snp_num < allow_snps_in_homs_maxNum and downstream_snp_num <allow_snps_in_homs_maxNum:
                exon_pairs[str(region_index)].append((H2_mip, H1_mip,))
                last_H2_end, last_H1_start, last_H1_end = H2_stop, H1_start, H1_stop

#Write the MIP
out_bed = open(output_bed, "w")
out_bed.write("track name={}\n".format(output_bed))

#blast database file
out_blast_db = open(H1_H2_fasta_database, "w")
#blast query file
out_blast_query = open(H1_H2_fasta_query, "w")


MIP_backbone = "TTCCAACCTTCGATCTGTGCUUUGCCGCTCCGAGAACTTG"
smMIP_tag = "NNNNNNNNNNNN"
Five_Phos = "/5Phos/"

out_design = open(output_design, "w")
#out_design.write("Mip_Name\tH2_positions\tgap_fill_positions\tH1_positions\tH2_seq\tH1_seq\tDesign_H2_seq\tDesign_H1_seq\tDesign_MIP\n")
out_design.write("Mip_Name\tupstream\tgap_fill_positions\tdownstream\tH2_positions\tgap_fill_positions\tH1_positions\tstrand\tupstream_seq\tdownstream_seq\tDesign_H2_seq\tDesign_H1_seq\tDesign_MIP\n")

info("Output result")
for region_index, (start, stop) in enumerate(region_list):
    out_bed.write("chr{}\t{}\t{}\t{}\n".format(chromosome, start, stop, region_index+1))
    exon_mips = exon_pairs[str(region_index)]
    for MIP_index, (H2_hom, H1_hom) in enumerate(exon_mips):
        H2_start, H2_stop = fetch_positions(region_index, H2_hom)
        H1_start, H1_stop = fetch_positions(region_index, H1_hom)
        H2_seq = fetch_sequence(H2_hom)
        H1_seq = fetch_sequence(H1_hom)

        upstream_start=H2_start
        upstream_stop=H2_stop
        downstream_start=H1_start
        downstream_stop=H1_stop

        out_bed.write("chr{}\t{}\t{}\t{}_{}-{}\n".format(chromosome, H2_start, H1_stop, GENE_NAME,region_index+1,MIP_index+1))

        gap_fill_start = max(H2_stop,H2_start) + 1
        gap_fill_stop = min(H1_start,H1_stop) - 1

        strand="+"
        if MIP_index % 2 == 0:
            design = "{}{}{}{}{}".format(Five_Phos, H1_seq, MIP_backbone, smMIP_tag, H2_seq)
            design_H2, design_H1 = H2_seq, H1_seq
        else:
            design = "{}{}{}{}{}".format(Five_Phos,revcomp(H2_seq),MIP_backbone, smMIP_tag, revcomp(H1_seq))
            design_H2, design_H1 = revcomp(H1_seq), revcomp(H2_seq)
            tmp1, tmp2 = H2_start, H2_stop
            H2_start, H2_stop = H1_start, H1_stop
            H1_start, H1_stop = tmp1, tmp2
            strand="-"

        #Self-blast analysis output
        out_blast_query.write(">{}.H1\n{}\n>{}.H2\n{}\n".format(H1_hom, H1_seq, H2_hom, H2_seq))
        for n_index in range(0,10):
            out_blast_db.write(">{}-{}\n{}{}{}\n".format(H2_hom, H1_hom, H2_seq, 'N' * n_index, H1_seq))

        out_design.write("{}_{}-{}\t{}-{}\t{}-{}\t{}-{}\t{}-{}\t{}-{}\t{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(GENE_NAME,region_index+1, MIP_index + 1,
            upstream_start,upstream_stop,gap_fill_start,gap_fill_stop,downstream_start,downstream_stop,
            H2_start, H2_stop, gap_fill_start, gap_fill_stop,
            H1_start, H1_stop,strand, H2_seq, H1_seq, design_H2, design_H1, design))

out_bed.close()
out_blast_query.close()
out_blast_db.close()
out_design.close()
info("Finished")
