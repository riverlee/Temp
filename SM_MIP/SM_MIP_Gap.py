#!/bin/env python
###################################
# Author: Jiang Li
# Email:  jiang_li@affymetrix.com
# Date:   Tue Jan 13 13:27:06 2015
###################################
import argparse
import re
## Custermized modules/functions
from mylog import info,error,warn,debug
from dust import score_dust
from mip_tm_clip import get_seq
from nib_frag import nibFragger
from local_snp_finder import local_snp_id

## Parser
parser = argparse.ArgumentParser(description="Design MIPs for somatic mutations (SM)")

parser.add_argument("--input_file", dest="input_file", help="SM input file",required=True)
parser.add_argument("--output_file",dest="output_file",help="Output file of designed MIPs for SM",required=True)

### Default options
parser.add_argument('--tm_min', dest="tm_min", help="MIPs with Tms lower than this are ignored (default 55)", type=float, default=55)
parser.add_argument('--tm_max', dest="tm_max", help="MIPs with Tms higher than this are ignored (default 62)", type=float, default=62)
parser.add_argument('--min_hom_length', dest="min_hom_length", help="Min hom length (default 17)", type=int, default=17)
parser.add_argument('--max_hom_length', dest="max_hom_length", help="Max hom length (default 28)", type=int, default=28)
parser.add_argument('--gc_threshold_min', dest="gc_threshold_min", help="Low GC content threshold (default .15)", type=float, default=.15)
parser.add_argument('--gc_threshold_max', dest="gc_threshold_max", help="High GC content threshold (default .85)", type=float, default=.85)
parser.add_argument('--gap_num', dest="gap_num",help="Gap fill bases", type=int, default=1)

args = parser.parse_args()

## Assign to variables
input_file, output_file,tm_min,tm_max = args.input_file,args.output_file,args.tm_min,args.tm_max
min_hom_length,max_hom_length,gc_threshold_min,gc_threshold_max = args.min_hom_length,args.max_hom_length,args.gc_threshold_min,args.gc_threshold_max

gap_num = args.gap_num

#### Functions
def revcomp(dna, reverse=True):
    bases = 'ATGCatgcTACGtacg'
    complement_dict = {bases[i]:bases[i+8] for i in range(8)}
    if reverse:
        dna = reversed(dna)
    result = [complement_dict[base] for base in dna]
    return ''.join(result)

def replaceString(s,to,first=True):
    n=len(to)
    if first:
        return to+s[n:]
    else:
        return s[:-n]+to
    
def make_Hom_pairs(upstream_list,downstream_list,hom_strand,MIP_Mut_Alignment,GapFillBase_M,GapFillBase_W):
    returnStr = ""

    if len(upstream_list)>0 and len(downstream_list) >0:
        for upstream_seq, upstream_tm, upstream_gc, upstream_len in upstream_list:
            for downstream_seq, downstream_tm, downstream_gc, downstream_len in downstream_list:
                # Do on forward strand
                h2,h1,h2_tm,h1_tm,h2_len,h1_len,h2_gc,h1_gc="","","","","","","",""
                
                if hom_strand =="+":
                    h2 = upstream_seq
                    h1 = downstream_seq
                    h2_tm = upstream_tm
                    h1_tm = downstream_tm
                    h2_len = upstream_len
                    h1_len = downstream_len
                    h2_gc = upstream_gc
                    h1_gc = downstream_gc 
                else:
                    h1 = revcomp(upstream_seq)
                    h2 = revcomp(downstream_seq)
                    h1_tm = upstream_tm
                    h2_tm = downstream_tm
                    h1_len = upstream_len
                    h2_len = downstream_len
                    h1_gc = upstream_gc
                    h2_gc = downstream_gc
# gene, mutation_AA,mutation_cDNA,chrome,start,end, ref, alt,gene_strand, cosmic, tumor_type, forwardSeq, reverseSeq 
# are outside  variables

                returnStr+=("\t".join((gene,mutation_AA,mutation_cDNA,chrome,str(start),str(end),ref,alt,gene_strand,cosmic,tumor_type,MIP_Mut_Alignment,GapFillBase_M,GapFillBase_W,upstream_seq,downstream_seq,hom_strand,h2,h1,str(h2_tm),str(h1_tm),str(h2_len),str(h1_len),str(h2_gc),str(h1_gc),forwardSeq,reverseSeq))+"\n")
    return returnStr


def getName(n,j):
    name = ""
    before =n-j-1 
    after = j
    
    if before==1:
        name="GF"
    elif before >1:
        name="GF"+str(before)

    name+="SM"

    if after ==1:
        name+="GF"
    elif after > 1:
        name+="GF"+str(after)

    return name


## Get MIP for SNP
def SNP_MIP_Gap(chrome,start,end,min_hom_length,max_hom_length,tm_min,tm_max,gc_threshold_min,gc_threshold_max,ref,alt):
    returnStr=""
    
    ## New code on 1/21/2015
    # SNP on the gap fill (use + for gap)
    # if gap fill is 2 bases, it will do -+, +-
    # if gap fill is 3 bases, it will do --+, -+-, +--
    ## Fetch sequences

    for n in range(1,gap_num+1):
        for j in range(0,n):
            upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length-(n-j-1),max_hom_length).lower()
            downstream_seq = nibFragger(chrome.replace("chr",""),end+1+j,max_hom_length).lower()
            upstream_list = get_seq(upstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
            downstream_list = get_seq(downstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)

            gapfill_W = nibFragger(chrome.replace("chr",""),start-(n-j-1),n)
            tmplist= list(gapfill_W)
            tmplist[n-j-1]=alt
            gapfill_M = "".join(tmplist)
            upstream_pos = chrome+":"+str(start-(n-j-1)-upstream_list[0][3])+"-"+str(start-(n-j-1)-1)
            downstream_pos = chrome+":"+str(end+1+j)+"-"+str(end+1+j+downstream_list[0][3]-1)

            returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="SNP_on_GF",GapFillBase_M=gapfill_M,GapFillBase_W=gapfill_W)
            returnStr = returnStr.rstrip("\n")
            returnStr+="\t"+upstream_pos+"\t"+downstream_pos+"\t"+getName(n,j)+"\n"
            returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="SNP_on_GF",GapFillBase_M=revcomp(gapfill_M),GapFillBase_W=revcomp(gapfill_W))
            returnStr = returnStr.rstrip("\n")
            returnStr+="\t"+upstream_pos+"\t"+downstream_pos+"\t"+getName(n,j)+"\n"
    return returnStr

""" # old codes (01/16/2015)
    if (re.search(r"[AT]",ref) and re.search(r"[CG]",alt)) or (re.search(r"[AT]",alt) and re.search(r"[CG]",ref)):
        upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length,max_hom_length).lower()
        downstream_seq = nibFragger(chrome.replace("chr",""),end+1,max_hom_length).lower()
        upstream_list = get_seq(upstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
        downstream_list = get_seq(downstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)
    
        returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="SNP_on_GF",GapFillBase_M=alt,GapFillBase_W=ref)
        returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="SNP_on_GF",GapFillBase_M=revcomp(alt),GapFillBase_W=revcomp(ref))
        

    # Gap fill will only cover the SM position
    ## SNP on the H2 forward strand and H1 reverse - Mutation Type
    upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length+1,max_hom_length).lower()
    downstream_seq = nibFragger(chrome.replace("chr",""),end+2,max_hom_length).lower()
    gapfill = nibFragger(chrome.replace("chr",""),end+1,1)
    ### replace the last chracter
    upstream_seq_replaced=replaceString(upstream_seq,alt,first=False)
    upstream_list = get_seq(upstream_seq_replaced, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
    downstream_list = get_seq(downstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)
    
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="SNP_on_H2_M",GapFillBase_M=gapfill,GapFillBase_W=gapfill)
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="SNP_on_H1_M",GapFillBase_M=revcomp(gapfill),GapFillBase_W=revcomp(gapfill))
    
    ## SNP on the H2 forward strand and H1 reverse - Wild Type
    upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length+1,max_hom_length).lower()
    downstream_seq = nibFragger(chrome.replace("chr",""),end+2,max_hom_length).lower()
    gapfill = nibFragger(chrome.replace("chr",""),end+1,1)
    ### replace the last chracter
    upstream_seq_replaced=replaceString(upstream_seq,ref,first=False)
    upstream_list = get_seq(upstream_seq_replaced, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
    downstream_list = get_seq(downstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)
    
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="SNP_on_H2_W",GapFillBase_M=gapfill,GapFillBase_W=gapfill)
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="SNP_on_H1_W",GapFillBase_M=revcomp(gapfill),GapFillBase_W=revcomp(gapfill))
    
    ## SNP on the H2 reverse Strand or H1 forward strand - Mutation Type
    upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length-1,max_hom_length).lower()
    downstream_seq = nibFragger(chrome.replace("chr",""),end,max_hom_length).lower()
    Mgapfill = nibFragger(chrome.replace("chr",""),end-1,1)
    ### replace the first chracter
    downstream_seq_replaced=replaceString(downstream_seq,alt,first=True)
    upstream_list = get_seq(upstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
    downstream_list = get_seq(downstream_seq_replaced, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)

    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="SNP_on_H2_M",GapFillBase_M=revcomp(gapfill),GapFillBase_W=revcomp(gapfill))
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="SNP_on_H1_M",GapFillBase_M=gapfill,GapFillBase_W=gapfill)

    ## SNP on the H2 reverse Strand or H1 forward strand - Wild Type
    upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length-1,max_hom_length).lower()
    downstream_seq = nibFragger(chrome.replace("chr",""),end,max_hom_length).lower()
    Mgapfill = nibFragger(chrome.replace("chr",""),end-1,1)
    ### replace the first chracter
    downstream_seq_replaced=replaceString(downstream_seq,ref,first=True)
    upstream_list = get_seq(upstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
    downstream_list = get_seq(downstream_seq_replaced, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)

    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="SNP_on_H2_W",GapFillBase_M=revcomp(gapfill),GapFillBase_W=revcomp(gapfill))
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="SNP_on_H1_W",GapFillBase_M=gapfill,GapFillBase_W=gapfill)
    return returnStr
"""

## Get MIP for MNP (Multi-nucleotide polymorphism) 
def MNP_MIP(chrome,start,end,min_hom_length,max_hom_length,tm_min,tm_max,gc_threshold_min,gc_threshold_max,ref,alt):    
    returnStr=""

    # MNP on the H2 forward and H1 reverse
    ## Fetch sequences
    upstream_seq = nibFragger(chrome.replace("chr",""),end-max_hom_length+1,max_hom_length).lower()
    downstream_seq = nibFragger(chrome.replace("chr",""),end+2,max_hom_length).lower()
    gapfill = nibFragger(chrome.replace("chr",""),end+1,1)
    upstream_seq = replaceString(upstream_seq,alt,first=False)

    upstream_list = get_seq(upstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
    downstream_list = get_seq(downstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)
    
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="MNP_on_H2",GapFillBase_M=gapfill,GapFillBase_W=gapfill)
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="MNP_on_H1",GapFillBase_M=revcomp(gapfill),GapFillBase_W=revcomp(gapfill))

    ## MNP on the H2 reverse Strand or H1 forward strand
    upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length-1,max_hom_length).lower()
    downstream_seq = nibFragger(chrome.replace("chr",""),start,max_hom_length).lower()
    gapfill = nibFragger(chrome.replace("chr",""),start-1,1)
    ### replace the first chracter
    downstream_seq_replaced=replaceString(downstream_seq,alt,first=True)
    upstream_list = get_seq(upstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
    downstream_list = get_seq(downstream_seq_replaced, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)

    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="MNP_on_H2",GapFillBase_M=revcomp(gapfill),GapFillBase_W=revcomp(gapfill))
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="MNP_on_H1",GapFillBase_M=gapfill,GapFillBase_W=gapfill)

    ## Do MNP_on_H2GF
    for i in range(1,len(ref)+1):
        upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length-1+i,max_hom_length).lower()
        downstream_seq = nibFragger(chrome.replace("chr",""),start+i,max_hom_length).lower()
        gapfillM=alt[i-1]
        gapfillW=ref[i-1]

        upOverlapWithMutation = i -1
        downOverlapWithMutation = len(ref)-i
        if upOverlapWithMutation>0:
            upstream_seq = replaceString(upstream_seq,alt[:i-1],first=False)
        if downOverlapWithMutation > 0:
            downstream_seq = replaceString(downstream_seq,alt[i:],first=True)

        upstream_list = get_seq(upstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
        downstream_list = get_seq(downstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)
        
        MIP_name = "MNP_on_"
        if upOverlapWithMutation>0:
            MIP_name+="H2"
        MIP_name+="Gap"

        if downOverlapWithMutation>0:
            MIP_name+="H1"
        returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment=MIP_name,GapFillBase_M=gapfillM,GapFillBase_W=gapfillW)
        

        MIP_name = "MNP_on_"
        if downOverlapWithMutation>0:
            MIP_name+="H2"
        MIP_name+="Gap"

        if upOverlapWithMutation>0:
            MIP_name+="H1"
        returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment=MIP_name,GapFillBase_M=revcomp(gapfillM),GapFillBase_W=revcomp(gapfillW))



    return returnStr
"""
    ## SNP on the H2 forward strand and H1 reverse
    upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length+1,max_hom_length).lower()
    downstream_seq = nibFragger(chrome.replace("chr",""),end+2,max_hom_length).lower()
    gapfill = nibFragger(chrome.replace("chr",""),end+1,1)
    ### replace the last chracter
    upstream_seq_replaced=replaceString(upstream_seq,alt,first=False)

    upstream_list = get_seq(upstream_seq_replaced, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
    downstream_list = get_seq(downstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)
    
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="SNP_on_H2",GapFillBase_M=gapfill,GapFillBase_W=gapfill)
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="SNP_on_H1",GapFillBase_M=revcomp(gapfill),GapFillBase_W=revcomp(gapfill))


    ## SNP on the H2 reverse Strand or H1 forward strand
    upstream_seq = nibFragger(chrome.replace("chr",""),start-max_hom_length-1,max_hom_length).lower()
    downstream_seq = nibFragger(chrome.replace("chr",""),end,max_hom_length).lower()
    gapfill = nibFragger(chrome.replace("chr",""),end-1,1)
    ### replace the first chracter
    downstream_seq_replaced=replaceString(downstream_seq,alt,first=True)
    upstream_list = get_seq(upstream_seq, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=False)
    downstream_list = get_seq(downstream_seq_replaced, MIN_LENGTH = min_hom_length, MIN_TM=tm_min, MAX_TM=tm_max, GC_MIN=gc_threshold_min, GC_MAX=gc_threshold_max,right2left=True)

    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="-",MIP_Mut_Alignment="SNP_on_H2",GapFillBase_M=revcomp(gapfill),GapFillBase_W=revcomp(gapfill))
    returnStr+=make_Hom_pairs(upstream_list,downstream_list,hom_strand="+",MIP_Mut_Alignment="SNP_on_H1",GapFillBase_M=gapfill,GapFillBase_W=gapfill)

    return returnStr
"""

def determineMutationType(ref,alt):
    mType="SNP"
    if len(ref) ==1 and len(alt) ==1 and re.search(r"[ATCG]",ref) and re.search(r"[ATCG]",alt):
        mType="SNP"
    elif len(ref)>1 and len(ref) == len(alt) and re.search(r"^[ATCG]+$",ref) and re.search(r"^[ATCG]+$",alt):
        mType="MNP"

    return mType

# Main Program
output_fh = open(output_file,"w")

header = "#Gene\tMutation AA\tMutation cDNA\tChr\tStart\tEnd\tRef\tAlt\tGene Strand\tCOSMIC\tTumor type\tMIP_Mut_Alignment\tGapFillBase_M\tGapFillBase_W\tUpstream\tDownstream\tHom strand\textracted Hom2\textracted Hom1\tHom2 tm\tHom1 tm\tHom2 len\tHom1 len\tH2 GC\tH1 GC\tFwd(50bp+Mutation+50bp)\tRev\tUpstream pos\tDownstream pos\tMIP_Mut_Alignment2\n"
output_fh.write(header)

for line in open(input_file):
    parts = line.rstrip().split("\t")
    if line.startswith("Gene"):
        continue

    info("Doing " + "|".join(parts[:3]))
    
    gene,mutation_AA,mutation_cDNA,chrome,start,end,ref,alt,gene_strand,cosmic,tumor_type=parts[0],parts[1],parts[2],parts[3],int(parts[4]),int(parts[5]),parts[6],parts[7],parts[8],parts[9],parts[10]

    
    ## Only works for SNP
    ## Situation 1: SNP_on_GF (snp in gap fill)
    
    flank5 = nibFragger(chrome.replace("chr",""),start-50,50)
    flank3 = nibFragger(chrome.replace("chr",""),end+1,50)
    forwardSeq=flank5.lower()+ "["+ref+"/"+alt+"]"+flank3.lower()
    reverseSeq = revcomp(flank3).lower()+"["+revcomp(ref)+"/"+revcomp(alt)+"]"+revcomp(flank5).lower()
    resultStr = ""
    
    mType = determineMutationType(ref,alt)
    if mType == "SNP":
        resultStr = SNP_MIP_Gap(chrome,start,end,min_hom_length,max_hom_length,tm_min,tm_max,gc_threshold_min,gc_threshold_max,ref,alt)
    else:
        print ("Right now, it only supports SNP")
        continue

    #elif mType == "MNP":
    #    resultStr = MNP_MIP(chrome,start,end,min_hom_length,max_hom_length,tm_min,tm_max,gc_threshold_min,gc_threshold_max,ref,alt)

    if resultStr != "" :
        output_fh.write(resultStr)
        output_fh.write("\n")
    else:
        debug("Failed " + "|".join(parts[:3]))
""" 
    result_snp_on_gf = SNP_on_GF(chrome,start,end,min_hom_length,max_hom_length,tm_min,tm_max,gc_threshold_min,gc_threshold_max,ref,alt)
    if result_snp_on_gf !="" :
        output_fh.write(result_snp_on_gf)

    result_snp_on_H2 = SNP_on_H2(chrome,start,end,min_hom_length,max_hom_length,tm_min,tm_max,gc_threshold_min,gc_threshold_max,ref,alt)
    if result_snp_on_H2 !="" :
        output_fh.write(result_snp_on_H2)

    result_snp_on_H1 = SNP_on_H1(chrome,start,end,min_hom_length,max_hom_length,tm_min,tm_max,gc_threshold_min,gc_threshold_max,ref,alt)
    if result_snp_on_H1 !="" :
        output_fh.write(result_snp_on_H1)
"""
output_fh.close()
                
