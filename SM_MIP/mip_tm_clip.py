#!/usr/bin/python
from math import log

deltaHDict = {'AA':8.0, 'CA':8.2, 'GA':8.8,  'TA':5.6,
              'AC':9.4, 'CC':10.9,'GC':10.5, 'TC':8.8,
              'AG':6.6, 'CG':11.8,'GG':10.9, 'TG':8.2,
              'AT':5.6, 'CT':6.6, 'GT':9.4,  'TT':8.0}
deltaSDict = {'AA':21.9, 'CA':21.0, 'GA':23.5, 'TA':15.2,
              'AC':25.5, 'CC':28.4, 'GC':26.4, 'TC':23.5,
              'AG':16.4, 'CG':29.0, 'GG':28.4, 'TG':21.0,
              'AT':15.2, 'CT':16.4, 'GT':25.5, 'TT':21.9}

def optimize_seq(seq, MIN_LENGTH = 17, MIN_TM = 55):
    #seq = seq[::-1]  #mip is at right hand side - optimize from there
    for seq_pos_index in range(MIN_LENGTH, len(seq)+1):
        test_seq =seq[:seq_pos_index]
        test_tm = seq_tm(test_seq)
        if test_tm > MIN_TM:
            return test_seq, test_tm
    return seq, seq_tm(seq)

def optimize_seq2(seq, MIN_LENGTH = 17, MIN_TM = 55):
    #seq = seq[::-1]  #mip is at right hand side - optimize from there

    for seq_pos_index in range(MIN_LENGTH, len(seq)+1):
        test_seq =seq[(len(seq)-i):len(seq)]
        test_tm = seq_tm(test_seq)
        if test_tm > MIN_TM:
            return test_seq, test_tm
    return seq, seq_tm(seq)
	
# Return a list of tuples (updated 01/16/2015, only output the optimize one)
def get_seq(seq, MIN_LENGTH = 17, MIN_TM=55, MAX_TM=62, GC_MIN=0.15, GC_MAX=0.85,right2left=True):
    result = list()
    for seq_pos_index in range(MIN_LENGTH, len(seq)+1):
        test_seq = ""
        if right2left:
            test_seq = seq[:seq_pos_index]
        else:
            test_seq = seq[(len(seq)-seq_pos_index):len(seq)]
        test_tm = seq_tm(test_seq)
        test_gc = seq_gc(test_seq)
        if test_tm > MIN_TM and test_tm < MAX_TM and test_gc > GC_MIN and test_gc < GC_MAX:
            result.append((test_seq,test_tm,test_gc,len(test_seq))) 
            return result
    return result

# Return a list of tuples
def get_seq_bk(seq, MIN_LENGTH = 17, MIN_TM=55, MAX_TM=62, GC_MIN=0.15, GC_MAX=0.85,right2left=True):
    result = list()
    for seq_pos_index in range(MIN_LENGTH, len(seq)+1):
        test_seq = ""
        if right2left:
            test_seq = seq[:seq_pos_index]
        else:
            test_seq = seq[(len(seq)-seq_pos_index):len(seq)]
        test_tm = seq_tm(test_seq)
        test_gc = seq_gc(test_seq)
        if test_tm > MIN_TM and test_tm < MAX_TM and test_gc > GC_MIN and test_gc < GC_MAX:
            result.append((test_seq,test_tm,test_gc,len(test_seq)))
    return result

def seq_gc(seq):
    gc_count = seq.count("C") + seq.count("G") + seq.count("c") + seq.count("g")
    gc_content = gc_count / float(len(seq))
    return gc_content
	
def seq_tm(seq):
    global deltaHDict
    global deltaSDict
    if len(seq) < 8:
        ATseq = [x for x in list(seq) if (x == 'A' or x == 'T')]
        CGseq = [x for x in list(seq) if (x == 'C' or x == 'G')]
        return 2 * len(ATseq) + 4 * len(CGseq)
    else:
        dDeltah, dDeltas = 0, 0

        for iIndex in range(len(seq) - 1):
            dDeltah += deltaHDict[seq[iIndex : iIndex + 2].upper()]
            dDeltas += deltaSDict[seq[iIndex : iIndex + 2].upper()]

        dRlnK = 1.987 * log(1.0 / (100 * 1.0e-9))
        dSaltadj = 7.21 * log(100.0 / 1000.0)

        return (1000.0 * (dDeltah - 3.4)) / (dDeltas + dRlnK) - 272.9 + dSaltadj

