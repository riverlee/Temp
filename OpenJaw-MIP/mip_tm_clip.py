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

    for seq_pos_index in range(MIN_LENGTH, len(seq)):
        test_seq =seq[:seq_pos_index]
        test_tm = seq_tm(test_seq)
        if test_tm > MIN_TM:
            return test_seq, test_tm
    return seq, seq_tm(seq)


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

