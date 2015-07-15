#!/usr/bin/python

import subprocess

def nibFragger(chromosome, position, seq_length=35):

    nibStart = int(position) - 1
    nibEnd = int(position) + int(seq_length)
    forward_total_sequence, reverse_total_sequence = "", ""

    nibLocation = "/nfs/gp/Hs/2009_02-hg19/nib/chr" + str(chromosome) + ".nib"
    proc = subprocess.Popen(["nibFrag %s %d %d + stdout" % (nibLocation,nibStart,nibEnd-1)], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    nib_data = out.decode('utf8').split("\n")
    for index in range(1, len(nib_data)):
        forward_total_sequence += nib_data[index].upper()
    return forward_total_sequence

