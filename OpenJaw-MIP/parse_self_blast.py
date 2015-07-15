import argparse
from collections import defaultdict
from os import listdir

parser = argparse.ArgumentParser(description='Parses self-blast results, returns confliciting homs. Append the result to your existing file')
parser.add_argument("--input_file", dest="input_file", help="File containing the results of self-blast", default="IDH1.out")
parser.add_argument("--conflict_count", dest="conflict_count", help="How many non-self matches are a significant problem (threshold) - default 8", type=int, default=8)
parser.add_argument("--min_align_length", dest="min_align_length", help="Minimum alignment length (default 14)", type=int, default=14)
args = parser.parse_args()

input_file, conflict_count, min_align_length = args.input_file, args.conflict_count, args.min_align_length
hom_hits = defaultdict(int)

for blast_result_line in open(input_file):
    if blast_result_line.startswith("#"): continue
    (hom_id, database_id, pct_id, align_length, mismatches, gaps, q_start, q_end, s_start, s_end, e_val, bit_score) = blast_result_line.rstrip().split("\t")
    hom_name = hom_id.split(".")[0]
    #Perhaps update the database name so that the split will not be kludgy like this
    database_parts = database_id.split("-")
    h2_name, h1_name = "{}-{}".format(database_parts[0], database_parts[1]), "{}-{}".format(database_parts[2], database_parts[3])

    #Toss self-hits
    if hom_name == h1_name or hom_name == h2_name:
        continue

    #Is this a reasonable hit
    if int(align_length) > min_align_length:
        hom_hits[hom_name] += 1

for hom_name in hom_hits:
    if hom_hits[hom_name] >= conflict_count:
        print (hom_name)



