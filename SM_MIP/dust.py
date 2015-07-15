import re
from itertools import permutations

def mask_sequence(sequence):
    #simple repeats
    two_mers = [''.join(x) for x in permutations('ACTG', 2)]
    three_mers = [''.join(x) for x in permutations('ACTG', 3)]
    for mer in two_mers + three_mers:
        if re.search('(%s){3,}' % (mer), sequence):
            sequence = re.sub('(%s){3,}' % (mer), 'N' * (re.search('(%s){3,}' % (mer), sequence).end() - re.search('(%s){3,}' % (mer), sequence).start()), sequence)
    return sequence

def score_dust(sequence):
    masked_sequence = mask_sequence(sequence)

    if re.search(r'AAAAAAA', sequence) or re.search(r'TTTTTTT', sequence): hp_run = 1
    elif re.search(r'GGGGGG', sequence) or re.search(r'CCCCCC', sequence): hp_run = 1
    else: hp_run = 0

    #Mimic Fan's script - penalties for sequences closer to the gap fill

    closest_5_chars, next_5_chars, rest_of_sequence =  masked_sequence[len(masked_sequence)-5:], \
                                                       masked_sequence[len(masked_sequence)-10:len(masked_sequence)-5], \
                                                       masked_sequence[:len(masked_sequence)-10]
    H2_score =  5 * closest_5_chars.count('N') + 3 * next_5_chars.count('N') + rest_of_sequence.count('N')
    dust_pct_H2 = (H2_score / float(len(masked_sequence) + 30))

    closest_5_chars, next_5_chars, rest_of_sequence =  masked_sequence[:5], \
                                                       masked_sequence[5:10], \
                                                       masked_sequence[10:]
    H1_score = 5 * closest_5_chars.count('N') + 3 * next_5_chars.count('N') + rest_of_sequence.count('N')
    dust_pct_H1 = (H2_score / float(len(masked_sequence) + 30))

    return hp_run, H1_score, H2_score, dust_pct_H1, dust_pct_H2



