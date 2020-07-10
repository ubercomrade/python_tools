'''
Copyright © 2018 Anton Tsukanov. Contacts: tsukanov@bionet.nsc.ru
License: http://www.gnu.org/licenses/gpl.txt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

import argparse
import sys
import os
import random
import bisect
from collections import Counter
from bisect import bisect, bisect_left
from operator import itemgetter


def read_fasta(path):

    fasta = list()
    append = fasta.append
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                line = line.strip().upper()
                #if not 'N' in line:
                append(line)
                append(complement(line))
    file.close()
    return(fasta)


def replace_n_to_randome_nucl(peak):
    container = []
    for n in peak:
        if n == 'N':
            container.append(random.choice(['A', 'C', 'G', 'T']))
    return(''.join(container))


def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
        for key in pwm:
            pwm[key] = tuple(pwm[key])
    file.close()
    return(pwm)


def score(seq, pwm):
    score = sum([pwm[letter][index] for index, letter in enumerate(seq)])
    return(score)


def to_score(norm_value, pwm):
    '''
    norm = (score - min) / (max - min) -> score = norm * (max - min) + min
    '''
    max_s = max_score(pwm)
    min_s = min_score(pwm)
    score = norm_value * (max_s - min_s) + min_s
    return(score)


def to_norm(score, pwm):
    max_s = max_score(pwm)
    min_s = min_score(pwm)
    norm_value = (score - min_s) / (max_s - min_s)
    return(norm_value)


def min_score(pwm):
    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += min(tmp)
    return(value)


def max_score(pwm):
    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += max(tmp)
    return(value)


# ALL
def calculate_scores(peaks, pwm, length_of_site):
    sites = (peak[i:length_of_site + i] for peak in peaks for i in range(len(peak) - length_of_site + 1))
    scores = []
    append = scores.append
    for site in sites:
        if check_nucleotides(site):
            append(score(site, pwm))
        else:
            continue
    return(scores)


# THREHOLD
# def calculate_scores(peaks, pwm, length_of_site, threshold):
#     scores = []
#     append = scores.append
#     for peak in peaks:
#         for i in range(len(peak) - length_of_site + 1):
#             site = peak[i:length_of_site + i]
#             if 'N' in site:
#                 continue
#             s = score(site, pwm)
#             if s > threshold:
#                 append(s)
#             else:
#                 continue
#     return(scores)


# NUMBER_OF_SITES
# def calculate_scores(peaks, pwm, length_of_site, number_of_sites):
#     list_of_scores = deque([], maxlen=number_of_sites)
#     for peak in peaks:
#         for i in range(len(peak) - length_of_site + 1):
#             site = peak[i:length_of_site + i]
#             s = score(site, pwm)
#             if len(list_of_scores) != number_of_sites:
#                 #bisect.insort(list_of_scores, s)
#                 list_of_scores.insert(bisect.bisect_left(list_of_scores, s), s)
#             else:
#                 maximum = list_of_scores[-1]
#                 minimum = list_of_scores[0]
#                 if s < minimum:
#                     continue
#                 elif s > maximum:
#                     list_of_scores.append(s)
#                 else:
#                     list_of_scores.popleft()
#                     list_of_scores.insert(bisect.bisect_left(list_of_scores, s), s)
#     return(list_of_scores)


def check_nucleotides(site):
    s = set(site)
    n = {'A', 'C', 'G', 'T'}
    if len(s - n) == 0:
        return(True)
    else:
        return(False)


def complement(seq):
    #return(seq[::-1].translate(seq.maketrans('ACGT', 'TGCA')))
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


# def get_threshold(scores, path_out, number_of_sites):
#     scores.sort(reverse=False) # sorted score from small to big
#     counts = Counter(scores)

#     found_sites = counts[list(counts.keys())[::-1][0]]
#     for index, k in enumerate(list(counts.keys())[::-1][1:]):
#         found_sites = counts[k] + found_sites
#         counts[k] = found_sites

#     fprs_table = []
#     append = fprs_table.append
#     for index, k in enumerate(list(counts.keys())[::-1]):
#         append(( k, counts[k] / number_of_sites ))


#     with open(path_out, "w") as file:
#         file.write("Scores\tFPR\n")
#         for (score, fpr) in fprs_table:
#             if fpr > 0.00051:
#                 break
#             else:
#                 file.write("{0}\t{1}\n".format(score, fpr))
#     file.close()
#     return(0)


def get_threshold(scores, number_of_sites, path_out):
    scores.sort(reverse=True) # big -> small
    with open(path_out, "w") as file:
        last_score = scores[0]
        for count, score in enumerate(scores[1:], 1):
            if score == last_score:
                continue
            elif count/number_of_sites > 0.0005:
                file.write("{0}\t{1}\n".format(last_score, count/number_of_sites))
                break
            elif score != last_score:
                file.write("{0}\t{1}\n".format(last_score, count/number_of_sites))
                last_score = score 
    file.close()
    return(0)
    

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('fasta', action='store', help='path to fasta')
    parser.add_argument('pwm', action='store', help='path to PWM file')
    parser.add_argument('out', action='store', help='path to write results')
    parser.add_argument('-p', '--false_positive', action='store', type=float, dest='false_positive',
                            required=False, help='value of FP (FP ~ P-VALUE)')


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())
    

def main():
    args = parse_args()

    pwm_path = args.pwm
    fasta_path = args.fasta
    path_out = args.out
    fp = args.false_positive

    peaks = read_fasta(fasta_path)
    pwm = read_pwm(pwm_path)
    length_of_site = len(pwm['A'])
    #number_of_sites = sum([len(range(len(peak) - length_of_site + 1)) for peak in peaks])
    #threshold = to_score(.7, pwm)
    scores = calculate_scores(peaks, pwm, length_of_site)#, threshold)
    number_of_sites = len(scores) 
    print(number_of_sites)
    get_threshold(scores, number_of_sites, path_out)

if __name__ == '__main__':
    main()
