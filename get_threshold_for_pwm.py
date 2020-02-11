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
from collections import Counter
from bisect import bisect, bisect_left
from operator import itemgetter
#import numpy as np


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


# def calculate_scores(peaks, pwm, length_of_site):
#     scores = [score(peak[i:length_of_site + i], pwm) for peak in peaks for i in range(len(peak) - length_of_site + 1)]
#     return(scores)


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


def create_fpr_table(scores, unique):
    container = []
    append = container.append
    #scores = np.array(scores)
    print('create_fpr_table')
    for u in unique:
        #append((u, len(scores[:np.searchsorted(scores, u)]))) 
        fpr = len(scores[:bisect_left(scores, u)]) / len(scores)
        if fpr > 0.0001:
            break
        else:
            append((u, fpr)) 
    return(container)


def get_threshold(scores, path_out, number_of_sites):
    scores.sort(reverse=False) # sorted score from small to big
    scores = [round(i, 3) for i in scores]
    counts = Counter(scores)

    found_sites = counts[list(counts.keys())[::-1][0]]
    for index, k in enumerate(list(counts.keys())[::-1][1:]):
        found_sites = counts[k] + found_sites
        counts[k] = found_sites

    fprs_table = []
    append = fprs_table.append
    for index, k in enumerate(list(counts.keys())[::-1]):
        append(( k, counts[k] / number_of_sites ))


    with open(path_out, "w") as file:
        file.write("Scores\tFPR\n")
        for (score, fpr) in fprs_table:
            file.write("{0}\t{1}\n".format(score, fpr))
    file.close()


# def get_threshold(scores, path_out):
#     scores.sort(reverse=True) # sorted score from big to small
#     fprs = [5*10**(-4), 3.33*10**(-4), 1.90*10**(-4), 1.02*10**(-4), 5.24*10**(-5)]
#     with open(path_out, "w") as file:
#         file.write("Scores\tFPR\n")
#         for fpr in fprs:
#             thr = scores[int(fpr * len(scores))]
#             file.write("{0}\t{1}\n".format(thr, fpr))
#     file.close()


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
    number_of_sites = sum([len(range(len(peak) - length_of_site + 1)) for peak in peaks])
    scores = calculate_scores(peaks, pwm, length_of_site)
    get_threshold(scores, path_out, number_of_sites)


if __name__ == '__main__':
    main()
