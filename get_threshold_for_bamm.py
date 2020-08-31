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
import itertools
import random
from math import log
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


def parse_bamm_and_bg_from_file(bamm_file, bg_file):

    # Read BaMM file
    if os.path.isfile(bamm_file):
        motif_order = 0
        with open(bamm_file) as file:
            for line in file:
                if line[0] != '\n':
                    motif_order = motif_order + 1
                else:
                    break

        # count the motif length
        motif_length = int(sum(1 for line in open(bamm_file)) / (motif_order + 1))

        # read in bamm model
        model = {}
        for k in range(motif_order):
            model[k] = []

        with open(bamm_file) as file:
            for j in range(motif_length):
                for k in range(motif_order):
                    model[k].append([float(p) for p in file.readline().split()])
                file.readline()

    else:
        print('File {0} does not exist'.format(bg_file))
        sys.exit()

    # Read BG file
    bg = {}
    order = 0
    if os.path.isfile(bg_file):
        with open(bg_file) as bgmodel_file:
            line = bgmodel_file.readline()  # skip the first line for K
            line = bgmodel_file.readline()  # skip the second line for Alpha
            while order < motif_order:
                line = bgmodel_file.readline()
                bg_freq = [float(p) for p in line.split()]
                bg[order] = bg_freq
                order += 1
    else:
        print('File {0} does not exist'.format(bg_file))
        sys.exit()
    return(model, bg, order-1)


def make_k_mers(order):
    #  make list with possible k-mer based on bHMM model
    tmp = itertools.product('ACGT', repeat=order + 1)
    k_mer = []
    for i in tmp:
        k_mer.append(''.join(i))
    k_mer_dict = dict()
    index = 0
    for i in k_mer:
        k_mer_dict[i] = index
        index += 1
    return(k_mer_dict)


def make_log_odds_bamm(bamm, bg):
    log_odds_bamm = dict()
    for order in bamm.keys():
        log_odds_bamm[order] = [list(map(lambda x: log(x[0] / x[1], 2), zip(bamm_col, bg[order]))) for bamm_col in bamm[order]]
    return(log_odds_bamm)


def bamm_to_dict(log_odds_bamm, order, k_mers):
    bamm_dict = {}
    for k_mer in k_mers:
        bamm_dict[k_mer] = list()
    for i in range(len(log_odds_bamm[order])):
        for index, k_mer in enumerate(k_mers):
            bamm_dict[k_mer].append(log_odds_bamm[order][i][index])
    return(bamm_dict)


def prepare_bamm(bamm_path, bg_path):
    bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    container = dict()
    log_odds_bamm = make_log_odds_bamm(bamm, bg)
    for i in range(order + 1):
        k_mers = make_k_mers(i)
        bamm_dict = bamm_to_dict(log_odds_bamm, i, k_mers)
        container.update(bamm_dict)
    return(container, order)


def min_score_bamm(bamm, order, length_of_site):
    scores = []
    min_score = 0
    for index in range(order):
        k_mers = itertools.product('ACGT', repeat=index + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for index in range(length_of_site - order):
        k_mers = itertools.product('ACGT', repeat=order + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for s in scores:
        min_score += min(s)
    return(min_score)


def max_score_bamm(bamm, order, length_of_site):
    scores = []
    max_score = 0
    for index in range(order):
        k_mers = itertools.product('ACGT', repeat=index + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for index in range(length_of_site - order):
        k_mers = itertools.product('ACGT', repeat=order + 1)
        scores.append([])
        for k in k_mers:
            k = ''.join(k)
            scores[-1].append(bamm[k][index])
    for s in scores:
        max_score += max(s)
    return(max_score)


def to_score(norm_value, bamm, order, length_of_site):
    '''
    norm = (score - min) / (max - min) -> score = norm * (max - min) + min
    '''
    max_s = max_score_bamm(bamm, order, length_of_site)
    min_s = min_score_bamm(bamm, order, length_of_site)
    score = norm_value * (max_s - min_s) + min_s
    return(score)


def to_norm(score, bamm, order, length_of_site):
    max_s = max_score_bamm(bamm, order, length_of_site)
    min_s = min_score_bamm(bamm, order, length_of_site)
    norm_value = (score - min_s) / (max_s - min_s)
    return(norm_value)


def score_bamm(site, bamm, order, length_of_site):
    score = float()
    for index in range(order):
        score += bamm[site[0:index + 1]][index]
    for index in range(length_of_site - order):
        score += bamm[site[index:index+order + 1]][index + order]
    return(score)


# def calculate_scores(peaks, bamm, order, length_of_site):
#     scores = [score_bamm(peak[i:length_of_site + i], bamm, order, length_of_site) for peak in peaks for i in range(len(peak) - length_of_site + 1)]
#     return(scores)


# def calculate_scores(peaks, bamm, order, length_of_site):
#     sites = (peak[i:length_of_site + i] for peak in peaks for i in range(len(peak) - length_of_site + 1))
#     scores = []
#     append = scores.append
#     for site in sites:
#         if check_nucleotides(site):
#             append(score_bamm(site, bamm, order, length_of_site))
#         else:
#             continue
#     return(scores)


def calculate_scores(peaks, bamm, order, length_of_site, threshold):
    scores = []
    append = scores.append
    for peak in peaks:
        for i in range(len(peak) - length_of_site + 1):
            site = peak[i:length_of_site + i]
            if 'N' in site:
                continue
            s = score_bamm(site, bamm, order, length_of_site)
            if s > threshold:
                append(s)
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

    
# def get_threshold(scores, path_out, number_of_sites):
#     scores.sort(reverse=False) # sorted score from small to big
#     #scores = [round(i, 3) for i in scores]
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
            elif count/number_of_sites > 0.01:
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
    parser.add_argument('bamm', action='store', help='path to .ihbcp file from BaMMmotif output')
    parser.add_argument('bg', action='store', help='path to .hbcp file from BaMMmotif output')
    parser.add_argument('out', action='store', help='path to write results')
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())
    

def main():
    args = parse_args()
    bamm_path = args.bamm
    bg_path = args.bg
    fasta_path = args.fasta
    path_out = args.out

    peaks = read_fasta(fasta_path)
    bamm, order = prepare_bamm(bamm_path, bg_path)
    length_of_site = len(bamm['A'])
    number_of_sites = sum([len(range(len(peak) - length_of_site + 1)) for peak in peaks])
    threshold = to_score(0.3, bamm, order, length_of_site)
    scores = calculate_scores(peaks, bamm, order, length_of_site, threshold)
    get_threshold(scores, number_of_sites, path_out)

if __name__ == '__main__':
    main()
