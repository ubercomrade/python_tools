import csv
import itertools
import os
import sys
import math
from math import log
import re
import random
import subprocess
import shutil
import argparse


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
        log_odds_bamm[order] = [list(map(lambda x: math.log(x[0] / x[1], 2), zip(bamm_col, bg[order]))) for bamm_col in bamm[order]]
    return(log_odds_bamm)


def bamm_to_dict(log_odds_bamm, order, k_mers):
    bamm_dict = {}
    for k_mer in k_mers:
        bamm_dict[k_mer] = list()
    for i in range(len(log_odds_bamm[order])):
        for index, k_mer in enumerate(k_mers):
            bamm_dict[k_mer].append(log_odds_bamm[order][i][index])
    return(bamm_dict)


def read_bamm(bamm_path, bg_path):
    bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    container = dict()
    log_odds_bamm = make_log_odds_bamm(bamm, bg)
    for i in range(order + 1):
        k_mers = make_k_mers(i)
        bamm_dict = bamm_to_dict(log_odds_bamm, i, k_mers)
        container.update(bamm_dict)
    return(container, order)


def write_fasta(peaks, path):
    with open(path, 'w') as file:
        for index, p in enumerate(peaks):
            file.write('>{}\n'.format(index))
            file.write(p + '\n')
    return(0)


def read_peaks(path):
    container = []
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                append(line.strip().upper())
    return(container)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def score_bamm(site, bamm, order, length_of_site):
    score = 0
    for index in range(order):
        score += bamm[site[0:index + 1]][index]
    for index in range(length_of_site - order):
        score += bamm[site[index:index+order + 1]][index + order]
    return(score)


def false_scores_bamm(peaks, bamm, order, length_of_site):
    false_scores = []
    append = false_scores.append
    for peak in peaks:
        complement_peak = complement(peak)
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if 'N' in site:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
            false_scores.append(score)
    return(false_scores)


def true_scores_bamm(peaks, bamm, order, length_of_site):
    true_scores = []
    for peak in peaks:
        complement_peak = complement(peak)
        best = -1000000
        full_peak = peak + 'N' * length_of_site + complement_peak
        n = len(full_peak) - length_of_site + 1
        for i in range(n):
            site = full_peak[i:length_of_site + i]
            if 'N' in site:
                continue
            score = score_bamm(site, bamm, order, length_of_site)
            if score >= best:
                best = score
        true_scores.append(best)
    return(true_scores)


def create_bamm_model(directory, order, meme):
    fasta_path = directory + '/train.fasta'
    args = ['BaMMmotif', directory, fasta_path, '--PWMFile', meme, '--EM', '--order', str(order), '--Order', str(order)]
    r = subprocess.run(args, capture_output=True)
    bamm_path = directory + '/train_motif_1.ihbcp'
    bg_path = directory + '/train.hbcp'
    bamm, order = read_bamm(bamm_path, bg_path)
    return(bamm, order)


# def creat_background(peaks, length_of_site, counter):
#     peaks = ''.join(peaks)
#     n = int((counter // ((len(peaks) - length_of_site + 1) * 2)) + 1)
#     print((counter // ((len(peaks) - length_of_site + 1) * 2)))
#     print(n)
#     shuffled_peaks = []
#     for i in range(n):
#         shuffled_peak = ''.join(random.sample(peaks, len(peaks)))
#         shuffled_peaks.append(shuffled_peak)
#     print((len(''.join(shuffled_peaks)) - length_of_site + 1) * 2)
#     return(shuffled_peaks)


def creat_background(peaks, length_of_site, counter):
    shuffled_peaks = []
    number_of_sites = 0
    while counter > number_of_sites:
        peak = random.choice(peaks)
        shuffled_peak = ''.join(random.sample(peak, len(peak)))
        shuffled_peaks.append(shuffled_peak)
        number_of_sites += (len(''.join(shuffled_peak)) - length_of_site + 1) * 2
    return(shuffled_peaks)


def bootstrap_bamm(peaks, length_of_site, counter, order, meme, tmp_dir):
    true_scores = []
    false_scores = []
    number_of_peaks = len(peaks)
    for i in range(10):
        train_peaks = random.choices(peaks, k=round(0.9 * number_of_peaks))
        test_peaks = [peak for peak in peaks if not peak in train_peaks]
        shuffled_peaks = creat_background(test_peaks, length_of_site, counter / 10)
        write_fasta(train_peaks, tmp_dir + '/train.fasta')
        bamm, order = create_bamm_model(tmp_dir, order, meme)
        for true_score in true_scores_bamm(test_peaks, bamm, order, length_of_site):
            true_scores.append(true_score)
        for false_score in false_scores_bamm(shuffled_peaks, bamm, order, length_of_site):
            false_scores.append(false_score)
    table = creat_table_bootstrap(true_scores, false_scores)
    return(table)


def creat_table_bootstrap(true_scores, false_scores):
    table = []
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    for tpr in [round(i * 0.01, 2) for i in range(5,105, 5)]:
        score = true_scores[round(true_length * tpr) - 1]
        actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
        fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
        table.append({'Scores': score, 'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
    return(table)


def write_table_bootstrap(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = data[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for line in data:
            writer.writerow(line)
    return(0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', action='store', help='path to file with peaks')
    parser.add_argument('results', action='store', help='path to write table with ROC')
    parser.add_argument('length', action='store', type=int, help='length of TFBS')
    parser.add_argument('meme', action='store', help='path to PWM in MEME format')
    parser.add_argument('-t', '--tmp', action='store', type=str, dest='tmp',
                        required=False, default='./bamm.tmp', help='tmp directory')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())
    
    
def main():
    args = parse_args()
    peaks_path = args.fasta
    results_path = args.results
    length_of_site = args.length
    meme = args.meme
    tmp_dir = args.tmp
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
        
    peaks = read_peaks(peaks_path)
    counter = 5000000
    order = 2
    table = bootstrap_bamm(peaks, length_of_site, counter, order, meme, tmp_dir)
    write_table_bootstrap(results_path, table)
    shutil.rmtree(tmp_dir)
    return(0)


def bootstrap_for_bamm(peaks_path, results_path, length_of_site, meme, tmp_dir, counter = 5000000, order=2):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    peaks = read_peaks(peaks_path)
    table = bootstrap_bamm(peaks, length_of_site, counter, order, meme, tmp_dir)
    write_table_bootstrap(results_path, table)
    shutil.rmtree(tmp_dir)
    return(0)


if __name__=="__main__":
    main()
