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
import time
import argparse


def make_pcm(motifs):
    matrix = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        matrix[i[0]] = []
    len_of_motif = len(motifs[0])
    for i in matrix.keys():
        matrix[i] = [0]*len_of_motif
    for i in range(len_of_motif):
        for l in motifs:
            matrix[l[i]][i] += 1
    return(matrix)


def make_pfm(pcm):
    number_of_sites = [0] * len(pcm['A'])
    for key in pcm.keys():
        for i in range(len(pcm[key])):
            number_of_sites[i] += pcm[key][i]
    pfm = dict()
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pfm[i[0]] = []
    first_key = list(pcm.keys())[0]
    nuc_pseudo = 1/len(pcm.keys())
    for i in range(len(pcm[first_key])):
        for nuc in pcm.keys():
            pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
    return(pfm)


def make_pwm(pfm):
    pwm = {}
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pwm[i[0]] = []
    first_key = list(pfm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def sites_to_pwm(sites):
    pcm = make_pcm(sites)
    pfm = make_pfm(pcm)
    pwm = make_pwm(pfm)
    return(pwm)


def parse_chipmunk(path):
    with open(path, 'r') as file:
        container = []
        for line in file:
            d = {'name': str(), 'start': int(), 'end': int(),
                 'seq': str(), 'strand': str()}
            if line.startswith('WORD|'):
                line = line[5:].strip()
                line = line.split()
                d['name'] = 'peaks_' + str(int(line[0]) - 1)
                d['start'] = int(line[1])
                d['end'] = int(line[1]) + len(line[2])
                d['seq'] = line[2]
                if line[4] == 'direct':
                    d['strand'] = '+'
                else:
                    d['strand'] = '-'
                container.append(d)
            else:
                continue
    seqs = [i['seq'] for i in container if not 'N' in i['seq']]
    return(seqs)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def run_chipmunk(path_to_java, path_to_chipmunk, fasta_path, path_out, motif_length_start, motif_length_end, cpu_count):
    args = [path_to_java, '-cp', path_to_chipmunk,
                   'ru.autosome.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', '1.0',
                   's:{}'.format(fasta_path),
                  '100', '10', '1', str(cpu_count), 'random']
    p = subprocess.run(args, shell=False, capture_output=True)
    out = p.stdout
    with open(path_out, 'wb') as file:
        file.write(out)
    return(0)


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


def score_pwm(seq, pwm):
    score = 0
    position = 0 
    length = len(seq)
    for index in range(length):
        score += pwm[seq[index]][position]
        position += 1
    return(score)


def false_scores_pwm(peaks, pwm):
    false_scores = []
    length_pwm = len(pwm['A'])
    append = false_scores.append
    for peak in peaks:
        complement_peak = complement(peak)
        length_pwm = len(pwm['A'])
        full_peak = peak + 'N' * length_pwm + complement_peak
        n = len(full_peak) - length_pwm + 1
        for i in range(n):
            site = peak[i:length_pwm + i]
            if 'N' in site:
                continue
            score = score_pwm(site, pwm)
            false_scores.append(score)
    return(false_scores)


def true_scores_pwm(peaks, pwm):
    true_scores = []
    for peak in peaks:
        complement_peak = complement(peak)
        length_pwm = len(pwm['A'])
        best = -1000000
        full_peak = peak + 'N' * length_pwm + complement_peak
        n = len(full_peak) - length_pwm + 1
        for i in range(n):
            site = full_peak[i:length_pwm + i]
            if 'N' in site:
                continue
            score = score_pwm(site, pwm)
            if score >= best:
                best = score
        true_scores.append(best)
    return(true_scores)


def creat_background(peaks, length_of_site, counter):
    shuffled_peaks = []
    number_of_sites = 0
    while counter > number_of_sites:
        peak = random.choice(peaks)
        shuffled_peak = ''.join(random.sample(peak, len(peak)))
        shuffled_peaks.append(shuffled_peak)
        number_of_sites += (len(''.join(shuffled_peak)) - length_of_site + 1) * 2
    return(shuffled_peaks)


def bootstrap_pwm(peaks, length_of_site, counter, path_to_java, path_to_chipmunk, tmp_dir, cpu_count):
    true_scores = []
    false_scores = []
    number_of_peaks = len(peaks)
    for i in range(5):
        train_peaks = random.choices(peaks, k=round(0.9 * number_of_peaks))
        test_peaks = [peak for peak in peaks if not peak in train_peaks]
        shuffled_peaks = creat_background(test_peaks, length_of_site, counter / 5)
        write_fasta(train_peaks, tmp_dir + '/train.fasta')
        run_chipmunk(path_to_java, path_to_chipmunk,
                     tmp_dir + '/train.fasta', tmp_dir + '/chipmunk_results.txt',
                     length_of_site, length_of_site, cpu_count)
        sites = parse_chipmunk(tmp_dir + '/chipmunk_results.txt')
        sites = list(set(sites))
        nsites = len(sites)
        pcm = make_pcm(sites)
        pfm = make_pfm(pcm)
        pwm = make_pwm(pfm)
        for true_score in true_scores_pwm(test_peaks, pwm):
            true_scores.append(true_score)
        for false_score in false_scores_pwm(shuffled_peaks, pwm):
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
    parser.add_argument('chipmunk', action='store', help='path to ChIPMunk source')
    parser.add_argument('-j', '--java', action='store', type=str, dest='java',
                        required=False, default='java', help='Path to java')
    parser.add_argument('-t', '--tmp', action='store', type=str, dest='tmp',
                        required=False, default='./pwm.tmp', help='tmp directory')
    parser.add_argument('-p', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=4, help='Number of processes to use, default: 4')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())
    
    
def main():
    args = parse_args()
    peaks_path = args.fasta
    results_path = args.results
    length_of_site = args.length
    path_to_chipmunk = args.chipmunk
    path_to_java = args.java
    tmp_dir = args.tmp
    cpu_count = args.cpu_count
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)        
    counter = 5000000
    peaks = read_peaks(peaks_path)
    table = bootstrap_pwm(peaks, length_of_site, counter, path_to_java, path_to_chipmunk, tmp_dir, cpu_count)
    write_table_bootstrap(results_path, table)
    shutil.rmtree(tmp_dir)
    return(0)


def bootstrap_for_pwm(peaks_path, results_path, length_of_site, path_to_java, path_to_chipmunk, tmp_dir, cpu_count, counter=5000000):
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    peaks = read_peaks(peaks_path)
    table = bootstrap_pwm(peaks, length_of_site, counter, path_to_java, path_to_chipmunk, tmp_dir, cpu_count)
    write_table_bootstrap(results_path, table)
    shutil.rmtree(tmp_dir)
    return(0)


if __name__=="__main__":
    main()

