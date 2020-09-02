import csv
import itertools
import os
import sys
import math
import re
import random
import subprocess
import shutil
import glob
from operator import itemgetter
import argparse


def read_peaks(path):
    container = []
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                append(line.strip().upper())
    return(container)


def read_inmode_bed(path):
    table = []
    with open(path) as file:
        for line in file:
            line = line.strip().split('\t')
            line[0] = int(line[0])
            line[1] = int(line[1])
            line[2] = int(line[2])
            line[4] = float(line[4])
            table.append(line)
    return(table)


def write_fasta(sites, tmp_dir, tag):
    with open('{0}/{1}.fa'.format(tmp_dir, tag), 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    return(0)


def true_scores_inmode(path_to_inmode, path_to_java, motif_length, tmp_dir, tag):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/Learned_DeNovo({1},2,2)_motif/XML_of_DeNovo({1},2,2)_motif.xml'.format(tmp_dir, motif_length),
            'id={0}/{1}.fa'.format(tmp_dir, tag), 'f=1.0', 'outdir={}'.format(tmp_dir), 'bs=false']
    r = subprocess.run(args, capture_output=True)
    scores = []
    table = read_inmode_bed('{0}/{1}'.format(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED"))
    table.sort(key=itemgetter(0, 4))
    last_index = 0
    for line in table:
        index = line[0]
        score = line[4]
        if last_index != index:
            scores.append(last_score)
        last_score = score
        last_index = index
    scores.append(score)
    scores = [math.log(float(i), 10) for i in scores]
    os.remove(tmp_dir + '/Binding_sites_from_SequenceScan(1.0).txt')
    os.remove(tmp_dir + '/Motif_hits_from_SequenceScan(1.0).BED')
    os.remove(tmp_dir + '/protocol_scan.txt')
    os.remove(tmp_dir + '/{}.fa'.format(tag))
    return(scores)


def false_scores_inmode(path_to_inmode, path_to_java, motif_length, tmp_dir, tag):
    scores = []
    args = [path_to_java, '-Xmx16G', '-Xms1G', 
            '-jar',
            path_to_inmode, 'scan',
            'i={0}/Learned_DeNovo({1},2,2)_motif/XML_of_DeNovo({1},2,2)_motif.xml'.format(tmp_dir, motif_length),
            'id={0}/{1}.fa'.format(tmp_dir, tag), 'f=1.0', 'outdir={}'.format(tmp_dir), 'bs=false']
    r = subprocess.run(args, capture_output=True)
    with open('{0}/{1}'.format(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED")) as file:
        for line in file:
            scores.append(math.log(float(line.split()[4]), 10))
    os.remove(tmp_dir + '/Binding_sites_from_SequenceScan(1.0).txt')
    os.remove(tmp_dir + '/Motif_hits_from_SequenceScan(1.0).BED')
    os.remove(tmp_dir + '/protocol_scan.txt')
    os.remove(tmp_dir + '/{}.fa'.format(tag))
    return(scores)


def make_inmode(path_to_inmode, path_to_java, motif_length, order, tmp_dir):
    args = [path_to_java, '-Xmx16G', '-Xms1G', '-jar', path_to_inmode,
    'denovo', 'i={}/train.fa'.format(tmp_dir), 'm={}'.format(motif_length), 'outdir={}'.format(tmp_dir),
    'mo={}'.format(order)]
    r = subprocess.run(args, capture_output=True)
    return(0)


def creat_background(peaks, length_of_site, counter):
    shuffled_peaks = []
    number_of_sites = 0
    while counter > number_of_sites:
        peak = random.choice(peaks)
        shuffled_peak = ''.join(random.sample(peak, len(peak)))
        shuffled_peaks.append(shuffled_peak)
        number_of_sites += (len(''.join(shuffled_peak)) - length_of_site + 1) * 2
    return(shuffled_peaks)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def bootstrap_inmode(peaks, length_of_site, counter, path_to_inmode, path_to_java, tmp_dir, order):
    true_scores = []
    false_scores = []
    number_of_peaks = len(peaks)
    for i in range(5):
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        train_peaks = random.choices(peaks, k=round(0.9 * number_of_peaks))
        test_peaks = [peak for peak in peaks  if not peak in train_peaks]
        shuffled_peaks = creat_background(test_peaks, length_of_site, counter / 5)
        write_fasta(train_peaks, tmp_dir, "train")
        write_fasta(test_peaks, tmp_dir, "test")
        write_fasta(shuffled_peaks, tmp_dir, "shuffled")
        make_inmode(path_to_inmode, path_to_java, length_of_site, order, tmp_dir)
        for true_score in true_scores_inmode(path_to_inmode, path_to_java, length_of_site, tmp_dir, "test"):
            true_scores.append(true_score)
        for false_score in false_scores_inmode(path_to_inmode, path_to_java, length_of_site, tmp_dir, "shuffled"):
            false_scores.append(false_score)
        shutil.rmtree(tmp_dir)
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
    parser.add_argument('inmode', action='store', help='path to InMoDe source')
    parser.add_argument('-j', '--java', action='store', type=str, dest='java',
                        required=False, default='java', help='Path to java')
    parser.add_argument('-t', '--tmp', action='store', type=str, dest='tmp',
                        required=False, default='./inmode.tmp', help='tmp directory')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())
    
    
def main():
    args = parse_args()
    peaks_path = args.fasta
    results_path = args.results
    path_to_inmode = args.inmode
    path_to_java = args.java
    length_of_site = args.length
    tmp_dir = args.tmp
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
        
    counter = 1000000
    order = 2
    peaks = read_peaks(peaks_path)
    table = bootstrap_inmode(peaks, length_of_site, counter, path_to_inmode, path_to_java, tmp_dir, order)
    write_table_bootstrap(results_path, table)
    return(0)
    
    
def bootstrap_for_inmode(peaks_path, results_path, length_of_site, path_to_inmode, path_to_java, tmp_dir, counter=5000000, order=2):
    peaks = read_peaks(peaks_path)
    table = bootstrap_inmode(peaks, length_of_site, counter, path_to_inmode, path_to_java, tmp_dir, order)
    write_table_bootstrap(results_path, table)
    return(0)


if __name__=="__main__":
    main()