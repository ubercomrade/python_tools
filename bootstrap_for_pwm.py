import csv
import itertools
import os
import sys
import argparse
import random

def read_sites(path):
    container = []
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                append(line.strip().upper())
    return(container)


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
    first_key = list(pcm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def motifs_to_pwm(motifs):
    pcm = make_pcm(motifs)
    pfm = make_pfm(pcm)
    pwm = make_pwm(pfm)
    return(pwm)


def score(seq, pwm):
    score = 0
    for position, letter in enumerate(seq):
        score += pwm[letter][position]
    return(score)


def calculate_scores(sites, pwm):
    scores = [score(pwm, site) for site in sites]
    return(scores)


def bootstrap_pwm(sites, size_of_random_sample):
    true_scores = []
    false_scores = []
    number_of_sites = len(sites)

    for i in range(10):
        train_sample = random.choices(sites, k=round(0.9 * number_of_sites))
        test_sample = [site for site in sites  if not site in train_sample]
        random_sample = [random.shuffle(random.choice(test_sample)) for i in range(len(test_sample) * size_of_random_sample)]
        pwm = make_pwm_from_sites(train_sample)
        for true_score in calculate_scores(test_sample, pwm):
            true_scores.append(true_score)
        for false_score in calculate_scores(randome_sample, pwm):
            false_scores.append(false_score)

    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)

    table = []
    for tpr in [round(i * 0.01, 2) for i in range(0,105, 5)]:
        score = true_scores[round(true_length * tpr)]
        actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
        fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
        table.append({'Scores': score, 'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
    return(table)


def write_table(path, data):
    with open('names.csv', 'w', newline='') as csvfile:
        fieldnames = data[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for line in data:
            writer.writerow(line)
    return(0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', action='store', help='path to write results')
    parser.add_argument('input', action='store', help='path to file with sites')
    parser.add_argument('-s', '--size', action='store', type=int, dest='false_positive',
                            default= 1000000, required=False, help='randome_sample times more than test_sample')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    path = args.input
    out = args.output
    size_of_random_sample = args.size
    sites = read_sites(path)
    table = bootstrap_pwm(sites, size_of_random_sample)
    write_table(out, table)
    return(0)

if __name__=='__main__':
    main()
