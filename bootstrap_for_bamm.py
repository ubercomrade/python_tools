import subprocess
import csv
import itertools
import os
import sys
import argparse
import random
import math
import shutil


def read_log_odds_zoops(path):
    container = []
    append = container.append
    with open(path) as file:
        file.readline()
        for line in file:
            site = line.split()[5]
            if not 'N' in site:
                append(site)
            else:
                continue
    return(container)


def read_sites(path):
    container = []
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                append(line.strip().upper())
            else:
                continue
    return(container)


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


def read_bamm(bamm_path, bg_path):
    bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    container = dict()
    log_odds_bamm = make_log_odds_bamm(bamm, bg)
    for i in range(order + 1):
        k_mers = make_k_mers(i)
        bamm_dict = bamm_to_dict(log_odds_bamm, i, k_mers)
        container.update(bamm_dict)
    return(container, order)


def score(site, bamm, order, length_of_site):
    score = float()
    for index in range(order):
        score += bamm[site[0:index + 1]][index]
    for index in range(length_of_site - order):
        score += bamm[site[index:index+order + 1]][index + order]
    return(score, site)


def calculate_scores(sites, bamm, order, length_of_site):
    scores = [score(site, bamm, oreder, length_of_site) for site in sites]
    return(scores)


def create_bamm_model(directory, order):
    fasta_path = directory + '/train.fasta'
    sites_path = directory + '/sites.txt'
    args = ['BaMMmotif', directory, fasta_path, '--bindingSiteFile', sites_path, '--EM', '--order', order, '--Order', order]
    subprocess.call(args)
    bamm_path = directory + '/train_motif_1.ihbcp'
    bg_path = directory + '/train.hbcp'
    bamm = read_bamm(bamm_path, bg_path)
    return(bamm)


def create_train_fasta_and_sites(sites, directory):
    with open(directory + '/train.fasta') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    with open(directory + '/sites.txt') as file:
        for site in sites:
            file.write('{}\n'.format(site))
    return(0)


def bootstrap_bamm(sites, size_of_random_sample, order, directory):
    true_scores = []
    false_scores = []
    number_of_sites = len(sites)
    len_of_site = len(sites[0])

    for i in range(10):
        train_sample = random.choices(sites, k=round(0.9 * number_of_sites))
        test_sample = [site for site in sites  if not site in train_sample]
        random_sample = [''.join(random.sample(list(random.choice(test_sample)), len_of_site)) for i in range(len(test_sample) * size_of_random_sample)]
        create_train_fasta_and_sites(sites[index_train], directory)
        bamm = create_bamm_model(directory)
        for true_score in calculate_scores(test_sample, bamm, order, length_of_site):
            true_scores.append(true_score)
        for false_score in calculate_scores(random_sample, bamm, order, length_of_site):
            false_scores.append(false_score)
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)

    table = []
    for tpr in [round(i * 0.01, 2) for i in range(0,105, 5)]:
        print(tpr)
        score = true_scores[round(true_length * tpr) - 1]
        actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
        fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
        table.append({'Scores': score, 'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
    print(table)
    return(table)


def write_table(path, data):
    with open(path, 'w') as csvfile:
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
    parser.add_argument('-s', '--size', action='store', type=int, dest='size',
                            default= 1000, required=False, help='randome_sample times more than test_sample')
    parser.add_argument('-o', '--order', action='store', type=int, dest='order',
                            default= 2, required=False, help='order of BaMM model')
    parser.add_argument('-t', '--tmp', action='store', dest='tmp_dir',
                            default= './tmp', required=False, help='dir for writing tmp files')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    path = args.input
    out = args.output
    size_of_random_sample = args.size
    order = args.order
    tmp_dir = args.tmp_dir
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    sites = read_log_odds_zoops(path)
    table = bootstrap_bamm(sites, size_of_random_sample, order, tmp_dir)
    write_table(out, table)
    shutil.rmtree(tmp_dir)
    return(0)


if __name__=='__main__':
    main()

