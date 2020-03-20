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

import math
import sys
import os
import random
import itertools
import argparse
import functools
import bisect
from operator import itemgetter
import re
import subprocess
from math import log


def parse_bamm_occurence(path):
    with open(path, 'r') as file:
        output = []
        file.readline()
        for line in file:
            d = {'name': int(), 'start': int(), 'end': int(),
                 'seq': str(), 'strand': str()}
            line = line.strip().split()

            d['name'] = int(line[0].split('::')[0].split('_')[1])
            d['start'] = int(line[3].split('..')[0])
            d['end'] = int(line[3].split('..')[1])
            d['seq'] = line[2]
            d['strand'] = line[2]
            output.append(d)
    return(output)


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
        print('File {0} does not exist'.format(bamm_file))
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


def create_bamm_model(tmp_dir, sites, background, model_order):
    #sites = remove_equalent_seq(sites, homology=0.95)
    pcm = make_pcm(sites)
    pfm = make_pfm_from_pcm(pcm)

    if not os.path.isdir(tmp_dir):
        os.mkdir(tmp_dir)

    nsites = len(sites)
    write_sites(tmp_dir, 'train', sites)
    write_meme(tmp_dir, 'train', pfm, background, nsites)

    # args = ['BaMMmotif', tmp_dir,
    #         tmp_dir + '/train.fasta',
    #        '--PWMFile', tmp_dir + '/train.meme',
    #         '--basename', 'train',
    #        '--EM',
    #        '--Order', str(model_order),
    #        '--order', str(model_order)]
    # r = subprocess.call(args)

    os.system('BaMMmotif {0} {1} --PWMFile {2} --basename train --EM --Order {3} --order {3} > {4}'.format(tmp_dir,
        tmp_dir + '/train.fasta', tmp_dir + '/train.meme', model_order, tmp_dir + '/logs.txt'))

    log_odds_bamm, order = prepare_bamm(tmp_dir + '/train_motif_1.ihbcp', tmp_dir + '/train.hbcp')
    return(log_odds_bamm, order)


def get_motifs_from_fasta(fasta, positions):
    motifs = []
    for line in positions:
        if line['strand'] == '+':
            start = line['start']
            end = line['end']
            if start < 0:
                continue
            elif end >= len(fasta[line['name']]):
                continue
            else:
                motifs.append(fasta[line['name']][start:end])
        else:
            start = len(fasta[line['name']]) - line['end']
            end = len(fasta[line['name']]) - line['start']
            if start < 0:
                continue
            elif end >= len(fasta[line['name']]):
                continue
            else:
                motifs.append(complement(fasta[line['name']])[start:end])
    motifs = [i for i in motifs if not 'N' in i]
    return(motifs)


def read_fasta(path):
    sequences = []
    with open(path, 'r') as file:
        sequences = [i.strip().upper() for i in file if i.strip()[0] != '>']
    return(sequences)


def complement(seq):
    '''
    Make reverse and compelent
    '''
    out_seq = str()
    for letter in seq:
        if letter == 'A':
            out_seq += 'T'
        elif letter == 'C':
            out_seq += 'G'
        elif letter == 'G':
            out_seq += 'C'
        elif letter == 'T':
            out_seq += 'A'
        elif letter == 'N':
            out_seq += 'N'
    out_seq = out_seq[::-1]
    return(out_seq)


# def complement(record):
#     output = dict(record)
#     strand = record['strand']
#     seq = str()
#     if strand == '+':
#         output['strand'] = '-'
#     else:
#         output['strand'] = '+'

#     seq = output['seq'].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
#     output['seq'] = seq
#     return(output)


def check_nucleotides(site):
    s = set(site)
    n = {'A', 'C', 'G', 'T'}
    if len(s - n) == 0:
        return(True)
    else:
        return(False)


def remove_equalent_seq(seq_list, homology=0.95):
    '''
    Удаление гомологичных последовательностей из списка (seq_list)
    Если кол-во совпадений при сравнении последовательности 1 и 2 >= длина последовательности * homology,
    то последовательность 1 удаляется из списка
    Функция возвращает новый список
    '''
    seq_list = list(seq_list)
    treshold = homology * len(seq_list[0])
    for seq1 in tuple(seq_list):
        sub_seq_list = list(seq_list)
        sub_seq_list.remove(seq1)
        for seq2 in sub_seq_list:
            score = len([i for i, j in zip(seq1, seq2) if i == j])
            if score >= treshold:
                seq_list.remove(seq1)
                break
    return(seq_list)


def score_bamm(site, bamm, order, length_of_site):
    score = float()
    for index in range(order):
        score += bamm[site[0:index + 1]][index]
    for index in range(length_of_site - order):
        score += bamm[site[index:index+order + 1]][index + order]
    return(score)


def make_pfm_from_pcm(pcm):
    '''
    Вычисление частотной матрицы на основе PCM.
    Для того чтобы избавиться от 0 значений частот используется pseudocount.
    Pseudocount может быть dict со стандартными значениями {'A':0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25} [1],
    либо pseudocount может быть str со значением sqroot [2].
    Подробнее о расчетах смотри:
    1)Wyeth W.Wasserman and Albin Sandelin
      APPLIED BIOINFORMATICS FOR THE IDENTIFICATION OF REGULATORY ELEMENTS
      doi:10.1038/nrg1315

    2)Victor G Levitsky
      Effective transcription factor binding site prediction using a
      combination of optimization, a genetic algorithm and discriminant
      analysis to capture distant interactions
      doi:10.1186/1471-2105-8-481

    В любых других условиях функция ничего не возвращает

    '''

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


def make_pwm_from_pcm(pcm, background, method='log-odds', pseudocount='1/N'):
    '''
    Функиця, которая считает PWM (position weight matrix) на основе PCM (position count matrix)
    с преобразованием log-odds (добавить новые)

    Ref:
    1)Wyeth W.Wasserman and Albin Sandelin
      APPLIED BIOINFORMATICS FOR THE IDENTIFICATION OF REGULATORY ELEMENTS
      doi:10.1038/nrg1315

    2)Victor G Levitsky
      Effective transcription factor binding site prediction using a
      combination of optimization, a genetic algorithm and discriminant
      analysis to capture distant interactions
      doi:10.1186/1471-2105-8-481

    3)Oliver D. King and Frederick P. Roth
      A non-parametric model for transcription factor binding sites
      doi: 10.1093/nar/gng117

    '''
    pwm = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pwm[i[0]] = []

    pfm = make_pfm_from_pcm(pcm, pseudocount)
    first_key = list(pcm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def make_pcm(motifs):
    '''
    input - список мотивов одинаковой длины
    output -  PCM
    Создает PCM на основе списка мотивов
    '''
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


def make_pwm(motifs, background):
    pcm = make_pcm(motifs)
    pwm = make_pwm_from_pcm(pcm, background, method='log-odds', pseudocount='1/N')
    return(pwm)


def random_seq(seq, k):
    '''
    Генерирует k случайных последовательностей с таким же нуклеотидным составом как и seq
    Возвращает list()
    '''
    out = list()
    for i in range(k):
        out.append(''.join(random.sample(seq, k=len(seq))))
    return(out)


def get_n_random_sites(sites, n):
    sites = random.choices(sites, k=n)
    out = [''.join(random.sample(site, k=len(site))) for site in sites]
    return(out)


# def get_site_scores(train_test, tmp_dir, background, model_order, k):
#     '''
#     Функция принимает три аргумента train - список сайтов для получения матрицы,
#     test - списов истинных сайтов для подсчета score матрицы, чтобы в дальнейшем считать TPR и
#     k - количество случайных сайтов на каждый истинный сайт для подсчета score матрицы,
#     чтобы в дальнейшем считать FPR
#                 score истинных сайтов            score случайных сайтов
#     output: [score_1, score_2, ... score_N], [score_1, score_2, ... score_L]
#     '''
#     train = train_test[0]
#     test = train_test[1]
#     log_odds_bamm, order = create_bamm_model(tmp_dir, train, background, model_order)
#     random_sites = get_n_random_sites(sites=train + [test], n=k)
#     length_of_site = len(train[0])
#     true_scores = score_bamm(test, log_odds_bamm, order, length_of_site)
#     false_scores = list(map(functools.partial(score_bamm, bamm=log_odds_bamm, order=order,
#         length_of_site=length_of_site), random_sites))
#     return(true_scores, false_scores)


# def calculate_score_for_train_test(seq, tmp_dir, model_order, background, k):
#     train_list = []
#     for index, test in enumerate(seq):
#         train = seq[:]
#         del train[index]
#         train_list.append(train)
#     test_list = seq[:]
#     train_test = zip(train_list, test_list)
#     res = list(map(functools.partial(get_site_scores,
#                                       tmp_dir=tmp_dir, model_order=model_order,
#                                       background=background, k=k), train_test))
#     return(res)


def get_site_scores(train_test, tmp_dir, background, model_order, times):
    '''
    Функция принимает три аргумента train - список сайтов для получения матрицы,
    test - списов истинных сайтов для подсчета score матрицы, чтобы в дальнейшем считать TPR и
    k - количество случайных сайтов на каждый истинный сайт для подсчета score матрицы,
    чтобы в дальнейшем считать FPR
                score истинных сайтов            score случайных сайтов
    output: [score_1, score_2, ... score_N], [score_1, score_2, ... score_L]
    '''
    train_sites = train_test[0]
    test_sites = train_test[1]
    #print(len(train_sites+test_sites))
    log_odds_bamm, order = create_bamm_model(tmp_dir, train_sites, background, model_order)
    random_sites = get_n_random_sites(sites=train_sites+test_sites, n=times)
    length_of_site = len(train_sites[0])
    print(times)
    true_scores = list(map(functools.partial(score_bamm, bamm=log_odds_bamm, order=order,
        length_of_site=length_of_site), test_sites)) 
    false_scores = list(map(functools.partial(score_bamm, bamm=log_odds_bamm, order=order,
        length_of_site=length_of_site), random_sites))
    return(true_scores, false_scores)


def calculate_score_for_train_test(seq, tmp_dir, model_order, background, times):
    train_list = []
    test_list = []
    step = round(len(seq) * 1 / 10)
    for k in range(10):
        start = k * step 
        stop = k * step + step
        train_list.append(seq[start:stop])
        test_list.append(seq[:start] + seq[stop:])
    train_test = zip(train_list, test_list)
    res = list(map(functools.partial(get_site_scores,
                                      tmp_dir=tmp_dir, model_order=model_order,
                                      background=background, times=times), train_test))
    return(res)


def roc(i, true_scores, false_scores):
    tpr = (len(true_scores) - bisect.bisect_left(true_scores, i)) / len(true_scores)
    fpr = (len(false_scores) - bisect.bisect_left(false_scores, i)) / len(false_scores)
    if fpr == 0.0:
        fpr = 1 / (2 * len(false_scores))
    tmp = {'tpr': tpr,
           'fpr': fpr,
           'score': i}
    return(tmp)


def calculate_fpr(tmp_dir, sites, background, times, model_order=2, tpr=0.5):

    #log_odds_bamm, order = create_bamm_model(tmp_dir, sites, background, model_order)
    scores = calculate_score_for_train_test(sites, tmp_dir, model_order, background, times=times)
    #true_scores = [i for i, j in scores]
    true_scores = [k for i, j in scores for k in i]
    false_scores = [k for i, j in scores for k in j]
    #print(len(true_scores), len(false_scores))

    true_scores.sort()
    false_scores.sort()

    results = list(map(functools.partial(roc, true_scores=true_scores,
                                          false_scores=false_scores), true_scores))
    
    results.sort(key=itemgetter('tpr'))
    tpr_list = [i['tpr'] for i in results]
    fpr = float(results[bisect.bisect(tpr_list, tpr)]['fpr'])
    return(fpr)


def write_meme(output, tag, pfm, background, nsites):
    with open(output + '/' + tag + '.meme', 'w') as file:
        file.write('MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies\n')
        file.write('A {0} C {1} G {2} T {3}\n\n'.format(background['A'], background['C'],
                                                        background['G'], background['T']))
        file.write('MOTIF {0}\n'.format(tag))
        file.write(
            'letter-probability matrix: alength= 4 w= {0} nsites= {1}\n'.format(len(pfm['A']), nsites))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def write_sites(output, tag, sites):
    with open(output + '/' + tag + '.fasta', 'w') as file:
        for index, site in enumerate(sites):
            file.write('>site_' + str(index) + '\n')
            file.write(site + '\n')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', action='store',
                        help='path to ChIPMunk results')
    parser.add_argument('fasta', action='store',
                        help='path to fasta file (peaks in fasta format)')
    parser.add_argument('output', action='store',
                        help='dir to write output files')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    parser.add_argument('-T', '--tmp', action='store', dest='tmp',
                        required=False, default='./tmp', help='tmp dir')
    parser.add_argument('-n', '--times', action='store', dest='times', type=int,
                        default=500000, required=False,
                        help='x times random sample will be larger than original sample')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():

    args = parse_args()
    bamm_sites_path = args.input
    fasta_path = args.fasta
    output_dir = args.output
    times = args.times
    tag = args.tag
    tmp_dir = args.tmp
    model_order = 3

    bamm_sites = parse_bamm_occurence(bamm_sites_path)
    fasta = read_fasta(fasta_path)
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
    fprs = []

    #log_odds_bamm, order = prepare_bamm(bamm_path, bg_path)

    optimal_motifs = get_motifs_from_fasta(fasta, bamm_sites)
    optimal_motifs = remove_equalent_seq(optimal_motifs, homology=0.95)
    optimal_fpr = calculate_fpr(tmp_dir, optimal_motifs, background, times, model_order=model_order, tpr=0.5)
    fprs.append(optimal_fpr)
    print(optimal_fpr)

    for i in bamm_sites:
        i['start'] -= 1
        i['end'] += 1
    motifs = get_motifs_from_fasta(fasta, bamm_sites)
    motifs = remove_equalent_seq(motifs, homology=0.95)
    print(len(motifs[0]))
    fpr = calculate_fpr(tmp_dir, motifs, background, times, model_order=model_order, tpr=0.5)
    fprs.append(fpr)
    print(fpr)

    while fpr <= optimal_fpr:
        optimal_fpr = fpr
        optimal_motifs = motifs.copy()
        for i in bamm_sites:
            i['start'] -= 1
            i['end'] += 1
        motifs = get_motifs_from_fasta(fasta, bamm_sites)
        motifs = remove_equalent_seq(motifs, homology=0.95)
        print(len(motifs[0]))
        fpr = calculate_fpr(tmp_dir, motifs, background, times, model_order=model_order, tpr=0.5)
        fprs.append(fpr)
        if fpr <= optimal_fpr and 1 - fprs[-1]/fprs[-2] < 0.1:
            print(fpr)
            break
        print(fpr)


    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    pcm = make_pcm(motifs)
    pfm = make_pfm_from_pcm(pcm)
    write_sites(output_dir, 'train', motifs)
    write_meme(output_dir, 'train', pfm, background, len(motifs))

    args = ['BaMMmotif', output_dir,
            output_dir + '/train.fasta',
           '--PWMFile', output_dir + '/train.meme',
            '--basename', 'oBaMM',
           '--EM',
           '--Order', str(model_order),
           '--order', str(model_order)]
    r = subprocess.call(args)

    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    nsites = len(motifs)
    write_sites(output=output_dir, tag=tag, sites=motifs)
    with open(output_dir + '/fprs.txt', 'w') as file:
        for i in fprs:
            file.write(str(i) + '\n')
    file.close()

if __name__=='__main__':
   main()
