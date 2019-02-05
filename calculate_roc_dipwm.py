import math
import csv
import sys
import os
import random
import itertools
import argparse
import functools
import multiprocessing as mp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st


def read_fasta(path):
    sequences = []
    with open(path, 'r') as file:
        sequences = [i.strip().upper() for i in file if i.strip()[0] != '>']
    return(sequences)


def background_freq_di(seq):
    s = ''.join(seq)
    background = {}
    di_nucleotides = itertools.product('ACGT', repeat=2)
    for i in di_nucleotides:
        background[''.join(i)] = s.count(''.join(i))
    sum_of_nuc = sum(background.values())
    for i in background.keys():
        background[i] = background[i]/sum_of_nuc
    return(background)


def make_dipfm_from_dipcm(pcm):
    '''
    Вычисление частотной матрицы на основе PCM.
    Для того чтобы избавиться от 0 значений частот используется pseudocount.
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

    number_of_sites = [0] * len(pcm['AA'])
    for key in pcm.keys():
        for i in range(len(pcm[key])):
            number_of_sites[i] += pcm[key][i]

    pfm = dict()
    di_nucleotides = itertools.product('ACGT', repeat=2)
    for i in di_nucleotides:
        pfm[''.join(i)] = []

    first_key = list(pcm.keys())[0]
    nuc_pseudo = 1/len(pcm.keys())
    for i in range(len(pcm[first_key])):
        for nuc in pcm.keys():
            pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
    return(pfm)


def make_dipwm_from_dipcm(pcm, background):
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
    di_nucleotides = itertools.product('ACGT', repeat=2)
    for i in di_nucleotides:
        pwm[''.join(i)] = []
    pfm = make_dipfm_from_dipcm(pcm)
    first_key = list(pcm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def make_dipcm(motifs):
    '''
    input - список мотивов одинаковой длины
    output -  PCM
    Создает PCM на основе списка мотивов
    '''
    matrix = {}
    di_nucleotides = itertools.product('ACGT', repeat=2)
    for i in di_nucleotides:
        matrix[''.join(i)] = []
    len_of_motif = len(motifs[0])
    for i in matrix.keys():
        matrix[i] = [0]*(len_of_motif - 1)

    for site in motifs:
        for i in range(len(site) - 1):
            matrix[site[i:i+2]][i] += 1
    return(matrix)


def make_dipwm(motifs, background):
    pcm = make_dipcm(motifs)
    pwm = make_dipwm_from_dipcm(pcm, background)
    return(pwm)


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


def score_dipwm(seq, dipwm):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и PWM
    тип diPWM
    '''
    length_of_seq = len(seq)
    score = 0
    for i in range(length_of_seq - 1):
        score += dipwm[seq[i:i+2]][i]
    return(score)


def to_norm(score, min_score, max_score):
    norm_value = (score - min_score) / (max_score - min_score)
    return(norm_value)


def get_min_score(pwm):
    '''
    Вичисляет минимальное значение score для матрицы
    '''

    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += min(tmp)
    return(value)


def get_max_score(pwm):
    '''
    Вичисляет минимальное значение score для матрицы
    '''

    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += max(tmp)
    return(value)


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


def get_site_scores(train_test, k):  # , background):
    '''
    Функция принимает три аргумента train - список сайтов для получения матрицы,
    test - списов истинных сайтов для подсчета score матрицы, чтобы в дальнейшем считать TPR и
    k - количество случайных сайтов на каждый истинный сайт для подсчета score матрицы,
    чтобы в дальнейшем считать FPR
                score истинных сайтов            score случайных сайтов
    output: [score_1, score_2, ... score_N], [score_1, score_2, ... score_L]
    '''
    train = train_test[0]
    test = train_test[1]
    background = background_freq_di(train)
    pwm = make_dipwm(train, background)
    #random_sites = [j for i in map(functools.partial(get_n_random_sites, n=k), train + test) for j in i]
    random_sites = get_n_random_sites(sites=train + [test], n=k)
    true_scores = score_dipwm(test, pwm)
    false_scores = list(map(functools.partial(score_dipwm, dipwm=pwm), random_sites))
    return(true_scores, false_scores)


def calculate_score_for_train_test(seq, k):
    train_list = []
    for index, test in enumerate(seq):
        train = seq[:]
        del train[index]
        train_list.append(train)
    test_list = seq[:]
    train_test = zip(train_list, test_list)

    with mp.Pool(mp.cpu_count()) as p:
        res = p.map(functools.partial(get_site_scores,
                                      k=k), train_test)
    return(res)


def roc(i, np_true_scores, np_false_scores):
    tpr = (len(np_true_scores) - np.searchsorted(np_true_scores, i, side='left')) / len(np_true_scores)
    fpr = (len(np_false_scores) - np.searchsorted(np_false_scores, i, side='left')) / len(np_false_scores)
    if fpr == 0.0:
        fpr = 1 / (2 * len(np_false_scores))
    tmp = {'tpr': tpr,
           'fpr': fpr,
           'score': i}
    return(tmp)


def closed_to_fpr(results_, fpr):
    df = pd.DataFrame(results_[np.logical_and(results_['fpr'] <=
                                              fpr + 0.0005, results_['fpr'] >= fpr - 0.0005)])
    df['delta'] = np.absolute(df['fpr'] - fpr)
    return(df.loc[df['delta'].idxmin()])


def closed_to_tpr(results_, tpr):
    df = pd.DataFrame(results_[np.logical_and(results_['tpr'] <=
                                              tpr + 0.01, results_['tpr'] >= tpr - 0.01)])
    # print(df)
    df['delta'] = np.absolute(df['tpr'] - tpr)
    return(df.loc[df['delta'].idxmin()])


def parse_args():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name')
    roc_parser = subparsers.add_parser('roc', help='return tsv file with roc [fpr, tpr, score]')
    roc_parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                            required=True, help='path to FASTA file with sites')
    roc_parser.add_argument('-n', '--tag', action='store', dest='tag',
                            required=True, help='file tag')
    roc_parser.add_argument('-t', '--times', action='store', dest='times', type=int,
                            default=100, required=False, help='x times random sample will be larger than original sample')
    roc_parser.add_argument('-o', '--output', action='store', dest='output',
                            required=True, help='dir for write output file')
    roc_parser.add_argument('-P', '--processes', action='store', type=int, dest='cpu_count',
    required=False, default=2, help='Number of processes to use, default: 2')

    tpr_parser = subparsers.add_parser('get_tpr', help='return tpr with fixed value of fpr')
    tpr_parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                            required=True, help='path to FASTA file with sites')
    tpr_parser.add_argument('-t', '--times', action='store', dest='times', type=int,
                            default=100, required=False, help='x times random sample will be larger than original sample')
    # tpr_parser.add_argument('-i', '--iterations', action='store', type=int, dest='iterations',
    #                        default=10, required=False, help='number of iterations')
    tpr_parser.add_argument('-v', '--fpr', action='store', type=float, dest='fpr',
                            required=True, help='value of FPR')
    tpr_parser.add_argument('-n', '--tag', action='store', dest='tag',
                            required=True, help='file tag')
    tpr_parser.add_argument('-o', '--output', action='store', dest='output',
                            required=True, help='dir for write output file')
    tpr_parser.add_argument('-P', '--processes', action='store', type=int, dest='cpu_count',
    required=False, default=2, help='Number of processes to use, default: 2')

    tpr_parser = subparsers.add_parser('get_fpr', help='return fpr with fixed value of tpr')
    tpr_parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                            required=True, help='path to FASTA file with sites')
    tpr_parser.add_argument('-t', '--times', action='store', dest='times', type=int,
                            default=100, required=False, help='x times random sample will be larger than original sample')
    # tpr_parser.add_argument('-i', '--iterations', action='store', type=int, dest='iterations',
    #                        default=10, required=False, help='number of iterations')
    tpr_parser.add_argument('-v', '--tpr', action='store', type=float, dest='tpr',
                            required=True, help='value of TPR')
    fpr_parser.add_argument('-n', '--tag', action='store', dest='tag',
                            required=True, help='file tag')
    fpr_parser.add_argument('-o', '--output', action='store', dest='output',
                            required=True, help='dir for write output file')
    fpr_parser.add_argument('-P', '--processes', action='store', type=int, dest='cpu_count',
    required=False, default=2, help='Number of processes to use, default: 2')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():

    args = parse_args()
    if args.subparser_name == 'get_tpr':
        #path = '/home/anton/DATA/TF/FOXA1_hg38_ENCSR819LGH/MOTIFS/FOXA1_14.fasta'
        path = args.input_fasta
        fpr = args.fpr
        #iterations = args.iterations
        times = args.times
        cpu_count = args.cpu_count

        seq = read_fasta(path)
        background = background_freq_di(seq)
        seq = remove_equalent_seq(seq, homology=0.95)
        pwm = make_dipwm(seq, background)
        min_score = get_min_score(pwm)
        max_score = get_max_score(pwm)

        #tpr_list = list()
        # for i in range(iterations):
        scores = calculate_score_for_train_test(seq, k=times)
        true_scores = [i for i, j in scores]
        false_scores = [k for i, j in scores for k in j]

        norm_true_scores = np.array([to_norm(score, min_score, max_score)
                                     for score in true_scores])
        norm_false_scores = np.array([to_norm(score, min_score, max_score)
                                      for score in false_scores])
        norm_true_scores.sort()
        norm_false_scores.sort()

        with mp.Pool(cpu_count) as p:
            results = p.map(functools.partial(roc, np_true_scores=norm_true_scores,
                                              np_false_scores=norm_false_scores), norm_true_scores)
        results = pd.DataFrame(results)
        results = results[['tpr', 'fpr', 'score']]

        if not os.path.isdir(output):
            os.mkdir(output)

        results.to_csv(output + '/' + tag +
                       '.tsv', sep='\t', header=True, index=False)

        tpr = closed_to_fpr(results, fpr)['tpr']
        # tpr_list.append(tpr)
        print(tpr)
        # for i in tpr_list:
        #    print(i)

    elif args.subparser_name == 'get_fpr':
        #path = '/home/anton/DATA/TF/FOXA1_hg38_ENCSR819LGH/MOTIFS/FOXA1_14.fasta'
        path = args.input_fasta
        tpr = args.tpr
        #iterations = args.iterations
        times = args.times
        cpu_count = args.cpu_count

        seq = read_fasta(path)
        background = background_freq_di(seq)
        seq = remove_equalent_seq(seq, homology=0.95)
        pwm = make_dipwm(seq, background)
        min_score = get_min_score(pwm)
        max_score = get_max_score(pwm)

        #fpr_list = list()
        # for i in range(iterations):
        scores = calculate_score_for_train_test(seq, k=times)
        true_scores = [i for i, j in scores]
        false_scores = [k for i, j in scores for k in j]

        norm_true_scores = np.array([to_norm(score, min_score, max_score)
                                     for score in true_scores])
        norm_false_scores = np.array([to_norm(score, min_score, max_score)
                                      for score in false_scores])
        norm_true_scores.sort()
        norm_false_scores.sort()

        with mp.Pool(cpu_count) as p:
            results = p.map(functools.partial(roc, np_true_scores=norm_true_scores,
                                              np_false_scores=norm_false_scores), norm_true_scores)
        results = pd.DataFrame(results)
        # print(results)
        results = results[['tpr', 'fpr', 'score']]

        if not os.path.isdir(output):
            os.mkdir(output)

        results.to_csv(output + '/' + tag +
                       '.tsv', sep='\t', header=True, index=False)

        fpr = closed_to_tpr(results, tpr)['fpr']
        # fpr_list.append(fpr)
        print(fpr)
        # for i in fpr_list:
        #    print(i)

    elif args.subparser_name == 'roc':
        #path = '/home/anton/DATA/TF/FOXA1_hg38_ENCSR819LGH/MOTIFS/FOXA1_14.fasta'
        path = args.input_fasta
        tag = args.tag
        output = args.output
        times = args.times
        cpu_count = args.cpu_count

        seq = read_fasta(path)
        background = background_freq_di(seq)
        seq = remove_equalent_seq(seq, homology=0.95)
        pwm = make_dipwm(seq, background)
        min_score = get_min_score(pwm)
        max_score = get_max_score(pwm)

        scores = calculate_score_for_train_test(seq, k=times)
        true_scores = [i for i, j in scores]
        false_scores = [k for i, j in scores for k in j]

        norm_true_scores = np.array([to_norm(score, min_score, max_score) for score in true_scores])
        norm_false_scores = np.array([to_norm(score, min_score, max_score)
                                      for score in false_scores])
        norm_true_scores.sort()
        norm_false_scores.sort()

        with mp.Pool(cpu_count) as p:
            results = p.map(functools.partial(roc, np_true_scores=norm_true_scores,
                                              np_false_scores=norm_false_scores), norm_true_scores)
        results = pd.DataFrame(results)
        results = results[['tpr', 'fpr', 'score']]

        if not os.path.isdir(output):
            os.mkdir(output)

        results.to_csv(output + '/' + tag + '.tsv', sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()
