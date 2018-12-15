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


def background_freq(seq):
    s = ''.join(seq)
    background = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        background[i[0]] = s.count(i[0])
    sum_of_nuc = sum(background.values())
    for i in background.keys():
        background[i] = background[i]/sum_of_nuc
    return(background)


def make_pfm_from_pcm(pcm, pseudocount='1/N'):
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

    if pseudocount == '1/N':
        first_key = list(pcm.keys())[0]
        nuc_pseudo = 1/len(pcm.keys())
        for i in range(len(pcm[first_key])):
            for nuc in pcm.keys():
                pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
        return(pfm)

    elif pseudocount == 'sqroot':
        total_sq_root = int()
        for i in pcm.keys():
            total_sq_root += pcm[i][0]
        total_sq_root = math.sqrt(total_sq_root)
        sq_root = total_sq_root/len(pcm.keys())

        first_key = list(pcm.keys())[0]
        for i in range(len(pcm[first_key])):
            for nuc in pcm.keys():
                pfm[nuc].append((pcm[nuc][i] + sq_root) / (number_of_sites[i] + total_sq_root))

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


def score(seq, pwm):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и PWM
    '''
    length_of_seq = len(seq)
    position = 0
    score = 0
    for letter in seq:
        score += pwm[letter][position]
        position += 1
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
    background = background_freq(train)
    pwm = make_pwm(train, background)
    #random_sites = [j for i in map(functools.partial(get_n_random_sites, n=k), train + test) for j in i]
    random_sites = get_n_random_sites(sites=train + [test], n=k)
    true_scores = score(test, pwm)
    false_scores = list(map(functools.partial(score, pwm=pwm), random_sites))
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


def main():
    path = '/home/anton/DATA/TF/FOXA1_hg38_ENCSR819LGH/MOTIFS/FOXA1_14.fasta'
    seq = read_fasta(path)
    background = background_freq(seq)
    seq = remove_equalent_seq(seq, homology=0.95)
    pwm = make_pwm(seq, background)
    min_score = get_min_score(pwm)
    max_score = get_max_score(pwm)
    fpr = 0.0005

    tpr_list = list()
    for s in range(20):
        scores = calculate_score_for_train_test(seq, k=400)
        true_scores = [i for i, j in scores]
        false_scores = [k for i, j in scores for k in j]

        norm_true_scores = np.array([to_norm(score, min_score, max_score) for score in true_scores])
        norm_false_scores = np.array([to_norm(score, min_score, max_score)
                                      for score in false_scores])
        norm_true_scores.sort()
        norm_false_scores.sort()

        with mp.Pool(mp.cpu_count()) as p:
            results = p.map(functools.partial(roc, np_true_scores=norm_true_scores,
                                              np_false_scores=norm_false_scores), norm_true_scores)
        results = pd.DataFrame(results)

        tpr = closed_to_fpr(results, fpr)['tpr']
        tpr_list.append(tpr)
    return(tpr_list)


if __name__ == '__main__':
    res = main()