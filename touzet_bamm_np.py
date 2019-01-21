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


import os
import numpy as np
import itertools
import random
import multiprocessing as mp
import functools
import pandas as pd
import argparse
import math
import sys
import argparse
import time


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
        file.close()

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
        file.close()

        # convert a bamm array to numpy array
        for k in range(motif_order):
            model[k] = np.array(model[k], dtype=float)
    else:
        print('File {0} does not exist'.format(bg_file))
        exit()

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
                bg[order] = np.array(bg_freq, dtype=float)
                order += 1
        file.close()
    else:
        print('File {0} does not exist'.format(bg_file))
        exit()
    return model, bg, order-1


def make_log_odds_bamm(bamm, bg):
    log_odds_bamm = dict()
    for order in bamm.keys():
        log_odds_bamm[order] = np.array(np.log2(bamm[order] / bg[order]))
    return(log_odds_bamm)


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


def bamm_to_dict(log_odds_bamm, order, k_mers):
    bamm_dict = {}
    for k_mer in k_mers:
        bamm_dict[k_mer] = list()

    for i in range(len(log_odds_bamm[order])):
        # print(i)
        for index, k_mer in enumerate(k_mers):
            # print(index)
            bamm_dict[k_mer].append(log_odds_bamm[order][i][index])
    return(bamm_dict)


def to_score(norm_value, min_score, max_score):
    '''
    norm = (score - min) / (max - min) -> score = norm * (max - min) + min
    '''
    score = norm_value * (max_score - min_score) + min_score
    return(score)


def to_norm(score, min_score, max_score):
    norm_value = (score - min_score) / (max_score - min_score)
    return(norm_value)


def min_score(pwm):
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


def max_score(pwm):
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


def slice_matrix(matrix, start, end):
    '''
    Возвращает срез матрицы
    '''
    sliced_matrix = dict()
    for key in matrix.keys():
        sliced_matrix[key] = matrix[key][start:end]
    return(sliced_matrix)


def possible_scores(pwm):
    scores = [[]]
    letters = list(pwm.keys())
    for letter in letters:
        scores[0].append(pwm[letter][0])
    for i in range(1, len(pwm[letters[0]])):
        tmp = []
        for letter in letters:
            for j in scores[i - 1]:
                tmp.append(j + pwm[letter][i])
        scores.append(tmp)
    return(scores[-1])


def score_distribution_np(matrix, alpha, beta, background, ndigits=20):
    matrix_length = len(matrix['AAA'])
    Q_old = [1]
    possible_scores_old = [0]
    for position in range(matrix_length):
        Q_new = []
        possible_scores_new = []
        matrix_slice = [matrix[i][position] for i in matrix.keys()]
        t = np.array([score + weight for score in possible_scores_old for weight in matrix_slice])
        t = np.around(t, decimals=ndigits)
        bs = max_score(slice_matrix(matrix, position + 1, matrix_length))
        ws = min_score(slice_matrix(matrix, position + 1, matrix_length))
        conditions = np.logical_and(alpha - bs <= t, t <= beta - ws)
        Q_new = np.array([Q * background[letter]
                          for Q in Q_old for letter in background.keys()])[conditions]
        possible_scores_new = t[conditions]
        Q_old = np.array(Q_new)
        possible_scores_old = np.array(possible_scores_new)
    return(Q_old, possible_scores_old)


def distribution_np_to_dict(Q, scores):
    keys = np.unique(scores)
    Q_out = dict()
    for key in keys:
        Q_out[key] = np.sum(Q[scores == key])
    return(Q_out)


def fast_pvalue_np(matrix, alpha, background, ndigits=20):
    matrix_length = len(matrix['AAA'])
    Q_old = [1]
    possible_scores_old = [0]
    P = 0.0
    for position in range(matrix_length):
        Q_new = []
        possible_scores_new = []
        matrix_slice = [matrix[i][position] for i in matrix.keys()]
        t = np.array([score + weight for score in possible_scores_old for weight in matrix_slice])
        t = np.around(t, decimals=ndigits)
        bs = max_score(slice_matrix(matrix, position + 1, matrix_length))
        ws = min_score(slice_matrix(matrix, position + 1, matrix_length))
        conditions_1 = (alpha - ws) <= t
        P += np.sum(np.array([Q * background[letter]
                              for Q in Q_old for letter in background.keys()])[conditions_1])
        conditions_2 = np.logical_and(np.logical_not(alpha - ws <= t), alpha - bs <= t)
        Q_new = np.array([Q * background[letter]
                          for Q in Q_old for letter in background.keys()])[conditions_2]
        possible_scores_new = t[conditions_2]
        Q_old = np.array(Q_new)
        possible_scores_old = np.array(possible_scores_new)
    return(P)


def pwm_shuffling(pwm):
    shuffled_pwm = dict()
    score_lists = list(zip(*pwm.values()))
    random.shuffle(score_lists)
    score_lists = list(zip(*score_lists))
    keys = list(pwm.keys())
    for i in range(len(score_lists)):
        shuffled_pwm[keys[i]] = score_lists[i]
    return(shuffled_pwm)


def round_pwm(pwm, ndigits):
    keys = list(pwm.keys())
    output_pwm = {}
    for key in keys:
        output_pwm[key] = [round(i, ndigits=ndigits) for i in pwm[key]]
    return(output_pwm)


def maximal_error_associated(pwm, rounded_pwm):
    E = 0
    letters = list(pwm.keys())
    length = len(pwm[letters[0]])
    for i in range(length):
        tmp = []
        for letter in letters:
            tmp.append(abs(abs(pwm[letter][i]) - abs(rounded_pwm[letter][i])))
        E += max(tmp)
    return(E)


def search_score(alpha, rounded_pwm, Q, error, ndigits, background):
    output = []
    tmp_pvalue = fast_pvalue_np(rounded_pwm, alpha+error, ndigits=ndigits, background=background)
    score_list = list(Q.keys())
    score_list.sort(reverse=False)
    if len(score_list) == 1:
        return(score_list[0])
    elif len(score_list) == 0:
        return(alpha)
    else:
        for delta in score_list:
            cumulative_Q = 0
            for i in score_list:
                if i < delta + error and i > delta:
                    cumulative_Q += Q[i]
            if cumulative_Q + fast_pvalue_np(rounded_pwm, delta + error, ndigits, background=background) >= tmp_pvalue:
                return(delta)
        return(alpha)


def from_pvalue_to_score(pwm, pvalue, ndigits, background, digit=1):
    rounded_pwm = round_pwm(pwm, ndigits)
    error = maximal_error_associated(pwm, rounded_pwm)
    Q, scores = score_distribution_np(rounded_pwm, alpha=min_score(
        rounded_pwm), beta=max_score(rounded_pwm), ndigits=ndigits, background=background)
    Q = distribution_np_to_dict(Q, scores)
    score_list = list(Q.keys())
    score_list.sort(reverse=True)
    alpha = score_list[0]
    cumulative_Q = Q[alpha]
    step = 1
    while not cumulative_Q >= pvalue:
        alpha = score_list[step]
        cumulative_Q += Q[alpha]
        step += 1
    pvalue_out = fast_pvalue_np(rounded_pwm, alpha, ndigits=ndigits, background=background)
    while fast_pvalue_np(rounded_pwm, alpha - error, ndigits=20, background=background) != fast_pvalue_np(rounded_pwm, alpha, ndigits=20, background=background):
        ndigits += digit
        rounded_pwm = round_pwm(pwm, ndigits)
        error = maximal_error_associated(pwm, rounded_pwm)
        Q, scores = score_distribution_np(rounded_pwm, alpha=(
            alpha - error), beta=(alpha + error), ndigits=ndigits, background=background)
        Q = distribution_np_to_dict(Q, scores)
        alpha = search_score(alpha, rounded_pwm, Q, error, ndigits, background)
        pvalue_out = fast_pvalue_np(rounded_pwm, alpha, ndigits=ndigits, background=background)
    return(alpha, pvalue_out)


def bg_to_dict(bg, order, k_mers):
    bg_dict = dict()
    for index, k_mer in enumerate(k_mers):
        bg_dict[k_mer] = bg[order][index]
    return(bg_dict)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--bamm', action='store', dest='input_bamm',
                        required=True, help='path to bamm file /format is .ihbcp/')
    parser.add_argument('-b', '--backgroud', action='store', dest='input_bg',
                        required=True, help='path to backgroud file /format is .hbcp/')
    # parser.add_argument('-f', '--fasta', action='store', dest='fasta',
    #                    required=False, help='path to BED file, needed to calculate backgroun frequences for nucleotieds. \
    #                    Without this parametr background frequances = {A: 0.25, C: 0.25, G: 0.25, T: 0.25}')
    parser.add_argument('-p', '--pvalue', action='store', dest='pvalue',
                        required=True, help='pvalue')
    # parser.add_argument('-pp', '--precisely', action='store_true', dest='precisely',
    #                    required=False, default=False,
    #                    help='Get more precisely value of score')
    parser.add_argument('-r', '--round', action='store', dest='ndigits',
                        required=False, default=int(1),
                        help='Initional matrix rounding in algorithm')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():

    args = parse_args()
    pvalue = float(args.pvalue)
    bamm_path = args.input_bamm
    bg_path = args.input_bg
    #fasta = args.fasta
    ndigits = int(args.ndigits)

    start_time = time.time()
    bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    log_odds_bamm = make_log_odds_bamm(bamm, bg)
    k_mers = make_k_mers(order)
    bamm_dict = bamm_to_dict(log_odds_bamm, order, k_mers)
    background = bg_to_dict(bg, order, k_mers)

    score_value, pvalue_out = from_pvalue_to_score(
        pwm=bamm_dict, pvalue=pvalue, ndigits=ndigits, background=background)
    with_out_round_pvalue = fast_pvalue_np(
        pwm=bamm_dict, score_value=score_value, background=background)
    end_time = time.time()
    print('\nSCORE = {0}\nPVALUE = {1}\nPVALUE_WITHOUT_ROUNDING = {2}'.format(
        score_value, pvalue_out, with_out_round_pvalue))
    print('TIME_OF_CALCULATION = {0} second'.format(end_time - start_time))


if __name__ == '__main__':
    main()
