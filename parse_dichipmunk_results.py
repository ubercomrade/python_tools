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
import csv
import sys
import os
import random
import itertools
import argparse
import numpy as np
import pandas as pd


def get_sequences(path):
    with open(path, 'r') as file:
        output = []
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
                output.append(d)
            else:
                continue
    df = pd.DataFrame(output)
    df = df[['name', 'start', 'end', 'strand', 'seq']]
    return(df)


def read_fasta(path):
    sequences = []
    with open(path, 'r') as file:
        sequences = [i.strip().upper() for i in file if i.strip()[0] != '>']
    return(sequences)


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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to ChIPMunk results')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='dir to write output files')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    parser.add_argument('-b', '--background', action='store', dest='background',
                        required=False, help='path to FASTA file, needed to calculate backgroun frequences for nucleotieds. \
                        Without this parametr background calculated based on sites')
    return(parser.parse_args())


def write_pwm(output, tag, pwm):
    matrix = [array for array in pwm.values()]
    with open(output + '/' + tag + '.pwm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in list(map(list, zip(*matrix))):
            i = [str(j) for j in i]
            file.write('\t'.join(i) + '\n')


def write_fasta(output, tag, sites):
    with open(output + '/' + tag + '.fasta', 'w') as file:
        for index, site in enumerate(sites):
            file.write('>site_' + str(index) + '\n')
            file.write(site + '\n')


def write_sites(output, tag, sites):
    with open(output + '/' + tag + '.sites', 'w') as file:
        for site in sites:
            file.write(site + '\n')


def main():
    args = parse_args()
    input_ = args.input
    output_ = args.output
    background_path = args.background
    tag = args.tag

    seq = get_sequences(input_)
    seq = remove_equalent_seq(seq_list=list(seq['seq']), homology=0.95)
    nsites = len(seq)
    if background_path is None:
        background = background_freq_di(seq)
    else:
        fasta = read_fasta(background_path)
        background = background_freq_di(fasta)

    pcm = make_dipcm(seq)
    pfm = make_dipfm_from_dipcm(pcm)
    pwm = make_dipwm_from_dipcm(pcm, background)

    if not os.path.isdir(output_):
        os.mkdir(output_)

    write_pwm(output=output_, tag=tag, pwm=pwm)
    write_fasta(output=output_, tag=tag, sites=seq)
    write_sites(output=output_, tag=tag, sites=seq)


if __name__ == '__main__':
    main()
