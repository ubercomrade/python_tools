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
import random
import itertools
import argparse
import numpy as np

def read_sites(path):
    '''
    Чтение мотивов из фаила, если это просто список строк с разделителем \n, то everyStr = True,
    если это fasta формат, то everyStr = False
    '''
    sequences = []
    with open(path, 'r') as file:
        sequences = [i.strip().upper() for i in file]
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

def write_pfm(output, tag, pfm):
    with open(output + '/' + tag + '.pfm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to sites')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='dir to PFM output file')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    return(parser.parse_args())


def main():
    args = parse_args()
    sites_path = args.input
    output_dir = args.output
    tag = args.tag

    seq = read_sites(sites_path)
    seq = remove_equalent_seq(seq_list=seq, homology=0.95)
    pcm = make_pcm(seq)
    pfm = make_pfm_from_pcm(pcm)
    write_pfm(output_dir, tag, pfm)

if __name__ == '__main__':
    main()
