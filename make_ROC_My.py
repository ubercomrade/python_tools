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
import random
import itertools
import argparse
import numpy as np


def read_fasta(path, everyStr=False):
    '''
    Чтение мотивов из фаила, если это просто список строк с разделителем \n, то everyStr = True,
    если это fasta формат, то everyStr = False
    '''
    sequences = []
    if everyStr:
        with open(path, 'r') as file:
            sequences = [i.strip().upper() for i in file]
    else:
        with open(path, 'r') as file:
            sequences = []
            seq = ''
            for line in file:
                if line[0] == '>':
                    sequences.append(seq)
                    seq = ''
                    continue
                else:
                    seq += line.strip().upper()
            sequences = sequences[1:]
    return(sequences)


def read_sites(path, every_str=False):
    '''
    Чтение мотивов из фаила, если это просто список строк с разделителем \n, то everyStr = True,
    если это fasta формат, то everyStr = False
    '''
    sequences = []
    if every_str:
        with open(path, 'r') as file:
            sequences = [i.strip().upper() for i in file]
    else:
        with open(path, 'r') as file:
            sequences = [i.strip().upper() for i in file if i.strip() != '>']
    return(sequences)


def remove_equalent_seq(seqList, homology=0.95):
    '''
    Удаление гомологичных последовательностей из списка (seqList)
    Если кол-во совпадений при сравнении последовательности 1 и 2 >= длина последовательности * homology,
    то последовательность 1 удаляется из списка
    Функция возвращает новый список
    '''
    seqList = list(seqList)
    treshold = homology * len(seqList[0])
    for seq1 in tuple(seqList):
        sub_seqList = list(seqList)
        sub_seqList.remove(seq1)
        for seq2 in sub_seqList:
            score = len([i for i, j in zip(seq1, seq2) if i == j])
            if score >= treshold:
                seqList.remove(seq1)
                break
    return(seqList)


def make_pfm_from_pcm(pcm, kind, pseudocount='1/N'):
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

    number_of_sites = int()
    for i in pcm.keys():
        number_of_sites += pcm[i][0]

    if kind == 'di':
        pfm = {}
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            pfm[''.join(i)] = []
    elif kind == 'mono':
        pfm = {}
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            pfm[i[0]] = []
    else:
        print('ALARM!')
        pass

    if pseudocount == '1/N':
        first_key = list(pcm.keys())[0]
        nuc_pseudo = 1/len(pcm.keys())
        for i in range(len(pcm[first_key])):
            for nuc in pcm.keys():
                pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites + 1))
        return(pfm)

    elif pseudocount == 'sqroot':
        total_sq_root = int()
        for i in pcm.keys():
            total_sq_root += pcm[i][0]
        total_sq_root = math.sqrt(total_sq_root)
        sq_root = totalSqRoot/len(pcm.keys())

        first_key = list(pcm.keys())[0]
        for i in range(len(PCM[first_key])):
            for nuc in pcm.keys():
                pfm[nuc].append((PCM[nuc][i] + sq_root) / (number_of_sites + total_sq_root))

        return(pfm)
    else:
        print('ALARM!')
        pass


def make_pwm_from_pcm(pcm, kind, background, method='log-odds', pseudocount='1/N'):
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
    if kind == 'di':
        pwm = {}
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            pwm[''.join(i)] = []
    elif kind == 'mono':
        pwm = {}
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            pwm[i[0]] = []
    else:
        print('ALARM!')
        pass

    pfm = make_pfm_from_pcm(pcm, kind, pseudocount)
    first_key = list(pcm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def make_pcm(motifs, kind):
    '''
    input - список мотивов одинаковой длины
    output -  PCM
    kind is type of matrix di or mono

    Создает PCM на основе списка мотивов
    '''
    if kind == 'di':
        matrix = {}
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            matrix[''.join(i)] = []
    elif kind == 'mono':
        matrix = {}
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            matrix[i[0]] = []

    len_of_motif = len(motifs[0])

    if kind == 'di':
        for i in matrix.keys():
            matrix[i] = [0]*(len_of_motif - 1)

        for i in range(len_of_motif - 1):
            for l in motifs:
                matrix[l[i:i+2]][i] += 1
    elif kind == 'mono':
        for i in matrix.keys():
            matrix[i] = [0]*len_of_motif

        for i in range(len_of_motif):
            for l in motifs:
                matrix[l[i]][i] += 1

    return(matrix)


def background_freq(seq, kind):

    s = ''.join(seq)
    if kind == 'mono':
        background = {}
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            background[i[0]] = s.count(i[0])

    elif kind == 'di':
        background = {}
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            background[''.join(i)] = s.count(i[0])

    sum_of_nuc = sum(background.values())
    for i in background.keys():
        background[i] = background[i]/sum_of_nuc
    return(background)


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


def score(seq, pwm, kind):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и PWM
    kind - тип PWM mono or di
    '''
    if kind == 'mono':
        length_of_seq = len(seq)
        position = 0
        score = 0
        for letter in seq:
            score += pwm[letter][position]
            position += 1
        return(score)

    elif kind == 'di':
        length_of_seq = len(seq)
        score = 0
        for i in range(len(seq) - 1):
            two_letters = seq[i:i+2]
            score += pwm[two_letters][i]
        return(score)
    else:
        pass


def scan_sequences_1(sequences, PWM, treshold, kind):
    '''
    sequences - список последовательностей (длина последовательностей совпадает с длиной матрицы)
    PWM - матрица
    Сканирование последовательностей матрицей PWM и вычесление score,
    Функция возвращает кол-во последовательностей для которыйх score >= treshold
    kind - тип PWM mono or di
    '''
    aboveTreshold = int()
    for i in sequences:
        s = score(i, PWM, kind)
        if s >= treshold:
            aboveTreshold += 1
    return(aboveTreshold)


def scan_sequences_2(sequences, PWM, kind):
    '''
    sequences - список последовательностей (длина последовательностей совпадает с длиной матрицы)
    PWM - матрица
    Сканирование последовательностей матрицей PWM и вычесление score,
    Функция возвращает scores для всех последовательностей
    kind - тип PWM mono or di
    '''
    scores = []
    for i in sequences:
        s = score(i, PWM, kind)
        scores.append(s)
    return(scores)


def random_seq(seq, k):
    '''
    Генерирует k случайных последовательностей с таким же нуклеотидным составом как и seq
    Возвращает list()
    '''
    out = list()
    for i in range(k):
        out.append(''.join(random.sample(seq, k=len(seq))))
    return(out)


def main(path, homology, kind):
    '''
    Calculating k-fold cross validation
    '''
    k = 10
    all_motifs = read_fasta(path, everyStr=False)
    print(len(all_motifs))
    #all_motifs = [i[50:-50] for i in all_motifs]  # for JASPAR
    #all_motifs = remove_equalent_seq(all_motifs, homology=homology)
    print(len(all_motifs))
    all_motifs = random.sample(all_motifs, len(all_motifs))
    all_motifs = all_motifs[:len(all_motifs) - len(all_motifs) % k]
    total_length = len(all_motifs)
    background = background_freq(all_motifs, kind=kind)

    results = {'FPR': [], 'TPR': [], 'Score': [], 'Norm_score': []}
    tp_scores = []
    fp_scores = []
    step = int(total_length / k)

    for i in range(step, total_length, step):
        train_set = all_motifs[i:] + all_motifs[:i - 30]
        test_set = all_motifs[i - step:i]
        random_set = []
        for seq in test_set:
            random_set += random_seq(seq, 10)

        PCM = make_pcm(train_set, kind=kind)
        PWM = make_pwm_from_pcm(PCM, background=background, kind=kind)

        tp_scores += scan_sequences_2(test_set, PWM, kind=kind)
        fp_scores += scan_sequences_2(random_set, PWM, kind=kind)

    all_scores = tp_scores + fp_scores
    all_scores.sort()
    uniq_scores = set(all_scores)
    uniq_scores = list(uniq_scores)
    uniq_scores.sort()
    PCM = make_pcm(all_motifs, kind=kind)
    PWM = make_pwm_from_pcm(PCM, background=background, kind=kind)
    max_score_value = max_score(PWM)
    min_score_value = min_score(PWM)

    for score in uniq_scores:
        norm_score = to_norm(score, min_score_value, max_score_value)
        fpr = len([i for i in fp_scores if i >= score])/len(fp_scores)
        tpr = len([i for i in tp_scores if i >= score])/len(tp_scores)
        if fpr == 0:
            fpr = 1/len(fp_scores)

        results['FPR'].append(fpr)
        results['TPR'].append(tpr)
        results['Score'].append(score)
        results['Norm_score'].append(norm_score)
    return(results)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='Path to file with binding sites')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='Path to file with results of calculation')
    parser.add_argument('-t', '--type', action='store', dest='kind', required=True,
                        help='Type of matrix di- or mononucliotide, after flag print di or mono')
    parser.add_argument('-H', '--removeHomologs', action='store', dest='homology',
                        required=True, help='portion of homology, for example 0.95')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    path = args.input  # Путь к фаилу с мотивами
    # Допустимый уровень гомологии между мотивами (если уровень гомологии выше, то последовательность выбрасывается)
    homology = float(args.homology)
    kind = args.kind  # Тип матрицы (mono or di)
    output_path = args.output  # Путь к фаилу для записи результатов

    #fileInput = '/Users/anton/Documents/Python/JASPAR/MA0491.1.sites'
    #everyStr = False
    #homology = 0.95
    #kind = 'mono'
    #fileOutput = 'test.csv'

    output = main(path, homology, kind)
    with open(output_path, 'w', newline='') as csvfile:
        fieldnames = ['FPR', 'TPR', 'Score', 'Norm_score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i in range(len(output['FPR'])):
            writer.writerow({'FPR': output['FPR'][i], 'TPR': output['TPR'][i],
                             'Score': output['Score'][i], 'Norm_score': output['Norm_score'][i]})
