
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



def make_pfm_from_pcm(pcm, kind, pseudocount= '1/N'):
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
        PFM = {}
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            pfm[''.join(i)] = []
    elif kind == 'mono':
        pfm = {}
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            PFM[i[0]] = []
    else:
        print('ALARM!')
        pass

    if pseudocount == '1/N':
        first_key = list(pcm.keys())[0]
        nuc_pseudo = 1/len(pcm.keys())
        for i in range(len(pcm[first_key])):
            for nuc in pcm.keys():
                pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites + 1))
        return(PFM)

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


def main(path, homology, kind, pseudocount):
    seq = read_fasta(fileInpute)
    seq = remove_equalent_seq(seq, homology=homology)
    background = background_freq(seq, kind=kind)
    PCM = make_pcm(seq, kind=kind)
    PWM = make_pwm_from_pcm(PCM, kind=kind, background=background, pseudocount=pseudocount)
    return(PWM)
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input', required=True , help='Path to file with motifs')
    parser.add_argument('-o', '--output', action='store', dest='output', required=True, help='Path to file with results of calculation')
    parser.add_argument('-t', '--type', action='store', dest='kind', required=True, help='Type of matrix di- or mononucliotide, after flag print di or mono')
    parser.add_argument('-H', '--removeHomologs', action='store', dest='homology', required=True, help='portion of homology, for example 0.95')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    
    path = args.input
    output = args.output 
    kind = args.kind
    pseudocount = '1/N'
    homology = float(args.homology)
    PWM = main(path, homology, kind, pseudocount)
    
    with open(output, 'w', newline='') as file:
        keys = sorted(list(PWM.keys()))
        file.write('>' + '\t'.join(keys) + '\n')
        for i in range(len(PWM[keys[0]])):
            for j in keys:
                if len(keys) - keys.index(j) != 1:
                    file.write(str(PWM[j][i]) + '\t')
                else:
                    file.write(str(PWM[j][i]) + '\n')

