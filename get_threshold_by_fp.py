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


import linecache
import argparse
import sys
import os
import itertools
import multiprocessing as mp
import functools
import pandas as pd
import numpy as np


def read_fasta(path):
    '''
    Чтение фаста фаила и запись каждых двух строчек в экземпляр класса BioRecord
    Все экземпляры хранятся в списке
    Функция возвращает список экземпляров класса BioRecord

    Шапка для FASTA: >uniq_id|chromosome|start-end|strand
    '''
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                line = line[1:].strip().split('|')
                record = dict()
                record['id'] = line[0]
                record['chromosome'] = line[1]
                record['start'] = line[2].split('-')[0]
                record['end'] = line[2].split('-')[1]
                try:
                    record['strand'] = line[3]
                except:
                    # print('Record with out strand. Strand is +')
                    record['strand'] = '+'
                continue
            record['seq'] = line.strip().upper()
            fasta.append(record)
    return(fasta)


def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        #inf = inf.strip().split('\t')
        #id_pfm = inf[0][1:].strip()
        #tf_name = inf[1].split('/')[0].strip()
        #tf_db = inf[1].split('/')[1].strip()
        #inf = {'id_pfm': id_pfm, 'tf_name': tf_name, 'tf_db': tf_db}
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(pwm)  # , inf)


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

        # convert a bamm array to numpy array
        for k in range(motif_order):
            model[k] = np.array(model[k], dtype=float)
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
                bg[order] = np.array(bg_freq, dtype=float)
                order += 1
    else:
        print('File {0} does not exist'.format(bg_file))
        sys.exit()
    return model, bg, order-1


def make_log_odds_bamm(bamm, bg):
    log_odds_bamm = dict()
    for order in bamm.keys():
        log_odds_bamm[order] = np.array(np.log2(bamm[order] / bg[order]))
    return(log_odds_bamm)


def read_fasta(path):
    '''
    Чтение фаста фаила и запись каждых двух строчек в экземпляр класса BioRecord
    Все экземпляры хранятся в списке
    Функция возвращает список экземпляров класса BioRecord

    Шапка для FASTA: >uniq_id|chromosome|start-end|strand
    '''
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                line = line[1:].strip().split('|')
                record = dict()
                record['id'] = line[0]
                record['chromosome'] = line[1]
                record['start'] = line[2].split('-')[0]
                record['end'] = line[2].split('-')[1]
                try:
                    record['strand'] = line[3]
                except:
                    # print('Record with out strand. Strand is +')
                    record['strand'] = '+'
                continue
            record['seq'] = line.strip().upper()
            fasta.append(record)
    return(fasta)


def score(seq, pwm):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и PWM
    kind - тип PWM mono or di
    '''
    length_of_seq = len(seq)
    position = 0
    score = 0
    for letter in seq:
        score += pwm[letter][position]
        position += 1
    return(score)


def score_bamm(seq, bamm, k_mers, order):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и bamm
    '''
    length_of_seq = len(seq)
    position = 0
    score = 0
    for position in range(length_of_seq - order):
        letter = seq[position:order + position + 1][::-1]
        loc = k_mers[letter]
        score += bamm[order][position][loc]
    return(score)


# def make_k_mers(order):
#    #  make list with possible k-mer based on bHMM model
#    tmp = itertools.product('ACGT', repeat=order + 1)
#    k_mer = []
#    for i in tmp:
#        k_mer.append(''.join(i[1:]) + i[0])
#    k_mer_dict = dict()
#    index = 0
#    for i in k_mer:
#        k_mer_dict[i] = index
#        index += 1
#    return(k_mer_dict)


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


def complement(record):
    '''
    Make reverse and compelent
    '''
    output = dict(record)
    strand = record['strand']
    seq = str()
    if strand == '+':
        output['strand'] = '-'
    else:
        output['strand'] = '+'
    for letter in output['seq']:
        if letter == 'A':
            seq += 'T'
        elif letter == 'C':
            seq += 'G'
        elif letter == 'G':
            seq += 'C'
        elif letter == 'T':
            seq += 'A'
    output['seq'] = seq[::-1]
    return(output)


def scan_seq_by_bamm(record, log_odds_bamm, order):

    k_mers = make_k_mers(order)
    motif_length = len(log_odds_bamm[0])
    reverse_record = complement(record)
    seq = record['seq']
    reverse_seq = reverse_record['seq']
    results = []

    # scan first strand
    for i in range(len(seq) - motif_length + 1):
        site_seq = seq[i:motif_length + i]
        s = score_bamm(site_seq, log_odds_bamm, k_mers, order)
        results.append(s)

    # scan second strand
    for i in range(len(seq) - motif_length + 1):
        site_seq = reverse_seq[i:motif_length + i]
        s = score_bamm(site_seq, log_odds_bamm, k_mers, order)
        results.append(s)
    return(results)


def scan_seq_by_pwm(record, pwm):
    results = []
    reverse_record = complement(record)
    length_pwm = len(pwm['A'])
    seq = record['seq']
    reverse_seq = reverse_record['seq']

    # first strand
    for i in range(len(seq) - length_pwm + 1):
        site_seq = seq[i:length_pwm + i]
        s = score(site_seq, pwm)
        results.append(s)

    # second strand
    for i in range(len(seq) - length_pwm + 1):
        site_seq = reverse_seq[i:length_pwm + i]
        s = score(site_seq, pwm)
        results.append(s)
    return(results)


def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name')
    pwm_parser = subparsers.add_parser('pwm', help='threshold for PWM')
    pwm_parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                            required=True, help='path to FASTA file with head: >uniq_id|chromosome|start-end|strand')
    pwm_parser.add_argument('-m', '--pwm', action='store', dest='input_pwm',
                            required=True, help='path to PWM file')
    pwm_parser.add_argument('-p', '--false_positive', action='store', type=float, dest='false_positive',
                            required=True, help='value of FP (FP ~ P-VALUE)')

    bamm_parser = subparsers.add_parser('bamm', help='threshold for bamm')
    bamm_parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                             required=True, help='path to FASTA file with head: >uniq_id|chromosome|start-end|strand')
    bamm_parser.add_argument('-m', '--bamm', action='store', dest='input_bamm',
                             required=True, help='path to bamm file /format is .ihbcp/')
    bamm_parser.add_argument('-b', '--backgroud', action='store', dest='input_bg',
                             required=True, help='path to backgroud file /format is .hbcp/')
    bamm_parser.add_argument('-p', '--false_positive', action='store', type=float, dest='false_positive',
                             required=True, help='value of FP (FP ~ P-VALUE)')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def get_threshold(results, fp):
    scores = np.array(results)
    scores.sort()
    scores = scores[::-1]  # sorted score from big to small
    i = int(len(scores) * fp) - 10  # position of score value
    step = 0
    while abs((scores >= scores[i + step]).sum() / len(scores) - fp) \
        >= abs((scores >= scores[i + step + 1]).sum() / len(scores) - fp) \
            or scores[i + step] == scores[i + step + 1]:
        step += 1
    thr = scores[i + step]  # possible threshold score
    calc_fp = (scores >= thr).sum() / len(scores)  # p-value with possible threshold score
    # delta = abs(calc_pval - pval)  # diff between p-value and calculated p-value
    return(thr, calc_fp)


def main():
    args = parse_args()

    if args.subparser_name == 'pwm':
        pwm_path = args.input_pwm
        fasta_path = args.input_fasta
        fp = args.false_positive

        fasta = read_fasta(fasta_path)
        pwm = read_pwm(pwm_path)

        # results = []
        # for record in fasta:
        #    results += scan_seq_by_pwm(record, pwm)

        with mp.Pool(4) as p:
            results = p.map(functools.partial(scan_seq_by_pwm, pwm=pwm), fasta)
        results = [i for i in results if i != []]
        results = [j for sub in results for j in sub]

        thr, calc_fp = get_threshold(results, fp)
        #print('Optimal score threshold = {0}\nCalculated FP = {1}'. format(thr, calc_fp))
        print(thr)

    elif args.subparser_name == 'bamm':
        bamm_path = args.input_bamm
        bg_path = args.input_bg
        fasta_path = args.input_fasta
        fp = args.false_positive

        fasta = read_fasta(fasta_path)
        bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
        log_odds_bamm = make_log_odds_bamm(bamm, bg)

        # results=[]
        # for record in fasta:
        #    results += scan_seq_by_bamm(record, log_odds_bamm, order)

        with mp.Pool(4) as p:
            results = p.map(functools.partial(scan_seq_by_bamm,
                                              log_odds_bamm=log_odds_bamm, order=order), fasta)
        results = [i for i in results if i != []]
        results = [j for sub in results for j in sub]

        thr, calc_fp = get_threshold(results, fp)
        #print('Optimal score threshold = {0}\nCalculated FP = {1}'. format(thr, calc_fp))
        print(thr)


if __name__ == '__main__':
    main()
