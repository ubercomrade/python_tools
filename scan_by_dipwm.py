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
                record['name'] = line[0]
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
    file.close()
    return(fasta)


def read_dipwm(path):
    dipwm = dict()
    with open(path, 'r') as file:
        inf = file.readline()
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            dipwm[''.join(i)] = []
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(dipwm.keys(), line):
                dipwm[letter].append(float(value))
    file.close()
    return(dipwm)  # , inf)


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


def scan_seq_by_dipwm(record, dipwm, threshold):
    results = []
    reverse_record = complement(record)
    length_pwm = len(dipwm['AA'])
    seq = record['seq']
    reverse_seq = reverse_record['seq']

    # first strand
    for i in range(len(seq) - length_pwm + 2):
        site_seq = seq[i:length_pwm + i]
        s = score_dipwm(site_seq, dipwm)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['start']) + i)
            site_dict['end'] = str(int(record['start']) + i + length_pwm)
            site_dict['site'] = site_seq
            site_dict['strand'] = record['strand']
            site_dict['score'] = s
            results.append(site_dict)

    # second strand
    for i in range(len(seq) - length_pwm + 2):
        site_seq = reverse_seq[i:length_pwm + i]
        s = score_dipwm(site_seq, dipwm)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['end']) - i - length_pwm)
            site_dict['end'] = str(int(record['end']) - i)
            site_dict['site'] = site_seq
            site_dict['strand'] = reverse_record['strand']
            site_dict['score'] = s
            results.append(site_dict)
    return(results)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file with head: >uniq_id|chromosome|start-end|strand')
    parser.add_argument('-m', '--pwm', action='store', dest='input_pwm',
                        required=True, help='path to diPWM file')
    parser.add_argument('-t', '--threshold', action='store', type=float, dest='threshold',
                        required=True, help='threshold for diPWM')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='path to BED like file for output')
    parser.add_argument('-P', '--processes', action='store', type=int, dest='cpu_count',
    required=False, default=2, help='Number of processes to use, default: 2')
    return(parser.parse_args())


def main():

    args = parse_args()
    pwm_path = args.input_pwm
    fasta_path = args.input_fasta
    threshold = args.threshold
    results_path = args.output
    cpu_count = args.cpu_count

    fasta = read_fasta(fasta_path)
    dipwm = read_dipwm(pwm_path)

    # results = []
    # for record in fasta:
    #    results += scan_seq_by_pwm(record, pwm, threshold)

    with mp.Pool(cpu_count) as p:
        results = p.map(functools.partial(scan_seq_by_dipwm,
                                          dipwm=dipwm, threshold=threshold), fasta)
    results = [i for i in results if i != []]
    results = [j for sub in results for j in sub]

    df = pd.DataFrame(results)
    # df['name'] = np.repeat('.', len(df))
    # df['score'] = np.repeat(0, len(df))
    df = df[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']]
    # print(df)
    df.to_csv(results_path, sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
