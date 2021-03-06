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
import re


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


def read_fasta(path):
    '''
    Чтение фаста фаила и запись каждых двух строчек в экземпляр класса BioRecord
    Все экземпляры хранятся в списке
    Функция возвращает список экземпляров класса BioRecord

    Шапка для FASTA: >uniq_id::chromosome:start-end(strand)
    '''
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            #print(line)
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                record = dict()
                record['name'] = line[0]
                record['chromosome'] = line[2]
                coordinates_strand = line[3]

                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                record['start'] = start
                record['end'] = end

                strand = re.findall(r'\(.\)', coordinates_strand[:-3])
                if not strand == []:
                    record['strand'] = strand[0].strip('()')
                else:
                    record['strand'] = '+'
            else:
                record['seq'] = line.strip().upper()
                fasta.append(record)
    file.close()
    return(fasta)


def score_bamm(seq, bamm, order):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и bamm
    '''
    length_of_seq = len(seq)
    score = 0
    #site = ''
    for position in range(len(seq) - order):
        score += bamm[seq[position:position + order + 1]][position]
    #    site += seq[position:position + order + 1][-1]
    return(score)#, site)


def bamm_to_dict(log_odds_bamm, order, k_mers):
    bamm_dict = {}
    for k_mer in k_mers:
        bamm_dict[k_mer] = list()
    for i in range(len(log_odds_bamm[order])):
        for index, k_mer in enumerate(k_mers):
            bamm_dict[k_mer].append(log_odds_bamm[order][i][index])
    return(bamm_dict)


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


def complement(line):
    '''
    Make reverse and compelent
    '''
    seq = ''
    for letter in line:
        if letter == 'A':
            seq += 'T'
        elif letter == 'C':
            seq += 'G'
        elif letter == 'G':
            seq += 'C'
        elif letter == 'T':
            seq += 'A'
    output = seq[::-1]
    return(output)


def scan_seq_by_bamm(record, log_odds_bamm, order, threshold):

    motif_length = len(log_odds_bamm[list(log_odds_bamm.keys())[0]])
    reverse_record = complement(record)
    seq = record['seq']
    reverse_seq = reverse_record['seq']
    results = []
    for i in range(len(seq) - motif_length - order + 1):
        site_seq = seq[i:motif_length + order + i]
        if 'N' in site_seq:
            continue
        s, site_seq = score_bamm(site_seq, log_odds_bamm, order)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['start']) + i + order)
            site_dict['end'] = str(int(record['start']) + i + motif_length + order)
            site_dict['site'] = site_seq
            site_dict['strand'] = record['strand']
            site_dict['score'] = s
            results.append(site_dict)

    # scan second strand
    for i in range(len(seq) - motif_length - order + 1):
        site_seq = reverse_seq[i:motif_length + order + i]
        if 'N' in site_seq:
            continue
        s, site_seq = score_bamm(site_seq, log_odds_bamm, order)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['end']) - i - motif_length - order)
            site_dict['end'] = str(int(record['end']) - i - order)
            site_dict['site'] = site_seq
            site_dict['strand'] = reverse_record['strand']
            site_dict['score'] = s
            results.append(site_dict)
    return(results)


def read_and_scan_genome(path, path_to_write, log_odds_bamm, order, threshold):
    end_of_line = ''
    results = open(path_to_write, 'w+')
    motif_length = len(log_odds_bamm[list(log_odds_bamm.keys())[0]])
    with open(path, 'r') as file:
        for line in file:
            #line = end_of_line + line.strip().upper()
            if line.startswith('>'):
                chr = line.split()[0][1:]
                position = 0
                end_of_line = ''
            else:
                line = end_of_line + line.strip().upper()
                #print(line)
                for i in range(len(line) - motif_length - order + 1):
                    seq = line[i:motif_length + order + i]
                    #print(seq)
                    if 'N' in seq:
                        position += 1
                        continue

                    # +
                    s = score_bamm(seq, log_odds_bamm, order)
                    if s >= threshold:
                        results.write('{0}\t{1}\t{2}\t.\t{3}\t+\t{4}\n'.format(chr, position, position + motif_length + order, s, seq))
                    # -
                    seq = complement(seq)
                    s = score_bamm(seq, log_odds_bamm, order)
                    if s >= threshold:
                        results.write('{0}\t{1}\t{2}\t.\t{3}\t-\t{4}\n'.format(chr, position, position + motif_length + order, s, seq))

                    position += 1
                    end_of_line = line[i + 1:]

                #end_of_line = seq = complement(seq)
                #position -= 1
    file.close()
    results.close()



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file')
    parser.add_argument('-m', '--bamm', action='store', dest='input_bamm',
                        required=True, help='path to bamm file /format is .ihbcp/')
    parser.add_argument('-b', '--backgroud', action='store', dest='input_bg',
                        required=True, help='path to backgroud file /format is .hbcp/')
    parser.add_argument('-t', '--threshold', action='store', type=float, dest='threshold',
                        required=True, help='threshold for model')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='path to BED like file for output')
    return(parser.parse_args())


def main():

    args = parse_args()
    bamm_path = args.input_bamm
    bg_path = args.input_bg
    fasta_path = args.input_fasta
    threshold = args.threshold
    results_path = args.output

    bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    log_odds_bamm = make_log_odds_bamm(bamm, bg)
    k_mers = make_k_mers(order=order)
    log_odds_bamm = bamm_to_dict(log_odds_bamm=log_odds_bamm, order=order, k_mers=k_mers)

    # results = []
    # for record in fasta:
    #    results += scan_seq_by_bamm(record, log_odds_bamm, order, threshold)

    read_and_scan_genome(fasta_path, results_path, log_odds_bamm, order, threshold)


if __name__ == '__main__':
    main()
