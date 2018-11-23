import os
import numpy as np
import itertools
import random
import pandas as pd
import argparse


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
                    #print('Record with out strand. Strand is +')
                    record['strand'] = '+'
                continue
            record['seq'] = line.strip().upper()
            fasta.append(record)
    return(fasta)


def score_bamm(seq, bamm, k_mers, order):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и bamm
    '''
    length_of_seq = len(seq)
    position = 0
    score = 0
    for position in range(length_of_seq - order):
        letter = seq[position:order + position + 1]
        loc = k_mers[letter]
        score += bamm[order][position][loc]
    return(score)


def make_k_mers(order):
    #  make list with possible k-mer based on bHMM model
    tmp = itertools.product('ACGT', repeat=order + 1)
    k_mer = []
    for i in tmp:
        k_mer.append(''.join(i[1:]) + i[0])
    k_mer_dict = dict()
    index = 0
    for i in k_mer:
        k_mer_dict[i] = index
        index += 1
    return(k_mer_dict)


def reverse_complement(record):
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


def scan_seq_by_bamm(record, log_odds_bamm, order, threshold):

    k_mers = make_k_mers(order)
    motif_length = len(log_odds_bamm[0])
    reverse_record = reverse_complement(record)
    seq = record['seq']
    reverse_seq = reverse_record['seq']
    results = []

    # scan first strand
    for i in range(len(seq) - motif_length + 1):
        site_seq = seq[i:motif_length + i]
        s = score_bamm(site_seq, log_odds_bamm, k_mers, order)
        if s >= threshold:
            site_dict = dict()
            site_dict['id'] = record['id']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['start']) + i)
            site_dict['end'] = str(int(record['start']) + i + motif_length)
            site_dict['site'] = site_seq
            site_dict['strand'] = record['strand']
            site_dict['score'] = s
            results.append(site_dict)

    # scan second strand
    for i in range(len(seq) - motif_length + 1):
        site_seq = reverse_seq[i:motif_length + i]
        s = score_bamm(site_seq, log_odds_bamm, k_mers, order)
        if s >= threshold:
            site_dict = dict()
            site_dict['id'] = record['id']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['end']) - i - motif_length)
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

    fasta = read_fasta(fasta_path)
    bamm, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    log_odds_bamm = make_log_odds_bamm(bamm, bg)
    results = []
    for record in fasta:
        results += scan_seq_by_bamm(record, log_odds_bamm, order, threshold)

    df = pd.DataFrame(results)
    df['name'] = np.repeat('.', len(df))
    #df['score'] = np.repeat(0, len(df))
    df = df[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'id', 'site']]
    df.to_csv(results_path, sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
