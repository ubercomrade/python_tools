import linecache
import argparse
import sys
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
                    #print('Record with out strand. Strand is +')
                    record['strand'] = '+'
                continue
            record['seq'] = line.strip().upper()
            fasta.append(record)
    return(fasta)


def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        inf = inf.strip().split('\t')
        id_pfm = inf[0][1:].strip()
        tf_name = inf[1].split('/')[0].strip()
        tf_db = inf[1].split('/')[1].strip()
        inf = {'id_pfm': id_pfm, 'tf_name': tf_name, 'tf_db': tf_db}
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(pwm, inf)


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


def scan_seq_by_pwm(pwm, record, threshold):
    results = []
    reverse_record = reverse_complement(record)
    length_pwm = len(pwm['A'])
    seq = record['seq']
    reverse_seq = reverse_record['seq']

    # first strand
    for i in range(len(seq) - length_pwm + 1):
        site_seq = seq[i:length_pwm + i]
        s = score(site_seq, pwm)
        if s >= threshold:
            site_dict = dict()
            site_dict['id'] = record['id']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['start']) + i)
            site_dict['end'] = str(int(record['start']) + i + length_pwm)
            site_dict['site'] = site_seq
            site_dict['strand'] = record['strand']
            site_dict['score'] = s
            results.append(site_dict)

    # second strand
    for i in range(len(seq) - length_pwm + 1):
        site_seq = reverse_seq[i:length_pwm + i]
        s = score(site_seq, pwm)
        if s >= threshold:
            site_dict = dict()
            site_dict['id'] = record['id']
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
                        required=True, help='path to PWM file')
    parser.add_argument('-t', '--threshold', action='store', type=float, dest='threshold',
                        required=True, help='threshold for PWM')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='path to BED like file for output')
    return(parser.parse_args())


def main():

    args = parse_args()
    pwm_path = args.input_pwm
    fasta_path = args.input_fasta
    threshold = args.threshold
    results_path = args.output

    fasta = read_fasta(fasta_path)
    pwm, inf = read_pwm(pwm_path)
    results = []
    for record in fasta:
        results += scan_seq_by_pwm(pwm, record, threshold)

    df = pd.DataFrame(results)
    df['name'] = np.repeat('.', len(df))
    #df['score'] = np.repeat(0, len(df))
    df = df[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'id', 'site']]
    df.to_csv(results_path, sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
