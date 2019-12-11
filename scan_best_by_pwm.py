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
import multiprocessing as mp
import functools
import re
import csv


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


def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(pwm)  # , inf)


def score(seq, pwm):
    score = sum([pwm[letter][index] for index, letter in enumerate(seq)])
    return(score)


def complement(record):
    output = dict(record)
    strand = record['strand']
    seq = str()
    if strand == '+':
        output['strand'] = '-'
    else:
        output['strand'] = '+'

    seq = output['seq'].replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
    output['seq'] = seq
    return(output)


def check_nucleotides(site):
    s = set(site)
    n = {'A', 'C', 'G', 'T'}
    if len(s - n) == 0:
        return(True)
    else:
        return(False)


def scan_seq_by_pwm(record, pwm):
    results = []
    reverse_record = complement(record)
    length_pwm = len(pwm['A'])
    seq = record['seq']
    reverse_seq = reverse_record['seq']
    threshold = -1000

    # first strand
    for i in range(len(seq) - length_pwm + 1):
        site_seq = seq[i:length_pwm + i]
        if not check_nucleotides(site_seq):
            continue
        s = score(site_seq, pwm)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['start']) + i)
            site_dict['end'] = str(int(record['start']) + i + length_pwm)
            site_dict['site'] = site_seq
            site_dict['strand'] = record['strand']
            site_dict['score'] = s
            threshold = s

    # second strand
    for i in range(len(seq) - length_pwm + 1):
        site_seq = reverse_seq[i:length_pwm + i]
        if not check_nucleotides(site_seq):
            continue
        s = score(site_seq, pwm)
        if s >= threshold:
            site_dict = dict()
            site_dict['name'] = record['name']
            site_dict['chromosome'] = record['chromosome']
            site_dict['start'] = str(int(record['end']) - i - length_pwm)
            site_dict['end'] = str(int(record['end']) - i)
            site_dict['site'] = site_seq
            site_dict['strand'] = reverse_record['strand']
            site_dict['score'] = s
            threshold = s

    results.append(site_dict)
    return(results)


def write_list(path, data):
    scores = [i['score'] for i in data]
    with open(path, "w") as file:
        for line in scores:
            file.write("{0}\n".format(line))
    file.close()
    pass


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file')
    parser.add_argument('-m', '--pwm', action='store', dest='input_pwm',
                        required=True, help='path to PWM file')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='path to write list with scores')

    return(parser.parse_args())


def main():

    args = parse_args()
    pwm_path = args.input_pwm
    fasta_path = args.input_fasta
    results_path = args.output

    fasta = read_fasta(fasta_path)
    pwm = read_pwm(pwm_path)

    results = []
    for record in fasta:
      results += scan_seq_by_pwm(record, pwm)

    write_list(results_path, results)


if __name__ == '__main__':
    main()
