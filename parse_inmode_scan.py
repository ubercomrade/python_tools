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
import os
import random
import itertools
import argparse
import functools
import multiprocessing as mp
import numpy as np
import pandas as pd


def read_sites(path):
    sequences = []
    with open(path, 'r') as file:
        sequences = [i.strip().upper() for i in file if i.strip()[0] != '>']
    return(sequences)


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
                record['name'] = line[0].split('_')[1]
                record['chr'] = line[1]
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


def complement(inseq):
    '''
    Make reverse and compelent
    '''
    seq = str()
    for letter in inseq:
        if letter == 'A':
            seq += 'T'
        elif letter == 'C':
            seq += 'G'
        elif letter == 'G':
            seq += 'C'
        elif letter == 'T':
            seq += 'A'
    seq = seq[::-1]
    return(seq)


def read_inmode_bed(path):
    bed = pd.read_csv(path,
                      sep='\t', header=None,
                      usecols=[0, 1, 2, 3, 4],
                      dtype= {'chromosome': str},
                      names=['id', 'start', 'end', 'strand', 'score'])
    bed['chr'] = '.'
    bed['site'] = '.'
    bed = bed[['chr', 'start', 'end', 'id',  'score', 'strand', 'site']]
    return(bed)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--inputFasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file')
    parser.add_argument('-bed', '--inputBED', action='store', dest='input_bed',
                        required=True, help='path to BED file with inmode scan output')
    parser.add_argument('-o', '--outputBED', action='store', dest='output',
                        required=True, help='path to output BED with sites')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    input_fasta = args.input_fasta
    input_bed = args.input_bed
    output = args.output

    fasta = read_fasta(input_fasta)
    bed = read_inmode_bed(input_bed)

    for index, line in bed.iterrows():
        if line['strand'] == '-':
            bed.loc[index, 'site'] = complement(fasta[line['id']]['seq'])[line['start']:line['end']]
        else:
            bed.loc[index, 'site'] = fasta[line['id']]['seq'][line['start']:line['end']]
        bed.loc[index, 'chr'] = fasta[line['id']]['chr']
        bed.loc[index, 'id'] = 'peaks_' + str(fasta[line['id']]['name'])
        bed.loc[index, 'start'] = int(line['start']) + int(fasta[line['id']]['start'])
        bed.loc[index, 'end'] = int(line['end']) + int(fasta[line['id']]['start'])

    bed.to_csv(output, sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
