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
import pickle
from itertools import islice
import pandas as pd
import numpy as np


# def read_bed(path):
#    bed = pd.read_csv(path,
#                      sep='\t', header=None,
#                      usecols=[0, 1, 2, 3, 4, 5],
#                      names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
#    return(bed)


def read_bed(path):

    try:
        bed = pd.read_csv(path,
                          sep='\t', header=None,
                          usecols=[0, 1, 2, 3, 4, 5],
                          dtype= {'chromosome': str, 'strand': str},
                          names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    except:
        bed = pd.read_csv(path,
                          sep='\t', header=None,
                          usecols=[0, 1, 2, 3, 4],
                          dtype= {'chromosome': str, 'strand': str},
                          names=['chromosome', 'start', 'end', 'name', 'score'])

        bed['strand'] = '.'
    #bed = bed.sort_values(by=['chromosome', 'start'])
    return(bed)


def stat_of_fasta_file(path):
    '''
    start line == 1
    '''
    inf_data = path + '.pickle'
    flag = os.path.isfile(inf_data)
    if flag:
        with open(inf_data, 'rb') as file:
            output = pickle.load(file)
        file.close()
        return output
    else:
        stat = {}
        length = int()
        with open(path, 'r') as file:
            count = 0
            for line in file:
                line = line.strip().split()[0]
                if line[0] == '>':
                    stat[line[1:]] = count + 2
                else:
                    if count == 2:
                        length = len(line)
                count += 1
        file.close()
        output = (stat, length)
        with open(inf_data, 'wb') as file:
            pickle.dump(output, file)
        file.close()
        return output


def complement(seq):
    '''
    Make reverse and compelent
    '''
    output = str()
    for letter in seq:
        if letter == 'A':
            output += 'T'
        elif letter == 'C':
            output += 'G'
        elif letter == 'G':
            output += 'C'
        elif letter == 'T':
            output += 'A'
    output = output[::-1]
    return(output)


def chek_nucleotides(line):
    flag = True
    for char in line:
        flag = char is 'A' or char is 'C' or char is 'G' or char is 'T'
        if not flag:
            break
    return flag


def bed_to_fasta(path_fasta, path_bed):
    bed_peaks = read_bed(path_bed)
    chromosomes, length = stat_of_fasta_file(path_fasta)
    results = []
    bed_peaks['seq'] = ''
    for index, row in bed_peaks.iterrows():
        seq = str()  # container for sequnce of peak
        try:
            chr_start_line = chromosomes[row['chromosome']]
        except:
            print('{0} was not found in FASTA'.format(row['chromosome']))
            continue

        if row['start'] < 0 or row['end'] < 0:
            continue


        peak_start_line = row['start'] // length + \
            chr_start_line  # number of line where seq started
        peak_end_line = row['end'] // length + \
            chr_start_line  # number of line where seq ended

        # start and end ofseq nucleotide
        position_start = row['start'] % length
        position_end = row['end'] - row['start'] + position_start


        file = open(path_fasta, 'r') #FASTA file
        lines = list(islice(file, peak_start_line, peak_end_line + 1))
        file.close()
        lines = [i.strip() for i in lines]    

        for k in lines:
            seq += k

        if chek_nucleotides(seq[position_start:position_end]):
            if not row['strand'] == '-':
                bed_peaks.loc[index, 'seq'] = seq[position_start:position_end]
            else:
                bed_peaks.loc[index, 'seq'] = complement(seq[position_start:position_end])
        else:
            continue
    return(bed_peaks)


def write_fasta(peaks_fasta, path_to_write):
    with open(path_to_write, 'w') as file:
        for index, record in peaks_fasta.iterrows():
            file.write('>{0}|{1}|{2}-{3}|{4}\n'.format(record['name'], record['chromosome'],
                record['start'], record['end'], record['strand']))
            file.write(record['seq'] + '\n')
    file.close()
    pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--inputFasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file')
    parser.add_argument('-bed', '--inputBED', action='store', dest='input_bed',
                        required=True, help='path to BED file')
    parser.add_argument('-of', '--outputFasta', action='store', dest='output_fasta',
                        required=True, help='path to output Fasta')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    input_fasta = args.input_fasta
    input_bed = args.input_bed
    output_fasta = args.output_fasta
    
    peaks_fasta = bed_to_fasta(input_fasta, input_bed)
    write_fasta(peaks_fasta, output_fasta)


if __name__ == '__main__':
    main()
