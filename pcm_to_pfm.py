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
import random
import itertools
import argparse


def read_matrix(path):
    with open(path, 'r') as file:
        header = file.readline()
        header = header.strip()
        header = header.strip('>')
        matrix = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split()
            for letter, value in zip(matrix.keys(), line):
                matrix[letter].append(float(value))
    file.close()
    return(matrix, header)


def make_pfm_from_pcm(pcm):
    number_of_sites = [0] * len(pcm['A'])
    for key in pcm.keys():
        for i in range(len(pcm[key])):
            number_of_sites[i] += pcm[key][i]

    pfm = dict()
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pfm[i[0]] = []

    first_key = list(pcm.keys())[0]
    nuc_pseudo = 1/len(pcm.keys())
    for i in range(len(pcm[first_key])):
        for nuc in pcm.keys():
            pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
    return(pfm)


def write_pfm(output, header, pfm):
    with open(output, 'w') as file:
        file.write('>{0}\n'.format(header))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))
    pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', action='store',
                        help='path to PCM file (HOMER/HOCOMOCO like format)')
    parser.add_argument('output', action='store',
                        help='path to write PFM')
    return(parser.parse_args())


def main():
    args = parse_args()
    matrix_path = args.input
    output_path = args.output

    pcm, header = read_matrix(matrix_path)
    pfm = make_pfm_from_pcm(pcm)
    write_pfm(output_path, header, pfm)
    pass

if __name__ == '__main__':
    main()
