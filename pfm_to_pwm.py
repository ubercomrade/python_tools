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


import itertools
import argparse
import sys
import math


def read_pfm(path):
    with open(path, 'r') as file:
        header = file.readline()
        pfm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split()
            for letter, value in zip(pfm.keys(), line):
                pfm[letter].append(float(value))
    file.close()
    return(pfm, header)


def make_pcm_from_pfm(pfm):
    matrix = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        matrix[i[0]] = []
    pseudo_count = 10**9
    for key in pfm.keys():
        for i in pfm[key]:
            matrix[key].append(round(i * pseudo_count))
    return(matrix)


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


def make_pwm_from_pcm(pcm, background):
    pwm = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pwm[i[0]] = []

    pfm = make_pfm_from_pcm(pcm)
    first_key = list(pcm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log10(pfm[j][i] / background[j]))
    return(pwm)


def background_freq(seq):
    s = ''.join(seq)
    background = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        background[i[0]] = s.count(i[0])
    sum_of_nuc = sum(background.values())
    for i in background.keys():
        background[i] = background[i]/sum_of_nuc
    return(background)


def read_fasta(path):
    fasta = []
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            else:
                fasta.append(line.strip().upper())
    file.close()
    return(fasta)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to PFM file')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='path to PWM output file')
    return(parser.parse_args())


def main():
    args = parse_args()
    input_homer = args.input
    output_pwm = args.output

    pfm, header = read_pfm(input_homer)
    pcm = make_pcm_from_pfm(pfm)
    background = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    pwm = make_pwm_from_pcm(pcm, background)
    with open(output_pwm, 'w') as file:
        file.write(header)
        for i in zip(pwm['A'], pwm['C'], pwm['G'], pwm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))
    file.close()


if __name__ == '__main__':
    main()
