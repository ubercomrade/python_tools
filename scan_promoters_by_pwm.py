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
import csv
import multiprocessing as mp
import functools
import re
import random


def read_fasta(path):

    fasta = list()
    append = fasta.append
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                line = line.strip().upper()
                #if 'N' in line:
                #    line = replace_n_to_randome_nucl(line)
                append(line)
                append(complement(line))
    file.close()
    return(fasta)


def replace_n_to_randome_nucl(peak):
    container = []
    for n in peak:
        if n == 'N':
            container.append(random.choice(['A', 'C', 'G', 'T']))
    return(''.join(container))



def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(pwm)


def score(seq, pwm):
    length_of_seq = len(seq)
    position = 0
    score = 0
    for letter in seq:
        score += pwm[letter][position]
        position += 1
    return(score)


def complement(seq):
    seq = seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
    return(seq)


def scan_seq_by_pwm(seq, pwm, threshold):
    results = []
    length_pwm = len(pwm['A'])
    # first strand
    counter = int()
    for i in range(len(seq) - length_pwm + 1):
        site = seq[i:length_pwm + i]
        if 'N' in site:
            continue
        s = score(site, pwm)
        counter += 1
        if s >= threshold:
            results.append(site)
    return(results, counter)


def write_sites(path, data):
    with open(path, 'w') as file:
        for line in data:
            file.write('{}\n'.format(line))
    pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA')
    parser.add_argument('-m', '--pwm', action='store', dest='input_pwm',
                        required=True, help='path to PWM file')
    parser.add_argument('-t', '--threshold', action='store', type=float, dest='threshold',
                        required=True, help='threshold for PWM')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='path to write sites')
    return(parser.parse_args())


def main():

    args = parse_args()
    pwm_path = args.input_pwm
    fasta_path = args.input_fasta
    threshold = args.threshold
    results_path = args.output

    fasta = read_fasta(fasta_path)
    pwm = read_pwm(pwm_path)
    results = []
    number_of_sites = 0
    for record in fasta:
      add_results, add_counter = scan_seq_by_pwm(record, pwm, threshold)
      results += add_results
      number_of_sites += add_counter
    length_of_site = len(pwm['A'])
    #number_of_sites = sum([len(range(len(peak) - length_of_site + 1)) for peak in fasta])
    print(len(results), number_of_sites, len(results) / number_of_sites)
    write_sites(results_path, results)

if __name__ == '__main__':
    main()
