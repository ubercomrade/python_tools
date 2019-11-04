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

import argparse
import sys
import os


def read_fasta(path):

    fasta = list()
    append = fasta.append
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                line = line.strip().upper()
                if not 'N' in line:
                    append(line)
                    append(complement(line))
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
        for key in pwm:
        	pwm[key] = tuple(pwm[key])
    file.close()
    return(pwm)


def score(seq, pwm):
    score = sum([pwm[letter][index] for index, letter in enumerate(seq)])
    return(score)


def calculate_scores(peaks, pwm, length_of_site):
    scores = [score(peak[i:length_of_site + i], pwm) for peak in peaks for i in range(len(peak) - length_of_site + 1)]
    return(scores)


def complement(seq):
    return(seq[::-1].translate(seq.maketrans('ACGT', 'TGCA')))
    #return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def get_threshold(scores, fp):
    scores.sort(reverse=True) # sorted score from big to small
    thr = scores[int(round(fp * len(scores)))]
    return(thr)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                            required=True, help='path to fasta')
    parser.add_argument('-m', '--pwm', action='store', dest='input_pwm',
                            required=True, help='path to PWM file')
    parser.add_argument('-p', '--false_positive', action='store', type=float, dest='false_positive',
                            required=True, help='value of FP (FP ~ P-VALUE)')
    parser.add_argument('-P', '--processes', action='store', type=int, dest='cpu_count',
    required=False, default=2, help='Number of processes to use, default: 2')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())
    

def main():
    #args = parse_args()

    #pwm_path = args.input_pwm
    #fasta_path = args.input_fasta
    #fp = args.false_positive
    #cpu_count = args.cpu_count

    fasta_path = "/home/anton/DATA/PROMOTERS/mm10.fa"
    pwm_path = "/home/anton/DATA/DB_DATA/TFS/CHOSEN/CHOSEN-TFS-SCAN-5000/AR_41593/MOTIFS/AR_41593_OPTIMAL_MOTIF.pwm"
    fp = 0.0001
    cpu_count = 2
    peaks = read_fasta(fasta_path)
    pwm = read_pwm(pwm_path)
    length_of_site = len(pwm['A'])
    scores = calculate_scores(peaks, pwm, length_of_site)

    thr = get_threshold(scores, fp)
    print(thr)


if __name__ == '__main__':
    main()
