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
    #return(seq[::-1].translate(seq.maketrans('ACGT', 'TGCA')))
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


# def get_threshold(scores, fpr):
#     scores.sort(reverse=True) # sorted score from big to small
#     thr = scores[int(round(fpr * len(scores)))]
#     return(thr)


def get_threshold(scores, path_out):
    scores.sort(reverse=True) # sorted score from big to small
    fprs = [5*10**(-4), 3.33*10**(-4), 1.90*10**(-4), 1.02*10**(-4), 5.24*10**(-5)]
    with open(path_out, "w") as file:
        file.write("Scores\tFPR\n")
        for fpr in fprs:
            thr = scores[int(fpr * len(scores))]
            file.write("{0}\t{1}\n".format(thr, fpr))
    file.close()


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('fasta', action='store', help='path to fasta')
    parser.add_argument('pwm', action='store', help='path to PWM file')
    parser.add_argument('out', action='store', help='path to write results')
    parser.add_argument('-p', '--false_positive', action='store', type=float, dest='false_positive',
                            required=False, help='value of FP (FP ~ P-VALUE)')


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())
    

def main():
    args = parse_args()

    pwm_path = args.pwm
    fasta_path = args.fasta
    path_out = args.out
    fp = args.false_positive

    peaks = read_fasta(fasta_path)
    pwm = read_pwm(pwm_path)
    length_of_site = len(pwm['A'])
    scores = calculate_scores(peaks, pwm, length_of_site)
    get_threshold(scores, path_out)


if __name__ == '__main__':
    main()
