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
from collections import Counter


def read_sites(path):
    sequences = []
    with open(path, 'r') as file:
        sequences = [i.strip().upper() for i in file if i.strip()[0] != '>']
        sequences = [i for i in sequences if not 'N' in i]
    return(sequences)


def remove_equalent_seq(seq_list, homology=0.95):
    seq_list = list(seq_list)
    treshold = homology * len(seq_list[0])
    for seq1 in tuple(seq_list):
        sub_seq_list = list(seq_list)
        sub_seq_list.remove(seq1)
        for seq2 in sub_seq_list:
            score = len([i for i, j in zip(seq1, seq2) if i == j])
            if score >= treshold:
                seq_list.remove(seq1)
                break
    return(seq_list)


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
            pwm[j].append(math.log2(pfm[j][i] / background[j]))
    return(pwm)


def make_pcm(motifs):
    matrix = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        matrix[i[0]] = []
    len_of_motif = len(motifs[0])
    for i in matrix.keys():
        matrix[i] = [0]*len_of_motif
    for i in range(len_of_motif):
        for l in motifs:
            matrix[l[i]][i] += 1
    return(matrix)


def check_legth_of_sites(seq):
    length = [len(i) for i in seq]
    c = Counter(length).most_common(1)[0][0]
    seq = [i for i in seq if len(i) == c]
    return(seq)


def write_meme(output, tag, pfm, background, nsites):
    with open(output + '/' + tag + '.meme', 'w') as file:
        file.write('MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies\n')
        file.write('A {0} C {1} G {2} T {3}\n\n'.format(background['A'], background['C'],
                                                        background['G'], background['T']))
        file.write('MOTIF {0}\n'.format(tag))
        file.write(
            'letter-probability matrix: alength= 4 w= {0} nsites= {1}\n'.format(len(pfm['A']), nsites))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0:.8f}\t{1:.8f}\t{2:.8f}\t{3:.8f}\n'.format(i[0], i[1], i[2], i[3]))


def write_pwm(output, tag, pwm):
    with open(output + '/' + tag + '.pwm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pwm['A'], pwm['C'], pwm['G'], pwm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def write_pfm(output, tag, pfm):
    with open(output + '/' + tag + '.pfm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0:.9f}\t{1:.9f}\t{2:.9f}\t{3:.9f}\n'.format(i[0], i[1], i[2], i[3]))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to fasta')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='dir to MODEL output file')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    parser.add_argument('-M', '--meme', action='store_true', dest='meme',
                        required=False, help='if the flag is active model will be writen in meme format')
    parser.add_argument('-P', '--pwm', action='store_true', dest='pwm',
                        required=False, help='if the flag is active model will be writen in pwm format')
    parser.add_argument('-p', '--pfm', action='store_true', dest='pfm',
                        required=False, help='if the flag is active model will be writen in pfm format')
    return(parser.parse_args())


def main():
    args = parse_args()
    fasta_path = args.input
    output_dir = args.output
    tag = args.tag
    pwm_flag = args.pwm
    pfm_flag = args.pfm
    meme_flag = args.meme

    seq = read_sites(fasta_path)
    seq = list(set(seq))
    seq = check_legth_of_sites(seq)
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
    nsites = len(seq)

    pcm = make_pcm(seq)
    pfm = make_pfm_from_pcm(pcm)
    pwm = make_pwm_from_pcm(pcm, background)

    if pwm_flag: write_pwm(output_dir, tag, pwm)
    if pfm_flag: write_pfm(output_dir, tag, pfm)
    if meme_flag: write_meme(output_dir, tag, pfm, background, nsites)

if __name__ == '__main__':
    main()
