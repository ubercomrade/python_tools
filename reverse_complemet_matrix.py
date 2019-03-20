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

import sys
import os
import argparse


def read_matrix(path):
    with open(path, 'r') as file:
        tag = file.readline().strip('>\n')
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(matrix, tag)


def read_meme(path):
    with open(path, 'r') as file:
        meme = file.readlines()
        inf_meme = meme[:9]
        pfm_meme = meme[9:]
        pfm_meme = [line.strip().split('\t') for line in pfm_meme]
        pfm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in pfm_meme:
            A,C,G,T = line
            pfm['A'].append(A)
            pfm['C'].append(C)
            pfm['G'].append(G)
            pfm['T'].append(T)
        tag = [i[6:] for i in inf_meme if i.startswith('MOTIF')][0].strip()
    file.close()
    return(pfm, inf_meme, tag)


def reverse_complemet_matrix(matrix):
    reverse_complement_matrix = {'A':matrix['T'][::-1],
                                'C':matrix['G'][::-1],
                                'G':matrix['C'][::-1],
                                'T':matrix['A'][::-1]}
    return(reverse_complement_matrix)


def write_meme(output, inf, tag, pfm):
    with open(output + '/' + tag + '.meme', 'w') as file:
        for line in inf:
            file.write(line)
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def write_pwm(output, tag, pwm):
    with open(output + '/' + tag + '.pwm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pwm['A'], pwm['C'], pwm['G'], pwm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def write_pfm(output, tag, pfm):
    with open(output + '/' + tag + '.pfm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to ChIPMunk results')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='dir to write output files')
    parser.add_argument('-m', '--meme', action='store_true', dest='flag_meme',
                        required=False, help='use this flag if matrix in meme format')
    parser.add_argument('-p', '--pwm', action='store_true', dest='flag_pwm',
                        required=False, help='use this flag if matrix is pwm')
    return(parser.parse_args())


def main():
    args = parse_args()
    path = args.input
    output = args.output
    meme_flag = args.flag_meme
    pwm_flag = args.flag_pwm

    if meme_flag:
        pfm, inf_meme, tag = read_meme(path)
        tag += '_RC'
        pfm = reverse_complemet_matrix(pfm)
        write_meme(output, inf_meme, tag, pfm)
    elif pwm_flag:
        pwm, tag = read_matrix(path)
        tag += '_RC'
        pwm = reverse_complemet_matrix(pwm)
        write_pwm(output, tag, pwm)
    else:
        pfm, tag = read_matrix(path)
        tag += '_RC'
        pfm = reverse_complemet_matrix(pfm)
        write_pfm(output, tag, pfm)

if __name__ == '__main__':
    main()
