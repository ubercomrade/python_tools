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


import os
import numpy as np
import itertools
import random
import multiprocessing as mp
import functools
import pandas as pd
import argparse


def parse_bamm_and_bg_from_file(bamm_file, bg_file):

    # Read BaMM file
    if os.path.isfile(bamm_file):
        motif_order = 0
        with open(bamm_file) as file:
            for line in file:
                if line[0] != '\n':
                    motif_order = motif_order + 1
                else:
                    break
        file.close()

        # count the motif length
        motif_length = int(sum(1 for line in open(bamm_file)) / (motif_order + 1))

        # read in bamm model
        model = {}
        for k in range(motif_order):
            model[k] = []

        with open(bamm_file) as file:
            for j in range(motif_length):
                for k in range(motif_order):
                    model[k].append([float(p) for p in file.readline().split()])
                file.readline()
        file.close()

        # convert a bamm array to numpy array
        for k in range(motif_order):
            model[k] = np.array(model[k], dtype=float)
    else:
        print('File {0} does not exist'.format(bg_file))
        exit()

    # Read BG file
    bg = {}
    order = 0
    if os.path.isfile(bg_file):
        with open(bg_file) as bgmodel_file:
            line = bgmodel_file.readline()  # skip the first line for K
            line = bgmodel_file.readline()  # skip the second line for Alpha
            while order < motif_order:
                line = bgmodel_file.readline()
                bg_freq = [float(p) for p in line.split()]
                bg[order] = np.array(bg_freq, dtype=float)
                order += 1
        file.close()
    else:
        print('File {0} does not exist'.format(bg_file))
        exit()
    return model, bg, order-1


def write_pfm(model, tag, output):
    pfm = model[0]
    with open(output + '/' + tag + '.pfm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in pfm:
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def write_meme(model, tag, bg, output):
    background = bg[0]
    pfm = model[0]
    with open(output + '/' + tag + '.meme', 'w') as file:
        file.write('MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies\n')
        file.write('A {0} C {1} G {2} T {3}\n\n'.format(background[0], background[1],
                                                        background[2], background[3]))
        file.write('MOTIF {0}\n'.format(tag))
        file.write(
            'letter-probability matrix: alength= 4 w= {0}\n'.format(len(pfm)))
        for i in pfm:
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--bamm', action='store', dest='input_bamm',
                        required=True, help='path to bamm file /format is .ihbcp/')
    parser.add_argument('-b', '--backgroud', action='store', dest='input_bg',
                        required=True, help='path to backgroud file /format is .hbcp/')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='dir to write output files')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    return(parser.parse_args())


def main():
    args = parse_args()
    bamm_path = args.input_bamm
    bg_path = args.input_bg
    output = args.output
    tag = args.tag

    model, bg, order = parse_bamm_and_bg_from_file(bamm_path, bg_path)
    write_pfm(model, tag, output)
    write_meme(model, tag, bg, output)

if __name__ == '__main__':
    main()
