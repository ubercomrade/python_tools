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
import argparse
import numpy as np
import pandas as pd


def statistic(fasta_path):
    stat = dict()
    with open(fasta_path, 'r') as fasta:
        for line in fasta:
            if line.startswith('>'):
                chr_ = line[1:].strip().split()[0]
                stat[chr_] = 0
                continue
            else:
                stat[chr_] += len(line.strip())
    fasta.close()
    return(stat)


def write_table(stat, out_path):
    with open(out_path, 'w')  as tsv:
        for chr, length in stat.items():
            tsv.write(str(chr) + '\t' + str(length) + '\n')
    tsv.close()




def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--fasta', action='store', dest='fasta',
                        required=True, help='path to genome in fasta format')
    parser.add_argument('-o', '--output', action='store', dest='tsv',
                                  required=True, help='path to read tsv file')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    fasta_path = args.fasta
    out_path = args.tsv

    stat = statistic(fasta_path)
    write_table(stat, out_path)

if __name__ == '__main__':
    main()
