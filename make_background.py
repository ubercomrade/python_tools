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


from random import sample
import argparse
import sys


def read_fasta(path):
    with open(path, 'r') as fasta:
        out = []
        for line in fasta:
            if line.startswith('>'):
                continue
            else:
                out.append(line.strip().upper())
    fasta.close()
    return(out)


def shuffle_fasta(fasta, times):
    out = []
    for seq in fasta:
        shuffled = [''.join(sample(seq, k=len(seq))) for i in range(times)]
        out += shuffled
    return(out)


def write(fasta, path):
    with open(path, 'w') as file:
        for index, seq in enumerate(fasta):
            file.write('>Background_' + str(index) + '\n')
            file.write(seq + '\n')
    file.close


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to FASTA file')
    parser.add_argument('-t', '--times', action='store', type=int, dest='times',
                        required=False, default=2, help='times bigger, default = 2')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='path to write output file with shuffled sequences')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    ifasta = args.input
    ofasta = args.output
    times = args.times

    fasta = read_fasta(ifasta)
    res = shuffle_fasta(fasta, times)
    write(res, ofasta)


if __name__ == '__main__':
    main()
