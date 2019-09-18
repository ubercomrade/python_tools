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
import pandas as pd
import numpy as np


def read_fasta(path):

    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                fasta.append(line.strip().upper())
    return(fasta)


def complement(seq):
    '''
    Make reverse and compelent
    '''
    output = str()
    for letter in seq:
        if letter == 'A':
            output += 'T'
        elif letter == 'C':
            output += 'G'
        elif letter == 'G':
            output += 'C'
        elif letter == 'T':
            output += 'A'
    output = output[::-1]
    return(output)


def write_results(path, data):
    with open(path, 'w') as file:
        for line in data:
            file.write(line + '\n')
    file.close()
    pass


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('input', action='store',
                        help='path to file with motifs in format: site\nsite\n...')
    parser.add_argument('output', action='store',
                        help='path to write results')
    parser.add_argument('-l', '--left', action='store', dest='left',
                        required=False, type=int, default=0, help='start position of motif in sites (def first = 0)')
    parser.add_argument('-r', '--right', action='store', dest='right',
                        required=False, type=int, default=-1, help='end position of motif in sites (def last = -1)')
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    input_path = args.input
    output_path = args.output
    left = args.left
    right = args.right
    
    container = list()
    sites = read_fasta(input_path)
    for site in sites:
        container.append(complement(site[left:right]))
    write_results(output_path, container)
    
if __name__=='__main__':
    main()
    