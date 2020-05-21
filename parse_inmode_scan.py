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
import os
import argparse
import re


def read_fasta(path):
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            #print(line)
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                record = dict()
                record['name'] = int(line[0].split('_')[1])
                #record['name'] = line[0]

                record['chr'] = line[2]
                coordinates_strand = line[3]

                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                record['start'] = start
                record['end'] = end

                strand = re.findall(r'\(.\)', coordinates_strand[:-3])
                if not strand == []:
                    record['strand'] = strand[0].strip('()')
                else:
                    record['strand'] = '+'
            else:
                record['seq'] = line.strip().upper()
                fasta.append(record)
    file.close()
    return(fasta)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').replace('N', 'n').upper()[::-1])


def read_inmode_bed(path):
    container = []
    with open(path) as file:
        for line in file:
            n, start, end, strand, score = line.split()
            container.append({
                'chr': '.',
                'start': int(start),
                'end': int(end),
                'id': int(n),
                'score': math.log(float(score), 10),
                'strand': strand,
                'site': '.'
            })
    return(container)


def write_bed(data, threshold, path):
    with open(path, 'w') as file:
        for line in data:
            if line['score'] >= threshold:
                line = list(line.values())
                line = [str(i) for i in line]
                file.write('\t'.join(line) + '\n')
            else:
                continue
    return(0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--inputFasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file')
    parser.add_argument('-bed', '--inputBED', action='store', dest='input_bed',
                        required=True, help='path to BED file with inmode scan output')
    parser.add_argument('-o', '--outputBED', action='store', dest='output',
                        required=True, help='path to output BED with sites')
    parser.add_argument('-t', '--threshold', action='store', dest='threshold',
                        required=True, type=float, help='Threshold for data filtering')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    input_fasta = args.input_fasta
    input_bed = args.input_bed
    output = args.output
    threshold = args.threshold

    fasta = read_fasta(input_fasta)
    bed = read_inmode_bed(input_bed)

    for index, line in enumerate(bed):
        record = fasta[int(line['id'])]
        if line['strand'] == '-':
            line['site'] = complement(record['seq'][line['start']:line['end']])
        else:
            line['site'] = record['seq'][line['start']:line['end']]
        line['chr'] = record['chr']
        line['id'] = 'peaks_' + str(record['name'])
        line['start'] = int(line['start']) + int(record['start'])
        line['end'] = int(line['end']) + int(record['start'])
    write_bed(bed, threshold, output)

if __name__ == '__main__':
    main()
