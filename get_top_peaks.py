#!/usr/bin/python3

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
import pandas


def read_file(path):
    peaks = []
    with open(path, 'r') as file:
        for line in file:
            if line.isspace():
                continue
            else:
                line = line.strip().split()
                peaks.append(line)
    file.close()
    return(peaks)


def clear_peaks(peaks):
    peaks = [i for i in peaks if len(i[0]) < 6 and i[0] != "chrMT"]
    return(peaks)


def get_top_peaks(peaks, amount, col):
    #peaks = [i for i in peaks if len(i[0]) <= 5]
    sorted_peaks = sorted(peaks, key=lambda i: float(i[col]), reverse=True)
    for index, line in enumerate(sorted_peaks):
        line[3] = 'peaks_' + str(index)
    return(sorted_peaks[:amount])


def get_legths(data):
    l = list()
    for line in data:
        l.append(int(line[2]) - int(line[1]))
    return(l)


def split_peaks(peaks):
    peaks1 = []  # for train
    peaks2 = []  # for test
    for i in range(len(peaks)):
        if i % 2 == 0:
            peaks1.append(peaks[i])
        else:
            peaks2.append(peaks[i])
    return(peaks1, peaks2)


def write_length(path, data):
    with open(path, "w") as file:
        for line in data:
            file.write("{0}\n".format(line))
    file.close()
    pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to BED file for reading')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='dir to write output file/files')
    parser.add_argument('-a', '--amount', action='store', dest='amount',
                        default=100, required=False,
                        help='Amount of top peaks')
    parser.add_argument('-c', '--column', action='store', type=int, dest='column',
                        default=4, required=False,
                        help='number of colmun [0, 1, 2, ...] for sorting default = 4')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='tag for output file/files name')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    path = args.input
    output = args.output
    col = args.column
    tag = args.tag
    amount = int(args.amount)
    file_extension = os.path.splitext(path)[1]

    if not os.path.isdir(output):
        os.mkdir(output)

    peaks = read_file(path)
    peaks = clear_peaks(peaks)
    peaks = get_top_peaks(peaks, amount, col)
    l = get_legths(peaks)
    write_length(output + '/' + tag + '.length.txt', l)
    with open(output + '/' + tag + '.bed', 'w') as file:
        for i in peaks:
            i[-1] = str(i[-1])
            file.write('\t'.join(i) + '\n')
    file.close()


if __name__ == '__main__':
    main()
