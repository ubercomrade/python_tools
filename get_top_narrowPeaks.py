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


def read_bed_cistrome(path):
    peaks = []
    with open(path, 'r') as file:
        for line in file:
            line = line.strip().split()
            peaks.append(line)
    file.close()
    return(peaks)


def get_top_peaks(peaks, amount):
    sorted_peaks = sorted(peaks, key=lambda i: float(i[6]), reverse=True)
    return(sorted_peaks[:amount])


def split_peaks(peaks):
    peaks1 = []  # for train
    peaks2 = []  # for test
    for i in range(len(peaks)):
        if i % 2 == 0:
            peaks1.append(peaks[i])
        else:
            peaks2.append(peaks[i])
    return(peaks1, peaks2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input_narrowPeaks',
                        required=True, help='path to BED file for reading')
    parser.add_argument('-o', '--output', action='store', dest='output_narrowPeaks',
                        required=True, help='dir to write output file/files')
    parser.add_argument('-a', '--amount', action='store', dest='amount',
                        default=100, required=False,
                        help='Amount of top peaks')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='tag for output file/files name')
    parser.add_argument('-s', '--split', action='store_true', dest='split',
                        default=False, required=False,
                        help='Get train and test peaks by spliting top peaks by 2')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    path = args.input_narrowPeaks
    output = args.output_narrowPeaks
    tag = args.tag
    amount = int(args.amount)
    split = args.split

    if not os.path.isdir(output):
        os.mkdir(output)

    if split:
        peaks = read_bed_cistrome(path)
        peaks = get_top_peaks(peaks, amount)
        peaks1, peaks2 = split_peaks(peaks)
        with open(output + '/' + tag + '_TRAIN' + '.narrowPeaks', 'w') as file:
            for i in peaks1:
                i[-1] = str(i[-1])
                file.write('\t'.join(i) + '\n')
        file.close()
        with open(output + '/' + tag + '_TEST' + '.narrowPeaks', 'w') as file:
            for i in peaks2:
                i[-1] = str(i[-1])
                file.write('\t'.join(i) + '\n')
        file.close()
    else:
        peaks = read_bed_cistrome(path)
        peaks = get_top_peaks(peaks, amount)
        with open(output + '/' + tag + '.narrowPeaks', 'w') as file:
            for i in peaks:
                i[-1] = str(i[-1])
                file.write('\t'.join(i) + '\n')
        file.close()
