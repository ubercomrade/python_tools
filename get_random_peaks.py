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
import random
import sys


def read_bed(path):
    peaks = []
    with open(path, 'r') as file:
        for line in file:
            peaks.append(line)
    file.close()
    return(peaks)


def get_random_peaks(peaks, amount):
    random_peaks = random.sample(peaks, k=amount)
    return(random_peaks)


def get_other_peaks(peaks, random_peaks):
    other_peaks = [peak for peak in peaks if peak not in random_peaks]
    retunr(other_peaks)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputBED', action='store', dest='input_bed',
                        required=True, help='path to BED file for reading')
    parser.add_argument('-o', '--outputBED', action='store', dest='output_bed',
                        required=True, help='path for writing BED file')
    parser.add_argument('-oo', '--output_other_BED', action='store', dest='output_other_bed',
                        required=False, help='path for writing BED file for ALL_PEAKS - RANDOM_PEAKS')
    parser.add_argument('-a', '--amount', action='store', dest='amount',
                        default=100, required=False,
                        help='Amount of top peaks')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    path = args.input_bed
    output = args.output_bed
    other = args.output_other_bed
    amount = int(args.amount)

    peaks = read_bed(path)
    random_peaks = get_random_peaks(peaks, amount)
    with open(output, 'w') as file:
        for i in random_peaks:
            file.write(i)
    file.close()

    if not other is None:
        other_peaks = get_other_peaks(peaks, random_peaks)
        with open(other, 'w') as file:
            for i in other_peaks:
                file.write(i)
        file.close()
