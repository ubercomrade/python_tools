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
import pandas as pd
import numpy as np


def read_bed_like_file(path):
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'name', 'score', 'strand', 'site'])
    return(df)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--peaks', action='store', dest='input_peaks',
                        required=True, help='path to peaks in BED format or like BED format \
                        that contain all scaned peaks coordinates (sites is 7 column)')

    parser.add_argument('-o', '--out', action='store', dest='out_path',
                        required=True, help='path to out file')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def write_sites(peaks, out_path):
    #sites = list(peaks['site'])
    sites = list(peaks['site'][np.logical_not(pd.isna(peaks['site']))])

    #sites = [i for i in sites if not 'nan']
    with open(out_path, 'w') as file:
        for line in sites:
            file.write(line + '\n')
    file.close()


def main():
    args = parse_args()
    peaks_path = args.input_peaks
    out_path = args.out_path
    peaks = read_bed_like_file(peaks_path)
    #print(list(df['site']))
    write_sites(peaks, out_path)


if __name__=='__main__':
    main()

