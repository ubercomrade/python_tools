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


import pandas as pd
import numpy as np
import argparse
import sys


def read_fantom_peaks(path, right=500, left=-500):

    fantom = pd.read_csv(path,
                         sep='\t',
                         comment='#',
                         dtype= {'chr': str},
                         names=['chr', 'chromStart', 'chromEnd', 'name', 'score',
                            'strand', 'thickStart', 'thickEnd','itemRgb',
                            'blockCount', 'blockSizes', 'blockStarts'])
    fantom['thickStart'] += left
    fantom['thickEnd'] += right
    return(fantom)


def slpit_attribute(record):
    record = record.split('; ')
    rec = []
    for i in record:
        if i == '':
            continue
        (a,b) = i.strip().split(maxsplit=1)
        rec.append((a.strip(),b.strip(';\" ')))
    return(dict(rec))


def read_peaks(path):
    df = pd.read_csv(path,
                     sep='\t', header=None,
                     usecols=[0, 1, 2, 3, 4, 5], dtype= {'chr': str},
                     names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    return(df)


def overlap(peak, enchancers):
    '''
    Does the range (start1, end1) overlap with (start2, end2)?
    Based on De Morgan's laws
    '''
    overlaps = enchancers[np.logical_and(np.less_equal(peak['start'], enchancers['thickEnd']),
                                 np.greater_equal(peak['end'], enchancers['thickStart']))]
    return(list(overlaps['name']))


def peaks_intersect_fantom(peaks, enchancers):
    chrs_of_enchancers = enchancers['chr'].unique()
    chrs_of_peaks = peaks['chr'].unique()
    chrs = np.intersect1d(chrs_of_enchancers, chrs_of_peaks)

    genes_id = []
    for chr_ in chrs:
        chr_peaks = pd.DataFrame(peaks[peaks['chr'] == chr_])
        chr_enchancers = pd.DataFrame(enchancers[enchancers['chr'] == chr_])
        for index, peak in chr_peaks.iterrows():
            genes_id += overlap(peak, chr_enchancers)

    genes_id = [i.split(';')[1:-2] for i in genes_id]
    IDs = [j for i in genes_id for j in i if not j.startswith('NM_')]
    return(set(IDs))


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fantom', action='store', dest='fantom',
                        required=True, help='path to fantom file in BED12 format with ecnhancer')
    parser.add_argument('-p', '--peaks', action='store', dest='peaks',
                        required=True, help='path to peaks file')
    parser.add_argument('-o', '--output', action='store', dest='genes_id',
                                  required=True, help='path to txt file to write genes_id')
    parser.add_argument('-l', '--left', action='store', type=int, dest='left',
                                  default=-500, required=False, help='left_tail + enchancer, default_value = 500')
    parser.add_argument('-r', '--right', action='store', type=int, dest='rigth',
                                  default=500, required=False, help='enchancer + rigth_tail, default_value = -500')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def write_results(out, res):
    #res = [i.capitalize() for i in res]
    res = [i.split('.')[0] for i in res]
    with open(out, 'w') as file:
        for gene_id in res:
            file.write(gene_id + '\n')


def main():
    args = parse_args()
    fantom = args.fantom
    peaks_path = args.peaks
    left = args.left
    rigth = args.rigth
    out = args.genes_id

    enchancers = read_fantom_peaks(fantom, rigth, left)
    peaks = read_peaks(peaks_path)


    peaks = read_peaks(peaks_path)
    peaks = peaks.sort_values(by=['chr', 'start'])
    res = peaks_intersect_fantom(peaks, enchancers)
    #res = np.unique(res)
    write_results(out, res)


if __name__ == '__main__':
    main()
