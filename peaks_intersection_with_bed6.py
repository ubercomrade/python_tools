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


def read_fantom_peaks(path):

    fantom = pd.read_csv(path,
                      sep='\t',comment='#', dtype= {'chr': str})
                     #names=['CAGE_peak_ID', 'short_description', 'description',
                     #       'chr', 'start', 'end', 'strand', 'ID'])
    fantom.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
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


def extend_border(record):

    right=5000
    left=-5000
    start = record[0]
    end = record[1]
    strand = record[2]

    if strand == '+':
        out_start = start + left
        out_end = start + right

    if strand == '-':
        out_start = end - right
        out_end = end - left
    return((out_start, out_end))


def overlap(peak, promoters):
    '''
    Does the range (start1, end1) overlap with (start2, end2)?
    Based on De Morgan's laws
    '''
    overlaps = promoters[np.logical_and(np.less_equal(peak['start'], promoters['end']),
                                 np.greater_equal(peak['end'], promoters['start']))]
    return(list(overlaps['name']))


def peaks_intersect_fantom(peaks, promoters):
    chrs_of_promoters = promoters['chr'].unique()
    chrs_of_peaks = peaks['chr'].unique()
    chrs = np.intersect1d(chrs_of_promoters, chrs_of_peaks)

    genes_id = []
    for chr_ in chrs:
        chr_peaks = pd.DataFrame(peaks[peaks['chr'] == chr_])
        chr_promoters = pd.DataFrame(promoters[promoters['chr'] == chr_])
        for index, peak in chr_peaks.iterrows():
            genes_id += overlap(peak, chr_promoters)
    return(genes_id)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fantom', action='store', dest='fantom',
                        required=True, help='path to fantom CAGE peaks')
    parser.add_argument('-p', '--peaks', action='store', dest='peaks',
                        required=True, help='path to peaks file')
    parser.add_argument('-o', '--output', action='store', dest='genes_id',
                                  required=True, help='path to txt file to write genes_id')
    parser.add_argument('-l', '--left', action='store', type=int, dest='left',
                                  default=-2000, required=False, help='left_tail + TSS, default_value = 2000')
    parser.add_argument('-r', '--right', action='store', type=int, dest='right',
                                  default=2000, required=False, help='TSS + right_tail, default_value = -2000')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def write_results(out, res):
    #res = [i.capitalize() for i in res]
    #res = [i.split('.')[0] for i in res]
    with open(out, 'w') as file:
        for gene_id in res:
            file.write(gene_id + '\n')


def main():

    args = parse_args()
    fantom_path = args.fantom
    peaks_path = args.peaks
    left = args.left
    right = args.right
    out = args.genes_id

    fantom = read_fantom_peaks(fantom_path)
    fantom = fantom.sort_values(by=['chr', 'start'])

    peaks = read_peaks(peaks_path)
    peaks = peaks.sort_values(by=['chr', 'start'])

    starts = list(fantom['start'])
    ends = list(fantom['end'])
    starts, ends = zip(*list(map(extend_border,
                                 zip(list(fantom['start']),
                                     list(fantom['end']),
                                     list(fantom['strand'])))))
    fantom['start'] = starts
    fantom['end'] = ends

    res = peaks_intersect_fantom(peaks, fantom)
    res = [i for i in res if isinstance(i, str)]
    res = [i.split(',') for i in res]
    res = [item for sublist in res for item in sublist]
    res = set(res)

    write_results(out, res)


if __name__ == '__main__':
    main()
