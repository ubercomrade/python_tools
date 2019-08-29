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


def read_bed_file(path):

    bed = pd.read_csv(path,
                      sep='\t',comment='#', header=None)
    bed = bed.rename(columns={0:'chr', 1:'start', 2:'end', 3:'name', 5:'strand'})
    return(bed)


def read_peaks(path):
    # df = pd.read_csv(path,
    #                  sep='\t', header=None,
    #                  usecols=[0, 1, 2, 3, 4, 5, 6], dtype= {'chr': str},
    #                  names=['chr', 'start', 'end', 'name', 'score', 'strand', 'sites'])
    df = pd.read_csv(path, sep='\t', header=None)
    
    df = df.rename(columns={0:'chr', 1:'start', 2:'end', 3:'name', 5:'strand'})
    return(df)


def get_site_near_tss(record):

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


def overlap(peak, sites):
    '''
    Does the range (start1, end1) overlap with (start2, end2)?
    Based on De Morgan's laws
    '''
    overlaps = sites[np.logical_and(np.less_equal(peak['start'], sites['end']),
                                 np.greater_equal(peak['end'], sites['start']))]
    if overlaps.empty:
        return('-')
    else:
        #print(str(overlaps['name']))
        return(overlaps.iloc[0,3])


def peaks_intersect_with_sites(bed, sites):
    chrs_of_sites = sites['chr'].unique()
    chrs_of_bed = bed['chr'].unique()
    chrs = np.intersect1d(chrs_of_sites, chrs_of_bed)

    genes_id = []
    for chr_ in chrs:
        chr_bed = pd.DataFrame(bed[bed['chr'] == chr_])
        chr_sites = pd.DataFrame(sites[sites['chr'] == chr_])
        for index, site in chr_sites.iterrows():
            genes_id.append(overlap(site, chr_bed))
    return(genes_id)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store',
                        help='path to bed file, 4th column (name) is nessessary. \
                        See https://www.ensembl.org/info/website/upload/bed.html')
    parser.add_argument('peaks', action='store',
                        help='path to peaks file')
    parser.add_argument('output', action='store',
                        help='path to txt file to write genes_id')
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
        for ID in res:
            if ID != '-':
                file.write(ID + '\n')


def main():

    args = parse_args()
    bed_path = args.bed
    peaks_path = args.peaks
    left = args.left
    right = args.right
    out = args.output

    bed = read_bed_file(bed_path)
    bed = bed.sort_values(by=['chr', 'start'])

    peaks = read_peaks(peaks_path)
    peaks = peaks.sort_values(by=['chr', 'start'])

    starts = list(bed['start'])
    ends = list(bed['end'])
    starts, ends = zip(*list(map(get_site_near_tss,
                                 zip(list(bed['start']),
                                     list(bed['end']),
                                     list(bed['strand'])))))
    bed['start'] = starts
    bed['end'] = ends
    res = peaks_intersect_with_sites(bed, peaks)
    #res = [i for i in res if isinstance(i, str)]
    #res = [i.split(',') for i in res]
    #res = [item for sublist in res for item in sublist]
    #res = set(res)
    peaks['name'] = res
    #peaks.to_csv(out, sep='\t', header=None, index=False)
    write_results(out, res)


if __name__ == '__main__':
    main()
