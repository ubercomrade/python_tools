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
from multiprocessing import Pool
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_bed_like_file(path):
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site'])
    return(df)


def read_peaks(path):
    df = pd.read_csv(path,
                     sep='\t', header=None,
                     usecols=[0, 1, 2, 3, 4, 5],
                     names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    df['strand'] = '.'
    if df['name'][0] == '.':
        name = 'peaks_'
        names = [name + str(i) for i in range(len(df))]
        df['name'] = names
    return(df)


def read_fasta_headers(path):
    '''
    Чтение фаста фаила и запись хедеров
    Шапка для FASTA: >uniq_id|chromosome|start-end|strand
    '''
    fasta_headers = list()
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                line = line[1:].strip().split('|')
                record = dict()
                record['name'] = line[0]
                record['chromosome'] = line[1]
                record['start'] = line[2].split('-')[0]
                record['end'] = line[2].split('-')[1]
                try:
                    record['strand'] = line[3]
                except:
                    #print('Record with out strand. Strand is +')
                    record['strand'] = '+'
                continue
            fasta_headers.append(record)
    file.close()
    df = pd.DataFrame(fasta_headers)
    df = df[['chromosome', 'start', 'end', 'name', 'strand']]
    return(df)


def overlap(site_1, site_2):
    '''
    Does the range (start1, end1) overlap with (start2, end2)?
    Based on De Morgan's laws
    '''
    return site_1['end'] >= site_2['start'] and site_2['end'] >= site_1['start']


def get_coord(site):
    return(np.array(site['start'], site['end']))



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--peaks', action='store', dest='input_peaks',
                        required=True, help='path to peaks in BED format or like BED format that contain all scaned peaks coordinates')
    parser.add_argument('-m', '--pwm', action='store', dest='pwm_path',
                        required=True, help='path to file that contain sites obtained by PWM in peaks file')
    parser.add_argument('-b', '--bamm', action='store', dest='bamm_path',
                        required=True, help='path to file that contain sites obtained by BAMM in peaks file')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    parser.add_argument('-o', '--out', action='store', dest='out_dir',
                        required=True, help='OUT_DIR')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():

    args = parse_args()
    bamm_path = args.bamm_path
    pwm_path = args.pwm_path
    peaks_path = args.input_peaks
    tag = args.tag
    out_dir = args.out_dir

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    data = read_peaks(peaks_path)
    names = [int(i.split('_')[1]) for i in data['name']]
    data['name'] = names

    bamm_sites = read_bed_like_file(bamm_path)
    bamm_sites['type'] = 'bamm'
    names = [int(i.split('_')[1]) for i in bamm_sites['name']]
    bamm_sites['name'] = names

    pwm_sites = read_bed_like_file(pwm_path)
    pwm_sites['type'] = 'pwm'
    names = [int(i.split('_')[1]) for i in pwm_sites['name']]
    pwm_sites['name'] = names


    id_pwm = np.unique(pwm_sites['name']) # ids of peaks where there are pwm_sites
    id_bamm = np.unique(bamm_sites['name']) # ids of peaks where there are bamm_sites

    only_pwm_sites_id = np.setdiff1d(id_pwm, id_bamm) # ids of peaks where there are only pwm_sites
    only_bamm_sites_id = np.setdiff1d(id_bamm, id_pwm) # ids of peaks where there are only bamm_sites

    only_pwm_sites = pwm_sites[np.in1d(pwm_sites['name'], only_pwm_sites_id)] # df with only pwm_sites
    only_bamm_sites = bamm_sites[np.in1d(bamm_sites['name'], only_bamm_sites_id)]  # df with only only bamm_sites

    only_pwm_peaks = data[np.in1d(data['name'], only_pwm_sites_id)] # df with only pwm_peaks
    only_bamm_peaks = data[np.in1d(data['name'], only_bamm_sites_id)]  # df with only only bamm_peaks

    id_no_sites = np.setdiff1d(np.setdiff1d(np.array(data['name']), id_bamm), id_pwm) # id of peaks where there are not sites
    no_sites_peaks = data[np.in1d(data['name'], id_no_sites)] # peaks where there are not sites


    #######
    #Get peaks with overlap and not overlaped sites
    perhaps_overlap_peaks = data[np.in1d(data['name'], np.intersect1d(id_pwm, id_bamm))]
    peaks_with_not_overlap_sites = pd.DataFrame()
    peaks_with_overlap_sites = pd.DataFrame()

    for index, peak in perhaps_overlap_peaks.iterrows():
        subset_bamm = pd.DataFrame(bamm_sites[np.logical_and(
            bamm_sites['name'] == peak['name'], bamm_sites['strand'] == '+')])
        subset_pwm = pd.DataFrame(pwm_sites[np.logical_and(
            pwm_sites['name'] == peak['name'], pwm_sites['strand'] == '+')])
        intersections_down = [overlap(row_pwm, row_bamm) for index_pwm, row_pwm in subset_pwm.iterrows(
            ) for index_bamm, row_bamm in subset_bamm.iterrows()]

        subset_bamm = pd.DataFrame(bamm_sites[np.logical_and(
            bamm_sites['name'] == peak['name'], bamm_sites['strand'] == '-')])
        subset_pwm = pd.DataFrame(pwm_sites[np.logical_and(
            pwm_sites['name'] == peak['name'], pwm_sites['strand'] == '-')])
        intersections_up = [overlap(row_pwm, row_bamm) for index_pwm, row_pwm in subset_pwm.iterrows(
            ) for index_bamm, row_bamm in subset_bamm.iterrows()]

        intersections = intersections_up + intersections_down
        if sum(intersections) == 0:
            peaks_with_not_overlap_sites = peaks_with_not_overlap_sites.append(peak)
        else:
            peaks_with_overlap_sites = peaks_with_overlap_sites.append(peak)

    peaks_with_not_overlap_sites['name'] =  [int(i) for i in peaks_with_not_overlap_sites['name']]
    peaks_with_overlap_sites['name'] =  [int(i) for i in peaks_with_overlap_sites['name']]
    ########

    ########
    #Make table with count of diff kind of peaks
    count = []
    top = [i * 1000 for i in range(1, len(data) // 1000 + 1)]

    for i in range(len(top)):
        subset_name_peaks = data.iloc[i*1000:(i+1)*1000]['name']
        count_pwm_sites = len(np.intersect1d(only_pwm_sites['name'], subset_name_peaks))
        count_bamm_sites = len(np.intersect1d(only_bamm_sites['name'], subset_name_peaks))
        count_overlap_sites = len(np.intersect1d(peaks_with_overlap_sites['name'], subset_name_peaks))
        count_not_overlap_sites = len(np.intersect1d(peaks_with_not_overlap_sites['name'], subset_name_peaks))
        count_no_sites = len(np.intersect1d(no_sites_peaks['name'], subset_name_peaks))

        count.append({'no_sites': count_no_sites,
         'pwm_sites': count_pwm_sites, 'bamm_sites': count_bamm_sites,
         'overlap_sites': count_overlap_sites,
         'not_overlap_sites': count_not_overlap_sites})

    count_ = pd.DataFrame(count)
    count = pd.DataFrame()
    count = count.append(count_.iloc[0])
    for i in range(1,len(count_)):
        count = count.append(count.iloc[i - 1] + count_.iloc[i], ignore_index=True)
    count['peaks'] = top
    count = count[['overlap_sites', 'pwm_sites', 'bamm_sites',
                       'not_overlap_sites', 'no_sites', 'peaks']]
    count.to_csv(out_dir + '/' + tag + '_COUNT.tsv', sep='\t', index=False)

    frequency = pd.DataFrame(count)
    for column in frequency:
        if column == 'peaks':
            continue
        frequency[column] = frequency[column] / frequency['peaks']
    frequency.to_csv(out_dir + '/' + tag + '_FREQUENCY.tsv', sep='\t', index=False)

    ax = frequency.plot.bar(stacked=True, x='peaks')
    ax.legend(loc=[1.05, 0.62])
    fig = ax.get_figure()
    fig.savefig(out_dir + '/' + tag + '_PIC.png', dpi=150,  bbox_inches='tight')


    only_bamm_sites = only_bamm_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    only_bamm_sites.to_csv(out_dir + '/' + tag + '_bamm.sites', sep='\t', index=False, header=False)

    only_pwm_sites = only_pwm_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    only_pwm_sites.to_csv(out_dir + '/' + tag + '_pwm.sites', sep='\t', index=False, header=False)

    only_bamm_peaks = only_bamm_peaks[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    only_bamm_peaks.to_csv(out_dir + '/' + tag + '_bamm.peaks', sep='\t', index=False, header=False)

    only_pwm_peaks = only_pwm_peaks[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    only_pwm_peaks.to_csv(out_dir + '/' + tag + '_pwm.peaks', sep='\t', index=False, header=False)

    peaks_with_overlap_sites = peaks_with_overlap_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    peaks_with_overlap_sites.to_csv(out_dir + '/' + tag + '_overlap_sites.peaks', sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()
