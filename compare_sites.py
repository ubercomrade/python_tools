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


# def read_peaks(path):
#    df = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2, 3, 4, 5],
#                    names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
#    if df['name'][0] == '.':
#        name = 'peaks_'
#        names = [name + str(i) for i in range(len(df))]
#        df['name'] = names
#    # print(df)
#    return(df)


def read_peaks(path):
    df = pd.read_csv(path,
                     sep='\t', header=None,
                     #usecols=[0, 1, 2, 3, 4, 5],
                     names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    #df.loc[np.isnan(df['strand']), 'strand'] = '.'
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


def calculate_fraction_of_sites(sub_data, sub_pwm_sites, sub_bamm_sites):

    # Varable for write results
    results = {'peaks': int(), 'no_sites': int(), 'pwm_sites': int(), 'bamm_sites': int(),
               'overlap_sites': int(), 'not_overlap_sites': int()}
    results['peaks'] = len(sub_data)
    results['no_sites'] = len(set(sub_data['name']) -
                              set(sub_bamm_sites['name']) - set(sub_pwm_sites['name']))
    results['bamm_sites'] = len(set(sub_bamm_sites['name']) - set(sub_pwm_sites['name']))
    results['pwm_sites'] = len(set(sub_pwm_sites['name']) - set(sub_bamm_sites['name']))

    # Caclculation fraction of sites in peaks
    both_method_peaks = set(sub_bamm_sites['name']) & set(sub_pwm_sites['name'])
    common_peaks_results = pd.DataFrame()

    for peak in both_method_peaks:
        subset_bamm = pd.DataFrame(sub_bamm_sites[np.logical_and(
            sub_bamm_sites['name'] == peak, sub_bamm_sites['strand'] == '+')])
        subset_pwm = pd.DataFrame(sub_pwm_sites[np.logical_and(
            sub_pwm_sites['name'] == peak, sub_pwm_sites['strand'] == '+')])
        intersections_down = [overlap(row_pwm, row_bamm) for index_pwm, row_pwm in subset_pwm.iterrows(
        ) for index_bamm, row_bamm in subset_bamm.iterrows()]

        subset_bamm = pd.DataFrame(sub_bamm_sites[np.logical_and(
            sub_bamm_sites['name'] == peak, sub_bamm_sites['strand'] == '-')])
        subset_pwm = pd.DataFrame(sub_pwm_sites[np.logical_and(
            sub_pwm_sites['name'] == peak, sub_pwm_sites['strand'] == '-')])
        intersections_up = [overlap(row_pwm, row_bamm) for index_pwm, row_pwm in subset_pwm.iterrows(
        ) for index_bamm, row_bamm in subset_bamm.iterrows()]

        intersections = intersections_up + intersections_down
        if sum(intersections) == 0:
            results['not_overlap_sites'] += 1
        else:
            results['overlap_sites'] += 1
    return(results)


def get_sub_data(top, peaks, sites):
    sub_data = peaks[:top]
    names = set(sub_data['name'])
    sub_sites = sites[[i in names for i in sites['name']]]
    return(sub_sites)


def get_top(top, peaks):
    sub_data = peaks[:top]
    return(sub_data)


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
    # MAIN
    #bamm_path = '~/DATA/TF/CEPBA_mm10_GSM2845732/SCAN/CEBPA_all_BAMM.tsv'
    #pwm_path = '~/DATA/TF/CEPBA_mm10_GSM2845732/SCAN/CEBPA_all_PWM.tsv'
    #peaks_path = '~/DATA/TF/CEPBA_mm10_GSM2845732/BED/CEBPA_all.bed'
    #tag = '123'
    #out_dir = '~/DATA/TEST'

    args = parse_args()
    bamm_path = args.bamm_path
    pwm_path = args.pwm_path
    peaks_path = args.input_peaks
    tag = args.tag
    out_dir = args.out_dir

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # Read bamm sites
    bamm_sites = read_bed_like_file(bamm_path)
    bamm_sites['type'] = 'bamm'
    bamm_sites['overlaps'] = False

    # Read pwm sites
    pwm_sites = read_bed_like_file(pwm_path)
    pwm_sites['type'] = 'pwm'
    pwm_sites['overlaps'] = False

    # Read peaks and give names to every peask
    data = read_peaks(peaks_path)

    # Get subdata by peak score (TOP 1000, 2000, ...)
    top = [i * 1000 for i in range(1, len(data) // 1000 + 1)]
    pwm_sub_sites = [get_sub_data(i, data, pwm_sites) for i in top]
    bamm_sub_sites = [get_sub_data(i, data, bamm_sites) for i in top]
    sub_data = [get_top(i, data) for i in top]

    data = zip(sub_data, pwm_sub_sites, bamm_sub_sites)
    with Pool(4) as p:
        total_results = p.starmap(calculate_fraction_of_sites, data)

    # Count data
    count = pd.DataFrame(total_results)
    count = count[['overlap_sites', 'pwm_sites', 'bamm_sites',
                   'not_overlap_sites', 'no_sites', 'peaks']]
    count.to_csv(out_dir + '/' + tag + '_COUNT.tsv', sep='\t', index=False)

    # Fraction data devided by fractiomal volume of peaks
    frequency = pd.DataFrame(total_results)
    frequency = frequency[['overlap_sites', 'pwm_sites',
                           'bamm_sites', 'not_overlap_sites', 'no_sites', 'peaks']]
    for column in frequency:
        if column == 'peaks':
            continue
        frequency[column] = frequency[column] / frequency['peaks']
    frequency.to_csv(out_dir + '/' + tag + '_FREQUENCY.tsv', sep='\t', index=False)

    # Picture
    ax = frequency.plot.bar(stacked=True, x='peaks')
    ax.legend(loc=[1.05, 0.62])
    fig = ax.get_figure()
    fig.savefig(out_dir + '/' + tag + '_PIC.png', dpi=150,  bbox_inches='tight')


if __name__ == '__main__':
    main()
