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

    try:
        bed = pd.read_csv(path,
                          sep='\t', header=None,
                          usecols=[0, 1, 2, 3, 4, 5],
                          dtype= {'chromosome': str},
                          names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    except:
        bed = pd.read_csv(path,
                          sep='\t', header=None,
                          usecols=[0, 1, 2, 3, 4],
                          dtype= {'chromosome': str},
                          names=['chromosome', 'start', 'end', 'name', 'score'])

        bed['strand'] = '.'

    if bed['name'][0] == '.':
        name = 'peaks_'
        names = [name + str(i) for i in range(len(df))]
        bed['name'] = names
    return(bed)


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
    parser.add_argument('-m', '--first_model', action='store', dest='first_model_path',
                        required=True, help='path to file that contain sites obtained by first_model in peaks file')
    parser.add_argument('-b', '--second_model', action='store', dest='second_model_path',
                        required=True, help='path to file that contain sites obtained by_second_model in peaks file')
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
    second_model_path = args.second_model_path
    first_model_path = args.first_model_path
    peaks_path = args.input_peaks
    tag = args.tag
    out_dir = args.out_dir

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    data = read_peaks(peaks_path)
    names = [int(i.split('_')[1]) for i in data['name']]
    data['name'] = names

    second_model_sites = read_bed_like_file(second_model_path)
    second_model_sites['type'] = 'second_model'
    names = [int(i.split('_')[1]) for i in second_model_sites['name']]
    second_model_sites['name'] = names

    first_model_sites = read_bed_like_file(first_model_path)
    first_model_sites['type'] = 'first_model'
    names = [int(i.split('_')[1]) for i in first_model_sites['name']]
    first_model_sites['name'] = names


    id_first_model = np.unique(first_model_sites['name']) # ids of peaks where there are first_model_sites
    id_second_model = np.unique(second_model_sites['name']) # ids of peaks where there are second_model_sites

    only_first_model_sites_id = np.setdiff1d(id_first_model, id_second_model) # ids of peaks where there are only first_model_sites
    only_second_model_sites_id = np.setdiff1d(id_second_model, id_first_model) # ids of peaks where there are only second_model_sites

    only_first_model_sites = first_model_sites[np.in1d(first_model_sites['name'], only_first_model_sites_id)] # df with only first_model_sites
    only_second_model_sites = second_model_sites[np.in1d(second_model_sites['name'], only_second_model_sites_id)]  # df with only only second_model_sites

    only_first_model_peaks = data[np.in1d(data['name'], only_first_model_sites_id)] # df with only first_model_peaks
    only_second_model_peaks = data[np.in1d(data['name'], only_second_model_sites_id)]  # df with only only_second_model_peaks

    id_no_sites = np.setdiff1d(np.setdiff1d(np.array(data['name']), id_second_model), id_first_model) # id of peaks where there are not sites
    no_sites_peaks = data[np.in1d(data['name'], id_no_sites)] # peaks where there are not sites


    #######
    #Get peaks with overlap and not overlaped sites
    perhaps_overlap_peaks = data[np.in1d(data['name'], np.intersect1d(id_first_model, id_second_model))]
    peaks_with_not_overlap_sites = pd.DataFrame()
    peaks_with_overlap_sites = pd.DataFrame()

    for index, peak in perhaps_overlap_peaks.iterrows():
        subset_second_model = pd.DataFrame(second_model_sites[second_model_sites['name'] == peak['name']])
        subset_first_model = pd.DataFrame(first_model_sites[first_model_sites['name'] == peak['name']])
        intersections = [overlap(row_first_model, row_second_model) for index_first_model, row_first_model in subset_first_model.iterrows(
            ) for index_second_model, row_second_model in subset_second_model.iterrows()]

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
        count_first_model_sites = len(np.intersect1d(only_first_model_sites['name'], subset_name_peaks))
        count_second_model_sites = len(np.intersect1d(only_second_model_sites['name'], subset_name_peaks))
        count_overlap_sites = len(np.intersect1d(peaks_with_overlap_sites['name'], subset_name_peaks))
        count_not_overlap_sites = len(np.intersect1d(peaks_with_not_overlap_sites['name'], subset_name_peaks))
        count_no_sites = len(np.intersect1d(no_sites_peaks['name'], subset_name_peaks))

        count.append({'no_sites': count_no_sites,
         'first_model_sites': count_first_model_sites, 'second_model_sites': count_second_model_sites,
         'overlap_sites': count_overlap_sites,
         'not_overlap_sites': count_not_overlap_sites})

    count_ = pd.DataFrame(count)
    count = pd.DataFrame()
    count = count.append(count_.iloc[0])
    for i in range(1,len(count_)):
        count = count.append(count.iloc[i - 1] + count_.iloc[i], ignore_index=True)
    count['peaks'] = top
    count = count[['overlap_sites', 'first_model_sites', 'second_model_sites',
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


    overlap_second_model_sites = second_model_sites[np.in1d(second_model_sites['name'], peaks_with_overlap_sites['name'])]
    overlap_first_model_sites = first_model_sites[np.in1d(first_model_sites['name'], peaks_with_overlap_sites['name'])]
    overlap_sites = pd.concat([overlap_second_model_sites, overlap_first_model_sites])
    overlap_sites = overlap_sites.sort_values(by=['name'])


    only_second_model_sites = only_second_model_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']]
    only_second_model_sites.to_csv(out_dir + '/' + tag + '_all_second_model.sites', sep='\t', index=False, header=False)
    only_second_model_sites[only_second_model_sites['name'] <= 10000 ].to_csv(out_dir + '/' + tag + '_10000_second_model.sites',
                                                             sep='\t', index=False, header=False)

    only_first_model_sites = only_first_model_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']]
    only_first_model_sites.to_csv(out_dir + '/' + tag + '_all_first_model.sites', sep='\t', index=False, header=False)
    only_first_model_sites[only_first_model_sites['name'] <= 10000 ].to_csv(out_dir + '/' + tag + '_10000_first_model.sites',
                                                             sep='\t', index=False, header=False)

    overlap_sites = overlap_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']]
    overlap_sites['start'] = [int(i) for i in overlap_sites['start']]
    overlap_sites['end'] = [int(i) for i in overlap_sites['end']]
    overlap_sites.to_csv(out_dir + '/' + tag + '_all_overlap.sites', sep='\t', index=False, header=False)
    overlap_sites[overlap_sites['name'] <= 10000 ].to_csv(out_dir + '/' + tag + '_10000_overlap.sites',
                                                             sep='\t', index=False, header=False)

    only_second_model_peaks = only_second_model_peaks[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    only_second_model_peaks.to_csv(out_dir + '/' + tag + '_all_second_model.peaks', sep='\t', index=False, header=False)
    only_second_model_peaks[only_second_model_peaks['name'] <= 10000 ].to_csv(out_dir + '/' + tag + '_10000_second_model.peaks',
                                                             sep='\t', index=False, header=False)

    only_first_model_peaks = only_first_model_peaks[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    only_first_model_peaks.to_csv(out_dir + '/' + tag + '_all_first_model.peaks', sep='\t', index=False, header=False)
    only_first_model_peaks[only_first_model_peaks['name'] <= 10000 ].to_csv(out_dir + '/' + tag + '_10000_first_model.peaks',
                                                             sep='\t', index=False, header=False)

    peaks_with_overlap_sites = peaks_with_overlap_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand']]
    peaks_with_overlap_sites['start'] = [int(i) for i in peaks_with_overlap_sites['start']]
    peaks_with_overlap_sites['end'] = [int(i) for i in peaks_with_overlap_sites['end']]
    peaks_with_overlap_sites.to_csv(out_dir + '/' + tag + '_all_overlap.peaks', sep='\t', index=False, header=False)
    peaks_with_overlap_sites[peaks_with_overlap_sites['name'] <= 10000 ].to_csv(out_dir + '/' + tag + '_10000_overlap.peaks',
                                                             sep='\t', index=False, header=False)
if __name__ == '__main__':
    main()
