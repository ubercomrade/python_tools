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
from matplotlib_venn import venn3


def read_bed_like_file(path):
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site'])
    return(df)


def read_peaks(path):

    try:
        bed = pd.read_csv(path,
                          sep='\t', header=None,
                          usecols=[0, 1, 2, 3, 4, 5],
                          dtype= {'chromosome': str, 'start': int, 'end': int,
                                 'score': float, 'strand': str},
                          names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    except:
        bed = pd.read_csv(path,
                          sep='\t', header=None,
                          usecols=[0, 1, 2, 3, 4],
                          dtype= {'chromosome': str, 'start': int, 'end': int,
                                 'score': float, 'strand': str},
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


def compare(vec1, vec2):
    return(np.array([i >= j for i in list(vec1) for j in list(vec2)]))

def overlap(sites_1, sites_2):
    '''
    Does the range (start1, end1) overlap with (start2, end2)?
    Based on De Morgan's laws
    '''
    if len(sites_1['start']) > 0 and len(sites_2['start']) > 0:
        #print(len(sites_1['start']), len(sites_1['end']), '1')
        #print(len(sites_2['start']), len(sites_2['end']), '2')
        return np.logical_and(compare(sites_1['end'], sites_2['start']),
                              compare(sites_2['end'], sites_1['start']))
    else:
        return(np.array(['False']))


def get_coord(site):
    return(np.array(site['start'], site['end']))



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--peaks', action='store', dest='input_peaks',
                        required=True, help='path to peaks in BED format or like BED format that contain all scaned peaks coordinates')
    parser.add_argument('-first', '--first_model', action='store', dest='first_model_path',
                        required=True, help='path to file that contain sites obtained by first_model in peaks file')
    parser.add_argument('-second', '--second_model', action='store', dest='second_model_path',
                        required=True, help='path to file that contain sites obtained by_second_model in peaks file')
    parser.add_argument('-third', '--third_model', action='store', dest='third_model_path',
                        required=True, help='path to file that contain sites obtained by_third_model in peaks file')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    parser.add_argument('-o', '--out', action='store', dest='out_dir',
                        required=True, help='OUT_DIR')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def peak_classification(peak, first_model_sites, second_model_sites, third_model_sites):

    first_model = first_model_sites[first_model_sites['name'] == peak['name']]
    second_model = second_model_sites[second_model_sites['name'] == peak['name']]
    third_model = third_model_sites[third_model_sites['name'] == peak['name']]

    if len(first_model) == 0 and len(second_model) == 0 and len(third_model) == 0:
        return('no_sites')
    elif len(first_model) != 0 and len(second_model) == 0 and len(third_model) == 0:
        return('first_model')
    elif len(first_model) == 0 and len(second_model) != 0 and len(third_model) == 0:
        return('second_model')
    elif len(first_model) == 0 and len(second_model) == 0 and len(third_model) != 0:
        return('third_model')
    elif len(first_model) != 0 and len(second_model) != 0 and len(third_model) == 0:
        if sum(overlap(first_model, second_model)) > 0:
            return('overlap_first_second_models')
        else:
            #return('not_overlap_first_second_models')
            #return('no_sites')
            return('not_overlap')
    elif len(first_model) != 0 and len(second_model) == 0 and len(third_model) != 0:
        if sum(overlap(first_model, third_model)) > 0:
            return('overlap_first_third_models')
        else:
            #return('not_overlap_first_third_models')
            #return('no_sites')
            return('not_overlap')
    elif len(first_model) == 0 and len(second_model) != 0 and len(third_model) != 0:
        if sum(overlap(second_model, third_model)) > 0:
            return('overlap_second_third_models')
        else:
            #return('not_overlap_second_third_models')
            #return('no_sites')
            return('not_overlap')
    elif len(first_model) != 0 and len(second_model) != 0 and len(third_model) != 0:
        if sum([i and j for i in overlap(first_model, second_model) for j in overlap(first_model, third_model)]) > 0:
            return('overlap_all_models')
        elif sum(overlap(first_model, second_model)) > 0 and sum(overlap(first_model, third_model)) == 0 and sum(overlap(second_model, third_model)) == 0:
            return('overlap_first_second_models')
        elif sum(overlap(first_model, second_model)) == 0 and sum(overlap(first_model, third_model)) > 0 and sum(overlap(second_model, third_model)) == 0:
            return('overlap_first_third_models')
        elif sum(overlap(first_model, second_model)) == 0 and sum(overlap(first_model, third_model)) == 0 and sum(overlap(second_model, third_model)) > 0:
            return('overlap_second_third_models')
        else:
            #return('no_sites')
            return('not_overlap')


def main():

    args = parse_args()
    first_model_path = args.first_model_path
    second_model_path = args.second_model_path
    third_model_path = args.third_model_path
    peaks_path = args.input_peaks
    tag = args.tag
    out_dir = args.out_dir

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    data = read_peaks(peaks_path)
    names = [int(i.split('_')[1]) for i in data['name']]
    data['name'] = names

    first_model_sites = read_bed_like_file(first_model_path)
    first_model_sites['type'] = 'first_model'
    names = [int(i.split('_')[1]) for i in first_model_sites['name']]
    first_model_sites['name'] = names

    second_model_sites = read_bed_like_file(second_model_path)
    second_model_sites['type'] = 'second_model'
    names = [int(i.split('_')[1]) for i in second_model_sites['name']]
    second_model_sites['name'] = names

    third_model_sites = read_bed_like_file(third_model_path)
    third_model_sites['type'] = 'third_model'
    names = [int(i.split('_')[1]) for i in third_model_sites['name']]
    third_model_sites['name'] = names

    classification = []
    for index, peak in data.iterrows():
        classification.append(peak_classification(peak, first_model_sites, second_model_sites, third_model_sites))


    #############################################
    #Make table with count of diff kind of peaks#
    #############################################
    count = []
    top = [i * 1000 for i in range(1, len(data) // 1000 + 1)]

    for i in range(len(top)):
        subset_classification = classification[i*1000:(i+1)*1000]
        count_first_model_sites = sum(['first_model' == i for i in subset_classification])
        count_second_model_sites = sum(['second_model' == i for i in subset_classification])
        count_third_model_sites = sum(['third_model' == i for i in subset_classification])
        count_no_sites = sum(['no_sites' == i for i in subset_classification])
        overlap_first_second_models = sum(['overlap_first_second_models' == i for i in subset_classification])
        overlap_first_third_models = sum(['overlap_first_third_models' == i for i in subset_classification])
        overlap_second_third_models = sum(['overlap_second_third_models' == i for i in subset_classification])
        overlap_all_models = sum(['overlap_all_models' == i for i in subset_classification])
        not_overlap = sum(['not_overlap' == i for i in subset_classification])


        count.append({'no_sites': count_no_sites,
                      'not_overlap': not_overlap,
                      'first_model_sites': count_first_model_sites,
                      'second_model_sites': count_second_model_sites,
                      'third_model_sites': count_third_model_sites,
                      'overlap_first_second_models': overlap_first_second_models,
                      'overlap_first_third_models': overlap_first_third_models,
                      'overlap_second_third_models': overlap_second_third_models,
                      'overlap_all_models': overlap_all_models})

    count_ = pd.DataFrame(count)
    count = pd.DataFrame()
    count = count.append(count_.iloc[0])
    for i in range(1,len(count_)):
        count = count.append(count.iloc[i - 1] + count_.iloc[i], ignore_index=True)
    count['peaks'] = top
    count = count[['first_model_sites', 'second_model_sites', 'overlap_first_second_models', 'third_model_sites',
                   'overlap_first_third_models', 'overlap_second_third_models',
                   'overlap_all_models', 'no_sites', 'not_overlap', 'peaks']]
    count.to_csv(out_dir + '/' + tag + '_COUNT.tsv', sep='\t', index=False)

    frequency = pd.DataFrame(count).copy()
    for column in frequency:
        if column == 'peaks':
            continue
        frequency[column] = frequency[column] / frequency['peaks']
    frequency.to_csv(out_dir + '/' + tag + '_FREQUENCY.tsv', sep='\t', index=False)


    venn3(subsets=np.around(np.array(frequency.iloc[4,:7]), 3), set_labels = ('PWM', 'BAMM', 'INMODE'))
    plt.savefig(out_dir + '/' + tag + '_PIC.png', dpi=150)

    ##################################
    only_first_model_sites = first_model_sites.loc[first_model_sites['name'].searchsorted(np.array([index for index, i in enumerate(classification) if i == 'first_model']))]
    only_second_model_sites = second_model_sites.loc[second_model_sites['name'].searchsorted(np.array([index for index, i in enumerate(classification) if i == 'second_model']))]
    only_third_model_sites = third_model_sites.loc[third_model_sites['name'].searchsorted(np.array([index for index, i in enumerate(classification) if i == 'third_model']))]


    only_first_model_sites = only_first_model_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']]
    only_first_model_sites.to_csv(out_dir + '/' + tag + '_all_first_model.sites', sep='\t', index=False, header=False)

    only_second_model_sites = only_second_model_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']]
    only_second_model_sites.to_csv(out_dir + '/' + tag + '_all_second_model.sites', sep='\t', index=False, header=False)

    only_third_model_sites = only_third_model_sites[['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']]
    only_third_model_sites.to_csv(out_dir + '/' + tag + '_all_third_model.sites', sep='\t', index=False, header=False)

if __name__ == '__main__':
    main()
