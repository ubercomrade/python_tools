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
import csv
from operator import itemgetter


def read_scan_file(path, scan_id):
    container = []
    with open(path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            container.append({
                'chromosome': row[0], 'start': int(row[1]), 'end': int(row[2]),
                'name': int(row[3].split('_')[1]), 'score': float(row[4]),
                'strand': row[5], 'site': row[6], 'type': scan_id
            })
    return(container)


def read_peaks(path):
    container = []
    with open(path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            container.append({
                'chromosome': row[0], 'start': int(row[1]),
                'end': int(row[2]), 'name': int(row[3].split('_')[1])
            })
    return(container)


def write_sites(path, data):
    with open(path, 'w') as file:
        for line in data:
            w = list(line.values())
            file.write('\t'.join(map(str, w[:-1])) + '\n')
    pass


def write_table(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = ['first_model_sites', 'second_model_sites', 'third_model_sites',
        'overlap_first_second_models', 'overlap_first_third_models', 'overlap_second_third_models',
        'overlap_all_models', 'overlap_first_second_and_first_third_models',
        'overlap_second_first_and_second_third_models', 
        'overlap_third_first_and_third_second_models', 'no_sites', 'not_overlap', 'peaks']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for line in data:
            writer.writerow(line)
    pass


def overlap(sites_1, sites_2):
    for site in sites_1:
        if get_intersect_index(site, sites_2):
            return(1)
        else:
            continue

    for site in sites_2:
        if get_intersect_index(site, sites_1):
            return(1)
        else:
            continue

    return(0)
    
    
def is_intersect(start1, end1, start2, end2):
    return(end1 > start2 and end2 > start1)


def get_intersect_index(interval, intervals):
    for index, i in enumerate(intervals):
        if is_intersect(interval['start'], interval['end'], i['start'], i['end']) and interval['type'] != i['type']:
            return(1)
        else:
            return(0)


def peak_classification(peak, first_model_sites, second_model_sites, third_model_sites):

    first_model = [i for i in first_model_sites if i['name'] == peak['name']]
    second_model = [i for i in second_model_sites if i['name'] == peak['name']]
    third_model = [i for i in third_model_sites if i['name'] == peak['name']]
    
    if len(first_model) == 0 and len(second_model) == 0 and len(third_model) == 0:
        return('no_sites')
    elif len(first_model) != 0 and len(second_model) == 0 and len(third_model) == 0:
        return('first_model')
    elif len(first_model) == 0 and len(second_model) != 0 and len(third_model) == 0:
        return('second_model')
    elif len(first_model) == 0 and len(second_model) == 0 and len(third_model) != 0:
        return('third_model')
    elif len(first_model) != 0 and len(second_model) != 0 and len(third_model) == 0:
        overlap_fs = overlap(first_model, second_model)
        if overlap_fs:
            return('overlap_first_second_models')
        else:
            return('not_overlap')
    elif len(first_model) != 0 and len(second_model) == 0 and len(third_model) != 0:
        overlap_ft = overlap(first_model, third_model)
        if overlap_ft:
            return('overlap_first_third_models')
        else:
            return('not_overlap')
    elif len(first_model) == 0 and len(second_model) != 0 and len(third_model) != 0:
        overlap_st = overlap(second_model, third_model)
        if overlap_st:
            return('overlap_second_third_models')
        else:
            return('not_overlap')
    elif len(first_model) != 0 and len(second_model) != 0 and len(third_model) != 0:
        overlap_fs = overlap(first_model, second_model)
        overlap_ft = overlap(first_model, third_model)
        overlap_st = overlap(second_model, third_model)
        if overlap_fs and overlap_ft and overlap_st:
            return('overlap_all_models')
        ###
        elif overlap_fs == 1 and overlap_ft == 1 and overlap_st == 0:
            return('overlap_first_second_and_first_third_models')
        elif overlap_fs == 1 and overlap_ft == 0 and overlap_st == 1:
            return('overlap_second_first_and_second_third_models')
        elif overlap_fs == 0 and overlap_ft == 1 and overlap_st == 1:
            return('overlap_third_first_and_third_second_models')
        ###    
        elif overlap_fs == 1 and overlap_ft == 0 and overlap_st == 0:
            return('overlap_first_second_models')
        elif overlap_fs == 0 and overlap_ft == 1 and overlap_st == 0:
            return('overlap_first_third_models')
        elif overlap_fs == 0 and overlap_ft == 0 and overlap_st == 1:
            return('overlap_second_third_models')

        ###
        else:
            return('not_overlap')


def creat_classification_of_peaks(data, first_model_sites, second_model_sites, third_model_sites):
    classification = []
    append = classification.append
    for peak in data:
        append(peak_classification(peak, first_model_sites, second_model_sites, third_model_sites))
    return(classification)


def creat_count_freq_tables(classification, data, step=100):
    count = []
    frequency = []
    top = [i * step for i in range(1, len(data) // step + 1)]

    for i in top:
        subset_classification = classification[:i]
        count_first_model_sites = sum(['first_model' == i for i in subset_classification])
        count_second_model_sites = sum(['second_model' == i for i in subset_classification])
        count_third_model_sites = sum(['third_model' == i for i in subset_classification])
        count_no_sites = sum(['no_sites' == i for i in subset_classification])
        overlap_first_second_models = sum(['overlap_first_second_models' == i for i in subset_classification])
        overlap_first_third_models = sum(['overlap_first_third_models' == i for i in subset_classification])
        overlap_second_third_models = sum(['overlap_second_third_models' == i for i in subset_classification])
        overlap_all_models = sum(['overlap_all_models' == i for i in subset_classification])
        not_overlap = sum(['not_overlap' == i for i in subset_classification])
        overlap_first_second_and_first_third_models = sum(['overlap_first_second_and_first_third_models' == i for i in subset_classification])
        overlap_second_first_and_second_third_models = sum(['overlap_second_first_and_second_third_models' == i for i in subset_classification])
        overlap_third_first_and_third_second_models = sum(['overlap_third_first_and_third_second_models' == i for i in subset_classification])

        count.append({'no_sites': count_no_sites,
                      'not_overlap': not_overlap,
                      'first_model_sites': count_first_model_sites,
                      'second_model_sites': count_second_model_sites,
                      'third_model_sites': count_third_model_sites,
                      'overlap_first_second_models': overlap_first_second_models,
                      'overlap_first_third_models': overlap_first_third_models,
                      'overlap_second_third_models': overlap_second_third_models,
                      'overlap_all_models': overlap_all_models,
                      'overlap_first_second_and_first_third_models': overlap_first_second_and_first_third_models,
                      'overlap_second_first_and_second_third_models': overlap_second_first_and_second_third_models,
                      'overlap_third_first_and_third_second_models': overlap_third_first_and_third_second_models,
                     'peaks': len(subset_classification)})

        frequency.append({'no_sites': count_no_sites / len(subset_classification),
                  'not_overlap': not_overlap / len(subset_classification),
                  'first_model_sites': count_first_model_sites / len(subset_classification),
                  'second_model_sites': count_second_model_sites / len(subset_classification),
                  'third_model_sites': count_third_model_sites / len(subset_classification),
                  'overlap_first_second_models': overlap_first_second_models / len(subset_classification),
                  'overlap_first_third_models': overlap_first_third_models / len(subset_classification),
                  'overlap_second_third_models': overlap_second_third_models / len(subset_classification),
                  'overlap_all_models': overlap_all_models / len(subset_classification),
                  'overlap_first_second_and_first_third_models': overlap_first_second_and_first_third_models / len(subset_classification),
                  'overlap_second_first_and_second_third_models': overlap_second_first_and_second_third_models / len(subset_classification),
                  'overlap_third_first_and_third_second_models': overlap_third_first_and_third_second_models / len(subset_classification),
                 'peaks': len(subset_classification)})

    return(count, frequency)



        
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
    parser.add_argument('-fname', '--first_name', action='store', dest='fname',
                        required=False, default='First', help='First data name for plot')
    parser.add_argument('-sname', '--second_name', action='store', dest='sname',
                        required=False, default='Second', help='Second data name for plot')
    parser.add_argument('-tname', '--third_name', action='store', dest='tname',
                        required=False, default='Third', help='Third data name for plot')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():

    args = parse_args()
    first_model_path = args.first_model_path
    second_model_path = args.second_model_path
    third_model_path = args.third_model_path
    peaks_path = args.input_peaks
    tag = args.tag
    out_dir = args.out_dir
    fname = args.fname
    sname = args.sname
    tname = args.tname

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    data = read_peaks(peaks_path)

    first_model_sites = read_scan_file(first_model_path, 'first_model')
    second_model_sites = read_scan_file(second_model_path, 'second_model')
    third_model_sites = read_scan_file(third_model_path, 'third_model')

    # Classification of peaks
    classification = creat_classification_of_peaks(data, first_model_sites, second_model_sites, third_model_sites)


    #############################################
    #Make table with count of diff kind of peaks#
    #############################################
    count, frequency = creat_count_freq_tables(classification, data, step=100)
    write_table(out_dir + '/' + tag + '_COUNT.tsv', count)
    write_table(out_dir + '/' + tag + '_FREQUENCY.tsv', frequency)


    #############
    #WRITE SITES#
    #############

    peak_ids = [index for index, i in enumerate(classification) if i == 'first_model']
    only_first_model_sites = [i for i in first_model_sites if i['name'] in peak_ids]
    peak_ids = [index for index, i in enumerate(classification) if i == 'second_model']
    only_second_model_sites = [i for i in second_model_sites if i['name'] in peak_ids]
    peak_ids = [index for index, i in enumerate(classification) if i == 'third_model']
    only_third_model_sites = [i for i in third_model_sites if i['name'] in peak_ids]


    write_sites(out_dir + '/' + tag + '_only_first_model.sites', only_first_model_sites)
    write_sites(out_dir + '/' + tag + '_only_second_model.sites', only_second_model_sites)
    write_sites(out_dir + '/' + tag + '_only_third_model.sites', only_third_model_sites)

if __name__ == '__main__':
    main()
