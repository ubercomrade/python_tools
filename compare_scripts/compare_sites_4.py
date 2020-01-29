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
from venn import generate_petal_labels, draw_venn, generate_colors
from numpy.random import choice


def read_bed(path):
    container = []
    with open(path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            container.append({
                'chr': row[0], 'start': int(row[1]), 'end': int(row[2])})
    container = sorted(container, key = itemgetter('chr', 'start'))
    return(container)


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

    
def is_intersect(interval, intervals):
    for i in intervals:
        if interval['start'] < i['end'] and interval['end'] > i['start']:
            return(1)
    return(0)
        
def get_indexes(peaks, sites):
    container = []
    append = container.append
    
    chrs = list(set([i['chr'] for i in peaks]))
    chrs.sort()
    
    index = 0
    for chr_ in chrs:
        sub_peaks = [i for i in peaks if i['chr'] == chr_]
        sub_sites = [i for i in sites if i['chr'] == chr_]
        if len(sub_sites) == 0:
            index += len(sub_peaks)
            continue
            
        for p in sub_peaks:
            if is_intersect(p, sub_sites):
                append(index)
            index += 1    
                
    return(container)


def write_table(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = list(data.keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerow(data)
    pass

        
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
    parser.add_argument('-fourth', '--fourth_model', action='store', dest='fourth_model_path',
                        required=True, help='path to file that contain sites obtained by_third_model in peaks file')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    parser.add_argument('-o', '--out', action='store', dest='out_dir',
                        required=True, help='OUT_DIR')
    parser.add_argument('-fname', '--first_name', action='store', dest='fname',
                        required=False, default='PWM', help='First data name for plot')
    parser.add_argument('-sname', '--second_name', action='store', dest='sname',
                        required=False, default='BAMM', help='Second data name for plot')
    parser.add_argument('-tname', '--third_name', action='store', dest='tname',
                        required=False, default='INMODE', help='Third data name for plot')
    parser.add_argument('-foname', '--fourth_name', action='store', dest='foname',
                    required=False, default='SITEGA', help='Fourth data name for plot')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def creat_petal(first_names, second_names, third_names, fourth_names):
    petal_labels = {'0001': str(len(fourth_names - (first_names | second_names | third_names))),
     '0010': str(len(third_names - (first_names | second_names | fourth_names))),
     '0011': str(len((fourth_names & third_names) - (first_names | second_names))),
     '0100': str(len(second_names - (fourth_names | first_names | third_names))),
     '0101': str(len((second_names & fourth_names) - (first_names | third_names))),
     '0110': str(len((second_names & third_names) - (first_names | fourth_names))),
     '0111': str(len((second_names & third_names & fourth_names) - first_names)),
     '1000': str(len(first_names - (fourth_names | second_names | third_names))),
     '1001': str(len((first_names & fourth_names) - (second_names | third_names))),
     '1010': str(len((first_names & third_names) - (second_names | fourth_names))),
     '1011': str(len((first_names & third_names & fourth_names - second_names))),
     '1100': str(len((first_names & second_names) - (third_names | fourth_names))),
     '1101': str(len((first_names & second_names & fourth_names) - third_names)),
     '1110': str(len((first_names & second_names & third_names) - fourth_names)),
     '1111': str(len(first_names & second_names & third_names & fourth_names ))}
    return(petal_labels)


def main():

    args = parse_args()
    first_model_path = args.first_model_path
    second_model_path = args.second_model_path
    third_model_path = args.third_model_path
    fourth_model_path = args.fourth_model_path
    peaks_path = args.input_peaks
    tag = args.tag
    out_dir = args.out_dir
    fname = args.fname
    sname = args.sname
    tname = args.tname
    foname = args.foname
    
    id_to_name = {'0001': 'SITEGA',
     '0010': 'INMODE',
     '0011': 'INMODExSITEGA',
     '0100': 'BAMM',
     '0101': 'BAMMxSITEGA',
     '0110': 'BAMMxINMODE',
     '0111': 'BAMMxINMODExSITEGA',
     '1000': 'PWM',
     '1001': 'PWMxSITEGA',
     '1010': 'PWMxINMODE',
     '1011': 'PWMxINMODExSITEGA',
     '1100': 'PWMxBAMM',
     '1101': 'PWMxBAMMxSITEGA',
     '1110': 'PWMxBAMMxINMODE',
     '1111': 'PWMxBAMMxINMODExSITEGA'}

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)


    peaks = read_bed(peaks_path)

    first_model_sites = read_bed(first_model_path)
    first_names = set(get_indexes(peaks, first_model_sites))

    second_model_sites = read_bed(second_model_path)
    second_names = set(get_indexes(peaks, second_model_sites))

    third_model_sites = read_bed(third_model_path) 
    third_names = set(get_indexes(peaks, third_model_sites))

    fourth_model_sites = read_bed(fourth_model_path)
    fourth_names = set(get_indexes(peaks, fourth_model_sites))
    
    petal_labels = creat_petal(first_names, second_names, third_names, fourth_names)

    
    ########################
    # WRITE RESULTS TO TSV #
    data = dict()
    for k in petal_labels.keys():
        data[id_to_name[k]] = petal_labels[k]
    data['PEAKS'] = len(peaks)
    write_table(out_dir + '/' + tag + '_COUNT.tsv', data)
    ########################
    
    #############
    # DRAW VENN #

    for k in petal_labels.keys():
        petal_labels[k] = '{:.2f}%'.format((int(petal_labels[k]) / len(peaks) * 100))
        
    ax = draw_venn(
        petal_labels=petal_labels, dataset_labels=['PWM', 'BAMM', 'INMODE', 'SITEGA'],
        hint_hidden=False, colors=generate_colors(n_colors=4),
        figsize=(8, 8), fontsize=14, legend_loc="best", ax=None
        )
    ax.figure.savefig(out_dir + '/' + tag + '_VENN.pdf', dpi=150, bbox_inches="tight")
    #############

if __name__ == '__main__':
    main()
