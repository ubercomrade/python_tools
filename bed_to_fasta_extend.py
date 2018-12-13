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


import linecache
import argparse
import sys
import os
import functools
from itertools import islice
import multiprocessing as mp
import pickle
import pandas as pd
import numpy as np


def read_bed(path):
    bed = pd.read_csv(path,
                      sep='\t', header=None,
                      usecols=[0, 1, 2, 3, 4, 5],
                      names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    return(bed)


def minimal_length(bed):
    '''
    Получить минимальную длину последовательности
    '''
    length = min(bed['end'] - bed['start'])
    return(length)


def maximal_length(bed):
    '''
    Получить максимальную длину последовательности
    '''
    length = max(bed['end'] - bed['start'])
    return(length)


def slice_to_size(bed, length):
    '''
    Привести все записи к единой длине
    '''
    start = bed['start']
    end = bed['end']
    mid = (end - start) // 2 + start
    new_start = mid - length // 2
    new_end = mid + (length - length // 2)

    bed['start'] = new_start
    bed['end'] = new_end
    return(bed)


def add_tail(bed, tail_length):
    '''
    Добавть хвосты ко всем записям
    '''
    start = bed['start']
    end = bed['end']
    new_start = start - tail_length
    new_end = end + tail_length
    bed['start'] = new_start
    bed['end'] = new_end
    return(bed)


def stat_of_fasta_file(path):
    '''
    start line == 1
    '''
    inf_data = path + '.pickle'
    flag = os.path.isfile(inf_data)
    if flag:
        with open(inf_data, 'rb') as file:
            output = pickle.load(file)
        file.close()
        return output
    else:
        stat = {}
        length = int()
        with open(path, 'r') as file:
            count = 0
            for line in file:
                line = line.strip().split()[0]
                if line[0] == '>':
                    stat[line[1:]] = count + 2
                else:
                    if count == 2:
                        length = len(line)
                count += 1
        file.close()
        output = (stat, length)
        with open(inf_data, 'wb') as file:
            pickle.dump(output, file)
        file.close()
        return output


def modify_bio_records(bed, to_min, to_max, to_size, tail):
    '''
    Модифицирование длины в зависимости от параметров
    '''
    if to_min:
        min_length = minimal_length(bed)
        bio_copy = slice_to_size(bed, min_length)
        return(bio_copy)
    elif to_max:
        max_length = maximal_length(bed)
        bio_copy = slice_to_size(bed, max_length)
        return(bio_copy)
    elif not (to_size is None):
        to_size = int(to_size)
        bio_copy = slice_to_size(bed, to_size)
        return(bio_copy)
    elif not (tail is None):
        tail = int(tail)
        bio_copy = add_tail(bed, tail)
        return(bio_copy)
    else:
        return(bed)


def complement(seq):
    '''
    Make reverse and compelent
    '''
    output = str()
    for letter in seq:
        if letter == 'A':
            output += 'T'
        elif letter == 'C':
            output += 'G'
        elif letter == 'G':
            output += 'C'
        elif letter == 'T':
            output += 'A'
    output = output[::-1]
    return(output)


def bed_to_fasta(path_fasta, path_bed, to_min, to_max, to_size, tail):
    bed_peaks = read_bed(path_bed)
    bed_peaks = modify_bio_records(bed_peaks, to_min, to_max, to_size, tail)
    chromosomes, length = stat_of_fasta_file(path_fasta)
    results = []
    for i in range(len(bed_peaks)):
        rec = {'name': str(), 'chr': str(), 'start': int(), 'end': int(), 'strand': str()}
        seq = str()  # container for sequnce of peak
        chr_name = bed_peaks.iloc[i]['chromosome']
        try:
            chr_start_line = chromosomes[chr_name]
        except:
            print('{0} was not found in FASTA'.format(chr_name))
            continue

        rec['chr'] = bed_peaks.iloc[i]['chromosome']
        if bed_peaks.iloc[i]['name'] != '.':
            rec['name'] = bed_peaks.iloc[i]['name']
        else:
            rec['name'] = 'peaks_' + str(i)

        if not bed_peaks.iloc[i]['start'] < 0:
            rec['start'] = bed_peaks.iloc[i]['start']
        else:
            continue

        if not bed_peaks.iloc[i]['end'] < 0:
            rec['end'] = bed_peaks.iloc[i]['end']
        else:
            continue

        if bed_peaks.iloc[i]['strand'] == '.' and isinstance(bed_peaks.iloc[i]['strand'], str):
            rec['strand'] = '+'
        elif not isinstance(bed_peaks.iloc[i]['strand'], str):
            if np.isnan(bed_peaks.iloc[i]['strand']):
                rec['strand'] = '+'
        else:
            rec['strand'] = bed_peaks.iloc[i]['strand']

        peak_start_line = bed_peaks.iloc[i]['start'] // length + \
            chr_start_line  # number of line where seq started
        peak_end_line = bed_peaks.iloc[i]['end'] // length + \
            chr_start_line  # number of line where seq ended
        # position number of nucleotide in first line
        position_start = bed_peaks.iloc[i]['start'] % length
        position_end = bed_peaks.iloc[i]['end'] - bed_peaks.iloc[i]['start'] + position_start
        for line_number in range(peak_start_line, peak_end_line + 1):
            seq += linecache.getline(path_fasta, line_number).strip()
        if not 'N' in seq[position_start:position_end]:
            if rec['strand'] == '+':
                rec['seq'] = seq[position_start:position_end]
            else:
                rec['seq'] = complement(seq[position_start:position_end])
        else:
            continue
        results.append(rec)
    linecache.clearcache()
    return(results)


def stat_of_wig_file(wig_path):
    stat = []
    with open(wig_path, 'r') as file:
        count = 0
        for line in file:
            if line.startswith('#'):
                line = line.strip().split()

                rec = {'line': count,
                       'chr': line[2].split(':')[0],
                       'start': int(line[2].split(':')[1].split('-')[0]),
                       'end': int(line[2].split(':')[1].split('-')[1])}
                stat.append(rec)
                count += 1
            else:
                count += 1
    file.close()
    return pd.DataFrame(stat)


def bed_to_wig_fasta(record, stat, wig_path):
    chr_ = record['chr']
    start = record['start']
    end = record['end']
    record['wig'] = str()

    line_number_start = int(stat[np.logical_and(np.logical_and(
        stat['chr'] == chr_, stat['start'] <= start), stat['end'] >= start)]['line'])
    #line_number_end = int(stat[np.logical_and(np.logical_and(stat['chr'] == chr_, stat['start'] <= end), stat['end'] >= end)]['line'])
    index_end = (stat[np.logical_and(np.logical_and(stat['chr'] == chr_,
                                                    stat['start'] <= end), stat['end'] >= end)].index + 1)[0]
    line_number_end = stat.loc[index_end]['line']

    if line_number_start == line_number_end:
        index_end = (stat[np.logical_and(np.logical_and(stat['chr'] == chr_,
                                                        stat['start'] <= end), stat['end'] >= end)].index + 1)[0]
        line_number_end = stat.loc[index_end]['line']

    section = {}
    with open(wig_path, 'r') as f:
        for index, wig_line in enumerate(f):
            if wig_line.startswith('#'):
                continue
            elif index >= line_number_start and index < line_number_end:
                wig_start = int(wig_line.strip().split()[1])
                wig_end = int(wig_line.strip().split()[2])
                wig_score = float(wig_line.strip().split()[3])
                for i in range(wig_start, wig_end):
                    section[i] = wig_score
            elif index == line_number_end:
                break
    f.close()
    # print(section)
    for i in range(start, end):
        record['wig'] += str(section[i]) + ' '
        record['wig'] = record['wig'][-1]

    return(record)


def bed_to_wig_fasta_2(record, stat, wig_path):
    chr_ = record['chr']
    start = record['start']
    end = record['end']
    record['wig'] = []

    line_number_start = int(stat[np.logical_and(np.logical_and(
        stat['chr'] == chr_, stat['start'] <= start), stat['end'] > start)]['line'])

    #line_number_end = int(stat[np.logical_and(np.logical_and(stat['chr'] == chr_, stat['start'] <= end), stat['end'] >= end)]['line'])
    index_end = (stat[np.logical_and(np.logical_and(stat['chr'] == chr_,
                                                    stat['start'] <= end), stat['end'] > end)].index + 1)[0]
    line_number_end = stat.loc[index_end]['line']

    if line_number_start == line_number_end:
        index_end = (stat[np.logical_and(np.logical_and(stat['chr'] == chr_,
                                                        stat['start'] <= end), stat['end'] > end)].index + 1)[0]
        line_number_end = stat.loc[index_end]['line']

    section = {}
    with open(wig_path, 'r') as f:
        res = []
        try:
            line = islice(f, int(line_number_start), int(line_number_end))
            for i in line:
                res.append(i)
        except StopIteration:
            print('Not lines in file')

        res = [i for i in res if not i.startswith('#')]
        for wig_line in res:
            wig_start = int(wig_line.strip().split()[1])
            wig_end = int(wig_line.strip().split()[2])
            wig_score = float(wig_line.strip().split()[3])
            for i in range(wig_start, wig_end):
                section[i] = wig_score
    f.close()
    # print(section)
    for i in range(start, end):
        record['wig'].append(str(section[i]))
    #record['wig'] = [str(i) for i in record['wig']]
    record['wig'] = ' '.join(record['wig'])
    # print(record['wig'])
    return(record)


def write_fasta(results, path_to_write):
    with open(path_to_write, 'w') as file:
        for record in results:
            file.write('>' + record['name'] + '|' + record['chr'] + '|' +
                       str(record['start']) + '-' + str(record['end']) + '|' +
                       record['strand'] + '\n')
            file.write(record['seq'] + '\n')
    file.close()
    pass


def write_wig_fata(results, path_to_write):
    with open(path_to_write, 'w') as file:
        for record in results:
            file.write('>' + record['wig'] + '\n')
            file.write(record['seq'] + '\n')
    file.close()
    pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--inputFasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file')
    parser.add_argument('-bed', '--inputBED', action='store', dest='input_bed',
                        required=True, help='path to BED file')
    parser.add_argument('-w', '--wig', action='store', dest='input_wig',
                        required=False, help='path to WIG file to make WIG_FASTA file')
    parser.add_argument('-of', '--outputFasta', action='store', dest='output_fasta',
                        required=True, help='path to output Fasta')
    parser.add_argument('-m', '--minimize', action='store_true', dest='to_min',
                        default=False, required=False,
                        help='With this parametr all peaks will change their length to the length of the shortest sequence in set')
    parser.add_argument('-M', '--maximaze', action='store_true', dest='to_max',
                        default=False, required=False,
                        help='With this parametr all peaks will change their length to the length of the longest sequence in set')
    parser.add_argument('-l', '--length', action='store', dest='length',
                        default=None, required=False,
                        help='With this parametr all peaks will change their length to the length equal to the given parameter')
    parser.add_argument('-t', '--tail', action='store', dest='tail',
                        default=None, required=False,
                        help='With this parametr all peaks will get tail in start and end of sequence')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    input_fasta = args.input_fasta
    input_bed = args.input_bed
    output_fasta = args.output_fasta
    to_min = args.to_min
    to_max = args.to_max
    to_size = args.length
    tail = args.tail
    wig_path = args.input_wig

    if os.path.isfile(wig_path):
        results = bed_to_fasta(input_fasta, input_bed, to_min, to_max, to_size, tail)
        stat = stat_of_wig_file(wig_path)
        with mp.Pool(40) as p:
            results = p.map(functools.partial(bed_to_wig_fasta_2,
                                              stat=stat, wig_path=wig_path), results)
        write_wig_fata(results, output_fasta)
    else:
        results = bed_to_fasta(input_fasta, input_bed, to_min, to_max, to_size, tail)
        write_fasta(results, output_fasta)


if __name__ == '__main__':
    main()
