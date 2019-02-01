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
import array
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
                      # usecols=[0, 1, 2, 3, 4, 5],
                      names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    #bed.loc[np.isnan(bed['strand']), 'strand'] = '.'
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
                rec['seq'] = seq[position_start:position_end].upper()
            else:
                rec['seq'] = complement(seq[position_start:position_end].upper())
        else:
            continue
        results.append(rec)
    linecache.clearcache()
    return(results)


def stat_of_wig_file(wig_path):
    stat = []
    with open(wig_path, 'r') as file:
        for line in file:
            if line.startswith('track'):
                continue
            elif line.startswith('fixedStep'):
                flag = 'fixedStep'
                break
            elif line.startswith('variableStep'):
                flag = 'variableStep'
                if 'span' in line:
                    span = True
                else:
                    span = False
                break
    file.close()

    if flag == 'fixedStep':
        with open(wig_path, 'r') as file:
            count = -1
            for line in file:
                count += 1
                if line.startswith('fixedStep'):
                    line = line.strip().split()
                    rec = {'line_start': count + 1,
                           'line_end': None,
                           'chr': line[1].split('=')[1],
                           'start': int(line[2].split('=')[1]),
                           'end': None,
                           'step': int(line[3].split('=')[1])}
                    stat.append(rec)
                if len(stat) >= 1:
                    stat[-1]['end'] = int((count - stat[-1]['line_start']) *
                                          stat[-1]['step'] + stat[-1]['start'])
                    stat[-1]['line_end'] = count
        file.close()

    elif flag == 'variableStep':
        with open(wig_path, 'r') as file:
            count = 0
            for line in file:
                count += 1
                if line.startswith('track'):
                    count -= 1
                    continue
                if line.startswith('variableStep'):
                    count -= 1
                    start = True
                    line = line.strip().split()
                    rec = {'line_start': count,
                           'line_end': None,
                           'chr': line[1].split('=')[1]}
                    if span:
                        rec['span'] = int(line[2].split('=')[1])
                    else:
                        rec['span'] = 1
                    stat.append(rec)
                if len(stat) >= 1:
                    stat[-1]['line_end'] = count

        file.close()
    return (pd.DataFrame(stat), flag)


def bed_to_wig_fasta_variable_step(record, wig, span):
    chr_ = record['chr']
    record_start = record['start']
    record_end = record['end']
    record['wig'] = []

    index_start = int(np.searchsorted(wig['position'], record_start)) - 3
    index_end = int(np.searchsorted(wig['position'], record_end)) + 3

    section = dict()
    for line in wig.iloc[index_start:index_end,].iterrows():
        #print(line)
        for i in range(span):
            section[int(line[1]['position']) + i] = int(line[1]['score'])
    #print(section)

    if not record_start in section: # Extend left-tail
        for i in range(record_start, min(section.keys())):
            section[i] = 1
        #print('left-tail was extended for {0} start-{1} end-{2}'.format(chr_, record_start, record_end))

    if not record_start in section: # Extend right-tail
        for i in range(max(section.keys()), record_end):
            section[i] = 1
        #print('right-tail was extended for {0} start-{1} end-{2}'.format(chr_, record_start, record_end))

        #print(section)
    for i in range(record_start, record_end):
        if not i in section:
            score = '1'
        else:
            score = str(section[i])
        record['wig'].append(score)
    record['wig'] = ' '.join(record['wig'])

    return(record)

def bed_to_wig_fasta_fixed_step(record, stat, wig_path):

    chr_ = record['chr']
    record_start = record['start']
    record_end = record['end']
    record['wig'] = []

    index_start = stat[np.logical_and(np.logical_and(stat['chr'] == chr_,
                                                         stat['start'] <= record_start),
                                          stat['end'] > record_start)].index

    line_number_start = int(stat.iloc[index_start]['line_start'])
    line_number_end = int(stat.iloc[index_start]['line_end'])
    step = int(stat.iloc[index_start]['step'])
    start = int(stat.iloc[index_start]['start'])
    end = int(stat.iloc[index_start]['end'])

    with open(wig_path, 'r') as f:
        res = []
        try:
            line = islice(f, int(line_number_start), int(line_number_end))
            for i in line:
                values = i.strip().split()
                values = [int(i) for i in values]
                res.append(array.array('i', values))
        except StopIteration:
                print('Not lines in file')
    f.close()

    res = [i for j in res for i in j] # unlist
    scores = np.array([i for i in res for j in range(step)])
    positions = np.array([i for i in range(start, end + 1)])
    for i in range(record_start, record_end):
        record['wig'].append(str(scores[np.searchsorted(positions, i)]))
    record['wig'] = ' '.join(record['wig'])

    return(record)


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
                        required=True, help='path to WIG file to make WIG_FASTA file')
    parser.add_argument('-of', '--outputFasta', action='store', dest='output_fasta',
                        required=True, help='path to output Fasta')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    input_fasta = args.input_fasta
    input_bed = args.input_bed
    output_fasta = args.output_fasta
    to_min = None
    to_max = None
    to_size = None
    tail = None
    wig_path = args.input_wig

    results = bed_to_fasta(input_fasta, input_bed, to_min, to_max, to_size, tail)
    stat, wig_flag = stat_of_wig_file(wig_path)
    if wig_flag == 'variableStep':
        chrs = {i['chr'] for i in results}
        wig = pd.read_csv(wig_path, sep='\t', names=['position', 'score'],
                         skiprows=1, header=None, dtype=int, comment='v')
        res = []
        for chr_ in chrs:
            sub_results = [record for record in results if record['chr'] == chr_]
            sub_stat = stat[stat['chr'] == chr_]
            sub_wig = wig.iloc[int(sub_stat['line_start']):int(sub_stat['line_end']),]
            span = int(sub_stat['span'])
            res += [bed_to_wig_fasta_variable_step(record, wig=sub_wig, span=span) for record in sub_results]

    elif wig_flag == 'fixedStep':
        with mp.Pool(mp.cpu_count()) as p:
            results = p.map(functools.partial(bed_to_wig_fasta_fixed_step, stat=stat,
                                          wig_path=wig_path), results)

    write_wig_fata(res, output_fasta)


if __name__ == '__main__':
    main()
