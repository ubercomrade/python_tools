import sys
import argparse
import re
from operator import itemgetter

def read_bed(path):
    container = {}
    min_score = 1000
    with open(path) as file:
        for line in file:
            chrom, start, end, peak, score, strand, site = line.strip().split()
            if float(score) < min_score:
                min_score = float(score)
            if not peak in container:
                container[peak] = [{'chr': chrom,
                                    'start': int(start),
                                    'end': int(end),
                                    'score': float(score),
                                    'strand': strand,
                                    'site': site.lower()}]
            else:
                container[peak].append({'chr': chrom,
                                    'start': int(start),
                                    'end': int(end),
                                    'score': float(score),
                                    'strand': strand,
                                    'site': site.lower()})
    file.close()
    return(container, min_score)


def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
        for key in pwm:
            pwm[key] = tuple(pwm[key])
    file.close()
    return(pwm)


def to_score(norm_value, pwm):
    '''
    norm = (score - min) / (max - min) -> score = norm * (max - min) + min
    '''
    max_s = max_score(pwm)
    min_s = min_score(pwm)
    score = norm_value * (max_s - min_s) + min_s
    return(score)


def to_norm(score, pwm):
    max_s = max_score(pwm)
    min_s = min_score(pwm)
    norm_value = (score - min_s) / (max_s - min_s)
    return(norm_value)


def min_score(pwm):
    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += min(tmp)
    return(value)


def max_score(pwm):
    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += max(tmp)
    return(value)


def read_fasta(path):
    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                record = dict()
                record['name'] = line[0]
                record['chr'] = line[2]
                coordinates_strand = line[3]

                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                record['start'] = int(start)
                record['end'] = int(end)

                strand = re.findall(r'\(.\)', coordinates_strand[:-3])
                if not strand == []:
                    record['strand'] = strand[0].strip('()')
                else:
                    record['strand'] = '+'
            else:
                record['seq'] = line.strip().upper()
                fasta.append(record)
    file.close()
    return(fasta)


# def write_mcot_format(path, bed, fasta, pwm, threshold):
#     with open(path, 'w') as file:
#         for index, line in enumerate(fasta):
#             file.write('>peaks_{0}::{1}:{2}-{3}({4})\tSEQ {5}\tTHR {6}\n'.format(index, 
#                                                                                  line['chr'],
#                                                                                  line['start'],
#                                                                                  line['end'],
#                                                                                  line['strand'],
#                                                                                  index + 1,
#                                                                                  to_norm(threshold, pwm)))
#             if line['name'] in bed:
#                 bed[line['name']] = sorted(bed[line['name']], key=itemgetter('start'))
#                 for site in bed[line['name']]:
#                     pos = site['start'] - line['start']
#                     file.write('{0}\t{1}\t{2}\t{3}\n'.format(pos, to_norm(site['score'], pwm), site['strand'], site['site']))
#     return(0)


def write_mcot_format(path, bed, fasta, threshold):
    with open(path, 'w') as file:
        for index, line in enumerate(fasta):
            file.write('>peaks_{0}::{1}:{2}-{3}({4})\tSEQ {5}\tTHR {6}\n'.format(index, 
                                                                                 line['chr'],
                                                                                 line['start'],
                                                                                 line['end'],
                                                                                 line['strand'],
                                                                                 index + 1,
                                                                                 threshold))
            if line['name'] in bed:
                bed[line['name']] = sorted(bed[line['name']], key=itemgetter('start'))
                for site in bed[line['name']]:
                    pos = site['start'] - line['start']
                    file.write('{0}\t{1}\t{2}\t{3}\n'.format(pos, site['score'], site['strand'], site['site']))
    return(0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store',
                        help='path to bed with sites')
    parser.add_argument('fasta', action='store',
                        help='path to fasta with peaks')
    # parser.add_argument('pwm', action='store',
    #                     help='path to pwm')
    parser.add_argument('output', action='store',
                        help='path to write results')
    parser.add_argument('-t', '--threshold', action='store', type=float, dest='threshold',
                        required=False, help='If not defined, then the minimum value is used from scores')    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    bed = args.bed
    fasta = args.fasta
    output = args.output
    threshold = args.threshold
    fasta = read_fasta(fasta)
    bed, min_score = read_bed(bed)
    if threshold is None:
        threshold = min_score
    write_mcot_format(output, bed, fasta, threshold) 
    return(0)


if __name__=='__main__':
    main()
