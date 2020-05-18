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
