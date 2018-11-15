import numpy as np
import pandas as pd
import csv
import argparse
import sys


def read_gtf(path):
    col_names = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase']
    data = dict()
    with open(path, 'r') as tsv:
        reader = csv.reader(tsv, delimiter='\t')
        counter = 0
        for line in reader:
            if line[0][0] == '#':
                continue
            for col, inf in zip(col_names, line[:-1]):
                if not col in data.keys():
                    data[col] = []
                    data[col].append(inf)
                else:
                    data[col].append(inf)
            for inf in line[-1].split(';')[:-1]:
                inf = inf.split()
                if inf[0] == 'tag':
                    continue
                if not inf[0] in data.keys():
                    data[inf[0]] = [None] * counter
                    data[inf[0]].append(inf[1].strip('/"'))
                else:
                    data[inf[0]].append(inf[1].strip('/"'))
            for key in data.keys():
                data[key] += (len(data['chr']) - len(data[key])) * [None]
            counter += 1
    tsv.close()
    gtf = pd.DataFrame(data)
    gtf['start'] = pd.to_numeric(gtf['start'])
    gtf['end'] = pd.to_numeric(gtf['end'])
    del data
    return(gtf)


def get_protein_coding_genes(gtf):
    return(gtf[np.logical_and(gtf['feature'] == 'gene', gtf['gene_type'] == 'protein_coding')])


def get_promoters(gtf, rigth=-2000, left=0):
    df = gtf[np.logical_and(gtf['feature'] == 'gene', gtf['gene_type'] == 'protein_coding')]
    promoters = {'chr': [], 'start': [], 'end': [],
                 'name': [], 'score': [], 'strand': [],
                 'signalValue': [], 'pValue': [], 'qValue': [],
                 'peak': []}

    for line, strand in enumerate(df['strand']):

        if strand == '+':
            promoters['start'].append(df['start'].iloc[line] + rigth)
            promoters['end'].append(df['start'].iloc[line] + left)

        if strand == '-':
            promoters['start'].append(df['end'].iloc[line] - left)
            promoters['end'].append(df['end'].iloc[line] - rigth)

        promoters['chr'].append(df['chr'].iloc[line])
        promoters['name'].append(df['gene_id'].iloc[line])
        promoters['score'].append(df['score'].iloc[line])
        promoters['strand'].append(df['strand'].iloc[line])
        promoters['signalValue'].append(0)
        promoters['pValue'].append(-1)
        promoters['qValue'].append(-1)
        promoters['peak'].append(-1)
    return(pd.DataFrame(promoters))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    parser_promoters = subparsers.add_parser('promoters', help='Get promoters from GTF')
    parser_promoters.add_argument('-i', '--input', action='store', dest='input_gtf',
                                  required=True, help='path to GTF file')
    parser_promoters.add_argument('-bed', '--output', action='store', dest='output_bed',
                                  required=True, help='path to BED file')
    parser_promoters.add_argument('-l', '--left', action='store', type=int, dest='left',
                                  default=-2000, required=False, help='left_tail + TSS, default_value = 0')
    parser_promoters.add_argument('-r', '--right', action='store', type=int, dest='rigth',
                                  default=-2000, required=False, help='TSS + rigth_tail, default_value = -2000')

    # if len(sys.argv) == 1:
    #    parser.print_help(sys.stderr)
    #    sys.exit(1)

    args = parser.parse_args()
    path = args.input_gtf
    bed = args.output_bed
    left = args.left
    rigth = args.rigth

    gtf = read_gtf(path)
    promoters = get_promoters(gtf, left=left, rigth=rigth)
    promoters.to_csv(bed,
                     header=False, index=False, sep='\t')
