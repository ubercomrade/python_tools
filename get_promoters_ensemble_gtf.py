import pandas as pd
import numpy as np
import argparse
import sys



def read_gtf(path):

    gtf = pd.read_csv(path,
                      sep='\t',comment='#', header=None, dtype= {'chr': str},
                     names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    attribute = gtf['attribute'].str.split('\";')
    gtf = gtf.drop(columns=['attribute'])

    res = []
    for record in attribute:
        rec = []
        for i in record:
            if i == '':
                continue
            (a,b) = i.strip().split(' \"', maxsplit=1)
            rec.append((a,b))
        res.append(dict(rec))
    attribute = pd.DataFrame(res)
    res = None
    gtf = pd.concat([gtf, attribute], axis=1, sort=False)
    return(gtf)


def get_promoters(gtf, rigth=-2000, left=0):
    df = gtf[np.logical_and(gtf['feature'] == 'gene', gtf['gene_biotype'] == 'protein_coding')]
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


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', action='store', dest='gtf',
                        required=True, help='path to GTF file')
    parser.add_argument('-o', '--output', action='store', dest='bed',
                                  required=True, help='path to BED file')
    parser.add_argument('-l', '--left', action='store', type=int, dest='left',
                                  default=0, required=False, help='left_tail + TSS, default_value = 0')
    parser.add_argument('-r', '--right', action='store', type=int, dest='rigth',
                                  default=-2000, required=False, help='TSS + rigth_tail, default_value = -2000')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    path = args.gtf
    bed = args.bed
    left = args.left
    rigth = args.rigth

    gtf = read_gtf(path)
    promoters = get_promoters(gtf, left=left, rigth=rigth)
    promoters.to_csv(bed,
                     header=False, index=False, sep='\t')


if __name__ == '__main__':
    main()
