import pandas as pd
import numpy as np
import argparse
import sys



def read_gtf(path):

    gtf = pd.read_csv(path,
                      sep='\t',comment='#', header=None, dtype= {'chr': str},
                     names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    attribute = gtf['attribute'].str.split('; ')
    gtf = gtf.drop(columns=['attribute'])

    res = []
    for record in attribute:
        rec = []
        for i in record:
            if i == '':
                continue
            (a,b) = i.strip().split(maxsplit=1)
            rec.append((a.strip(),b.strip(';\" ')))
        res.append(dict(rec))
    attribute = pd.DataFrame(res)
    res = None
    gtf = pd.concat([gtf, attribute], axis=1, sort=False)
    return(gtf)


def read_peaks(path):
    df = pd.read_csv(path,
                     sep='\t', header=None,
                     usecols=[0, 1, 2, 3, 4, 5], dtype= {'chr': str},
                     names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    return(df)


def get_promoters(gtf, rigth=-5000, left=5000):
    df = gtf[np.logical_and(gtf['feature'] == 'gene', gtf['gene_biotype'] == 'protein_coding')]
    promoters = {'chr': [], 'start': [], 'end': [],
                 'name': [], 'score': [], 'strand': [],
                 'signalValue': [], 'pValue': [], 'qValue': [],
                 'peak': [], 'gene_id': []}

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
        promoters['gene_id'].append(df['gene_id'].iloc[line])
    return(pd.DataFrame(promoters))


def overlap(peak, promoters):
    '''
    Does the range (start1, end1) overlap with (start2, end2)?
    Based on De Morgan's laws
    '''
    overlaps = promoters[np.logical_and(np.less_equal(peak['start'], promoters['end']),
                                 np.greater_equal(peak['end'], promoters['start']))]
    return(list(overlaps['gene_id']))


def peaks_intersect_genes(peaks, promoters):
    chrs_of_promoters = promoters['chr'].unique()
    chrs_of_peaks = peaks['chr'].unique()
    chrs = np.intersect1d(chrs_of_promoters, chrs_of_peaks)

    genes_id = []
    for chr_ in chrs:
        chr_peaks = pd.DataFrame(peaks[peaks['chr'] == chr_])
        chr_promoters = pd.DataFrame(promoters[promoters['chr'] == chr_])
        for index, peak in chr_peaks.iterrows():
            genes_id += overlap(peak, chr_promoters)
    return(genes_id)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', action='store', dest='gtf',
                        required=True, help='path to GTF file')
    parser.add_argument('-p', '--peaks', action='store', dest='peaks',
                        required=True, help='path to peaks file')
    parser.add_argument('-o', '--output', action='store', dest='genes_id',
                                  required=True, help='path to txt file to write genes_id')
    parser.add_argument('-l', '--left', action='store', type=int, dest='left',
                                  default=5000, required=False, help='left_tail + TSS, default_value = 5000')
    parser.add_argument('-r', '--right', action='store', type=int, dest='rigth',
                                  default=-5000, required=False, help='TSS + rigth_tail, default_value = -5000')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def write_results(out, res):
    res = [i.capitalize() for i in res]
    with open(out, 'w') as file:
        for gene_id in res:
            file.write(gene_id + '\n')


def main():
    args = parse_args()
    gtf_path = args.gtf
    peaks_path = args.peaks
    left = args.left
    rigth = args.rigth
    out = args.genes_id

    gtf = read_gtf(gtf_path)
    promoters = get_promoters(gtf, left=left, rigth=rigth)
    promoters = promoters.sort_values(by=['chr', 'start'])

    peaks = read_peaks(peaks_path)
    peaks = peaks.sort_values(by=['chr', 'start'])
    res = peaks_intersect_genes(peaks, promoters)
    res = np.unique(res)
    write_results(out, res)


if __name__ == '__main__':
    main()
