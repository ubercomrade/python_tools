import pandas as pd
import numpy as np
import os
import argparse
import sys
import subprocess


def get_bed_wig(bed3, bigwig, bedwig):
        args = ['bwtool', 'extract', 'bed',
               bed3,
               bigwig,
               bedwig]
        r = subprocess.call(args)


def get_bed3_form_bed6(bed6, bed3):
    bed = pd.read_csv(bed6, sep='\t', header=None)
    bed = bed.loc[:,0:2]
    bed.to_csv(bed3, sep='\t', header=None, index=False)


def read_fasta(path):

    fasta = list()
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                fasta.append(line.strip().upper())
    return(fasta)


def read_bedwig(path):
    bed = pd.read_csv(path, sep='\t', header=None)
    scores = list(bed[4])
    return(scores)


def write_fastawig(path, fasta, bedwig):
    with open(path, "w") as file:
        for f, w in zip(fasta, bedwig):
            file.write(">{0}\n{1}\n".format(' '.join(w.split(',')), f))
    file.close()
    pass


def get_summit_regions(bedwig, peaks, shoulder=50):
    df = pd.read_csv(bedwig, sep='\t', header=None,
                names=['chr', 'start', 'end', 'name', 'length', 'scores'])

    summits = []
    for i in range(len(df)):
        scores = df.loc[i, 'scores']
        scores = scores.split(',')
        scores = [float(i) if i != 'NA' else 1.0 for i in scores]
        start, end = scores.index(max(scores)), len(scores) - scores[::-1].index(max(scores)) - 1
        summit = (start + end) // 2
        summits.append(summit)
        start = summit - shoulder
        end = summit + shoulder

        begin = df.loc[i, 'start']
        df.loc[i, 'start'] = begin + start
        df.loc[i, 'end'] = begin + end

    df_write = df.loc[:,['chr', 'start', 'end', 'name']]
    df_write['scores'] = '.'
    df_write['strand'] = '.'
    df_write.to_csv(peaks, sep='\t', header=None, index=False)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--inputBED', action='store', dest='bed6',
                        required=True, help='path to BED file')
    parser.add_argument('-w', '--inputBigWig', action='store', dest='bigwig',
                        required=True, help='path to BigWig file')
    parser.add_argument('-f', '--inputFasta', action='store', dest='fasta',
                        required=True, help='path to BigWig file')
    parser.add_argument('-o', '--outputDir', action='store', dest='odir',
                        required=True, help='dir for write file')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output file')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    odir = args.odir
    bed6 = args.bed6
    bigwig = args.bigwig
    tag = args.tag
    fasta_in = args.fasta
    bed3 = odir + '/' + 'tmp.bed3'
    bedwig = odir + '/' + 'tmp.bedwig'
    fasta_out = odir + '/' + tag + '.faWig'

    get_bed3_form_bed6(bed6, bed3)
    get_bed_wig(bed3, bigwig, bedwig)
    f = read_fasta(fasta_in)
    bdw = read_bedwig(bedwig)
    for i in range(len(bdw)):
        bdw[i] = bdw[i].replace("NA","0.0")
    write_fastawig(fasta_out, f, bdw)
    os.remove(bed3)
    os.remove(bedwig)


if __name__=='__main__':
    main()
