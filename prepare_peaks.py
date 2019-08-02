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
    bed = bed.loc[:,0:3]
    bed.to_csv(bed3, sep='\t', header=None, index=False)


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
    parser.add_argument('-o', '--outputDir', action='store', dest='odir',
                        required=True, help='dir for write file')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output file')
    parser.add_argument('-s', '--shoulder', action='store', dest='shoulder', default=50,
                        required=False, type=int, help='summit +/- shoulder (extend peak) default=50')

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
    shoulder = args.shoulder
    bed3 = odir + '/' + 'tmp.bed3'
    bedwig = odir + '/' + 'tmp.bedwig'
    peaks = odir + '/' + tag + '.bed'

    get_bed3_form_bed6(bed6, bed3)
    get_bed_wig(bed3, bigwig, bedwig)
    get_summit_regions(bedwig, peaks, shoulder=shoulder)
    os.remove(bed3)
    os.remove(bedwig)


if __name__=='__main__':
    main()
