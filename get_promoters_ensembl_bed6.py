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

import pandas as pd
import numpy as np
import argparse
import sys


def read_bed6(path):

    bed = pd.read_csv(path,
                      sep='\t',comment='#', dtype= {'chr': str})
    bed.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    return(bed)


def get_promoters(bed, rigth=-2000, left=0):

    df = bed

    promoters = {'chr': [], 'start': [], 'end': [],
                 'name': [], 'score': [], 'strand': []}

    for line, strand in enumerate(df['strand']):

        if strand == '+':
            promoters['start'].append(df['start'].iloc[line] + rigth)
            promoters['end'].append(df['start'].iloc[line] + left)

        if strand == '-':
            promoters['start'].append(df['end'].iloc[line] - left)
            promoters['end'].append(df['end'].iloc[line] - rigth)

        promoters['chr'].append(df['chr'].iloc[line])
        promoters['name'].append(df['name'].iloc[line])
        promoters['score'].append(df['score'].iloc[line])
        promoters['strand'].append(df['strand'].iloc[line])

    return(pd.DataFrame(promoters))


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', action='store', dest='bed',
                        required=True, help='path to bed ensembl file')
    parser.add_argument('-o', '--output', action='store', dest='write',
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
    path = args.bed
    write = args.write
    left = args.left
    rigth = args.rigth

    bed = read_bed6(path)
    promoters = get_promoters(bed, left=left, rigth=rigth)
    promoters.to_csv(write,
                     header=False, index=False, sep='\t')


if __name__ == '__main__':
    main()
