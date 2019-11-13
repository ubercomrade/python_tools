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

import os
import sys
import shlex
import subprocess
import argparse
import glob
import numpy as np
import pandas as pd
import math


def inmode_scan(path_to_inmode, path_java, input_data, input_model, tmp_dir,
                     fpr_for_thr=1):

    args = [path_java, '-Xmx6G', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode, 'scan',
            'i={}'.format(input_model),
            'id={}'.format(input_data),
           'f={}'.format(fpr_for_thr),
           'outdir={}'.format(tmp_dir)]
    r = subprocess.call(args)
    pass


def parse_inmode_results(inmode_bed, out):

    container = list()
    append = container.append
    df = pd.read_csv(inmode_bed, sep='\t', header=None)
    df = df.sort_values(by=[0,4])
    index = 0
    for number, line in df.iterrows():
        if line[0] != index:
            append(last_score)
        index = line[0]
        last_score = line[4]
    append(line[4])

    with open(out, 'w') as file:
        for i in container:
            file.write('{}\n'.format(math.log2(float(i))))
    file.close()
    pass
    


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', action='store', help='path to FASTA file')
    parser.add_argument('model', action='store', help='path to xml inmode model file')
    parser.add_argument('output', action='store', help='prth to write results')
    parser.add_argument('inmode', action='store', help='path to InMoDe')
    parser.add_argument('-j', '--java', action='store', dest='java',
                            required=False, default='java', help='path to java')
    parser.add_argument('-t', '--tmp', action='store', dest='tmp',
                            required=False, default='./tmp', help='tmp dir')
    return(parser.parse_args())


def main():

    args = parse_args()
    path_to_model = args.model
    path_to_inmode = args.inmode
    fasta_path = args.fasta
    out = args.output
    tmp_dir = args.tmp
    path_to_java = args.java

    inmode_scan(path_to_inmode, path_to_java, fasta_path, path_to_model, tmp_dir)
    inmode_bed = glob.glob(tmp_dir + '/*.BED')[0]
    parse_inmode_results(inmode_bed, out)
    os.system("rm -r {}".format(tmp_dir))

if __name__ == '__main__':
    main()
