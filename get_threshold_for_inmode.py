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

import argparse
import sys
import os
import subprocess
from collections import Counter
from bisect import bisect, bisect_left
from operator import itemgetter
import math


def read_fasta(path):

    fasta = list()
    append = fasta.append
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                line = line.strip().upper()
                if not 'N' in line:
                    append(line)
                    append(complement(line))
    file.close()
    return(fasta)


def complement(seq):
    #return(seq[::-1].translate(seq.maketrans('ACGT', 'TGCA')))
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').replace('N', 'n').upper()[::-1])


def calculate_scores(path_to_inmode, path_to_model, path_to_fasta, path_to_java, tmp_dir):
    container = list()
    append = container.append

    # command = "{0} -Xmx4G -Xms1024m --add-modules java.xml.bind -jar {1} scan i={2} id={3} b='From file' d={3} f=0.01 outdir={4} bs=true".format(path_to_java,
    # 	path_to_inmode,
    # 	path_to_model,
    # 	path_to_fasta,
    # 	tmp_dir)
    args = [path_to_java, '-Xmx4096m', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode, 'scan',
        'i={}'.format(path_to_model),
        'id={}'.format(path_to_fasta),
        'b={}'.format('From file'),
        'd={}'.format(path_to_fasta),
       'f={}'.format(0.005),
       'outdir={}'.format(tmp_dir)]
    r = subprocess.call(args)
    # print(command)
    # os.system(command)

    with open(tmp_dir + "/Motif_hits_from_SequenceScan(0.005).BED") as file:
    	for line in file:
    		append(float(line.strip().split()[4]))

    return(container)


def get_threshold(scores, number_of_sites, path_out):
    scores.sort(reverse=False) # sorted score from small to big
    scores = [round(math.log(float(i), 10), 3) for i in scores]
    counts = Counter(scores)

    found_sites = counts[list(counts.keys())[::-1][0]]
    for index, k in enumerate(list(counts.keys())[::-1][1:]):
        found_sites = counts[k] + found_sites
        counts[k] = found_sites

    fprs_table = []
    append = fprs_table.append
    for index, k in enumerate(list(counts.keys())[::-1]):
        append(( k, counts[k] / number_of_sites ))


    with open(path_out, "w") as file:
        file.write("Scores\tFPR\n")
        for (score, fpr) in fprs_table:
            file.write("{0}\t{1}\n".format(score, fpr))
    file.close()


# def get_threshold(scores, number_of_sites, path_out):
#     scores.sort(reverse=True) # sorted score from big to small
#     fprs = [5*10**(-4), 3.33*10**(-4), 1.90*10**(-4), 1.02*10**(-4), 5.24*10**(-5)]
#     with open(path_out, "w") as file:
#         file.write("Scores\tFPR\n")
#         for fpr in fprs:
#         	thr = scores[int(fpr * number_of_sites)]
#         	file.write("{0}\t{1}\n".format(thr, fpr))
#     file.close()



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('fasta', action='store', help='path to fasta')
    parser.add_argument('model', action='store', help='path to .xml file (Inmode model)')
    parser.add_argument('inmode', action='store', help='path to inmode program')
    parser.add_argument('len', action='store', type=int, help='len of TF site')
    parser.add_argument('out', action='store', help='path to write results')
    parser.add_argument('-j', '--java', action='store', dest='java',
                            required=False, default="java", help='path to java')
    parser.add_argument('-t', '--tmp', action='store', dest='tmp_dir',
                            required=False, default="./tmp", help='tmp dir')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()

    path_to_model = args.model
    path_to_fasta = args.fasta
    path_out = args.out
    path_to_java = args.java
    path_to_inmode = args.inmode
    length_of_site = args.len
    tmp_dir = args.tmp_dir

    if not os.path.isdir(tmp_dir):
    	os.mkdir(tmp_dir)

    peaks = read_fasta(path_to_fasta)
    number_of_sites = sum([len(range(len(peak) - length_of_site + 1)) for peak in peaks])
    print(path_to_fasta)
    print(path_to_model)
    print(tmp_dir)

    scores = calculate_scores(path_to_inmode, path_to_model, path_to_fasta, path_to_java, tmp_dir)
    get_threshold(scores, number_of_sites, path_out)

    os.system("rm -r {}".format(tmp_dir))


if __name__ == '__main__':
    main()