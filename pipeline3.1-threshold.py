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
import math
import bisect
from operator import itemgetter
import numpy as np
import pandas as pd


def inmode_denovo(path_to_java, path_to_inmode, fasta_path, motif_length,
                 model_order, outdir):

    args = [path_to_java, '--add-modules', 'java.xml.bind', '-jar' ,path_to_inmode, 'denovo',
            'i={}'.format(fasta_path),
            'm={}'.format(motif_length),
           'mo={}'.format(model_order),
           'outdir={}'.format(outdir)]
    r = subprocess.call(args)
    pass


def bed_to_fasta(path_to_fa, path_to_bed, out):

    args = ['bedtools', 'getfasta' , '-s', '-name+',
            '-fi', path_to_fa,
            '-bed', path_to_bed,
            '-fo', out]
    r = subprocess.call(args)
    pass


def get_threshold(path, fpr_for_thr):
    conteiner = list()
    append = conteiner.append
    
    with open(path, 'r') as file:
        file.readline()
        for line in file:
            append(tuple(map(float, line.strip().split())))
    file.close()
    
    conteiner = sorted(conteiner, key=itemgetter(1))
    getcount = itemgetter(1)
    score = conteiner[bisect.bisect_left(list(map(getcount, conteiner)), fpr_for_thr)][0]
    return(score)



def get_top_peaks(path_to_python_tools, bed_in, bed_out, size, tag):
    args = ['python3', path_to_python_tools + '/get_top_peaks.py',
           '-i', bed_in,
           '-o', bed_out,
           '-a', str(size),
           '-c', '4',
           '-t', tag]
    r = subprocess.call(args)
    pass


def bootstrap_pwm(path_to_python_tools, out_path, sites):
    args = ['julia', path_to_python_tools + '/bootstrap_for_pwm.jl',
            out_path,
            sites, '-s', '2000000']
    r = subprocess.call(args)
    pass

def bootstrap_inmode(path_to_python_tools, path_to_java, out_path, sites, path_to_inmode):
    args = ['julia', path_to_python_tools + '/bootstrap_for_inmode.jl',
            out_path,
            sites,
            path_to_inmode, '-s', '2000000', '-j', path_to_java]
    r = subprocess.call(args)
    pass

def bootstrap_bamm(path_to_python_tools, out_path, sites):
    args = ['julia', path_to_python_tools + '/bootstrap_for_bamm.jl',
            out_path,
            sites, '-s', '2000000']
    r = subprocess.call(args)
    pass


def run_chipmunk_fasta(path_to_java, path_to_chipmunk, fasta_path, path_out, motif_length_start, motif_length_end, try_size, cpu_count, zoops):
    args = [path_to_java, '-cp', path_to_chipmunk,
                   'ru.autosome.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', zoops,
                   's:' + fasta_path,
                  try_size, '10', '1', cpu_count, 'random']
    p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    out = p.communicate()
    with open(path_out, 'wb') as file:
        file.write(out[0])
    pass


def inmode_scan(path_to_java, path_to_inmode, input_data, input_model, backgroud_path,
                     fpr_for_thr, outdir):

    args = [path_to_java, '-Xmx6G', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode, 'scan',
            'i={}'.format(input_model),
            'id={}'.format(input_data),
            'b={}'.format('From file'),
            'd={}'.format(backgroud_path),
           'f={}'.format(fpr_for_thr),
           'outdir={}'.format(outdir)]
    r = subprocess.call(args)
    pass


# def scan_best_by_inmode(path_to_python_tools, output, input_model, fasta_in, path_to_inmode):
#     args = ['julia', path_to_python_tools + '/scan_best_by_inmode.jl',
#            output,
#            input_model,
#            fasta_in,
#            path_to_inmode]
#     r = subprocess.call(args)
#     pass

def scan_best_by_inmode(path_to_python_tools, output, input_model, fasta_in, path_to_inmode, path_to_java):
    args = ['pypy', path_to_python_tools + '/scan_best_by_inmode_alt.py',
           fasta_in,
           input_model,
           output,
           path_to_inmode,
           '-j', path_to_java]
    r = subprocess.call(args)
    pass


def plot_best_score(path_to_python_tools, model1, model2, thr1, thr2, length, out, name1, name2):
    args = ['python3', path_to_python_tools + '/plot_best_score.py',
           model1,
           model2,
           thr1,
           thr2,
           length,
           out,
           '-n1', name1,
           '-n2', name2]
    r = subprocess.call(args)
    pass


def scan_best_by_pwm(path_to_python_tools, output, input_model, fasta_in, cpu_count):
    args = ['python3', path_to_python_tools + 'scan_best_by_pwm.py',
            '-f', fasta_in,
            '-m', input_model,
            '-o', output,
            '-P', cpu_count]
    r = subprocess.call(args)
    pass


def scan_best_by_bamm(path_to_python_tools, output, input_bamm_model, bg_model, fasta_in, cpu_count):
    args = ['python3', path_to_python_tools + '/scan_best_by_bamm.py',
            '-f', fasta_in,
            '-m', input_bamm_model,
            '-b', bg_model,
            '-o', output,
            '-P', cpu_count]
    r = subprocess.call(args)
    pass


def run_tomtom(query, model, outdir):
    args = ['tomtom', query, model, '-oc', outdir]
    r = subprocess.call(args)
    pass


def montecarlo(path_to_python_tools, scores1, scores2, thr1, thr2, length, results, name1, name2, fpr1, fpr2):
    with open(results, 'a') as file:
        file.write('{0}:thr={1},fpr={2};{3}:thr={4},fpr={5}\n'.format(name1, thr1, fpr1, name2, thr2, fpr2))
    file.close()
    args = ['monteCarlo', '{}'.format(scores1), '{}'.format(scores2), '{}'.format(thr1), '{}'.format(thr2), '{}'.format(length), '{}'.format(results)]
    r = subprocess.call(args)
    pass


# def montecarlo(path_to_python_tools, scores1, scores2, thr1, thr2, length, results, name1, name2, fpr1, fpr2):
#     with open(results, 'a') as file:
#         file.write('{0}:thr={1},fpr={2};{3}:thr={4},fpr={5}\n'.format(name1, thr1, fpr1, name2, thr2, fpr2))
#     file.close()
#     args = ['pypy', path_to_python_tools + 'montecarlo.py', '{}'.format(scores1), '{}'.format(scores2), '{}'.format(thr1), '{}'.format(thr2), '{}'.format(length), '{}'.format(results)]
#     r = subprocess.call(args)
#     pass


def corr_test(path_to_python_tools, scores1, scores2, results, name1, name2):
    with open(results, 'a') as file:
        file.write('{0};{1}\n'.format(name1, name2))
    file.close()
    args = ['python3', path_to_python_tools + 'corr_test.py', '{}'.format(scores1), '{}'.format(scores2), '{}'.format(results)]
    r = subprocess.call(args)
    pass


def binome_test(path_to_python_tools, scores1, scores2, thr1, thr2, results, name1, name2, fpr1, fpr2):
    with open(results, 'a') as file:
        file.write('{0}:thr={1},fpr={2};{3}:thr={4},fpr={5}\n'.format(name1, thr1, fpr1, name2, thr2, fpr2))
    file.close()
    args = ['python3', path_to_python_tools + 'binome_test.py', '{}'.format(scores1), '{}'.format(scores2), '{}'.format(thr1), '{}'.format(thr2), '{}'.format(results)]
    r = subprocess.call(args)
    pass


def make_model(path_to_python_tools, path_in, dir_out, tag):
    args = ['python3', path_to_python_tools + '/make_model.py',
            '-i', path_in,
            '-o', dir_out,
            '-t', tag,
            '-M']
    r = subprocess.call(args)
    pass  



def pipeline_inmode_bamm(bed_path, training_sample_size, testing_sample_size,
                      path_to_out, path_to_python_tools, path_to_java, path_to_inmode, path_to_imd, path_to_chipmunk,
                      path_to_promoters, path_to_genome, path_to_tss, path_to_hocomoco, cpu_count,
                      zoops, try_size, model_order):

    main_out = path_to_out + '/' + os.path.basename(bed_path).split('.')[0]
    zoops = str(zoops)
    model_order = str(model_order)
    try_size=str(try_size)
    cpu_count = str(cpu_count)
    motif_length_start = str(8)
    motif_length_end = str(14)
    path_to_tss = str(path_to_tss)


    if not path_to_python_tools[-1] == '/':
        path_to_python_tools += '/'


    if not os.path.isdir(main_out):
        os.mkdir(main_out)

    chipmunk = main_out + '/CHIPMUNK'
    scan = main_out + '/SCAN'
    motifs = main_out + '/MOTIFS'
    fasta = main_out + '/FASTA'
    bed = main_out + '/BED'
    bootstrap = main_out + '/BOOTSTRAP'
    scan_best = main_out + '/SCAN-BEST'
    compare_sites = main_out + '/COMPARE_SITES'
    gene_ids =  main_out + '/IDs_COMPARE'
    tomtom = main_out + '/TOMTOM'
    tag = os.path.basename(bed_path).split('.')[0]

    fname = 'PWM'
    sname = 'BAMM'
    tname = 'INMODE'

    ########################
    #      CREATE DIRS     #
    ########################

    if not os.path.isdir(main_out + '/CHIPMUNK'):
        os.mkdir(main_out + '/CHIPMUNK')
    if not os.path.isdir(main_out + '/SCAN'):
        os.mkdir(main_out + '/SCAN')
    if not os.path.isdir(main_out + '/BOOTSTRAP'):
        os.mkdir(main_out + '/BOOTSTRAP')
    if not os.path.isdir(main_out + '/SCAN-BEST'):
        os.mkdir(main_out + '/SCAN-BEST')
    if not os.path.isdir(main_out + '/MOTIFS'):
        os.mkdir(main_out + '/MOTIFS')
    if not os.path.isdir(main_out + '/FASTA'):
        os.mkdir(main_out + '/FASTA')
    if not os.path.isdir(main_out + '/BED'):
        os.mkdir(main_out + '/BED')
    if not os.path.isdir(main_out + '/COMPARE_SITES'):
        os.mkdir(main_out + '/COMPARE_SITES')
    if not os.path.isdir(main_out + '/IDs_COMPARE'):
        os.mkdir(main_out + '/IDs_COMPARE')
    if not os.path.isdir(main_out + '/TOMTOM'):
        os.mkdir(main_out + '/TOMTOM')



    ##########################
    #  CALCULATE THRESHOLDS  #
    ##########################

    if not os.path.isfile(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt'):
        print('Calculate threshold for BAMM based on promoters and FPR')
        args = ['pypy', path_to_python_tools + '/get_threshold_for_bamm.py',
                path_to_promoters,
                motifs + '/' + tag + '_motif_1.ihbcp',
                motifs + '/' + tag + '.hbcp',
                motifs + '/' + tag + '_BAMM_THRESHOLDS.txt']
        r = subprocess.call(args)
    else:
        print('Thresholds for bamm already calculated')

    if not os.path.isfile(motifs + '/' + tag + '_PWM_THRESHOLDS.txt'):
        print('Calculate threshold for PWM based on promoters and FPR')
        args = ['pypy', path_to_python_tools + '/get_threshold_for_pwm.py',
                path_to_promoters,
                motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
                motifs + '/' + tag + '_PWM_THRESHOLDS.txt']
        r = subprocess.call(args)
    else:
        print('Thresholds for pwm already calculated')

    if not os.path.isfile(motifs + '/' + tag + '_PWM_THRESHOLDS.txt'):
        args = ['python3', path_to_python_tools + '/get_threshold_for_inmode.py',
                path_to_promoters,
                glob.glob(motifs + '/Learned_DeNovo*/*.xml')[0],
                path_to_inmode,
                str(motif_length),
                motifs + '/' + tag + '_INMODE_THRESHOLDS.txt',
                '-j', path_to_java]
        r = subprocess.call(args)
    else:
        print('Thresholds for inmode already calculated')



def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store', help='path to BED file')
    parser.add_argument('promoters', action='store', help='path to promoters fasta file')
    parser.add_argument('genome', action='store', help='path to genome fasta file')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=True, help='size of training sample')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=True, help='size of testing sample')
    parser.add_argument('-p', '--python', action='store', dest='python_tools',
                        required=True, help='dir with python tools')
    parser.add_argument('-I', '--inmode', action='store', dest='inmode',
                        required=True, help='path to inmode')
    parser.add_argument('-J', '--java', action='store', dest='java',
                    required=False, default="java", help='path to Java')
    parser.add_argument('-c', '--chipmunk', action='store', dest='chipmunk',
                        required=True, help='path to chipmunk')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='output dir')
    parser.add_argument('-z', '--zoops', action='store', type=float, dest='zoops',
                        default=1.0, required=False,
                        help='zero-or-one-occurrence-per-sequence (ZOOPS). You should specify the \
                        zoops factor parameter, a value between 0 and 1.0. Default value = 1.0')
    parser.add_argument('-l', '--try_limit', action='store', type=int, dest='try_limit',
                        default=100, required=False,
                        help=' This is an internal number of motif optimization runs. \
                        For a random seeding, this would be simply equal to the number of seeds. \
                        It can be as high as your computational power \
                        (100-1000 seems to be generally enough depending on your dataset). Default value = 100')
    parser.add_argument('-m', '--order_model', action='store', type=int, dest='model_order',
                        default=2, required=False,
                        help='Order of model. Default value = 2')
    parser.add_argument('-C', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=2, help='Number of processes to use, default: 2')
    parser.add_argument('-tss', action='store', dest='path_to_tss',
                        required=True, help='path to BED file with transcripts')
    parser.add_argument('-i', '--imd', action='store', dest='path_to_imd',
                        required=True, help='path to DisentanglerCLI to run imd')
    parser.add_argument('-H', '--hocomoco', action='store', dest='path_to_hocomoco',
                        required=True, help='path to HOCOMOCO database in meme format for TOMTOM')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():

    args = parse_args()

    bed_path = args.bed
    path_to_out = args.output
    training_sample_size = args.train_size
    testing_sample_size = args.test_size


    path_to_python_tools = args.python_tools
    path_to_java = args.java
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    path_to_promoters = args.promoters
    path_to_genome = args.genome
    path_to_imd = args.path_to_imd
    path_to_tss = args.path_to_tss
    path_to_hocomoco = args.path_to_hocomoco


    zoops=args.zoops
    cpu_count = args.cpu_count
    try_size=args.try_limit
    model_order=args.model_order

    pipeline_inmode_bamm(bed_path, training_sample_size, testing_sample_size,
                          path_to_out, path_to_python_tools, path_to_java, path_to_inmode, path_to_imd, path_to_chipmunk,
                          path_to_promoters, path_to_genome, path_to_tss, path_to_hocomoco, cpu_count,
                          zoops, try_size, model_order)

if __name__ == '__main__':
    main()
