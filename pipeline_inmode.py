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


def InMoDeCLI_denovo(path_to_inmode, fasta_path, motif_length,
                 model_order, outdir):

    args = ['/home/anton/Programs/jdk-9/bin/java', '--add-modules', 'java.xml.bind', '-jar' ,path_to_inmode + 'InMoDeCLI-1.1.jar', 'denovo',
            'i={}'.format(fasta_path),
            'm={}'.format(motif_length),
           'mo={}'.format(model_order),
           'outdir={}'.format(outdir)]
    r = subprocess.call(args)
    pass


def bedToFasta(path_to_fa, path_to_bed, out):

    args = ['bedtools', 'getfasta' , '-s', '-name+',
            '-fi', path_to_fa,
            '-bed', path_to_bed,
            '-fo', out]
    r = subprocess.call(args)
    pass



def InMoDeCLI_scan(path_to_inmode, input_data, input_model, backgroud_path,
                     fpr_for_thr, outdir):

    args = ['/home/anton/Programs/jdk-9/bin/java', '-Xmx4096m', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode + 'InMoDeCLI-1.1.jar', 'scan',
    #args = ['/Users/anton/Documents/Programs/jre-9.0.4.jre/Contents/Home/bin/java', '-Xmx3072m', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode + 'InMoDeCLI-1.1.jar', 'scan',
            'i={}'.format(input_model),
            'id={}'.format(input_data),
            'b={}'.format('From file'),
            'd={}'.format(backgroud_path),
           'f={}'.format(fpr_for_thr),
           'outdir={}'.format(outdir)]
    r = subprocess.call(args)
    pass


def pipeline_inmode_imd(bed_path, bigwig_path, training_sample_size, testing_sample_size,
                      fpr_for_thr, path_to_out, path_to_python_tools, path_to_inmode, path_to_imd,
                      path_to_promoters, path_to_genome, path_to_tss, cpu_count, shoulder, model_order, motif_length):

    main_out = path_to_out + '/' + os.path.basename(bed_path).split('.')[0]
    zoops = str(zoops)
    model_order = str(model_order)
    try_size=str(try_size)
    cpu_count = str(cpu_count)
    shoulder = str(shoulder)
    motif_length = str(8)
    path_to_tss = str(path_to_tss)

    if not path_to_python_tools[-1] == '/':
        path_to_python_tools += '/'

    if not dir_with_chipmunk[-1] == '/':
        dir_with_chipmunk += '/'

    if not os.path.isdir(main_out):
        os.mkdir(main_out)

    scan = main_out + '/SCAN'
    motifs = main_out + '/MOTIFS'
    fasta = main_out + '/FASTA'
    bed = main_out + '/BED'
    compare_sites = main_out + '/COMPARE_SITES'
    gene_ids =  main_out + '/IDs_COMPARE'
    tag = os.path.basename(bed_path).split('.')[0]

    name = 'INMODE'

    if not os.path.isdir(main_out + '/CHIPMUNK'):
        os.mkdir(main_out + '/CHIPMUNK')
    if not os.path.isdir(main_out + '/SCAN'):
        os.mkdir(main_out + '/SCAN')
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

    if not os.path.isfile(bed + '/' + tag + '_' + str(training_sample_size) + '.bed'):
        #Get top training_sample_size bed peaks
        print('Get top {0} bed peaks for {1}'.format(training_sample_size, tag))
        args = ['python3', path_to_python_tools + 'get_top_peaks.py',
               '-i', bed_path,
               '-o', bed,
               '-a', str(training_sample_size),
               '-c', '4',
               '-t', tag + '_' + str(training_sample_size)]
        r = subprocess.call(args)

        if shoulder != '-1':
            args = ['python3', path_to_python_tools + 'prepare_peaks.py',
                   '-b', bed + '/' + tag + '_' + str(training_sample_size) + '.bed',
                    '-w', bigwig_path,
                   '-o', bed,
                    '-s', shoulder,
                   '-t', tag + '_' + str(training_sample_size)]
            r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_' + str(training_sample_size) + '.bed'))

    if not os.path.isfile(bed + '/' + tag + '_' + str(testing_sample_size) + '.bed'):
        #Get top testing_sample_size bed peaks
        print('Get top {1} bed peaks for {0}'.format(tag, testing_sample_size))
        args = ['python3', path_to_python_tools + 'get_top_peaks.py',
               '-i', bed_path,
               '-o', bed,
               '-a', str(testing_sample_size),
               '-c', '4',
               '-t', tag + '_' + str(testing_sample_size)]
        r = subprocess.call(args)

        if shoulder != '-1':
            args = ['python3', path_to_python_tools + 'prepare_peaks.py',
                   '-b', bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
                    '-w', bigwig_path,
                   '-o', bed,
                    '-s', shoulder,
                   '-t', tag + '_' + str(testing_sample_size)]
            r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_' + str(testing_sample_size) + '.bed'))

    if not os.path.isfile(fasta + '/' + tag + '_' + str(training_sample_size) +'.fa'):
        #Bed peaks to fasta
        print('Bed peaks to fasta for {0}'.format(tag))
        bedToFasta(path_to_genome,
            bed + '/' + tag + '_' + str(training_sample_size) +'.bed',
            fasta + '/' + tag + '_' + str(training_sample_size) +'.fa')

    else:
        print('File {0} already exists'.format(tag + '_' + str(training_sample_size) +'.fa'))

    if not os.path.isfile(fasta + '/' + tag + '_' + str(testing_sample_size) +'.fa'):
        bedToFasta(path_to_genome,
            bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
            fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa')

    else:
        print('File {0} already exists'.format(tag + '_' + str(testing_sample_size) +'.fa'))


    ####################################
    #CALCULATE INMODE MODEL WITH EM ALG#
    ####################################
    if glob.glob(motifs + '/Learned_DeNovo*') == []:
        print('Calculate InMoDe model')
        InMoDeCLI_denovo(path_to_inmode,
                         fasta_path=fasta + '/' + tag + '_'+ str(training_sample_size) + '.fa',
                         motif_length=motif_length,
                         model_order=model_order,
                         outdir=motifs)
    else:
        print('InMoDe model already exists')


    ################################################
    #CALCULATE THRESHOLDS FOR INMODE MODEL AND SCAN#
    ################################################

    if not os.path.isfile(scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed'):
        print('Scan by inmode model')
        InMoDeCLI_scan(path_to_inmode,
                       input_data=fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                       input_model=glob.glob(motifs + '/Learned_DeNovo*/*.xml')[0],
                       backgroud_path=path_to_promoters,
                       fpr_for_thr=fpr_for_thr,
                       outdir=scan)

        args = ['python3', path_to_python_tools + 'parse_inmode_scan.py',
                '-if', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                '-bed', glob.glob(scan + '/*.BED')[0],
                '-o', scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed']
        r = subprocess.call(args)
    else:
        print(tag + '_INMODE_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed', '- EXISTS')


    #####################
    #APPLY IMD ALGORITHM#
    #####################

    args = ['/home/anton/Programs/jdk-9/bin/java', '--add-modules', 'java.xml.bind',
    '-jar', path_to_imd, 'imd',
    'i=' + scan + '/' + 'inmode.all.sites.txt',
    'outdir=' + scan + '/inmode']
    r = subprocess.call(args)


    ### COMPARE IDS ###

    print('Compare sites ({0})'.format(tag))
    args = ['python3', path_to_python_tools + 'compare_gene_ids.py',
            '-first', gene_ids + '/' + 'pwm.ids.txt',
            '-second', gene_ids + '/' + 'bamm.ids.txt',
            '-third', gene_ids + '/' + 'inmode.ids.txt',
            '-o', gene_ids,
            '-fname', fname,
            '-sname', sname,
            '-tname', tname]
    r = subprocess.call(args)

    ### IMD FOR SITES ###

    print('EXTRACT SITES ({0})'.format(tag))

    args = ['python3', path_to_python_tools + 'extract_sites.py',
    '-p', scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
    '-o', scan + '/' + 'inmode.all.sites.txt']
    r = subprocess.call(args)

    args = ['/home/anton/Programs/jdk-9/bin/java', '--add-modules', 'java.xml.bind',
    '-jar', path_to_imd, 'imd',
    'i=' + scan + '/' + 'inmode.all.sites.txt',
    'outdir=' + scan + '/inmode']
    r = subprocess.call(args)



def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', action='store', dest='bed_path',
                        required=True, help='path to BED file')
    parser.add_argument('-w', '--bigwig', action='store', dest='wig_path',
                        required=True, help='path to BIGWIG file')
    parser.add_argument('-P', '--promoters', action='store', dest='promoters',
                        required=True, help='path to promoters fasta file')
    parser.add_argument('-g', '--genome', action='store', dest='genome',
                        required=True, help='path to genome fasta file')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=True, help='size of training sample')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=True, help='size of testing sample')
    parser.add_argument('-f', '--fpr', action='store', dest='fpr',
                        required=False, default=0.0001, type=float,
                        help='FPR value required to calculate threshold values \
                        default=0.0001')
    parser.add_argument('-l', '--motif_length', action='store', dest='motif_length',
                        required=False, default=15, type=float,
                        help='Length of motif (def = 15)')
    parser.add_argument('-p', '--python', action='store', dest='python_tools',
                        required=True, help='dir with python tools')
    parser.add_argument('-I', '--InMoDe', action='store', dest='inmode',
                        required=True, help='path to InMoDe')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='output dir')
    parser.add_argument('-m', '--order_model', action='store', type=int, dest='model_order',
                        default=2, required=False,
                        help='Order of model. Default value = 2')
    parser.add_argument('-C', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=2, help='Number of processes to use, default: 2')
    parser.add_argument('-s', '--shoulder', action='store', dest='shoulder', default=50,
                        required=False, type=int, help='summit +/- shoulder (extend peak) default=50')
    parser.add_argument('-tss', action='store', dest='path_to_tss',
                        required=True, help='path to BED file with transcripts')
    parser.add_argument('-i', '--imd', action='store', dest='path_to_imd',
                        required=True, help='path to DisentanglerCLI to run imd')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():

    args = parse_args()

    dir_with_chipmunk = args.chipmunk
    bed_path = args.bed_path
    bigwig_path = args.wig_path
    path_to_out = args.output
    training_sample_size = args.train_size
    testing_sample_size = args.test_size
    fpr_for_thr = args.fpr
    shoulder = args.shoulder
    motif_length = args.motif_length


    path_to_python_tools = args.python_tools
    path_to_inmode = args.inmode
    path_to_promoters = args.promoters
    path_to_genome = args.genome
    path_to_imd = args.path_to_imd
    path_to_tss = args.path_to_tss

    cpu_count = args.cpu_count
    model_order=args.model_order
    
    def pipeline_inmode_imd(bed_path, bigwig_path, training_sample_size, testing_sample_size,
                      fpr_for_thr, path_to_out, path_to_python_tools, path_to_inmode, path_to_imd,
                      path_to_promoters, path_to_genome, path_to_tss, cpu_count, shoulder, model_order,motif_length)

if __name__ == '__main__':
    main()
