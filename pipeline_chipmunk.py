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
import numpy as np


def pipeline_chipmunk(bed_path, training_sample_size, testing_sample_size, peak_length,
                      fpr_for_thr, path_to_out, path_to_python_tools, dir_with_chipmunk,
                      path_to_promoters, path_to_genome, cpu_count,
                      zoops, try_size, bamm_order, recalculate_model):

    main_out = path_to_out + '/' + os.path.basename(bed_path).split('.')[0]
    zoops = str(zoops)
    bamm_order = str(bamm_order)
    try_size=str(try_size)
    cpu_count = str(cpu_count)
    peak_length = str(peak_length)
    motif_length_start = str(8)
    motif_length_end = str(12)

    if not path_to_python_tools[-1] == '/':
        path_to_python_tools += '/'

    if not dir_with_chipmunk[-1] == '/':
        dir_with_chipmunk += '/'

    if not os.path.isdir(main_out):
        os.mkdir(main_out)

    chipmunk = main_out + '/CHIPMUNK'
    scan = main_out + '/SCAN'
    motifs = main_out + '/MOTIFS'
    fasta = main_out + '/FASTA'
    bed = main_out + '/BED'
    compare_sites = main_out + '/COMPARE_SITES'
    tag = os.path.basename(bed_path).split('.')[0]

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
    else:
        print('File {0} already exists'.format(tag + '_' + str(testing_sample_size) + '.bed'))

    if not os.path.isfile(fasta + '/' + tag + '_' + str(training_sample_size) +'.fa'):
        #Bed peaks to fasta
        print('Bed peaks to fasta for {0}'.format(tag))
        args = ['python3', path_to_python_tools + 'bed_to_fasta.py',
               '-if', path_to_genome,
               '-bed', bed + '/' + tag + '_' + str(training_sample_size) +'.bed',
               '-of', fasta + '/' + tag + '_' + str(training_sample_size) +'.fa',
               '-l', peak_length]
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_' + str(training_sample_size) +'.fa'))

    if not os.path.isfile(fasta + '/' + tag + '_' + str(testing_sample_size) +'.fa'):
        #Bed peaks to fasta
        args = ['python3', path_to_python_tools + 'bed_to_fasta.py',
               '-if', path_to_genome,
               '-bed', bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
               '-of', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa']
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_' + str(testing_sample_size) +'.fa'))

    ########################
    #FIND MODEL BY ChIPMunk#
    ########################
    if not os.path.isfile(chipmunk + '/CHIPMUNK_MOTIF.txt'):
        print('ChIPMunk find motifs for {0}'.format(tag))
        args = ['java', '-cp', dir_with_chipmunk + 'chipmunk.jar',
                       'ru.autosome.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', zoops,
                       's:' + fasta + '/' + tag + '_'+ str(training_sample_size) + '.fa',
                      try_size, '10', '1', cpu_count, 'random']
        path_out = chipmunk + '/CHIPMUNK_MOTIF.txt'
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        out = p.communicate()
        with open(path_out, 'wb') as file:
            file.write(out[0])
    else:
        print('File {0} already exists'.format(chipmunk + '/CHIPMUNK_MOTIF.txt'))

    ###########################################################################
    #Parse results of ChIPMunk into files .meme, .pwm and .fasta (multi fasta)#
    ###########################################################################
    args = ['python3', path_to_python_tools + 'parse_chipmunk_results.py',
           '-i', chipmunk + '/CHIPMUNK_MOTIF.txt',
           '-o', chipmunk,
           '-t', tag + '_' + 'CHIPMUNK_MOTIF']
    r = subprocess.call(args)


    ##############################################################################
    #Get oPWM from ChIPMunk results. OUTPUT: .meme, .pwm and .fasta (multi fasta)#
    ##############################################################################

    args = ['python3', path_to_python_tools + 'make_oPWM.py',
            '-c', chipmunk + '/CHIPMUNK_MOTIF.txt',
            '-f', fasta + '/' + tag + '_'+ str(training_sample_size) + '.fa',
            '-n', '5000',
            '-P', cpu_count,
            '-o', chipmunk,
            '-t', tag + '_' + 'OPTIMAL_MOTIF']
    r = subprocess.call(args)





    ##################################
    #CALCULATE BAMM MODEL WITH EM ALG#
    ##################################

    if not os.path.isfile(motifs + '/' + tag + '_' + bamm_order + '_motif_1.ihbcp') or True:

        #Get BaMM motif
        print('Get Bamm motifs for {0}'.format(tag))
        args = ['BaMMmotif', motifs,
                fasta + '/' + tag + '_' + str(training_sample_size) + '.fa',
               '--PWMFile', motifs + '/' + tag + '_OPTIMAL_MOTIF.meme',
                '--basename', tag + '_' + bamm_order,
               '--EM',
               #'--CGS',
               #'--extend', '2',
               '--Order', bamm_order,
               '--order', bamm_order]
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_' + bamm_order + '_motif_1.ihbcp'))



    #####################################################
    #IF FILES WITH SCANED MOTIFS EXIST -> COMPARE MOTIFS#
    #####################################################

    if os.path.isfile(scan + '/' + tag + '_' + bamm_order + '_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed') and\
    os.path.isfile(scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed'):
        print(str(fpr_for_thr) + ' scan files already exist')
        #Compare sites (IF FILES WITH SCANED MOTIFS EXIST)
        print('Compare sites ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'compare_sites.py',
                '-p', bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
                '-m', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-b', scan + '/' + tag + '_' + bamm_order + '_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-t', tag + '_' + str(fpr_for_thr) + '_' + bamm_order,
                '-o', compare_sites]
        r = subprocess.call(args)

    else:

        #################################
        #CALCULATE THRESHOLDS FOR MODELS#
        #################################

        #Calculate threshold for PWM based on promoters and FPR = fpr_for_thr
        print('Calculate threshold for PWM based on promoters and FPR = {0} ({1})'.format(fpr_for_thr, tag))
        args = ['python3', path_to_python_tools + 'get_threshold_by_fp_numpy.py', 'pwm',
                '-f', path_to_promoters,
                '-m', motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
                '-p', str(fpr_for_thr),
                '-P', cpu_count]
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        thr_pwm = p.communicate()[0].decode('utf-8').strip()

        #Calculate threshold for BAMM based on promoters and FPR = fpr_for_thr
        print('Calculate threshold for BAMM based on promoters and FPR = {0} ({1})'.format(fpr_for_thr, tag))
        args = ['python3', path_to_python_tools + 'get_threshold_by_fp_numpy.py', 'bamm',
                '-f', path_to_promoters,
                '-m', motifs + '/' + tag + '_' + bamm_order + '_motif_1.ihbcp',
                '-b', motifs + '/' + tag + '_' + bamm_order + '.hbcp',
                '-p', str(fpr_for_thr),
                '-P', cpu_count]
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        thr_bamm = p.communicate()[0].decode('utf-8').strip()

        print('PWM = ',thr_pwm, 'BAMM = ', thr_bamm)

        ######################
        #SCAN SEQ-S BY MODELS#
        ######################

        #Scan peaks by PWM with thr_pwm
        print('Scan peaks by PWM with thr_pwm ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'scan_by_pwm.py',
                '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                '-m', motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
                '-t', thr_pwm,
                '-o', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-P', cpu_count]
        r = subprocess.call(args)

        #Scan peaks by BAMM with thr_bamm
        print('Scan peaks by BAMM with thr_pwm ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'scan_by_bamm.py',
                '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                '-m', motifs + '/' + tag + '_' + bamm_order + '_motif_1.ihbcp',
                '-b', motifs + '/' + tag + '_' + bamm_order + '.hbcp',
                '-t', thr_bamm,
                '-o', scan + '/' + tag + '_' + bamm_order + '_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-P', cpu_count]
        r = subprocess.call(args)



        ##############################
        #COMPARE SITES OF DIFF MODELS#
        ##############################

        #Compare sites
        print('Compare sites ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'compare_sites.py',
                '-p', bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
                '-m', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-b', scan + '/' + tag + '_' + bamm_order + '_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-t', tag + '_' + str(fpr_for_thr) + '_' + bamm_order,
                '-o', compare_sites]
        r = subprocess.call(args)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', action='store', dest='bed_path',
                        required=True, help='path to BED file')
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
    parser.add_argument('-p', '--python', action='store', dest='python_tools',
                        required=True, help='dir with python tools')
    parser.add_argument('-c', '--chipmunk', action='store', dest='chipmunk',
                        required=True, help='dir with chipmunk')
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
    parser.add_argument('-m', '--bamm_order_model', action='store', type=int, dest='bamm_order',
                        default=2, required=False,
                        help='Order of BaMM model. Default value = 2')
    parser.add_argument('-C', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=2, help='Number of processes to use, default: 2')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():

    args = parse_args()
    path_to_python_tools = args.python_tools
    dir_with_chipmunk = args.chipmunk
    path_to_promoters = args.promoters
    path_to_genome = args.genome
    bed_path = args.bed_path
    path_to_out = args.output
    training_sample_size = args.train_size
    testing_sample_size = args.test_size
    fpr_for_thr = args.fpr

    zoops=args.zoops
    cpu_count = args.cpu_count
    try_size=args.try_limit
    bamm_order=args.bamm_order
    recalculate_model=False
    peak_length = 100

    pipeline_chipmunk(bed_path, training_sample_size, testing_sample_size, peak_length,
                      fpr_for_thr, path_to_out, path_to_python_tools, dir_with_chipmunk,
                      path_to_promoters, path_to_genome, cpu_count,
                      zoops, try_size, bamm_order, recalculate_model)

if __name__ == '__main__':
    main()
