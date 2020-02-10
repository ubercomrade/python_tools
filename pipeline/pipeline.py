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
from shutil import copyfile


def prepare_data(path_to_python_tools, path_to_genome, bed_path, bed, fasta, train_sample_size, test_sample_size):

    ########################
    #     GET TOP PEAKS    #
    ########################

    if not os.path.isfile(bed + '/' + 'train_sample.bed'):
        #Get top training_sample_size bed peaks
        print('Get top {0} bed peaks'.format(train_sample_size))
        bed_out = bed + '/'
        get_top_peaks(path_to_python_tools, bed_path, bed_out, train_sample_size, 'train_sample')      
    else:
        print('File {0} already exists'.format('train_sample.bed'))

    if not os.path.isfile(bed + '/' + 'test_sample.bed'):
        #Get top testing_sample_size bed peaks
        print('Get top {0} bed peaks'.format(test_sample_size))
        bed_out = bed + '/'
        get_top_peaks(path_to_python_tools, bed_path, bed_out, test_sample_size, 'test_sample')
    else:
        print('File {0} already exists'.format('test_sample.bed'))

    ########################
    #     BED TO FASTA     #
    ########################

    if not os.path.isfile(fasta + '/' + 'train_sample.fa'):
        #Bed peaks to fasta
        print('Get fasta from bed: {}'.format('train_sample.bed'))
        bed_to_fasta(path_to_genome,
            bed + '/train_sample.bed',
            fasta + '/train_sample.fa')
    else:
        print('File {0} already exists'.format('train_sample.fa'))

    if not os.path.isfile(fasta + '/' + 'test_sample.fa'):
        print('Get fasta from bed: {}'.format('test_sample.bed'))
        bed_to_fasta(path_to_genome,
            bed + '/test_sample.bed',
            fasta + '/test_sample.fa')
    else:
        print('File {0} already exists'.format('test_sample.fa'))
    pass


def get_inmode_model(models_path, fasta_path, path_to_java, path_to_inmode, motif_length, model_order):
    
    inmode_model_path = models_path + '/inmode_model'

    if not os.path.isdir(inmode_model_path):
        os.mkdir(inmode_model_path)

    if glob.glob(inmode_model_path + '/Learned_DeNovo*') == []:
        print('calculate inmode model')

        args = [path_to_java, '--add-modules', 'java.xml.bind', '-jar', path_to_inmode, 'denovo',
                'i={}'.format(fasta_path),
                'm={}'.format(motif_length),
               'mo={}'.format(model_order),
               'outdir={}'.format(inmode_model_path)]
        r = subprocess.call(args)
    else:
        print('inmode model already exists')

    copyfile(glob.glob(inmode_model_path + '/Learned_DeNovo*/XML*')[0], inmode_model_path + 'inmode_model.xml')
    copyfile(glob.glob(inmode_model_path + '/Learned_DeNovo*/Binding_sites*')[0], inmode_model_path + 'inmode_sites.txt')
    pass


def run_chipmunk_fasta(path_to_java, path_to_chipmunk, fasta_path, path_out, motif_length_start, motif_length_end, cpu_count):
    args = [path_to_java, '-cp', path_to_chipmunk,
                   'ru.autosome.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', '1.0',
                   's:{}'.format(fasta_path),
                  '100', '10', '1', cpu_count, 'random']
    p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    out = p.communicate()
    with open(path_out, 'wb') as file:
        file.write(out[0])
    pass


def get_chipmunk_model(models_path, fasta_path, path_to_python_tools, path_to_java, path_to_chipmunk, motif_length_start, motif_length_end, cpu_count):


    chipmunk_model_path = models_path + '/chipmunk_model'

    if not os.path.isdir(chipmunk_model_path):
        os.mkdir(chipmunk_model_path)


    ########################
    #FIND MODEL BY CHIPMUNK#
    ########################

    
    if not os.path.isfile(chipmunk_model_path + '/chipmunk_motif.txt'):
        #Create fastaWig for chipmunk
        print('Create pmw model by ChIPMunk')
        run_chipmunk_fasta(path_to_java, path_to_chipmunk,
        fasta_path,
        chipmunk_model_path + '/chipmunk_motif.txt',
        motif_length_start, motif_length_end, cpu_count)
    else:
        print('File {0} already exists (model exists)'.format(chipmunk_model_path + '/chipmunk_motif.txt'))

    ###########################################################################
    #Parse results of chipmunk into files .meme, .pwm and .fasta (multi fasta)#
    ###########################################################################

    args = ['python3', path_to_python_tools + 'parse_chipmunk_results.py',
           '-i', chipmunk_model_path + '/chipmunk_motif.txt',
           '-o', chipmunk_model_path,
           '-t', 'initial_chipmunk_motif']
    r = subprocess.call(args)


    ##############################################################################
    #Get oPWM from chipmunk results. OUTPUT: .meme, .pwm and .fasta (multi fasta)#
    ##############################################################################
    if not os.path.isfile(chipmunk_model_path + '/optimazed_chipmunk_motif.meme'):
        args = ['python3', path_to_python_tools + 'make_oPWM.py',
                '-c', chipmunk_model_path + '/chipmunk_motif.txt',
                '-f', fasta_path,
                '-n', '5000',
                '-P', cpu_count,
                '-o', chipmunk_model_path,
                '-t', 'optimazed_chipmunk_motif']
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(chipmunk_model_path + '/optimazed_chipmunk_motif.meme'))
    pass


def get_bamm_model(models_path, fasta_train, meme_model, model_order):
        #Get BaMM motif
    bamm_model_path = models_path + '/bamm'
    if not os.path.isdir(bamm_model_path):
        os.mkdir(bamm_model_path)

    args = ['BaMMmotif', bamm_model_path,
            fasta_train,
           '--PWMFile', meme_model,
            '--basename', 'bamm',
           '--EM',
           '--Order', model_order,
           '--order', model_order,
           '--scoreSeqset',
           '--saveLogOdds']
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


def bootstrap_inmode(path_to_python_tools, path_to_java, out_path, sites, path_to_inmode, tmp_dir):
    args = ['julia', path_to_python_tools + '/bootstrap_for_inmode.jl',
            out_path,
            sites,
            path_to_inmode, '-s', '2000000', '-j', path_to_java, '-t', tmp_dir]
    r = subprocess.call(args)
    pass


def bootstrap_bamm(path_to_python_tools, out_path, sites):
    args = ['julia', path_to_python_tools + '/bootstrap_for_bamm.jl',
            out_path,
            sites, '-s', '2000000']
    r = subprocess.call(args)
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


def scan_best_by_pwm(path_to_python_tools, output, input_model, fasta_in):
    args = ['pypy', path_to_python_tools + 'scan_best_by_pwm.py',
            '-f', fasta_in,
            '-m', input_model,
            '-o', output]
    r = subprocess.call(args)
    pass


def scan_best_by_bamm(path_to_python_tools, output, input_bamm_model, bg_model, fasta_in):
    args = ['pypy', path_to_python_tools + '/scan_best_by_bamm.py',
            '-f', fasta_in,
            '-m', input_bamm_model,
            '-b', bg_model,
            '-o', output]
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


def make_model(path_to_python_tools, path_in, dir_out, tag):
    args = ['python3', path_to_python_tools + '/make_model.py',
            '-i', path_in,
            '-o', dir_out,
            '-t', tag,
            '-M']
    r = subprocess.call(args)
    pass


def get_motif_length(models):
    with open(models + '/chipmunk_model/optimazed_chipmunk_motif.fasta', 'r') as file:
        for i in file:
            if i.startswith('>'):
                continue
            else:
                motif_length = len(i.strip())
                break
    file.close()
    return(motif_length)



def pipeline(bed_path, train_sample_size, test_sample_size,
                      path_to_out, path_to_python_tools, path_to_java, path_to_inmode, path_to_chipmunk,
                      path_to_promoters, path_to_genome, path_to_hocomoco, cpu_count):

    main_out = path_to_out
    model_order = str(2)
    cpu_count = str(cpu_count)
    motif_length_start = str(8)
    motif_length_end = str(14)


    if not path_to_python_tools[-1] == '/':
        path_to_python_tools += '/'


    if not os.path.isdir(main_out):
        os.mkdir(main_out)


    models = main_out + '/models'
    fasta = main_out + '/fasta'
    bed = main_out + '/bed'
    scan = main_out + '/scan'
    scan_best = main_out + '/scan-best'
    compare_sites = main_out + '/compare'
    tomtom = main_out + '/tomtom'
    #tag = os.path.basename(bed_path).split('.bed')[0]


    ########################
    #      CREATE DIRS     #
    ########################

    if not os.path.isdir(models):
        os.mkdir(models)
    if not os.path.isdir(fasta):
        os.mkdir(fasta)
    if not os.path.isdir(bed):
        os.mkdir(bed)
    if not os.path.isdir(scan):
        os.mkdir(scan)
    if not os.path.isdir(scan_best):
        os.mkdir(scan_best)
    if not os.path.isdir(compare_sites):
        os.mkdir(compare_sites)
    if not os.path.isdir(tomtom):
        os.mkdir(tomtom)



    # PREPARE BED AND FASTA FILES #
    prepare_data(path_to_python_tools, path_to_genome, bed_path, bed, fasta, train_sample_size, test_sample_size)

    fasta_train = fasta + '/train_sample.fa'
    fasta_test = fasta + '/test_sample.fa'
    bed_test = bed + '/test_sample.bed'
    bed_train = bed + '/train_sample.bed'

    # CALCULATE CHIPMUNK MODEL #
    get_chipmunk_model(models, fasta_train, path_to_python_tools, 
        path_to_java, path_to_chipmunk,
        motif_length_start, motif_length_end,
        cpu_count)

    # GET MOTIF LENGTH #
    motif_length = get_motif_length(models)

    # CALCULATE INMODE MODEL WITH EM ALG #
    get_inmode_model(models, fasta_train, path_to_java, path_to_inmode, motif_length,
                 model_order)

    # CALCULATE BAMM MODEL WITH EM ALG #
    get_bamm_model(models, fasta_train, meme_model, model_order)

    ###################
    #    BOOTSTRAP    #
    ###################

    # if not os.path.isfile(bootstrap + "/pwm.tsv"):
    #     print('RUNNIN BOOTSTRAP FOR PWM')
    #     bootstrap_pwm(path_to_python_tools, bootstrap + "/pwm.tsv",
    #     motifs + '/' + tag + '_OPTIMAL_MOTIF.fasta')
    # else:
    #     print('Bootstrap for pwm already calculated')

    # if not os.path.isfile(bootstrap + "/bamm.tsv"):
    #     print('RUNNIN BOOTSTRAP FOR BAMM')
    #     bootstrap_bamm(path_to_python_tools, bootstrap + "/bamm.tsv", motifs + '/' + tag + '_motif_1.logOddsZoops')
    # else:
    #     print('Bootstrap for bamm already calculated')
        
    # if not os.path.isfile(bootstrap + "/inmode.tsv"):
    #     print('RUNNIN BOOTSTRAP FOR INMODE')
    #     bootstrap_inmode(path_to_python_tools, path_to_java, bootstrap + "/inmode.tsv", glob.glob(motifs + '/Learned_DeNovo*/Binding_sites_of_DeNovo*motif.txt')[0], path_to_inmode, motifs + '/tmp')
    # else:
    #     print('Bootstrap for inmode already calculated')

    ##########################
    #  CALCULATE THRESHOLDS  #
    ##########################

    # if not os.path.isfile(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt'):
    #     print('Calculate threshold for BAMM based on promoters and FPR')
    #     args = ['pypy', path_to_python_tools + '/get_threshold_for_bamm.py',
    #             path_to_promoters,
    #             motifs + '/' + tag + '_motif_1.ihbcp',
    #             motifs + '/' + tag + '.hbcp',
    #             motifs + '/' + tag + '_BAMM_THRESHOLDS.txt']
    #     r = subprocess.call(args)
    # else:
    #     print('Thresholds for bamm already calculated')

    # if not os.path.isfile(motifs + '/' + tag + '_PWM_THRESHOLDS.txt'):
    #     print('Calculate threshold for PWM based on promoters and FPR')
    #     args = ['pypy', path_to_python_tools + '/get_threshold_for_pwm.py',
    #             path_to_promoters,
    #             motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
    #             motifs + '/' + tag + '_PWM_THRESHOLDS.txt']
    #     r = subprocess.call(args)
    # else:
    #     print('Thresholds for pwm already calculated')

    # if not os.path.isfile(motifs + '/' + tag + '_INMODE_THRESHOLDS.txt'):
    #     args = ['python3', path_to_python_tools + '/get_threshold_for_inmode.py',
    #             path_to_promoters,
    #             glob.glob(motifs + '/Learned_DeNovo*/*.xml')[0],
    #             path_to_inmode,
    #             str(motif_length),
    #             motifs + '/' + tag + '_INMODE_THRESHOLDS.txt',
    #             '-j', path_to_java,
    #             '-t', motifs + '/tmp']
    #     r = subprocess.call(args)
    # else:
    #     print('Thresholds for inmode already calculated')

    ############################################
    #  RUN LOOP THRUE SEVERAL FPR (THRESHOLD)  #
    ############################################

    # fprs = [5*10**(-4), 1.90*10**(-4), 5.24*10**(-5)]
    # for fpr_for_thr in fprs:

    #     #############################################
    #     #CALCULATE THRESHOLDS FOR BAMM MODEL AND SCAN#
    #     #############################################

    #     if not os.path.isfile(scan + '/' + tag + '_BAMM_' + str(testing_sample_size) +'_' + '{:.2e}'.format(fpr_for_thr) + '.bed'):
            
    #         thr_bamm = get_threshold(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt', fpr_for_thr)
    #         print('BAMM = ',thr_bamm)

    #         #Scan peaks by BAMM with thr_bamm
    #         print('Scan peaks by BAMM with thr_pwm ({0})'.format(tag))
    #         args = ['pypy', path_to_python_tools + 'scan_by_bamm.py',
    #                 '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
    #                 '-m', motifs + '/' + tag + '_motif_1.ihbcp',
    #                 '-b', motifs + '/' + tag + '.hbcp',
    #                 '-t', str(thr_bamm),
    #                 '-o', scan + '/' + tag + '_BAMM_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed']
    #         r = subprocess.call(args)
    #     else:
    #         print(tag + '_BAMM_' + str(testing_sample_size) +'_' + '{:.2e}'.format(fpr_for_thr) + '.bed', '- EXISTS')

    #     #############################################
    #     #CALCULATE THRESHOLDS FOR PWM MODEL AND SCAN#
    #     #############################################


    #     if not os.path.isfile(scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + '{:.2e}'.format(fpr_for_thr) + '.bed'):
            
    #         thr_pwm = get_threshold(motifs + '/' + tag + '_PWM_THRESHOLDS.txt', fpr_for_thr)
    #         print('PWM = ',thr_pwm)

    #         #Scan peaks by PWM with thr_pwm
    #         print('Scan peaks by PWM with thr_pwm ({0})'.format(tag))
    #         args = ['pypy', path_to_python_tools + 'scan_by_pwm.py',
    #                 '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
    #                 '-m', motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
    #                 '-t', str(thr_pwm),
    #                 '-o', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + '{:.2e}'.format(fpr_for_thr) + '.bed']
    #         r = subprocess.call(args)

    #     else:
    #         print(tag + '_PWM_' + str(testing_sample_size) +'_' + '{:.2e}'.format(fpr_for_thr) + '.bed', '- EXISTS')


    #     ################################################
    #     #CALCULATE THRESHOLDS FOR INMODE MODEL AND SCAN#
    #     ################################################

    #     if not os.path.isfile(scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed'):
    #         print('Scan by inmode model')
    #         inmode_scan(path_to_java, path_to_inmode,
    #                        input_data=fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
    #                        input_model=glob.glob(motifs + '/Learned_DeNovo*/*.xml')[0],
    #                        backgroud_path=path_to_promoters,
    #                        fpr_for_thr=fpr_for_thr,
    #                        outdir=scan + '/tmp')

    #         args = ['python3', path_to_python_tools + 'parse_inmode_scan.py',
    #                 '-if', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
    #                 '-bed', glob.glob(scan + '/tmp' + '/*.BED')[0],
    #                 '-o', scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed']
    #         r = subprocess.call(args)
    #         os.system("rm -r {}".format(scan + '/tmp'))

    #     else:
    #         print(tag + '_INMODE_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed', '- EXISTS')

    #     ##############################
    #     #COMPARE SITES OF DIFF MODELS#
    #     ##############################

    #     # if not os.path.isfile(compare_sites + '/' + tag + '_' + '{:.2e}'.format(fpr_for_thr) + '_COUNT.tsv'):
    #     # print('Compare sites ({0})'.format(tag))
    #     args = ['pypy', path_to_python_tools + 'compare_sites3-pypy.py',
    #             '-p', bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
    #             '-first', scan + '/' + tag + '_PWM_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed',
    #             '-second', scan + '/' + tag + '_BAMM_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed',
    #             '-third', scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed',
    #             '-t', tag + '_' + '{:.2e}'.format(fpr_for_thr),
    #             '-o', compare_sites,
    #             '-fname', fname,
    #             '-sname', sname,
    #             '-tname', tname]
    #     r = subprocess.call(args)
    #     # else:
    #     #     print('Sites already compared')


    #     ###################
    #     #  EXTRACT SITES  #
    #     ###################

    #     print('EXTRACT SITES ({0})'.format(tag))

    #     if not os.path.isfile(scan + '/' + 'pwm.sites.{:.2e}.txt'.format(fpr_for_thr)):
    #         args = ['python3', path_to_python_tools + 'extract_sites.py',
    #         '-p', scan + '/' + tag + '_PWM_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed',
    #         '-o', scan + '/' + 'pwm.sites.{:.2e}.txt'.format(fpr_for_thr)]
    #         r = subprocess.call(args)
    #     else:
    #         print('pwm.sites.{:.2e} already extracted'.format(fpr_for_thr))

    #     if not os.path.isfile(scan + '/' + 'bamm.sites.{:.2e}.txt'.format(fpr_for_thr)):
    #         args = ['python3', path_to_python_tools + 'extract_sites.py',
    #         '-p', scan + '/' + tag + '_BAMM_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed',
    #         '-o', scan + '/' + 'bamm.sites.{:.2e}.txt'.format(fpr_for_thr)]
    #         r = subprocess.call(args)
    #     else:
    #         print('bamm.sites.{:.2e} already extracted'.format(fpr_for_thr))

    #     if not os.path.isfile(scan + '/' + 'inmode.sites.{:.2e}.txt'.format(fpr_for_thr)):
    #         args = ['python3', path_to_python_tools + 'extract_sites.py',
    #         '-p', scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + '{:.2e}'.format(fpr_for_thr) + '.bed',
    #         '-o', scan + '/' + 'inmode.sites.{:.2e}.txt'.format(fpr_for_thr)]
    #         r = subprocess.call(args)
    #     else:
    #         print('inmode.sites.{:.2e} already extracted'.format(fpr_for_thr))

    #     ############
    #     # END LOOP #
    #     ############


    ########################################
    # CREATE MODELS FROM SITES AND COMPARE #
    ########################################

    # print('Run tomtom')

    # make_model(path_to_python_tools, scan + '/' + 'inmode.sites.{:.2e}.txt'.format(fpr_for_thr),
    #     tomtom, 'inmode')
    # make_model(path_to_python_tools, scan + '/' + 'pwm.sites.{:.2e}.txt'.format(fpr_for_thr),
    #     tomtom, 'pwm')
    # make_model(path_to_python_tools, scan + '/' + 'bamm.sites.{:.2e}.txt'.format(fpr_for_thr),
    #     tomtom, 'bamm')

    # run_tomtom(tomtom + '/pwm.meme', tomtom + '/bamm.meme', tomtom + '/pwm.bamm')
    # run_tomtom(tomtom + '/pwm.meme', tomtom + '/inmode.meme', tomtom + '/pwm.inmode')
    # run_tomtom(tomtom + '/inmode.meme', tomtom + '/bamm.meme', tomtom + '/inmode.bamm')

    # run_tomtom(path_to_hocomoco, tomtom + '/pwm.meme', tomtom + '/pwm.hocomoco')
    # run_tomtom(path_to_hocomoco, tomtom + '/inmode.meme', tomtom + '/inmode.hocomoco')
    # run_tomtom(path_to_hocomoco, tomtom + '/bamm.meme', tomtom + '/bamm.hocomoco')

    # ###########
    # #SCAN BEST#
    # ###########

    # if not os.path.isfile(scan_best + '/inmode.scores.txt'):
    #     print("Scan best inmode")
    #     scan_best_by_inmode(path_to_python_tools, scan_best + '/inmode.scores.txt',
    #     glob.glob(motifs + '/Learned_DeNovo*/*.xml')[0],
    #     fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
    #     path_to_inmode, path_to_java)
    # else:
    #     print('best scores of inmode already exists')

    # if not os.path.isfile(scan_best + '/pwm.scores.txt'):
    #     print("Scan best pwm")
    #     scan_best_by_pwm(path_to_python_tools, scan_best + '/pwm.scores.txt',
    #     motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
    #     fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa')
    # else:
    #     print('best scores of pwm already exists')

    # if not os.path.isfile(scan_best + '/bamm.scores.txt'):
    #     print("Scan best bamm")
    #     scan_best_by_bamm(path_to_python_tools, scan_best + '/bamm.scores.txt',
    #     motifs + '/' + tag + '_motif_1.ihbcp',
    #     motifs + '/' + tag + '.hbcp',
    #     fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa')
    # else:
    #     print('best scores of bamm already exists')


    # ###########
    # #PLOT BEST#
    # ###########

    # length = bed + '/' + tag + '_' + str(testing_sample_size) + '.length.txt'
    # fpr_for_thr = 5.24*10**(-5)

    # pwm_scores = scan_best + '/pwm.scores.txt'
    # bamm_scores = scan_best + '/bamm.scores.txt'
    # inmode_scores = scan_best + '/inmode.scores.txt'
    # results_corr = scan_best + '/corr.results.txt'

    # thr_bamm = str(get_threshold(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt', fpr_for_thr))
    # thr_pwm = str(get_threshold(motifs + '/' + tag + '_PWM_THRESHOLDS.txt', fpr_for_thr))
    # thr_inmode = str(math.log2(get_threshold(motifs + '/' + tag + '_INMODE_THRESHOLDS.txt', fpr_for_thr)))

    # if not os.path.isfile(scan_best + '/pwm-inmode-scores.pdf'):
    #     plot_best_score(path_to_python_tools, inmode_scores, pwm_scores,
    #         thr_inmode, thr_pwm, length, scan_best + '/pwm-inmode-scores.pdf', 'inmode scores', 'pwm scores')
    # else:
    #     print('pwm-inmode-scores.pdf exists')

    # if not os.path.isfile(scan_best + '/pwm-bamm-scores.pdf'):
    #     plot_best_score(path_to_python_tools, pwm_scores, bamm_scores,
    #         thr_pwm, thr_bamm, length, scan_best + '/pwm-bamm-scores.pdf', 'pwm scores', 'bamm scores')
    # else:
    #     print('pwm-bamm-scores.pdf exists')

    # if not os.path.isfile(scan_best + '/bamm-inmode-scores.pdf'):
    #     plot_best_score(path_to_python_tools, bamm_scores, inmode_scores,
    #         thr_bamm, thr_inmode, length, scan_best + '/bamm-inmode-scores.pdf', 'bamm scores', 'inmode scores')
    # else:    
    #     print('bamm-inmode-scores.pdf exists')


    # ############
    # #MONTECARLO#
    # ############

    # pwm_scores = scan_best + '/pwm.scores.txt'
    # bamm_scores = scan_best + '/bamm.scores.txt'
    # inmode_scores = scan_best + '/inmode.scores.txt'
    # results_montecarlo = scan_best + '/montecarlo.results.txt'
    # results_binome = scan_best + '/binome.results.txt'

    # print('Run montecarlo')

    # fprs = [5*10**(-4), 1.90*10**(-4), 5.24*10**(-5)]
    # for fpr in fprs:
    #     thr1 = str(get_threshold(motifs + '/' + tag + '_PWM_THRESHOLDS.txt', fpr))
    #     thr2 = str(get_threshold(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt', fpr))
    #     scores1 = pwm_scores
    #     scores2 = bamm_scores
    #     name1, name2 = 'PWM', 'BAMM'
    #     montecarlo(path_to_python_tools, scores1, scores2, thr1, thr2, length, results_montecarlo, name1, name2, fpr, fpr)


    #     thr1 = str(get_threshold(motifs + '/' + tag + '_PWM_THRESHOLDS.txt', fpr))
    #     thr2 = str(math.log2(get_threshold(motifs + '/' + tag + '_INMODE_THRESHOLDS.txt', fpr)))
    #     scores1 = pwm_scores
    #     scores2 = inmode_scores
    #     name1, name2 = 'PWM', 'INMODE'
    #     montecarlo(path_to_python_tools, scores1, scores2, thr1, thr2, length, results_montecarlo, name1, name2, fpr, fpr)


    #     thr1 = str(get_threshold(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt', fpr))
    #     thr2 = str(math.log2(get_threshold(motifs + '/' + tag + '_INMODE_THRESHOLDS.txt', fpr)))
    #     scores1 = bamm_scores
    #     scores2 = inmode_scores
    #     name1, name2 = 'BAMM', 'INMODE'
    #     montecarlo(path_to_python_tools, scores1, scores2, thr1, thr2, length, results_montecarlo, name1, name2, fpr, fpr)

    # print('Finish!')


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store', help='path to BED file')
    parser.add_argument('promoters', action='store', help='path to promoters fasta file')
    parser.add_argument('genome', action='store', help='path to genome fasta file')
    parser.add_argument('output', action='store', help='output dir')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=False, default=1000, help='size of training sample, by default size is equal to 1000')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=False, default=10000, help='size of testing sample, by default size is equal to 10000')
    parser.add_argument('-p', '--python', action='store', dest='python_tools',
                        required=True, help='dir with python tools')
    parser.add_argument('-I', '--inmode', action='store', dest='inmode',
                        required=True, help='path to inmode')
    parser.add_argument('-J', '--java', action='store', dest='java',
                    required=False, default="java", help='path to Java')
    parser.add_argument('-c', '--chipmunk', action='store', dest='chipmunk',
                        required=True, help='path to chipmunk')
    parser.add_argument('-C', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=2, help='Number of processes to use, default: 2')
    parser.add_argument('-H', '--hocomoco', action='store', dest='path_to_hocomoco',
                        required=False, help='path to HOCOMOCO database in meme format for TOMTOM')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():

    args = parse_args()

    bed_path = args.bed
    path_to_out = args.output
    train_sample_size = args.train_size
    test_sample_size = args.test_size


    path_to_python_tools = args.python_tools
    path_to_java = args.java
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    path_to_promoters = args.promoters
    path_to_genome = args.genome
    path_to_hocomoco = args.path_to_hocomoco
    cpu_count = args.cpu_count


    pipeline(bed_path, train_sample_size, test_sample_size,
                          path_to_out, path_to_python_tools, path_to_java, path_to_inmode, path_to_chipmunk,
                          path_to_promoters, path_to_genome, path_to_hocomoco, cpu_count)

if __name__ == '__main__':
    main()