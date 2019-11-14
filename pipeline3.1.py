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


# def get_threshold(path, fpr_for_thr):
#     df = pd.read_csv(path, sep='\t')
#     df = df.sort_values(by=["FPR"])
#     threshold = df['Scores'].iloc[np.searchsorted(df["FPR"], fpr_for_thr, side='right') - 1]
#     return(threshold)


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
            sites]
    r = subprocess.call(args)
    pass

def bootstrap_inmode(path_to_python_tools, out_path, sites, path_to_inmode):
    args = ['julia', path_to_python_tools + '/bootstrap_for_inmode.jl',
            out_path,
            sites,
            path_to_inmode]
    r = subprocess.call(args)
    pass

def bootstrap_bamm(path_to_python_tools, out_path, sites):
    args = ['julia', path_to_python_tools + '/bootstrap_for_bamm.jl',
            out_path,
            sites]
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


def run_chipmunk_fasta_wig(path_to_java, path_to_chipmunk, fasta_path, path_out, motif_length_start, motif_length_end, try_size, cpu_count, zoops):
    args = [path_to_java, '-cp', path_to_chipmunk,
                   'ru.autosome.ChIPMunk', str(motif_length_start), str(motif_length_end), 'yes', zoops,
                   'p:' + fasta_path,
                  try_size, '10', '1', cpu_count, 'random']
    p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    out = p.communicate()
    with open(path_out, 'wb') as file:
        file.write(out[0])
    pass


def create_fasta_wig(path_to_python_tools, bed_in, bw_in, fasta_in, dir_out, tag):
    args = ['python3', path_to_python_tools + '/prepare_fasta_for_chipmunk.py',
           '-b', bed_in,
           '-w', bw_in,
           '-f', fasta_in,
           '-o', dir_out,
           '-t', tag]
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


# def scan_best_by_inmode(path_to_python_tools, output, input_model, fasta_in, path_to_inmode):
#     args = ['julia', path_to_python_tools + '/scan_best_by_inmode.jl',
#            output,
#            input_model,
#            fasta_in,
#            path_to_inmode]
#     r = subprocess.call(args)
#     pass

def scan_best_by_inmode(path_to_python_tools, output, input_model, fasta_in, path_to_inmode, path_to_java):
    args = ['python3', path_to_python_tools + '/scan_best_by_inmode.py',
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


def montecarlo(scores1, scores2, thr1, thr2, length, results, name1, name2, fpr1, fpr2):
    with open(results, 'a') as file:
        file.write('{0}:thr={1},fpr={2};{3}:thr={4},fpr={5}\n'.format(name1, thr1, fpr1, name2, thr2, fpr2))
    file.close()
    args = ['monteCarlo', '{}'.format(scores1), '{}'.format(scores2), '{}'.format(thr1), '{}'.format(thr2), '{}'.format(length), '{}'.format(results),]
    r = subprocess.call(args)
    pass



def pipeline_inmode_bamm(bed_path, bigwig_path, training_sample_size, testing_sample_size,
                      fpr_for_thr, path_to_out, path_to_python_tools, path_to_java, path_to_inmode, path_to_imd, path_to_chipmunk,
                      path_to_promoters, path_to_genome, path_to_tss, cpu_count,
                      zoops, try_size, model_order):

    main_out = path_to_out + '/' + os.path.basename(bed_path).split('.')[0]
    zoops = str(zoops)
    model_order = str(model_order)
    try_size=str(try_size)
    cpu_count = str(cpu_count)
    motif_length_start = str(8)
    motif_length_end = str(20)
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



    ########################
    #     GET TOP PEAKS    #
    ########################

    if not os.path.isfile(bed + '/' + tag + '_' + str(training_sample_size) + '.bed'):
        #Get top training_sample_size bed peaks
        print('Get top {0} bed peaks for {1}'.format(training_sample_size, tag))
        bed_out = bed + '/'
        get_top_peaks(path_to_python_tools, bed_path, bed_out, training_sample_size, tag + '_' + str(training_sample_size))


    else:
        print('File {0} already exists'.format(tag + '_' + str(training_sample_size) + '.bed'))

    if not os.path.isfile(bed + '/' + tag + '_' + str(testing_sample_size) + '.bed'):
        #Get top testing_sample_size bed peaks
        print('Get top {1} bed peaks for {0}'.format(tag, testing_sample_size))
        bed_out = bed + '/'
        get_top_peaks(path_to_python_tools, bed_path, bed_out, testing_sample_size, tag + '_' + str(testing_sample_size))
    else:
        print('File {0} already exists'.format(tag + '_' + str(testing_sample_size) + '.bed'))

    ########################
    #     BED TO FASTA     #
    ########################

    if not os.path.isfile(fasta + '/' + tag + '_' + str(training_sample_size) +'.fa'):
        #Bed peaks to fasta
        print('Bed peaks to fasta for {0}'.format(tag))
        bed_to_fasta(path_to_genome,
            bed + '/' + tag + '_' + str(training_sample_size) +'.bed',
            fasta + '/' + tag + '_' + str(training_sample_size) +'.fa')

    else:
        print('File {0} already exists'.format(tag + '_' + str(training_sample_size) +'.fa'))

    if not os.path.isfile(fasta + '/' + tag + '_' + str(testing_sample_size) +'.fa'):
        bed_to_fasta(path_to_genome,
            bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
            fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa')
    else:
        print('File {0} already exists'.format(tag + '_' + str(testing_sample_size) +'.fa'))


    ########################
    #FIND MODEL BY CHIPMUNK#
    ########################
    if not os.path.isfile(chipmunk + '/CHIPMUNK_MOTIF.txt'):
        if bigwig_path != "PASS":
            #Create fastaWig for chipmunk
            create_fasta_wig(path_to_python_tools, bed + '/' + tag + '_' + str(training_sample_size) +'.bed',
            bigwig_path, fasta + '/' + tag + '_' + str(training_sample_size) +'.fa', fasta, tag + '_' + str(training_sample_size))

            # Run chipmunk
            run_chipmunk_fasta_wig(path_to_java, path_to_chipmunk,
            fasta + '/' + tag + '_'+ str(training_sample_size) + '.faWig',
            chipmunk + '/CHIPMUNK_MOTIF.txt',
            motif_length_start, motif_length_end,
            try_size, cpu_count, zoops)

        else:
            print('chipmunk find motifs for {0}'.format(tag))
            run_chipmunk_fasta(path_to_java, path_to_chipmunk,
            fasta + '/' + tag + '_'+ str(training_sample_size) + '.fa',
            chipmunk + '/CHIPMUNK_MOTIF.txt',
            motif_length_start, motif_length_end,
            try_size, cpu_count, zoops)
    else:
        print('File {0} already exists'.format(chipmunk + '/chipmunk_MOTIF.txt'))

    ###########################################################################
    #Parse results of chipmunk into files .meme, .pwm and .fasta (multi fasta)#
    ###########################################################################

    args = ['python3', path_to_python_tools + 'parse_chipmunk_results.py',
           '-i', chipmunk + '/CHIPMUNK_MOTIF.txt',
           '-o', chipmunk,
           '-t', tag + '_' + 'CHIPMUNK_MOTIF']
    r = subprocess.call(args)


    ##############################################################################
    #Get oPWM from chipmunk results. OUTPUT: .meme, .pwm and .fasta (multi fasta)#
    ##############################################################################
    if not os.path.isfile(motifs + '/' + tag + '_' + 'OPTIMAL_MOTIF.meme'):
        args = ['python3', path_to_python_tools + 'make_oPWM.py',
                '-c', chipmunk + '/CHIPMUNK_MOTIF.txt',
                '-f', fasta + '/' + tag + '_'+ str(training_sample_size) + '.fa',
                '-n', '5000',
                '-P', cpu_count,
                '-o', motifs,
                '-t', tag + '_' + 'OPTIMAL_MOTIF']
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(motifs + '/PEAKS039334_OPTIMAL_MOTIF.meme'))


    ##################
    #GET MOTIF LENGTH#
    ##################

    with open(motifs + '/' + tag + '_' + 'OPTIMAL_MOTIF.fasta', 'r') as file:
        for i in file:
            if i.startswith('>'):
                continue
            else:
                motif_length = len(i.strip())
                break
    file.close()


    ####################################
    #CALCULATE INMODE MODEL WITH EM ALG#
    ####################################
    if glob.glob(motifs + '/Learned_DeNovo*') == []:
        print('Calculate inmode model')
        inmode_denovo(path_to_java, path_to_inmode,
                         fasta_path=fasta + '/' + tag + '_'+ str(training_sample_size) + '.fa',
                         motif_length=motif_length,
                         model_order=model_order,
                         outdir=motifs)
    else:
        print('inmode model already exists')

    ##################################
    #CALCULATE BAMM MODEL WITH EM ALG#
    ##################################

    #Get BaMM motif
    if not os.path.isfile(motifs + '/' + tag + '_' + 'motif_1.ihbcp'):
        print('Get Bamm motifs for {0}'.format(tag))
        args = ['BaMMmotif', motifs,
                fasta + '/' + tag + '_' + str(training_sample_size) + '.fa',
               '--PWMFile', motifs + '/' + tag + '_OPTIMAL_MOTIF.meme',
                '--basename', tag,
               '--EM',
               '--Order', model_order,
               '--order', model_order,
               '--scoreSeqset',
               '--saveLogOdds']
        r = subprocess.call(args)
    else:
        print('BaMM model already exists')


    ###################
    #    BOOTSTRAP    #
    ###################

    # print('RUNNIN BOOTSTRAP FOR PWM')
    # bootstrap_pwm(path_to_python_tools, bootstrap + "/pwm.tsv",
    # motifs + '/' + tag + '_OPTIMAL_MOTIF.fasta')
    
    # print('RUNNIN BOOTSTRAP FOR BAMM')
    # bootstrap_bamm(path_to_python_tools, bootstrap + "/bamm.tsv", motifs + '/' + tag + '_motif_1.logOddsZoops')
    
    # print('RUNNIN BOOTSTRAP FOR INMODE')
    # bootstrap_inmode(path_to_python_tools, bootstrap + "/inmode.tsv", glob.glob(motifs + '/Learned_DeNovo*/Binding_sites_of_DeNovo*motif.txt')[0], path_to_inmode)


    ##########################
    #  CALCULATE THRESHOLDS  #
    ##########################


    # print('Calculate threshold for BAMM based on promoters and FPR')
    # args = ['pypy', path_to_python_tools + '/get_threshold_for_bamm.py',
    #         path_to_promoters,
    #         motifs + '/' + tag + '_motif_1.ihbcp',
    #         motifs + '/' + tag + '.hbcp',
    #         motifs + '/' + tag + '_BAMM_THRESHOLDS.txt']
    # r = subprocess.call(args)


    # print('Calculate threshold for PWM based on promoters and FPR')
    # args = ['pypy', path_to_python_tools + '/get_threshold_for_pwm.py',
    #         path_to_promoters,
    #         motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
    #         motifs + '/' + tag + '_PWM_THRESHOLDS.txt']
    # r = subprocess.call(args)


    # print('Calculate threshold for Inmode based on promoters and FPR')
    # args = ['python3', path_to_python_tools + '/get_threshold_for_inmode.py',
    #         path_to_promoters,
    #         glob.glob(motifs + '/Learned_DeNovo*/*.xml')[0],
    #         path_to_inmode,
    #         str(motif_length),
    #         motifs + '/' + tag + '_INMODE_THRESHOLDS.txt',
    #         '-j', path_to_java]
    # r = subprocess.call(args)

    ############################################
    #  RUN LOOP THRUE SEVERAL FPR (THRESHOLD)  #
    ############################################

    fprs = [5*10**(-4), 1.90*10**(-4), 5.24*10**(-5)]
    for fpr_for_thr in fprs:
        break

        #############################################
        #CALCULATE THRESHOLDS FOR BAMM MODEL AND SCAN#
        #############################################

        if not os.path.isfile(scan + '/' + tag + '_BAMM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed'):
            
            thr_bamm = get_threshold(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt', fpr_for_thr)
            print('BAMM = ',thr_bamm)

            #Scan peaks by BAMM with thr_bamm
            print('Scan peaks by BAMM with thr_pwm ({0})'.format(tag))
            args = ['python3', path_to_python_tools + 'scan_by_bamm.py',
                    '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                    '-m', motifs + '/' + tag + '_motif_1.ihbcp',
                    '-b', motifs + '/' + tag + '.hbcp',
                    '-t', str(thr_bamm),
                    '-o', scan + '/' + tag + '_BAMM_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
                    '-P', cpu_count]
            r = subprocess.call(args)
        else:
            print(tag + '_BAMM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed', '- EXISTS')

        #############################################
        #CALCULATE THRESHOLDS FOR PWM MODEL AND SCAN#
        #############################################


        if not os.path.isfile(scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed'):
            
            thr_pwm = get_threshold(motifs + '/' + tag + '_PWM_THRESHOLDS.txt', fpr_for_thr)
            print('PWM = ',thr_pwm)

            #Scan peaks by PWM with thr_pwm
            print('Scan peaks by PWM with thr_pwm ({0})'.format(tag))
            args = ['python3', path_to_python_tools + 'scan_by_pwm.py',
                    '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                    '-m', motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
                    '-t', str(thr_pwm),
                    '-o', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                    '-P', cpu_count]
            r = subprocess.call(args)

        else:
            print(tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed', '- EXISTS')


        ################################################
        #CALCULATE THRESHOLDS FOR INMODE MODEL AND SCAN#
        ################################################

        if not os.path.isfile(scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed'):
            print('Scan by inmode model')
            inmode_scan(path_to_java, path_to_inmode,
                           input_data=fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                           input_model=glob.glob(motifs + '/Learned_DeNovo*/*.xml')[0],
                           backgroud_path=path_to_promoters,
                           fpr_for_thr=fpr_for_thr,
                           outdir=scan + '/tmp')

            args = ['python3', path_to_python_tools + 'parse_inmode_scan.py',
                    '-if', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                    '-bed', glob.glob(scan + '/tmp' + '/*.BED')[0],
                    '-o', scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed']
            r = subprocess.call(args)
            os.system("rm -r {}".format(scan + '/tmp'))

        else:
            print(tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed', '- EXISTS')

        ##############################
        #COMPARE SITES OF DIFF MODELS#
        ##############################

        print('Compare sites ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'compare_sites3.py',
                '-p', bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
                '-first', scan + '/' + tag + '_PWM_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
                '-second', scan + '/' + tag + '_BAMM_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
                '-third', scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
                '-t', tag + '_' + str(fpr_for_thr),
                '-o', compare_sites,
                '-fname', fname,
                '-sname', sname,
                '-tname', tname]
        r = subprocess.call(args)


        ####################
        #WORK WITH GENE IDS#
        ####################

        ## GET GEN IDS ###
        print('GET GEN IDS ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'peaks_intersection_with_bed.py',
        path_to_tss,
        scan + '/' + tag + '_PWM_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
        gene_ids + '/' + 'pwm.ids.{}.txt'.format(fpr_for_thr)]
        r = subprocess.call(args)

        print('GET GEN IDS ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'peaks_intersection_with_bed.py',
        path_to_tss,
        scan + '/' + tag + '_BAMM_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
        gene_ids + '/' + 'bamm.ids.{}.txt'.format(fpr_for_thr)]
        r = subprocess.call(args)

        print('GET GEN IDS ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'peaks_intersection_with_bed.py',
        path_to_tss,
        scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
        gene_ids + '/' + 'inmode.ids.{}.txt'.format(fpr_for_thr)]
        r = subprocess.call(args)

        ### COMPARE IDS ###

        print('Compare sites ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'compare_gene_ids.py',
                '-first', gene_ids + '/' + 'pwm.ids.{}.txt'.format(fpr_for_thr),
                '-second', gene_ids + '/' + 'bamm.ids.{}.txt'.format(fpr_for_thr),
                '-third', gene_ids + '/' + 'inmode.ids.{}.txt'.format(fpr_for_thr),
                '-o', gene_ids + '/' + 'compare.ids.{}.pdf'.format(fpr_for_thr),
                '-fname', fname,
                '-sname', sname,
                '-tname', tname]
        r = subprocess.call(args)


        ###################
        #  IMD FOR SITES  #
        ###################

        print('EXTRACT SITES ({0})'.format(tag))

        args = ['python3', path_to_python_tools + 'extract_sites.py',
        '-p', scan + '/' + tag + '_PWM_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
        '-o', scan + '/' + 'pwm.sites.{}.txt'.format(fpr_for_thr)]
        r = subprocess.call(args)

        args = ['python3', path_to_python_tools + 'extract_sites.py',
        '-p', scan + '/' + tag + '_BAMM_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
        '-o', scan + '/' + 'bamm.sites.{}.txt'.format(fpr_for_thr)]
        r = subprocess.call(args)

        args = ['python3', path_to_python_tools + 'extract_sites.py',
        '-p', scan + '/' + tag + '_INMODE_' + str(testing_sample_size) + '_' + str(fpr_for_thr) + '.bed',
        '-o', scan + '/' + 'inmode.sites.{}.txt'.format(fpr_for_thr)]
        r = subprocess.call(args)


        # args = ['/home/anton/Programs/jdk-9/bin/java', '--add-modules', 'java.xml.bind',
        # '-jar', path_to_imd, 'imd',
        # 'i=' + scan + '/' + 'pwm.sites.{}.txt'.format(fpr_for_thr),
        # 'outdir=' + scan + '/pwm_{}'.format(fpr_for_thr)]
        # r = subprocess.call(args)
        
        # args = ['/home/anton/Programs/jdk-9/bin/java', '--add-modules', 'java.xml.bind',
        # '-jar', path_to_imd, 'imd',
        # 'i=' + scan + '/' + 'bamm.sites.{}.txt'.format(fpr_for_thr),
        # 'outdir=' + scan + '/bamm_{}'.format(fpr_for_thr)]
        # r = subprocess.call(args)
        
        # args = ['/home/anton/Programs/jdk-9/bin/java', '--add-modules', 'java.xml.bind',
        # '-jar', path_to_imd, 'imd',
        # 'i=' + scan + '/' + 'inmode.sites.{}.txt'.format(fpr_for_thr),
        # 'outdir=' + scan + '/inmode_{}'.format(fpr_for_thr)]
        # r = subprocess.call(args)

        ############
        # END LOOP #
        ############

    ###########
    #SCAN BEST#
    ###########

    scan_best_by_inmode(path_to_python_tools, scan_best + '/inmode.scores.txt',
    glob.glob(motifs + '/Learned_DeNovo*/*.xml')[0],
    fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
    path_to_inmode, path_to_java)

    scan_best_by_pwm(path_to_python_tools, scan_best + '/pwm.scores.txt',
    motifs + '/' + tag + '_OPTIMAL_MOTIF.pwm',
    fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
    cpu_count)

    scan_best_by_bamm(path_to_python_tools, scan_best + '/bamm.scores.txt',
    motifs + '/' + tag + '_motif_1.ihbcp',
    motifs + '/' + tag + '.hbcp',
    fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
    cpu_count)


    ###########
    #PLOT BEST#
    ###########

    length = bed + '/' + tag + '_' + str(training_sample_size) + '.length.txt'
    fpr_for_thr = 1.90*10**(-4)

    thr_bamm = str(get_threshold(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt', fpr_for_thr))
    thr_pwm = str(get_threshold(motifs + '/' + tag + '_PWM_THRESHOLDS.txt', fpr_for_thr))
    thr_inmode = str(math.log2(get_threshold(motifs + '/' + tag + '_INMODE_THRESHOLDS.txt', fpr_for_thr)))

    plot_best_score(path_to_python_tools, scan_best + '/inmode.scores.txt', scan_best + '/pwm.scores.txt',
        thr_inmode, thr_pwm, length, scan_best + '/pwm-inmode-scores.pdf', 'inmode scores', 'pwm scores')


    plot_best_score(path_to_python_tools, scan_best + '/pwm.scores.txt', scan_best + '/bamm.scores.txt',
        thr_pwm, thr_bamm, length, scan_best + '/pwm-bamm-scores.pdf', 'pwm scores', 'bamm scores')


    plot_best_score(path_to_python_tools, scan_best + '/bamm.scores.txt', scan_best + '/inmode.scores.txt',
        thr_bamm, thr_inmode, length, scan_best + '/bamm-inmode-scores.pdf', 'bamm scores', 'inmode scores')


    ############
    #MONTECARLO#
    ############

    pwm_scores = scan_best + '/pwm.scores.txt'
    bamm_scores = scan_best + '/bamm.scores.txt'
    inmode_scores = scan_best + '/inmode.scores.txt'
    results = scan_best + '/montecarlo.results.txt'

    fprs = [5*10**(-4), 1.90*10**(-4), 5.24*10**(-5)]
    for fpr1 in fprs:
        for fpr2 in fprs:
            thr1 = str(get_threshold(motifs + '/' + tag + '_PWM_THRESHOLDS.txt', fpr1))
            thr2 = str(get_threshold(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt', fpr2))
            scores1 = pwm_scores
            scores2 = bamm_scores
            name1, name2 = 'PWM', 'BAMM'
            montecarlo(scores1, scores2, thr1, thr2, length, results, name1, name2, fpr1, fpr2)

            thr1 = str(get_threshold(motifs + '/' + tag + '_PWM_THRESHOLDS.txt', fpr1))
            thr2 = str(math.log2(get_threshold(motifs + '/' + tag + '_INMODE_THRESHOLDS.txt', fpr2)))
            scores1 = pwm_scores
            scores2 = inmode_scores
            name1, name2 = 'PWM', 'INMODE'
            montecarlo(scores1, scores2, thr1, thr2, length, results, name1, name2, fpr1, fpr2)

            thr1 = str(get_threshold(motifs + '/' + tag + '_BAMM_THRESHOLDS.txt', fpr1))
            thr2 = str(math.log2(get_threshold(motifs + '/' + tag + '_INMODE_THRESHOLDS.txt', fpr2)))
            scores1 = bamm_scores
            scores2 = inmode_scores
            name1, name2 = 'BAMM', 'INMODE'
            montecarlo(scores1, scores2, thr1, thr2, length, results, name1, name2, fpr1, fpr2)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', action='store', dest='bed_path',
                        required=True, help='path to BED file')
    parser.add_argument('-w', '--bigwig', action='store', dest='wig_path',
                        required=False, default="pass", help='path to BIGWIG file')
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

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():

    args = parse_args()

    bed_path = args.bed_path
    bigwig_path = args.wig_path
    path_to_out = args.output
    training_sample_size = args.train_size
    testing_sample_size = args.test_size
    fpr_for_thr = args.fpr


    path_to_python_tools = args.python_tools
    path_to_java = args.java
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    path_to_promoters = args.promoters
    path_to_genome = args.genome
    path_to_imd = args.path_to_imd
    path_to_tss = args.path_to_tss


    zoops=args.zoops
    cpu_count = args.cpu_count
    try_size=args.try_limit
    model_order=args.model_order

    pipeline_inmode_bamm(bed_path, bigwig_path, training_sample_size, testing_sample_size,
                          fpr_for_thr, path_to_out, path_to_python_tools, path_to_java, path_to_inmode, path_to_imd, path_to_chipmunk,
                          path_to_promoters, path_to_genome, path_to_tss, cpu_count,
                          zoops, try_size, model_order)

if __name__ == '__main__':
    main()
