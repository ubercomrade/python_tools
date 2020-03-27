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
import itertools


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
        print('{0} already exists'.format('train_sample.bed'))

    if not os.path.isfile(bed + '/' + 'test_sample.bed'):
        #Get top testing_sample_size bed peaks
        print('Get top {0} bed peaks'.format(test_sample_size))
        bed_out = bed + '/'
        get_top_peaks(path_to_python_tools, bed_path, bed_out, test_sample_size, 'test_sample')
    else:
        print('{0} already exists'.format('test_sample.bed'))

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
        print('{0} already exists'.format('train_sample.fa'))

    if not os.path.isfile(fasta + '/' + 'test_sample.fa'):
        print('Get fasta from bed: {}'.format('test_sample.bed'))
        bed_to_fasta(path_to_genome,
            bed + '/test_sample.bed',
            fasta + '/test_sample.fa')
    else:
        print('{0} already exists'.format('test_sample.fa'))
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

    copyfile(glob.glob(inmode_model_path + '/Learned_DeNovo*/XML*')[0], inmode_model_path + '/inmode_model.xml')
    copyfile(glob.glob(inmode_model_path + '/Learned_DeNovo*/Binding_sites*')[0], inmode_model_path + '/inmode_sites.txt')
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


def get_pwm_model(models_path, fasta_path, path_to_python_tools, path_to_java, path_to_chipmunk, motif_length_start, motif_length_end, cpu_count):


    chipmunk_model_path = models_path + '/pwm_model'

    if not os.path.isdir(chipmunk_model_path):
        os.mkdir(chipmunk_model_path)

    # FIND MODEL BY CHIPMUNK #
    
    if not os.path.isfile(chipmunk_model_path + '/initial_pwm_model.pwm'):
        #Create fastaWig for chipmunk
        print('Create pmw model by ChIPMunk')
        run_chipmunk_fasta(path_to_java, path_to_chipmunk,
        fasta_path,
        chipmunk_model_path + '/pwm_model.txt',
        motif_length_start, motif_length_end, cpu_count)
    else:
        print('{0} already exists (model exists)'.format(chipmunk_model_path + '/initial_pwm_model.pwm'))

    # Parse results of chipmunk into files .meme, .pwm and .fasta (multi fasta) #
    
    args = ['python3', path_to_python_tools + 'parse_chipmunk_results.py',
           '-i', chipmunk_model_path + '/pwm_model.txt',
           '-o', chipmunk_model_path,
           '-t', 'initial_pwm_model']
    r = subprocess.call(args)

    # Get oPWM from chipmunk results. OUTPUT: .meme, .pwm and .fasta (multi fasta) #
    
    if not os.path.isfile(chipmunk_model_path + '/optimazed_pwm_model.meme'):
        args = ['python3', path_to_python_tools + 'make_oPWM.py',
                '-c', chipmunk_model_path + '/pwm_model.txt',
                '-f', fasta_path,
                '-n', '5000',
                '-P', cpu_count,
                '-o', chipmunk_model_path,
                '-t', 'optimazed_pwm_model']
        r = subprocess.call(args)
    else:
        print('{0} already exists'.format(chipmunk_model_path + '/optimazed_pwm_model.meme'))
    pass


def get_bamm_model(models_path, fasta_train, meme_model, model_order):

    #Get BaMM motif
    bamm_model_path = models_path + '/bamm_model'
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


def calculate_thresholds_for_bamm(path_to_python_tools, path_to_promoters, bamm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/bamm_model_thresholds.txt'):
        print('Calculate threshold for bamm based on promoters and fpr')
        args = ['pypy3', path_to_python_tools + '/get_threshold_for_bamm.py',
                path_to_promoters,
                bamm_model_dir + '/bamm_motif_1.ihbcp',
                bamm_model_dir + '/bamm.hbcp',
                thresholds_dir + '/bamm_model_thresholds.txt']
        r = subprocess.call(args)
    else:
        print('Thresholds for bamm already calculated')
    pass


def calculate_thresholds_for_pwm(path_to_python_tools, path_to_promoters, pwm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/pwm_model_thresholds.txt'):
        print('Calculate threshold for pwm based on promoters and fpr')
        args = ['pypy3', path_to_python_tools + '/get_threshold_for_pwm.py',
                path_to_promoters,
                pwm_model_dir + '/optimazed_pwm_model.pwm',
                thresholds_dir + '/pwm_model_thresholds.txt']
        r = subprocess.call(args)
    else:
        print('Thresholds for pwm already calculated')
    pass


def calculate_thresholds_for_inmode(path_to_python_tools, path_to_promoters, inmode_model_dir, thresholds_dir, motif_length, path_to_inmode, path_to_java):
    if not os.path.isfile(thresholds_dir + '/inmode_model_thresholds.txt'):
        args = ['python3', path_to_python_tools + '/get_threshold_for_inmode.py',
                path_to_promoters,
                inmode_model_dir + '/inmode_model.xml',
                path_to_inmode,
                str(motif_length),
                thresholds_dir + '/inmode_model_thresholds.txt',
                '-j', path_to_java,
                '-t', thresholds_dir + '/tmp']
        r = subprocess.call(args)
    else:
        print('Thresholds for inmode already calculated')
    pass


def scan_by_pwm(path_to_python_tools, fasta_test, model_path, scan, threshold_table_path, fpr):

    thr_pwm = get_threshold(threshold_table_path, fpr)
    pwm_scan_path = scan + '/pwm_{:.2e}.bed'.format(fpr)
    print('Scan peaks by PWM with FPR: {0} THR: {1}'.format(fpr, thr_pwm))
    args = ['pypy', path_to_python_tools + 'scan_by_pwm.py',
            '-f', fasta_test,
            '-m', model_path,
            '-t', str(thr_pwm),
            '-o', pwm_scan_path]
    r = subprocess.call(args)
    pass


def scan_by_bamm(path_to_python_tools, fasta_test, model_path, bg_model_path, scan, threshold_table_path, fpr):

    thr_bamm = get_threshold(threshold_table_path, fpr)
    bamm_scan_path = scan + '/bamm_{:.2e}.bed'.format(fpr)
    print('Scan peaks by BAMM with FPR: {0} THR: {1}'.format(fpr, thr_bamm))
    args = ['pypy', path_to_python_tools + 'scan_by_bamm.py',
            '-f', fasta_test,
            '-m', model_path,
            '-b', bg_model_path,
            '-t', str(thr_bamm),
            '-o', bamm_scan_path]
    r = subprocess.call(args)
    pass


def scan_by_inmode(path_to_python_tools, fasta_test, model_path, scan, threshold_table_path, fpr, path_to_java, path_to_inmode, path_to_promoters):

    inmode_scan_dir = scan + '/tmp'
    inmode_scan_path = scan + '/inmode_{:.2e}.bed'.format(fpr)

    thr_inmode = get_threshold(threshold_table_path, fpr)
    print('Scan peaks by INMODE with FPR: {0} THR: {1}'.format(fpr, thr_inmode))
    args = [path_to_java, '-Xmx6G', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode, 'scan',
            'i={}'.format(model_path),
            'id={}'.format(fasta_test),
            'b={}'.format('From file'),
            'd={}'.format(path_to_promoters),
           'f={}'.format(fpr),
           'outdir={}'.format(inmode_scan_dir)]
    r = subprocess.call(args)

    args = ['python3', path_to_python_tools + 'parse_inmode_scan.py',
            '-if', fasta_test,
            '-bed', glob.glob(inmode_scan_dir + '/*.BED')[0],
            '-o', inmode_scan_path]
    r = subprocess.call(args)

    os.system("rm -r {}".format(inmode_scan_dir))
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
    args = ['pypy3', path_to_python_tools + '/scan_best_by_inmode_alt.py',
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
    args = ['pypy3', path_to_python_tools + 'scan_best_by_pwm.py',
            '-f', fasta_in,
            '-m', input_model,
            '-o', output]
    r = subprocess.call(args)
    pass


def scan_best_by_bamm(path_to_python_tools, output, input_bamm_model, bg_model, fasta_in):
    args = ['pypy3', path_to_python_tools + '/scan_best_by_bamm.py',
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
            '-M', '-p']
    r = subprocess.call(args)
    pass


def get_motif_length(models):
    with open(models + '/pwm_model/optimazed_pwm_model.fasta', 'r') as file:
        for i in file:
            if i.startswith('>'):
                continue
            else:
                motif_length = len(i.strip())
                break
    file.close()
    return(motif_length)


def compare_by_pair(bed, tool_names, scan, compare, path_to_python_tools):
    list(itertools.combinations(a, 2))

    args = ['pypy3', path_to_python_tools + '/compare_scripts/compare_sites_2.py',
                    '-p', bed,
                    '-first', first,
                    '-second', second,
                    '-t', tag, '-fname', name1, '-sname', name2,
                    '-o', compare_sites]
    r = subprocess.call(args)
    pass


def extracting_sites(path_to_python_tools, scan_file, output):
        args = ['python3', path_to_python_tools + 'extract_sites.py',
            '-p', scan_file,
            '-o', output]
        r = subprocess.call(args)


def compare_2(bed, first, second, tag, name1, name2, compare_sites, path_to_python_tools):
    args = ['pypy3', path_to_python_tools + '/compare_scripts/compare_sites_2.py',
                    '-p', bed,
                    '-first', first,
                    '-second', second,
                    '-t', tag, '-fname', name1, '-sname', name2,
                    '-o', compare_sites]
    r = subprocess.call(args)
    pass


def compare_4(bed, first, second, third, fourth, tag, compare_sites, path_to_python_tools):
    args = ['python3', path_to_python_tools + '/compare_scripts/compare_sites_4.py',
                    '-p', bed,
                    '-first', first,
                    '-second', second,
                    '-third', third,
                    '-fourth', fourth,
                    '-t', tag,
                    '-o', compare_sites]
    r = subprocess.call(args)
    pass


def montecarlo(path_to_python_tools, scores1, scores2, thr1, thr2, length, results):
        args = ['monteCarlo', '{}'.format(scores1), '{}'.format(scores2), '{}'.format(thr1), '{}'.format(thr2), '{}'.format(length), '{}'.format(results)]
        r = subprocess.call(args)
        pass


def pipeline(tools, bed_path, fpr, train_sample_size, test_sample_size,
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
    bootstrap = models + '/bootstrap'
    thresholds = models + '/thresholds'
    fasta = main_out + '/fasta'
    bed = main_out + '/bed'
    scan = main_out + '/scan'
    scan_best = main_out + '/scan-best'
    compare_sites = main_out + '/compare'
    tomtom = main_out + '/tomtom'
    montecarlo = main_out + '/montecarlo'
    #tools = ['pwm', 'bamm', 'inmode', 'sitega']
    #tag = os.path.basename(bed_path).split('.bed')[0]


    ########################
    #      CREATE DIRS     #
    ########################

    if not os.path.isdir(models):
        os.mkdir(models)
    if not os.path.isdir(bootstrap):
        os.mkdir(bootstrap)
    if not os.path.isdir(thresholds):
        os.mkdir(thresholds)
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
    if not os.path.isdir(montecarlo):
        os.mkdir(montecarlo)



    # PREPARE BED AND FASTA FILES #
    prepare_data(path_to_python_tools, path_to_genome, bed_path, bed, fasta, train_sample_size, test_sample_size)

    fasta_train = fasta + '/train_sample.fa'
    fasta_test = fasta + '/test_sample.fa'
    bed_test = bed + '/test_sample.bed'
    bed_train = bed + '/train_sample.bed'

    # CALCULATE CHIPMUNK MODEL #
    pwm_model = models + '/pwm_model/optimazed_pwm_model.pwm'
    pwm_threshold_table = thresholds + '/pwm_model_thresholds.txt'
    get_pwm_model(models, fasta_train, path_to_python_tools, 
        path_to_java, path_to_chipmunk,
        motif_length_start, motif_length_end,
        cpu_count)
    calculate_thresholds_for_pwm(path_to_python_tools, path_to_promoters, models + '/pwm_model', thresholds)
    # bootstrap_pwm(path_to_python_tools,
    #     bootstrap + "/pwm_model.tsv",
    #     models + '/chipmunk_model/optimazed_pwm_model.fasta')
    scan_by_pwm(path_to_python_tools, fasta_test, pwm_model, scan, pwm_threshold_table, fpr)
    scan_best_by_pwm(path_to_python_tools, scan_best + '/pwm.scores.txt',
         pwm_model,
         fasta_test)
    extracting_sites(path_to_python_tools, scan + '/pwm_{:.2e}.bed'.format(fpr), tomtom + '/pwm.sites.txt')
    make_model(path_to_python_tools, tomtom + '/pwm.sites.txt', tomtom, 'pwm')
    run_tomtom(path_to_hocomoco, tomtom + '/pwm.meme', tomtom + '/pwm')

    # GET MOTIF LENGTH #
    motif_length = get_motif_length(models)

    # CALCULATE INMODE MODEL WITH EM ALG #
    inmode_model = models + '/inmode_model/inmode_model.xml'
    inmode_threshold_table = thresholds + '/inmode_model_thresholds.txt'

    get_inmode_model(models, fasta_train, path_to_java,
        path_to_inmode, motif_length, model_order)
    calculate_thresholds_for_inmode(path_to_python_tools, path_to_promoters, models + '/inmode_model', thresholds, motif_length, path_to_inmode, path_to_java)
    # bootstrap_inmode(path_to_python_tools, path_to_java, 
    #     bootstrap + "/inmode.tsv", 
    #     models + '/inmode_model/inmode_sites.txt',
    #     path_to_inmode, models + '/tmp')
    scan_by_inmode(path_to_python_tools, fasta_test, inmode_model, scan, inmode_threshold_table, fpr, path_to_java, path_to_inmode, path_to_promoters)
    scan_best_by_inmode(path_to_python_tools, scan_best + '/inmode.scores.txt',
        inmode_model,
        fasta_test,
        path_to_inmode, path_to_java)
    extracting_sites(path_to_python_tools, scan + '/inmode_{:.2e}.bed'.format(fpr), tomtom + '/inmode.sites.txt')
    make_model(path_to_python_tools, tomtom + '/inmode.sites.txt', tomtom, 'inmode')
    run_tomtom(path_to_hocomoco, tomtom + '/inmode.meme', tomtom + '/inmode')

    # CALCULATE BAMM MODEL WITH EM ALG #
    meme_model = models + '/pwm_model/optimazed_pwm_model.meme'
    bamm_threshold_table = thresholds + '/bamm_model_thresholds.txt'
    bamm_model = models + '/bamm_model/bamm_motif_1.ihbcp'
    bg_bamm_model = models + '/bamm_model/bamm.hbcp'

    get_bamm_model(models, fasta_train, meme_model, model_order)   
    calculate_thresholds_for_bamm(path_to_python_tools, path_to_promoters, models + '/bamm_model', thresholds)
    # bootstrap_bamm(path_to_python_tools,
    #     bootstrap + "/bamm.tsv",
    #     models + '/bamm_model/bamm_motif_1.logOddsZoops')
    scan_by_bamm(path_to_python_tools, fasta_test, bamm_model, bg_bamm_model, scan, bamm_threshold_table, fpr)
    scan_best_by_bamm(path_to_python_tools, scan_best + '/bamm.scores.txt',
        bamm_model,
        bg_bamm_model,
        fasta_test)
    extracting_sites(path_to_python_tools, scan + '/bamm_{:.2e}.bed'.format(fpr), tomtom + '/bamm.sites.txt')
    make_model(path_to_python_tools, tomtom + '/bamm.sites.txt', tomtom, 'bamm')
    run_tomtom(path_to_hocomoco, tomtom + '/bamm.meme', tomtom + '/bamm')
   

    # COMPARE SITES #
    print('COMPARE SITES')
    pair_tools = list(itertools.combinations(tools, 2))
    for tool1, tool2 in pair_tools:
        tag = 'compare'
        scan1 = scan + '/{0}_{1:.2e}.bed'.format(tool1, fpr)
        scan2 = scan + '/{0}_{1:.2e}.bed'.format(tool2, fpr)
        compare_2(bed_test, scan1, scan2, tag, tool1, tool2, compare_sites, path_to_python_tools)

    if len(tools) == 4:
        compare_4(bed_test,
            scan + '/{0}_{1:.2e}.bed'.format('pwm', fpr),
            scan + '/{0}_{1:.2e}.bed'.format('bamm', fpr),
            scan + '/{0}_{1:.2e}.bed'.format('inmode', fpr),
            scan + '/{0}_{1:.2e}.bed'.format('sitega', fpr),
            tag,
            compare_sites,
            path_to_python_tools)


    # MONTECARLO #
    # for tool1, tool2 in pair_tools:
    #     thr1 = str(get_threshold(thresholds + '/{}_model_thresholds.txt'.format(tool1), fpr))
    #     thr2 = str(get_threshold(thresholds + '/{}_model_thresholds.txt'.format(tool2), fpr))
    #     scores1 = scan_best + '/{}.scores.txt'.format(tool1)
    #     scores2 = scan_best + '/{}.scores.txt'.format(tool2)
    #     name1, name2 = tool1, tool2
    #     results_montecarlo = montecarlo + '/montecarlo.results.{0}.{1}.{2:.2e}.txt'.format(tool1, tool2, fpr)
    #     montecarlo(path_to_python_tools, scores1, scores2, thr1, thr2, length, results_montecarlo)



def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store', help='path to BED file')
    parser.add_argument('promoters', action='store', help='path to promoters fasta file')
    parser.add_argument('genome', action='store', help='path to genome fasta file')
    parser.add_argument('output', action='store', help='output dir')
    parser.add_argument('models', action='store', choices=['pwm', 'bamm', 'inmode', 'sitega'], metavar='N', nargs='+',
         help='list of models to use (pwm, bamm, inmode, sitega)')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=False, default=500, help='size of training sample, by default size is equal to 500')
    parser.add_argument('-f', '--FPR', action='store', type=float, dest='fpr',
                        required=False, default=1.9*10**(-4), help='FPR, def=1.9*10^(-4)')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=False, default=4000, help='size of testing sample, by default size is equal to 4000')
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
    fpr = args.fpr
    tools = args.models

    path_to_python_tools = args.python_tools
    path_to_java = args.java
    path_to_chipmunk = args.chipmunk
    path_to_inmode = args.inmode
    path_to_promoters = args.promoters
    path_to_genome = args.genome
    path_to_hocomoco = args.path_to_hocomoco
    cpu_count = args.cpu_count


    pipeline(tools, bed_path, fpr, train_sample_size, test_sample_size,
                          path_to_out, path_to_python_tools, path_to_java, path_to_inmode, path_to_chipmunk,
                          path_to_promoters, path_to_genome, path_to_hocomoco, cpu_count)

if __name__ == '__main__':
    main()
