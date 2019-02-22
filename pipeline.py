import os
import sys
import shlex
import subprocess
import argparse
import numpy as np


def work_with_tf_mono_version(bed_path, wig_path, training_sample_size, testing_sample_size,
                              list_fpr_for_thr, path_to_out, path_to_python_tools, dir_with_chipmunk,
                              path_to_promoters, path_to_genome, cpu_count,
                              wig_flag='wiggle', zoops=1.0, try_size=100,
                              bamm_order=2, recalculate_model=True):

    main_out = path_to_out + '/' + os.path.basename(bed_path).split('.')[0]
    zoops = str(zoops)
    bamm_order = str(bamm_order)
    try_size=str(try_size)
    cpu_count = str(cpu_count)

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
        #Get top 1000 bed peaks
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
        #Get top 8000 bed peaks
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
               '-of', fasta + '/' + tag + '_' + str(training_sample_size) +'.fa']
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

    if not os.path.isfile(fasta + '/' + tag + '_'+ str(training_sample_size) + '_WIG.fa'):
        #Bed peaks to wig_fasta
        print('Bed peaks to wig_fasta for {0}'.format(tag))
        if wig_flag == 'bedgraph':
            args = ['python3', path_to_python_tools + 'bed_to_bedgraph_fasta.py',
                   '-if', path_to_genome,
                   '-bed', bed + '/' + tag + '_' + str(training_sample_size) + '.bed',
                    '-w', wig_path,
                   '-of', fasta + '/' + tag + '_'+ str(training_sample_size) + '_WIG.fa']
            r = subprocess.call(args)
        elif wig_flag == 'wiggle':
            args = ['python3', path_to_python_tools + 'bed_to_wig_fasta.py',
                   '-if', path_to_genome,
                   '-bed', bed + '/' + tag + '_' + str(training_sample_size) + '.bed',
                    '-w', wig_path,
                   '-of', fasta + '/' + tag + '_'+ str(training_sample_size) + '_WIG.fa']
            r = subprocess.call(args)
        else:
            print('WRONG FLAG!')
            pass
    else:
        print('File {0} already exists'.format(tag + '_'+ str(training_sample_size) + '_WIG.fa'))

    #Calculate model
    if not os.path.isfile(motifs + '/' + tag + '_OPTIMAL.pwm') \
    or recalculate_model:

        #initial motif length
        ####################
        motif_length = 10
        ###################

        #ChIPMunk (find motif with current motif length)
        print('ChIPMunk find motifs for {0}'.format(tag))
        args = ['java', '-cp', dir_with_chipmunk + 'chipmunk.jar',
                       'ru.autosome.ChIPMunk', str(motif_length), str(motif_length), 'yes', zoops,
                       'p:' + fasta + '/' + tag + '_'+ str(training_sample_size) + '_WIG.fa',
                      try_size, '10', '1', cpu_count, 'random']
        path_out = chipmunk + '/RESULTS_' + str(motif_length) + '.txt'
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        out = p.communicate()
        with open(path_out, 'wb') as file:
            file.write(out[0])

        #Parse results of ChIPMunk into files .meme, .pwm and .fasta (multi fasta)
        args = ['python3', path_to_python_tools + 'parse_chipmunk_results.py',
               '-i', chipmunk + '/RESULTS_' + str(motif_length) + '.txt',
               '-o', motifs,
               '-t', tag + '_' + str(motif_length)]
        r = subprocess.call(args)

        #Calculate fpr while tpr = 0.5
        args = ['python3', path_to_python_tools + 'calculate_roc_pwm.py', 'get_fpr',
                '-f', motifs + '/' + tag + '_' + str(motif_length) + '.fasta',
                    '-t', '30000',
                    '-v', '0.5',
                   '-n', tag + '_' + str(motif_length),
                   '-o', motifs,
                   '-P', cpu_count]

        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        optimal = float(p.communicate()[0].decode('utf-8'))

        ####################
        optimal_length = 10
        #####################

        for motif_length in range(11, 30):
            #ChIPMunk (find motif with current motif length)
            args = ['java', '-cp', dir_with_chipmunk + 'chipmunk.jar',
                       'ru.autosome.ChIPMunk', str(motif_length), str(motif_length), 'yes', zoops,
                       'p:' + fasta + '/' + tag + '_' + str(training_sample_size) + '_WIG.fa',
                      try_size, '10', '1', cpu_count, 'random']
            path_out = chipmunk + '/RESULTS_' + str(motif_length) + '.txt'
            p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
            out = p.communicate()
            with open(path_out, 'wb') as file:
                file.write(out[0])

            #Parse results of ChIPMunk into files .meme, .pwm and .fasta (multi fasta)
            args = ['python3', path_to_python_tools + 'parse_chipmunk_results.py',
                    '-i', chipmunk + '/RESULTS_' + str(motif_length) + '.txt',
                    '-o', motifs,
                    '-t', tag + '_' + str(motif_length)]
            r = subprocess.call(args)

            #Calculate fpr while tpr = 0.5
            args = ['python3', path_to_python_tools + 'calculate_roc_pwm.py', 'get_fpr',
                    '-f', motifs + '/' + tag + '_' + str(motif_length) + '.fasta',
                    '-t', '30000',
                    '-v', '0.5',
                   '-n', tag + '_' + str(motif_length),
                   '-o', motifs,
                   '-P', cpu_count]

            p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
            tmp = float(p.communicate()[0].decode('utf-8'))
            #print(st.ttest_ind(optimal, tmp))
            print('FPR = {0} for matrix with length {1}, FPR = {2} for matrix with length {3}'.format(tmp, motif_length,
                                                                                                      optimal, optimal_length))

            if tmp <= optimal:
                optimal = tmp
                optimal_length = motif_length
                continue
            else:
                break

        #for file in os.listdir(motifs):
        #    os.remove(motifs + '/' + file)

        args = ['python3', path_to_python_tools + 'parse_chipmunk_results.py',
                '-i', chipmunk + '/RESULTS_' + str(optimal_length) + '.txt',
                '-o', motifs,
                '-t', tag + '_OPTIMAL']
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_OPTIMAL.pwm'))

################################################################################
    print('Calculate ROC for optimal matrix {0}'.format(tag))
    args = ['python3', path_to_python_tools + 'calculate_roc_pwm.py', 'roc',
                    '-f', motifs + '/' + tag + '_OPTIMAL.fasta',
                    '-t', '30000',
                   '-n', tag + '_' + 'OPTIMAL',
                   '-o', motifs,
                   '-P', cpu_count]
    p = subprocess.call(args)
################################################################################

    if not os.path.isfile(motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '_motif_1.ihbcp') or True:

        #Get BaMM motif
        print('Get Bamm motifs for {0}'.format(tag))
        args = ['BaMMmotif', motifs,
                fasta + '/' + tag + '_' + str(training_sample_size) + '_WIG.fa',
               '--PWMFile', motifs + '/' + tag + '_OPTIMAL.meme',
                '--basename', tag + '_OPTIMAL_ORDER_' + bamm_order,
               #'--EM',
               #'--CGS',
               '--Order', bamm_order,
               '--order', bamm_order]
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_OPTIMAL_ORDER_' + bamm_order + '_motif_1.ihbcp'))

    for fpr_for_thr in list_fpr_for_thr:

        if os.path.isfile(compare_sites + '/' + tag + '_' + str(fpr_for_thr) + '_' + '_FREQUENCY.tsv'):
            continue
        else:

            #Calculate threshold for PWM based on promoters and FPR = fpr_for_thr
            print('Calculate threshold for PWM based on promoters and FPR = {0} ({1})'.format(fpr_for_thr, tag))
            args = ['python3', path_to_python_tools + 'get_threshold_by_fp_numpy.py', 'pwm',
                    '-f', path_to_promoters,
                    '-m', motifs + '/' + tag + '_OPTIMAL.pwm',
                    '-p', str(fpr_for_thr),
                    '-P', cpu_count]
            p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
            thr_pwm = p.communicate()[0].decode('utf-8').strip()

            #Calculate threshold for BAMM based on promoters and FPR = fpr_for_thr
            print('Calculate threshold for BAMM based on promoters and FPR = {0} ({1})'.format(fpr_for_thr, tag))
            args = ['python3', path_to_python_tools + 'get_threshold_by_fp_numpy.py', 'bamm',
                    '-f', path_to_promoters,
                    '-m', motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '_motif_1.ihbcp',
                    '-b', motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '.hbcp',
                    '-p', str(fpr_for_thr),
                    '-P', cpu_count]
            p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
            thr_bamm = p.communicate()[0].decode('utf-8').strip()

            print('PWM = ',thr_pwm, 'BAMM = ', thr_bamm)
            #Scan peaks by PWM with thr_pwm
            print('Scan peaks by PWM with thr_pwm ({0})'.format(tag))
            args = ['python3', path_to_python_tools + 'scan_by_pwm.py',
                    '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                    '-m', motifs + '/' + tag + '_OPTIMAL.pwm',
                    '-t', thr_pwm,
                    '-o', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                    '-P', cpu_count]
            r = subprocess.call(args)


            #Scan peaks by BAMM with thr_bamm
            print('Scan peaks by BAMM with thr_pwm ({0})'.format(tag))
            args = ['python3', path_to_python_tools + 'scan_by_bamm.py',
                    '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                    '-m', motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '_motif_1.ihbcp',
                    '-b', motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '.hbcp',
                    '-t', thr_bamm,
                    '-o', scan + '/' + tag + '_BAMM_ORDER_' + bamm_order + '_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                    '-P', cpu_count]
            r = subprocess.call(args)


            #Compare sites
            print('Compare sites ({0})'.format(tag))
            args = ['python3', path_to_python_tools + 'compare_sites.py',
                    '-p', bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
                    '-m', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                    '-b', scan + '/' + tag + '_BAMM_ORDER_' + bamm_order + '_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                    '-t', tag + '_' + str(fpr_for_thr) + '_BAMM_ORDER_' + bamm_order,
                    '-o', compare_sites]
            r = subprocess.call(args)



def work_with_tf_di_version(bed_path, wig_path, training_sample_size, testing_sample_size,
                              list_fpr_for_thr, path_to_out, path_to_python_tools, dir_with_chipmunk,
                              path_to_promoters, path_to_genome,
                              wig_flag='wiggle', zoops=1.0, try_size=100,
                              bamm_order=2, recalculate_model=True):

    main_out = path_to_out + '/' + os.path.basename(bed_path).split('.')[0]
    zoops = str(zoops)
    bamm_order = str(bamm_order)
    try_size=str(try_size)
    cpu_count = str(cpu_count)

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
        #Get top 1000 bed peaks
        print('Get top {0} bed peaks for {1}'.format(training_sample_size, tag))
        args = ['python3', path_to_python_tools + 'get_top_peaks.py',
               '-i', bed_path,
               '-o', bed,
               '-a', str(training_sample_size),
               '-c', '4',
               '-t', tag + '_' + str(training_sample_size)]
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(bed + '/' + tag + '_' + str(training_sample_size) + '.bed'))

    if not os.path.isfile(bed + '/' + tag + '_' + str(testing_sample_size) + '.bed'):
        #Get top 8000 bed peaks
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
               '-of', fasta + '/' + tag + '_' + str(training_sample_size) +'.fa']
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

    if not os.path.isfile(fasta + '/' + tag + '_'+ str(training_sample_size) + '_WIG.fa'):
        #Bed peaks to wig_fasta
        print('Bed peaks to wig_fasta for {0}'.format(tag))
        if wig_flag == 'bedgraph':
            args = ['python3', path_to_python_tools + 'bed_to_bedgraph_fasta.py',
                   '-if', path_to_genome,
                   '-bed', bed + '/' + tag + '_' + str(training_sample_size) + '.bed',
                    '-w', wig_path,
                   '-of', fasta + '/' + tag + '_'+ str(training_sample_size) + '_WIG.fa']
            r = subprocess.call(args)
        elif wig_flag == 'wiggle':
            args = ['python3', path_to_python_tools + 'bed_to_wig_fasta.py',
                   '-if', path_to_genome,
                   '-bed', bed + '/' + tag + '_' + str(training_sample_size) + '.bed',
                    '-w', wig_path,
                   '-of', fasta + '/' + tag + '_'+ str(training_sample_size) + '_WIG.fa']
            r = subprocess.call(args)
        else:
            print('WRONG FLAG!')
            pass
    else:
        print('File {0} already exists'.format(tag + '_'+ str(training_sample_size) + '_WIG.fa'))


    #Calculate model
    if not os.path.isfile(motifs + '/' + tag + '_OPTIMAL.pwm') \
    or recalculate_model:

        #initial motif length
        ####################
        motif_length = 10
        ###################

        #ChIPMunk (find motif with current motif length)
        print('ChIPMunk find motifs for {0}'.format(tag))
        args = ['java', '-cp', dir_with_chipmunk + 'chipmunk.jar',
                       'ru.autosome.di.ChIPMunk', str(motif_length), str(motif_length), 'yes', zoops,
                       'p:' + fasta + '/' + tag + '_'+ str(training_sample_size) + '_WIG.fa',
                      try_size, '10', '1', cpu_count, 'random']
        path_out = chipmunk + '/RESULTS_' + str(motif_length) + '.txt'
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        out = p.communicate()
        with open(path_out, 'wb') as file:
            file.write(out[0])

        #Parse results of ChIPMunk into files .meme, .pwm and .fasta (multi fasta)
        args = ['python3', path_to_python_tools + 'parse_dichipmunk_results.py',
               '-i', chipmunk + '/RESULTS_' + str(motif_length) + '.txt',
               '-o', motifs,
               '-t', tag + '_' + str(motif_length)]
        r = subprocess.call(args)

        #Calculate fpr while tpr = 0.5
        args = ['python3', path_to_python_tools + 'calculate_roc_dipwm.py', 'get_fpr',
                '-f', motifs + '/' + tag + '_' + str(motif_length) + '.fasta',
                    '-t', '30000',
                    '-v', '0.5',
                   '-n', tag + '_' + str(motif_length),
                   '-o', motifs,
                   '-P', cpu_count]

        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        optimal = float(p.communicate()[0].decode('utf-8'))

        ####################
        optimal_length = 10
        #####################

        for motif_length in range(11, 30):
            #ChIPMunk (find motif with current motif length)
            args = ['java', '-cp', dir_with_chipmunk + 'chipmunk.jar',
                       'ru.autosome.di.ChIPMunk', str(motif_length), str(motif_length), 'yes', zoops,
                       'p:' + fasta + '/' + tag + '_' + str(training_sample_size) + '_WIG.fa',
                      try_size, '10', '1', cpu_count, 'random']
            path_out = chipmunk + '/RESULTS_' + str(motif_length) + '.txt'
            p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
            out = p.communicate()
            with open(path_out, 'wb') as file:
                file.write(out[0])

            #Parse results of ChIPMunk into files .meme, .pwm and .fasta (multi fasta)
            args = ['python3', path_to_python_tools + 'parse_dichipmunk_results.py',
                    '-i', chipmunk + '/RESULTS_' + str(motif_length) + '.txt',
                    '-o', motifs,
                    '-t', tag + '_' + str(motif_length)]
            r = subprocess.call(args)


            #Calculate fpr while tpr = 0.5
            args = ['python3', path_to_python_tools + 'calculate_roc_dipwm.py', 'get_fpr',
                    '-f', motifs + '/' + tag + '_' + str(motif_length) + '.fasta',
                    '-t', '30000',
                    '-v', '0.5',
                   '-n', tag + '_' + str(motif_length),
                   '-o', motifs,
                   '-P', cpu_count]

            p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
            tmp = float(p.communicate()[0].decode('utf-8'))
            #print(st.ttest_ind(optimal, tmp))
            print('FPR = {0} for matrix with length {1}, FPR = {2} for matrix with length {3}'.format(tmp, motif_length,
                                                                                                      optimal, optimal_length))


            if tmp <= optimal:
                optimal = tmp
                optimal_length = motif_length
                continue
            else:
                break

        #for file in os.listdir(motifs):
        #    os.remove(motifs + '/' + file)

        args = ['python3', path_to_python_tools + 'parse_dichipmunk_results.py',
                '-i', chipmunk + '/RESULTS_' + str(optimal_length) + '.txt',
                '-o', motifs,
                '-t', tag + '_OPTIMAL']
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_OPTIMAL.pwm'))

    if not os.path.isfile(motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '_motif_1.ihbcp') \
    or recalculate_model:

        #Get BaMM motif
        print('Get Bamm motifs for {0}'.format(tag))
        args = ['BaMMmotif', motifs,
                fasta + '/' + tag + '_' + str(training_sample_size) + '_WIG.fa',
               '--PWMFile', motifs + '/' + tag + '_OPTIMAL.meme',
                '--basename', tag + '_OPTIMAL_ORDER_' + bamm_order,
               #'--EM',
               #'--CGS',
               '--Order', bamm_order,
               '--order', bamm_order]
        r = subprocess.call(args)
    else:
        print('File {0} already exists'.format(tag + '_OPTIMAL_ORDER_' + bamm_order + '_motif_1.ihbcp'))

    for fpr_for_thr in list_fpr_for_thr:

        #Calculate threshold for PWM based on promoters and FPR = fpr_for_thr
        print('Calculate threshold for PWM based on promoters and FPR = {0} ({1})'.format(fpr_for_thr, tag))
        args = ['python3', path_to_python_tools + 'get_threshold_by_fp_numpy.py', 'pwm',
                '-f', path_to_promoters,
                '-m', motifs + '/' + tag + '_OPTIMAL.pwm',
                '-p', str(fpr_for_thr),
                '-P', cpu_count]
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        thr_pwm = p.communicate()[0].decode('utf-8').strip()

        #Calculate threshold for BAMM based on promoters and FPR = fpr_for_thr
        print('Calculate threshold for BAMM based on promoters and FPR = {0} ({1})'.format(fpr_for_thr, tag))
        args = ['python3', path_to_python_tools + 'get_threshold_by_fp_numpy.py', 'bamm',
                '-f', path_to_promoters,
                '-m', motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '_motif_1.ihbcp',
                '-b', motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '.hbcp',
                '-p', str(fpr_for_thr),
                '-P', cpu_count]
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        thr_bamm = p.communicate()[0].decode('utf-8').strip()

        print('PWM = ',thr_pwm, 'BAMM = ', thr_bamm)
        #Scan peaks by PWM with thr_pwm
        print('Scan peaks by PWM with thr_pwm ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'scan_by_pwm.py',
                '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                '-m', motifs + '/' + tag + '_OPTIMAL.pwm',
                '-t', thr_pwm,
                '-o', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-P', cpu_count]
        r = subprocess.call(args)


        #Scan peaks by BAMM with thr_bamm
        print('Scan peaks by BAMM with thr_pwm ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'scan_by_bamm.py',
                '-f', fasta + '/' + tag + '_' + str(testing_sample_size) + '.fa',
                '-m', motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '_motif_1.ihbcp',
                '-b', motifs + '/' + tag + '_OPTIMAL_ORDER_' + bamm_order + '.hbcp',
                '-t', thr_bamm,
                '-o', scan + '/' + tag + '_BAMM_ORDER_' + bamm_order + '_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-P', cpu_count]
        r = subprocess.call(args)


        #Compare sites
        print('Compare sites ({0})'.format(tag))
        args = ['python3', path_to_python_tools + 'compare_sites.py',
                '-p', bed + '/' + tag + '_' + str(testing_sample_size) + '.bed',
                '-m', scan + '/' + tag + '_PWM_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-b', scan + '/' + tag + '_BAMM_ORDER_' + bamm_order + '_' + str(testing_sample_size) +'_' + str(fpr_for_thr) + '.bed',
                '-t', tag + '_' + str(fpr_for_thr) + '_BAMM_ORDER_' + bamm_order,
                '-o', compare_sites]
        r = subprocess.call(args)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', action='store', dest='bed_path',
                        required=True, help='path to BED file')
    parser.add_argument('-w', '--wiggle', action='store', dest='wig_path',
                        required=True, help='path to WIG file or BedGraph file (if BedPraph use -B)')
    parser.add_argument('-P', '--promoters', action='store', dest='promoters',
                        required=True, help='path to promoters fasta file')
    parser.add_argument('-g', '--genome', action='store', dest='genome',
                        required=True, help='path to genome fasta file')
    parser.add_argument('-t', '--train', action='store', type=int, dest='train_size',
                        required=True, help='size of training sample')
    parser.add_argument('-T', '--test', action='store', type=int, dest='test_size',
                        required=True, help='size of testing sample')
    parser.add_argument('-f', '--fpr', action='store', dest='fpr', nargs='*',
                        required=False, default=[0.00005, 0.0001, 0.00025, 0.0005],
                        help='list of FPR values required to calculate threshold values \
                        Example: val1 val2 val3... , default: 0.00005, 0.0001, 0.00025, 0.0005')
    parser.add_argument('-p', '--python', action='store', dest='python_tools',
                        required=True, help='dir with python tools')
    parser.add_argument('-c', '--chipmunk', action='store', dest='chipmunk',
                        required=True, help='dir with chipmunk')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='output dir')
    parser.add_argument('-B', '--bedGraph_format', action='store_true', dest='bedGraph_flag',
                        required=False, help='if you use BedGraph in input you shuld use this flag')
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
    wig_path = args.wig_path

    path_to_out = args.output
    training_sample_size = args.train_size
    testing_sample_size = args.test_size

    list_fpr_for_thr = args.fpr
    list_fpr_for_thr = [float(i) for i in list_fpr_for_thr]
    #print(list_fpr_for_thr)

    bed_graph_flag = args.bedGraph_flag
    if bed_graph_flag:
        wig_flag='bedgraph'
    else:
        wig_flag='wiggle'

    zoops=args.zoops
    cpu_count = args.cpu_count
    try_size=args.try_limit
    bamm_order=args.bamm_order
    recalculate_model=False

    work_with_tf_mono_version(bed_path, wig_path, training_sample_size, testing_sample_size,
                              list_fpr_for_thr, path_to_out, path_to_python_tools, dir_with_chipmunk,
                              path_to_promoters, path_to_genome, cpu_count,
                              wig_flag=wig_flag, zoops=zoops, try_size=try_size,
                              bamm_order=2, recalculate_model=False)

if __name__ == '__main__':
    main()
