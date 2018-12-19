import os
import sys
import shlex
import subprocess
import numpy as np
import scipy.stats as st

input_bed = 'path_to_file_with_bed_full_names'
input_wig = 'path_to_file_with_wig_full_names'
fpr_for_thr = 'thr'
path_to_promoters = 'path_to_promoters'
path_to_genome = 'path_to_genome'
path_to_out = 'path'


def work_with_tf(bed_path, wig_path, fpr_for_thr, path_to_out):

    main_out = path_to_out + '/' + os.path.basename(bed_path).split('.')[0]

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

    # Get top 1000 bed peaks
    print('Get top 1000 bed peaks for {0}'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/get_top_peaks.py',
            '-i', bed_path,
            '-o', bed,
            '-a', '1000',
            '-c', '4',
            '-t', tag + '_1000']
    r = subprocess.call(args)

    # Get top 8000 bed peaks
    print('Get top 8000 bed peaks for {0}'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/get_top_peaks.py',
            '-i', bed_path,
            '-o', bed,
            '-a', '8000',
            '-c', '4',
            '-t', tag + '_8000']
    r = subprocess.call(args)

    # Bed peaks to fasta
    print('Bed peaks to fasta for {0}'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/bed_to_fasta.py',
            '-if', path_to_genome,
            '-bed', bed + '/' + tag + '_1000.bed',
            '-of', fasta + '/' + tag + '_1000.fa']
    r = subprocess.call(args)

    # Bed peaks to fasta
    args = ['python3', '/home/anton/Programs/python_tools/bed_to_fasta.py',
            '-if', path_to_genome,
            '-bed', bed + '/' + tag + '_8000.bed',
            '-of', fasta + '/' + tag + '_8000.fa']
    r = subprocess.call(args)

    # Bed peaks to wig_fasta
    print('Bed peaks to wig_fasta for {0}'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/bed_to_wig_fasta.py',
            '-if', path_to_genome,
            '-bed', bed + '/' + tag + '_1000.bed',
            '-w', wig_path,
            '-of', fasta + '/' + tag + '_1000_WIG.fa']
    r = subprocess.call(args)

    # initial motif length
    motif_length = 10
    # ChIPMunk (find motif with current motif length)
    print('ChIPMunk find motifs for {0}'.format(tag))
    args = ['java', '-cp', '/home/anton/Programs/ChIPMunk/chipmunk.jar',
            'ru.autosome.ChIPMunk', str(motif_length), str(motif_length), 'yes', '1.0',
            'p:' + fasta + '/' + tag + '_1000_WIG.fa',
            '100', '10', '1', '4', 'random']
    path_out = chipmunk + '/RESULTS_' + str(motif_length) + '.txt'
    p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    out = p.communicate()
    with open(path_out, 'wb') as file:
        file.write(out[0])

    # Parse results of ChIPMunk into files .meme, .pwm and .fasta (multi fasta)
    args = ['python3', '/home/anton/Programs/python_tools/parse_chipmunk_results.py',
            '-i', chipmunk + '/RESULTS_' + str(motif_length) + '.txt',
            '-o', motifs,
            '-t', tag + '_' + str(motif_length)]
    r = subprocess.call(args)

    # Calculate tpr while fpr = 0.0005
    args = ['python3', '/home/anton/Programs/python_tools/calculate_roc.py', 'tpr',
            '-f', motifs + '/' + tag + '_' + str(motif_length) + '.fasta',
            '-i', '10',
            '-t', '300',
            '-v', '0.0005']

    p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    optimal = np.array([float(i) for i in p.communicate()[0].decode('utf-8').split()])

    optimal_length = 10
    for motif_length in range(11, 21):
        # ChIPMunk (find motif with current motif length)
        args = ['java', '-cp', '/home/anton/Programs/ChIPMunk/chipmunk.jar',
                'ru.autosome.ChIPMunk', str(motif_length), str(motif_length), 'yes', '1.0',
                'p:' + fasta + '/' + tag + '_1000_WIG.fa',
                '100', '10', '1', '4', 'random']
        path_out = chipmunk + '/RESULTS_' + str(motif_length) + '.txt'
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        out = p.communicate()
        with open(path_out, 'wb') as file:
            file.write(out[0])

        # Parse results of ChIPMunk into files .meme, .pwm and .fasta (multi fasta)
        args = ['python3', '/home/anton/Programs/python_tools/parse_chipmunk_results.py',
                '-i', chipmunk + '/RESULTS_' + str(motif_length) + '.txt',
                '-o', motifs,
                '-t', tag + '_' + str(motif_length)]
        r = subprocess.call(args)

        # Calculate tpr while fpr = 0.0005
        args = ['python3', '/home/anton/Programs/python_tools/calculate_roc.py', 'tpr',
                '-f', motifs + '/' + tag + '_' + str(motif_length) + '.fasta',
                '-i', '10',
                '-t', '300',
                '-v', '0.0005']
        p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
        tmp = np.array([float(i) for i in p.communicate()[0].decode('utf-8').split()])
        print(st.ttest_ind(optimal, tmp))
        print(np.mean(tmp), np.mean(optimal))
        if st.ttest_ind(optimal, tmp)[1] < 0.01 and np.mean(tmp) > np.mean(optimal):
            optimal = tmp[:]
            optimal_length = motif_length
        elif np.mean(tmp) > np.mean(optimal):
            optimal = tmp[:]
            optimal_length = motif_length
            break
        else:
            break

    for file in os.listdir(motifs):
        os.remove(motifs + '/' + file)

    args = ['python3', '/home/anton/Programs/python_tools/parse_chipmunk_results.py',
            '-i', chipmunk + '/RESULTS_' + str(optimal_length) + '.txt',
            '-o', motifs,
            '-t', tag + '_OPTIMAL']
    r = subprocess.call(args)

    # Get BaMM motif
    print('Get Bamm motifs for {0}'.format(tag))
    args = ['BaMMmotif', motifs,
            fasta + '/' + tag + '_1000_WIG.fa',
            '--PWMFile', motifs + '/' + tag + '_OPTIMAL.meme',
            '--basename', tag + '_OPTIMAL'
            '--EM']
    r = subprocess.call(args)

    # Calculate threshold for PWM based on promoters and FPR = 0.0001
    print('Calculate threshold for PWM based on promoters and FPR = 0.0001 ({0})'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/get_threshold_by_fp.py', 'pwm',
            '-f', path_to_promoters,
            '-m', motifs + '/' + tag + '_OPTIMAL.pwm',
            '-p', '0.0001']
    p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    thr_pwm = p.communicate()[0].decode('utf-8').strip()

    # Calculate threshold for BAMM based on promoters and FPR = 0.0001
    print('Calculate threshold for BAMM based on promoters and FPR = 0.0001 ({0})'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/get_threshold_by_fp.py', 'bamm',
            '-f', '/home/anton/DATA/PROMOTERS/hg38.fa',
            '-m', motifs + '/' + tag + '_OPTIMAL--EM_motif_1.ihbcp',
            '-b', motifs + '/' + tag + '_OPTIMAL--EM.hbcp',
            '-p', '0.0001']
    p = subprocess.Popen(args, shell=False, stdout=subprocess.PIPE)
    thr_bamm = p.communicate()[0].decode('utf-8').strip()

    # Scan peaks by PWM with thr_pwm
    print('Scan peaks by PWM with thr_pwm ({0})'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/scan_by_pwm.py',
            '-f', fasta + '/' + tag + '_8000.fa',
            '-m', motifs + '/' + tag + '_OPTIMAL.pwm',
            '-t', thr_pwm,
            '-o', scan + '/' + tag + '_PWM_8000.bed']
    r = subprocess.call(args)

    # Scan peaks by BAMM with thr_bamm
    print('Scan peaks by BAMM with thr_pwm ({0})'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/scan_by_bamm.py',
            '-f', fasta + '/' + tag + '_8000.fa',
            '-m', motifs + '/' + tag + '_OPTIMAL--EM_motif_1.ihbcp',
            '-b', motifs + '/' + tag + '_OPTIMAL--EM.hbcp',
            '-t', thr_bamm,
            '-o', scan + '/' + tag + '_BAMM_8000.bed']
    r = subprocess.call(args)

    # Compare sites
    print('Compare sites ({0})'.format(tag))
    args = ['python3', '/home/anton/Programs/python_tools/compare_sites.py',
            '-p', bed + '/' + tag + '_8000.bed',
            '-m', scan + '/' + tag + '_PWM_8000.bed',
            '-b', scan + '/' + tag + '_BAMM_8000.bed',
            '-t', tag,
            '-o', compare_sites]
    r = subprocess.call(args)
