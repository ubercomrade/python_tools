import os
import sys
import shlex
import subprocess
import argparse
import glob
import numpy as np


def InMoDeCLI_scan(path_to_inmode, input_data, input_model, outdir,
                     fpr_for_thr=1):
    #backgroud_path
    #args = ['/home/anton/Programs/jdk-9/bin/java', '-Xmx4096m', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode, 'scan',
    args = ['/Users/anton/Documents/Programs/jre-9.0.4.jre/Contents/Home/bin/java', '-Xmx3072m', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode, 'scan',
            'i={}'.format(input_model),
            'id={}'.format(input_data),
           'f={}'.format(fpr_for_thr),
           'outdir={}'.format(outdir)]
    r = subprocess.call(args)
    pass

def parse_inmode_results(python_path, fasta_path, results_dir, tag):

    args = ['python3', python_path + 'parse_inmode_scan.py',
                    '-if', fasta_path,
                    '-bed', glob.glob(results_dir + '/*.BED')[0],
                    '-o', results_dir + '/' + tag + '_BEST_INMODE.bed']
    r = subprocess.call(args)
    pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file')
    parser.add_argument('-m', '--model', action='store', dest='input_model',
                        required=True, help='path to xml inmode model file')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='dir to write results')
    parser.add_argument('-I', '--InMoDe', action='store', dest='inmode',
                        required=True, help='path to InMoDe')
    parser.add_argument('-p', '--python', action='store', dest='python_tools',
                        required=True, help='dir with python tools')
    parser.add_argument('-t', '--tag', action='store', dest='tag',
                        required=True, help='TAG for output files')
    
    return(parser.parse_args())


def main():

    args = parse_args()
    model_path = args.input_model
    inmode_path = args.inmode
    fasta_path = args.input_fasta
    results_dir = args.output
    python_path = args.python_tools
    tag = args.tag

    InMoDeCLI_scan(inmode_path, fasta_path, model_path, results_dir)
    parse_inmode_results(python_path, fasta_path, results_dir, tag)


if __name__ == '__main__':
    main()







