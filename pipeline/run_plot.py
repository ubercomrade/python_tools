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


def run_tomtom(query, model, outdir):
    args = ['tomtom', query, model, '-oc', outdir]
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


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('tag', action='store', help='tag')
    parser.add_argument('output', action='store', help='output dir')
    parser.add_argument('fpr', action='store', type=float, help='fasle positive rates')
    parser.add_argument('-p', '--python', action='store', dest='python_tools',
                        required=True, help='dir with python tools')
    parser.add_argument('-H', '--hocomoco', action='store', dest='path_to_hocomoco',
                        required=True, help='path to HOCOMOCO database in meme format for TOMTOM')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():

    args = parse_args()

    tag = args.tag
    path_to_out = args.output
    fpr = args.fpr
    path_to_python_tools = args.python_tools
    path_to_hocomoco = args.path_to_hocomoco

    main_out = path_to_out + '/' + tag
    if not path_to_python_tools[-1] == '/':
        path_to_python_tools += '/'
    
    tomtom = main_out + '/TOMTOM'
    scan = main_out + '/SCAN'
    

    ########################
    #      CREATE DIRS     #
    ########################

    if not os.path.isdir(main_out + '/TOMTOM'):
        os.mkdir(main_out + '/TOMTOM')

    ########################################
    # CREATE MODELS FROM SITES AND COMPARE #
    ########################################

    print('Run tomtom')

    make_model(path_to_python_tools, scan + '/' + 'inmode.sites.{:.2e}.txt'.format(fpr),
        tomtom, 'inmode')
    make_model(path_to_python_tools, scan + '/' + 'pwm.sites.{:.2e}.txt'.format(fpr),
        tomtom, 'pwm')
    make_model(path_to_python_tools, scan + '/' + 'bamm.sites.{:.2e}.txt'.format(fpr),
        tomtom, 'bamm')

    run_tomtom(tomtom + '/pwm.meme', tomtom + '/bamm.meme', tomtom + '/pwm.bamm')
    run_tomtom(tomtom + '/pwm.meme', tomtom + '/inmode.meme', tomtom + '/pwm.inmode')
    run_tomtom(tomtom + '/inmode.meme', tomtom + '/bamm.meme', tomtom + '/inmode.bamm')

    run_tomtom(path_to_hocomoco, tomtom + '/pwm.meme', tomtom + '/pwm.hocomoco')
    run_tomtom(path_to_hocomoco, tomtom + '/inmode.meme', tomtom + '/inmode.hocomoco')
    run_tomtom(path_to_hocomoco, tomtom + '/bamm.meme', tomtom + '/bamm.hocomoco')

    
if __name__ == '__main__':
    main()
