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


import argparse
import sys
import seaborn as sns
import matplotlib.pyplot as plt


def read_list(path):
    container = list()
    append = container.append
    with open(path, "r") as file:
        for line in file:
            append(float(line.strip()))
    file.close()
    return(container)


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('scores1', action='store', help='path to file with scores for first model')
    parser.add_argument('scores2', action='store', help='path to file with scores for second model')
    parser.add_argument('thr1', action='store', type=float, help='threshold for first model')
    parser.add_argument('thr2', action='store', type=float, help='threshold for second model')
    parser.add_argument('length', action='store', help='path to file with length of peaks')
    parser.add_argument('out', action='store', help='path to write figure with (.pdf, .png format)')
    parser.add_argument('-n1', '--first_model_name', action='store', dest='first_model_name',
                        default="First model", required=False, help='Name for first model (X - axis)')
    parser.add_argument('-n2', '--second_model_name', action='store', dest='second_model_name',
                        default="Second model", required=False, help='Name for second model (Y - axis)')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())

def main():
    args = parse_args()

    first_model_path = args.scores1
    second_model_path = args.scores2
    first_model_thr = args.thr1
    second_model_thr = args.thr2
    length = args.length
    out = args.out
    xlabel = args.first_model_name
    ylabel = args.second_model_name
    
    
    first_model_scores = read_list(first_model_path)
    second_model_scores = read_list(second_model_path)
    
    sns.scatterplot(pwm_score, bamm_score)
    plt.axvline(first_model_thr, color="r")
    plt.axhline(second_model_thr, color="r")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(out, dpi=250)
    
if __name__=="__main__":
    main()
    
