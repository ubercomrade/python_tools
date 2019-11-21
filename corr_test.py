import random
import operator
import argparse
import math
import sys
from bisect import bisect
from scipy import stats


def read_scores(path):
    with open(path, 'r') as file:
        values = [float(line.strip()) for line in file]
    file.close()
    return(values)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('score1', action='store', help='path to file with score of first model')
    parser.add_argument('score2', action='store', help='path to file with score of second model')
    parser.add_argument('out', action='store', help='path to write results')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def write_results(corr, p_val, wpath):
    #Real Rand SD Zscore Log10Pval
    with open(wpath, 'a') as file:
        file.write('{0:.4f}\t{1:.4f}\n'.format(corr, p_val))
    file.close()


def main():
    args = parse_args()
    path1 = args.score1
    path2 = args.score2
    wpath = args.out
    scores1 = read_scores(path1)
    scores2 = read_scores(path2)
    corr, p_val = stats.pearsonr(scores1, scores2)
    write_results(corr, p_val, wpath)
    
if __name__=='__main__':
    main()