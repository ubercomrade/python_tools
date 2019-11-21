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
    parser.add_argument('thr1', action='store', type=float, help='threshold for first model')
    parser.add_argument('thr2', action='store', type=float, help='threshold for second model')
    parser.add_argument('out', action='store', help='path to write results')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def write_results(p_val, wpath):
    #Real Rand SD Zscore Log10Pval
    with open(wpath, 'a') as file:
        file.write('{:.4f}\n'.format(p_val))
    file.close()


def main():
    args = parse_args()
    path1 = args.score1
    path2 = args.score2
    thr1 = args.thr1
    thr2 = args.thr2
    wpath = args.out
    scores1 = read_scores(path1)
    scores2 = read_scores(path2)

    meet = sum(1 if s1 > thr1 and s2 > thr2 else 0 for s1, s2 in zip(scores1, scores2))   
    upper1 = sum(1 if s1 > thr1 else 0 for s1 in scores1)
    upper2 = sum(1 if s2 > thr2 else 0 for s2 in scores2)
    p_val = stats.binom_test(x=meet, n=len(scores1), p=(upper1 / len(scores1)) * (upper2 / len(scores2)), alternative='greater')
    write_results(p_val, wpath)
    
if __name__=='__main__':
    main()