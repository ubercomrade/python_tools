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

def read_length(path):
    with open(path, 'r') as file:
        values = [int(line.strip()) for line in file]
    file.close()
    return(values)


def sort_scores(scores, length):
    getcount = operator.itemgetter(0)
    sorted_score_len = sorted(zip(scores,length), key=operator.itemgetter(1))
    return(list(map(getcount, sorted_score_len)))


def split_scores(scores, length, shift=0, step=100):
    last_index = 0
    container = list()
    append = container.append
    for i in range(min(length) + shift, max(length) + 50, step):
        index = bisect(length, i)
        s = scores[last_index:index]
        if len(s) < 10:
            continue
        else:
            append(s)
            last_index = index
    container[-1] = container[-1] + scores[last_index:]  
    return(container)
    
    
# def montecarlo(scores1, scores2, length, thr1, thr2, shift=50, step=100):
#     container = list()
#     append = container.append
#     scores1 = split_scores(scores1, length, shift=0, step=step)
#     scores2 = split_scores(scores2, length, shift=0, step=step)
#     for i in range(10000):

#         scores1r = [random.sample(subscores, len(subscores)) for subscores in scores1]
#         scores1r = [score for subscores in scores1r for score in subscores]
#         scores2r = [random.sample(subscores, len(subscores)) for subscores in scores2]
#         scores2r = [score for subscores in scores2r for score in subscores]
#         scores1r = split_scores(scores1r, length, shift=shift, step=step)
#         scores2r = split_scores(scores2r, length, shift=shift, step=step)
#         scores1r = [random.sample(subscores, len(subscores)) for subscores in scores1r]
#         scores1r = [score for subscores in scores1r for score in subscores]
#         scores2r = [random.sample(subscores, len(subscores)) for subscores in scores2r]
#         scores2r = [score for subscores in scores2r for score in subscores]
#         append(sum(1 if s1 >= thr1 and s2 >= thr2 else 0 for s1, s2 in zip(scores1r, scores2r)))
#     return(container)


def montecarlo(scores1, scores2, length, thr1, thr2, shift=50, step=100):

    container = list()
    append = container.append
    scores1 = split_scores(scores1, length, shift=0, step=step)

    for i in range(20000):

        scores1r = [random.sample(subscores, len(subscores)) for subscores in scores1]
        scores1r = [score for subscores in scores1r for score in subscores]
        scores1r = split_scores(scores1r, length, shift=shift, step=step)
        scores1r = [random.sample(subscores, len(subscores)) for subscores in scores1r]
        scores1r = [score for subscores in scores1r for score in subscores]
        append(sum(1 if s1 > thr1 and s2 > thr2 else 0 for s1, s2 in zip(scores1r, scores2)))

    return(container)


def mean(l):
    return(float(sum(l)) / float(len(l)))


def stdev(l):
    m = mean(l)
    D = sum((i - m)**2 for i in l) / (len(l) - 1)
    return(math.sqrt(D))


def write_results(meet, mean_r_meet, std_r_meet, z_score, p_val, wpath):
    #Real Rand SD Zscore Log10Pval
    with open(wpath, 'a') as file:
        file.write('{0}\t{1:.4f}\t{2:.4f}\t{3:.4f}\t{4:.4f}\n'.format(meet, mean_r_meet, std_r_meet, z_score, p_val))
    file.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('score1', action='store', help='path to file with score of first model')
    parser.add_argument('score2', action='store', help='path to file with score of second model')
    parser.add_argument('thr1', action='store', type=float, help='threshold for first model')
    parser.add_argument('thr2', action='store', type=float, help='threshold for second model')
    parser.add_argument('length', action='store', help='path to file with length of peaks')
    parser.add_argument('out', action='store', help='path to write results')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    path1 = args.score1
    path2 = args.score2
    thr1 = args.thr1
    thr2 = args.thr2
    length_path = args.length
    wpath = args.out
    scores1 = read_scores(path1)
    scores2 = read_scores(path2)

    meet = sum(1 if s1 > thr1 and s2 > thr2 else 0 for s1, s2 in zip(scores1, scores2))
    length = read_length(length_path)
    scores1 = sort_scores(scores1, length)
    scores2 = sort_scores(scores2, length)
    r_meets = montecarlo(scores1, scores2, length, thr1, thr2, shift=25, step=50)
    mean_r_meet, std_r_meet = mean(r_meets), stdev(r_meets)
    z_score = (meet - mean_r_meet) / std_r_meet
    p_val = stats.norm.sf(abs(z_score))
    #print(meet, mean_r_meet, std_r_meet, z_score, p_val)
    write_results(meet, mean_r_meet, std_r_meet, z_score, p_val, wpath)
    
if __name__=='__main__':
    main()