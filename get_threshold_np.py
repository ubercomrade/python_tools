import os
import math
import random
import itertools
import sys
import argparse
import time
import numpy as np



'''
Чтение PCM из фаила загруженного из БД Hocomoco
[A,C,G,T]
'''

def read_matrix(path, file_format):
    '''
    Чтение PCM из фаила из разных источников
    v 0.1 (чтение из фаила HOCOMOCO)
    '''
    if file_format == 'HOCOMOCO':
        with open(path) as file:
            file.readline()
            matrix = [list(map(float, line.strip().split())) for line in file]
            matrix_dict = {'A':list(), 'C':list(), 'G':list(), 'T':list()}
            for i in matrix:
                matrix_dict['A'].append(i[0])
                matrix_dict['C'].append(i[1])
                matrix_dict['G'].append(i[2])
                matrix_dict['T'].append(i[3])
    else:
        pass
    return (matrix_dict)



def read_sites(path, every_str=False):
    '''
    Чтение мотивов из фаила, если это просто список строк с разделителем \n, то everyStr = True,
    если это fasta формат, то everyStr = False
    '''
    sequences = []
    if every_str:
        with open(path, 'r') as file:
            sequences = [i.strip().upper() for i in file]
    else:
        with open(path, 'r') as file:
            sequences = [i.strip().upper() for i in file if i.strip() != '>']
    return(sequences)


def remove_equalent_seq(seq_list, homology=0.95):
    '''
    Удаление гомологичных последовательностей из списка (seq_list)
    Если кол-во совпадений при сравнении последовательности 1 и 2 >= длина последовательности * homology,
    то последовательность 1 удаляется из списка
    Функция возвращает новый список
    '''
    seq_list = list(seq_list)
    treshold = homology * len(seq_list[0])
    for seq1 in tuple(seq_list):
        sub_seq_list = list(seq_list)
        sub_seq_list.remove(seq1)
        for seq2 in sub_seq_list:
            score = len([i for i, j in zip(seq1, seq2) if i == j])
            if score >= treshold:
                seq_list.remove(seq1)
                break
    return(seq_list)



def make_pfm_from_pcm(pcm, kind, pseudocount= '1/N'):
    '''
    Вычисление частотной матрицы на основе PCM.
    Для того чтобы избавиться от 0 значений частот используется pseudocount.
    Pseudocount может быть dict со стандартными значениями {'A':0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25} [1],
    либо pseudocount может быть str со значением sqroot [2].
    Подробнее о расчетах смотри:
    1)Wyeth W.Wasserman and Albin Sandelin
      APPLIED BIOINFORMATICS FOR THE IDENTIFICATION OF REGULATORY ELEMENTS
      doi:10.1038/nrg1315

    2)Victor G Levitsky
      Effective transcription factor binding site prediction using a
      combination of optimization, a genetic algorithm and discriminant
      analysis to capture distant interactions
      doi:10.1186/1471-2105-8-481

    В любых других условиях функция ничего не возвращает

    '''

    number_of_sites = int()
    for i in pcm.keys():
        number_of_sites += pcm[i][0]

    if kind == 'di':
        pfm = dict()
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            pfm[''.join(i)] = []
    elif kind == 'mono':
        pfm = dict()
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            pfm[i[0]] = []
    else:
        print('ALARM!')
        pass

    if pseudocount == '1/N':
        first_key = list(pcm.keys())[0]
        nuc_pseudo = 1/len(pcm.keys())
        for i in range(len(pcm[first_key])):
            for nuc in pcm.keys():
                pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites + 1))
        return(pfm)

    elif pseudocount == 'sqroot':
        total_sq_root = int()
        for i in pcm.keys():
            total_sq_root += pcm[i][0]
        total_sq_root = math.sqrt(total_sq_root)
        sq_root = total_sq_root/len(pcm.keys())

        first_key = list(pcm.keys())[0]
        for i in range(len(pcm[first_key])):
            for nuc in pcm.keys():
                pfm[nuc].append((pcm[nuc][i] + sq_root) / (number_of_sites + total_sq_root))

        return(pfm)
    else:
        print('ALARM!')
        pass


def make_pwm_from_pcm(pcm, kind, background, method='log-odds', pseudocount='1/N'):
    '''
    Функиця, которая считает PWM (position weight matrix) на основе PCM (position count matrix)
    с преобразованием log-odds (добавить новые)

    Ref:
    1)Wyeth W.Wasserman and Albin Sandelin
      APPLIED BIOINFORMATICS FOR THE IDENTIFICATION OF REGULATORY ELEMENTS
      doi:10.1038/nrg1315

    2)Victor G Levitsky
      Effective transcription factor binding site prediction using a
      combination of optimization, a genetic algorithm and discriminant
      analysis to capture distant interactions
      doi:10.1186/1471-2105-8-481

    3)Oliver D. King and Frederick P. Roth
      A non-parametric model for transcription factor binding sites
      doi: 10.1093/nar/gng117

    '''
    if kind == 'di':
        pwm = {}
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            pwm[''.join(i)] = []
    elif kind == 'mono':
        pwm = {}
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            pwm[i[0]] = []
    else:
        print('ALARM!')
        pass

    pfm = make_pfm_from_pcm(pcm, kind, pseudocount)
    first_key = list(pcm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def make_pcm(motifs, kind):
    '''
    input - список мотивов одинаковой длины
    output -  PCM
    kind is type of matrix di or mono

    Создает PCM на основе списка мотивов
    '''
    if kind == 'di':
        matrix = {}
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            matrix[''.join(i)] = []
    elif kind == 'mono':
        matrix = {}
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            matrix[i[0]] = []

    len_of_motif = len(motifs[0])

    if kind == 'di':
        for i in matrix.keys():
            matrix[i] = [0]*(len_of_motif - 1)

        for i in range(len_of_motif - 1):
            for l in motifs:
                matrix[l[i:i+2]][i] += 1
    elif kind == 'mono':
        for i in matrix.keys():
            matrix[i] = [0]*len_of_motif

        for i in range(len_of_motif):
            for l in motifs:
                matrix[l[i]][i] += 1

    return(matrix)



def background_freq(seq, kind):

    s = ''.join(seq)
    if kind == 'mono':
        background = {}
        mono_nucleotides = itertools.product('ACGT', repeat=1)
        for i in mono_nucleotides:
            background[i[0]] = s.count(i[0])

    elif kind == 'di':
        background = {}
        di_nucleotides = itertools.product('ACGT', repeat=2)
        for i in di_nucleotides:
            background[''.join(i)] = s.count(i[0])

    sum_of_nuc = sum(background.values())
    for i in background.keys():
        background[i] = background[i]/sum_of_nuc
    return(background)


def to_score(norm_value, min_score, max_score):
    '''
    norm = (score - min) / (max - min) -> score = norm * (max - min) + min
    '''
    score = norm_value * (max_score - min_score) + min_score
    return(score)


def to_norm(score, min_score, max_score):
    norm_value = (score - min_score) / (max_score - min_score)
    return(norm_value)


def min_score(pwm):
    '''
    Вичисляет минимальное значение score для матрицы
    '''
    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += min(tmp)
    return(value)


def max_score(pwm):
    '''
    Вичисляет минимальное значение score для матрицы
    '''
    value = int()
    keys = list(pwm.keys())
    length_pwm = len(pwm[keys[0]])
    for i in range(length_pwm):
        tmp = []
        for j in keys:
            tmp.append(pwm[j][i])
        value += max(tmp)
    return(value)


def score(seq, pwm, kind):
    '''
    Вспомагательная функция, считает score для строки с такой же длиной как и PWM
    kind - тип PWM mono or di
    '''
    if kind == 'mono':
        length_of_seq = len(seq)
        position = 0
        score = 0
        for letter in seq:
            score += pwm[letter][position]
            position += 1
        return(score)

    elif kind == 'di':
        length_of_seq = len(seq)
        score = 0
        for i in range(len(seq) - 1):
            two_letters = seq[i:i+2]
            score += pwm[two_letters][i]
        return(score)
    else:
        pass


def slice_matrix(matrix, start, end):
    '''
    Возвращает срез матрицы
    '''
    sliced_matrix = dict()
    for key in matrix.keys():
        sliced_matrix[key] = matrix[key][start:end]
    return(sliced_matrix)


def possible_scores(pwm):
    scores = [[]]
    letters = list(pwm.keys())
    for letter in letters:
        scores[0].append(pwm[letter][0])
    for i in range(1, len(pwm[letters[0]])):
        tmp = []
        for letter in letters:
            for j in scores[i - 1]:
                tmp.append(j + pwm[letter][i])
        scores.append(tmp)
    return(scores[-1])


def score_distribution_np(matrix, alpha, beta, ndigits=20,  background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
    matrix_length = len(matrix['A'])
    Q_old = [1]
    possible_scores_old = [0]
    for position in range(matrix_length):
        Q_new = []
        possible_scores_new = []
        matrix_slice = [matrix[i][position] for i in matrix.keys()]
        t = np.array([score + weight for score in possible_scores_old for weight in matrix_slice])
        t = np.around(t, decimals=ndigits)
        bs = max_score(slice_matrix(matrix, position + 1, matrix_length))
        ws = min_score(slice_matrix(matrix, position + 1, matrix_length))
        conditions = np.logical_and(alpha - bs <= t, t <= beta - ws)
        Q_new = np.array([Q * background[letter] for Q in Q_old for letter in background.keys()])[conditions]
        possible_scores_new = t[conditions]
        Q_old = np.array(Q_new)
        possible_scores_old = np.array(possible_scores_new)
    return(Q_old, possible_scores_old)


def distribution_np_to_dict(Q, scores):
    keys = np.unique(scores)
    Q_out = dict()
    for key in keys:
        Q_out[key] = np.sum(Q[scores == key])
    return(Q_out)


def fast_pvalue_np(matrix, alpha, ndigits=20, background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
    matrix_length = len(matrix['A'])
    Q_old = [1]
    possible_scores_old = [0]
    P = 0.0
    for position in range(matrix_length):
        Q_new = []
        possible_scores_new = []
        matrix_slice = [matrix[i][position] for i in matrix.keys()]
        t = np.array([score + weight for score in possible_scores_old for weight in matrix_slice])
        t = np.around(t, decimals=ndigits)
        bs = max_score(slice_matrix(matrix, position + 1, matrix_length))
        ws = min_score(slice_matrix(matrix, position + 1, matrix_length))
        conditions_1 = (alpha - ws) <= t
        P += np.sum(np.array([Q * background[letter] for Q in Q_old for letter in background.keys()])[conditions_1])
        conditions_2 = np.logical_and(np.logical_not(alpha - ws <= t), alpha - bs <= t)
        Q_new = np.array([Q * background[letter] for Q in Q_old for letter in background.keys()])[conditions_2]
        possible_scores_new = t[conditions_2]
        Q_old = np.array(Q_new)
        possible_scores_old = np.array(possible_scores_new)
    return(P)


def pwm_shuffling(pwm):
    shuffled_pwm = dict()
    score_lists = list(zip(*pwm.values()))
    random.shuffle(score_lists)
    score_lists =  list(zip(*score_lists))
    keys = list(pwm.keys())
    for i in range(len(score_lists)):
        shuffled_pwm[keys[i]] = score_lists[i]
    return(shuffled_pwm)


def round_pwm(pwm, ndigits):
    keys = list(pwm.keys())
    output_pwm = {}
    for key in keys:
        output_pwm[key] = [round(i, ndigits=ndigits) for i in pwm[key]]
    return(output_pwm)


def maximal_error_associated(pwm, rounded_pwm):
    E = 0
    letters = list(pwm.keys())
    length = len(pwm[letters[0]])
    for i in range(length):
        tmp = []
        for letter in letters:
            tmp.append(abs(abs(pwm[letter][i]) - abs(rounded_pwm[letter][i])))
        E += max(tmp)
    return(E)


def search_score(alpha, rounded_pwm, Q, error, ndigits, background):
    output = []
    tmp_pvalue = fast_pvalue_np(rounded_pwm, alpha+error, ndigits=ndigits, background=background)
    score_list = list(Q.keys())
    score_list.sort(reverse=False)
    if len(score_list) == 1:
        return(score_list[0])
    elif len(score_list) == 0:
        return(alpha)
    else:
        for delta in score_list:
            cumulative_Q = 0
            for i in score_list:
                if i < delta + error and i > delta:
                    cumulative_Q += Q[i]
            if cumulative_Q + fast_pvalue_np(rounded_pwm, delta + error, ndigits, background=background) >= tmp_pvalue:
                return(delta)
        return(alpha)


def from_pvalue_to_score(pwm, pvalue, ndigits, background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}, digit=1):
    rounded_pwm = round_pwm(pwm, ndigits)
    error = maximal_error_associated(pwm, rounded_pwm)
    Q, scores = score_distribution_np(rounded_pwm, alpha=min_score(rounded_pwm), beta=max_score(rounded_pwm), ndigits=ndigits, background=background)
    Q = distribution_np_to_dict(Q, scores)
    score_list = list(Q.keys())
    score_list.sort(reverse=True)
    alpha = score_list[0]
    cumulative_Q = Q[alpha]
    step = 1
    while not cumulative_Q >= pvalue:
        alpha = score_list[step]
        cumulative_Q += Q[alpha]
        step += 1
    pvalue_out = fast_pvalue_np(rounded_pwm, alpha, ndigits=ndigits, background=background)
    while fast_pvalue_np(rounded_pwm, alpha - error, ndigits=20, background=background) != fast_pvalue_np(rounded_pwm, alpha, ndigits=20, background=background):
        ndigits += digit
        rounded_pwm = round_pwm(pwm, ndigits)
        error = maximal_error_associated(pwm, rounded_pwm)
        Q, scores = score_distribution_np(rounded_pwm, alpha=(alpha - error), beta=(alpha + error), ndigits=ndigits, background=background)
        Q = distribution_np_to_dict(Q, scores)
        alpha = search_score(alpha, rounded_pwm, Q, error, ndigits, background)
        pvalue_out = fast_pvalue_np(rounded_pwm, alpha, ndigits=ndigits, background=background)
    return(alpha, pvalue_out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to HOCOMOCO matrix file')
    parser.add_argument('-f', '--fasta', action='store', dest='fasta',
                        required=False, help='path to BED file, needed to calculate backgroun frequences for nucleotieds. \
                        Without this parametr background frequances = {A: 0.25, C: 0.25, G: 0.25, T: 0.25}')
    parser.add_argument('-p', '--pvalue', action='store', dest='pvalue',
                        required=True, help='pvalue')
    parser.add_argument('-r', '--round', action='store', dest='ndigits',
                        required=False, default=int(1),
                        help='Initional matrix rounding in algorithm')


    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    pvalue = float(args.pvalue)
    path = args.input
    fasta = args.fasta
    kind = 'mono'
    ndigits = int(args.ndigits)
    #if ndigits is None:
        #ndigits = 1

    #if os.path.isfile(fasta):
    if not (fasta is None):
        start_time = time.time()
        all_sites = read_sites(fasta, every_str=False)
        background = background_freq(all_sites, kind=kind)
        pwm = read_matrix(path, 'HOCOMOCO')
        score_value, pvalue_out = from_pvalue_to_score(pwm=pwm, pvalue=pvalue, ndigits=ndigits, background=background)
        with_out_round_pvalue = fast_pvalue_np(pwm, score_value, background=background)
        end_time = time.time()
        print('\nSCORE = {0}\nPVALUE = {1}\nPVALUE_WITHOUT_ROUNDING = {2}'.format(score_value, pvalue_out, with_out_round_pvalue))
        print('TIME_OF_CALCULATION = {0} second'.format(end_time - start_time))

    else:
        start_time = time.time()
        pwm = read_matrix(path, 'HOCOMOCO')
        score_value, pvalue_out = from_pvalue_to_score(pwm=pwm, pvalue=pvalue, ndigits=ndigits)
        with_out_round_pvalue = fast_pvalue_np(pwm, score_value)
        end_time= time.time()
        print('\nSCORE = {0}\nPVALUE = {1}\nPVALUE_WITHOUT_ROUNDING = {2}'.format(score_value, pvalue_out, with_out_round_pvalue))
        print('TIME_OF_CALCULATION = {0} second\n'.format(end_time - start_time))

