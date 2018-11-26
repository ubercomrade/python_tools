import itertools
import argparse
import sys
import math


def read_homer(path):
    with open(path, 'r') as file:
        inf = file.readline()
        inf = inf.strip().split()
        id_pfm = inf[0][1:].strip()
        tf_name = inf[1].split(':')[1].split('/')[0]
        tf_db = inf[1].split(':')[1].split('/')[1]
        inf = {'id_pfm': id_pfm, 'tf_name': tf_name, 'tf_db': tf_db, 'alength': 4}

        pfm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pfm.keys(), line):
                pfm[letter].append(float(value))
        inf['w'] = len(pfm['A'])
    file.close()
    return(pfm, inf)


def make_pcm_from_homer(pfm_homer):
    '''
    input - PFM (HOMER)
    output -  PCM
    Создает PCM на основе PFM
    '''
    matrix = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        matrix[i[0]] = []

    pseudo_count = 10**9
    for key in pfm_homer.keys():
        for i in pfm_homer[key]:
            matrix[key].append(int(i * pseudo_count))
    return(matrix)


def make_pfm_from_pcm(pcm, pseudocount='1/N'):
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

    number_of_sites = [0] * len(pcm['A'])
    for key in pcm.keys():
        for i in range(len(pcm[key])):
            number_of_sites[i] += pcm[key][i]

    pfm = dict()
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pfm[i[0]] = []

    if pseudocount == '1/N':
        first_key = list(pcm.keys())[0]
        nuc_pseudo = 1/len(pcm.keys())
        for i in range(len(pcm[first_key])):
            for nuc in pcm.keys():
                pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
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
                pfm[nuc].append((pcm[nuc][i] + sq_root) / (number_of_sites[i] + total_sq_root))

    return(pfm)


def make_pwm_from_pcm(pcm, background, method='log-odds', pseudocount='1/N'):
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
    pwm = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pwm[i[0]] = []

    pfm = make_pfm_from_pcm(pcm, pseudocount)
    first_key = list(pcm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def make_pwm_from_homer(pfm, background={'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}):
    pwm = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pwm[i[0]] = []
    first_key = list(pfm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def background_freq(seq):
    s = ''.join(seq)
    background = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        background[i[0]] = s.count(i[0])
    sum_of_nuc = sum(background.values())
    for i in background.keys():
        background[i] = background[i]/sum_of_nuc
    return(background)


def read_fasta(path):
    fasta = []
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            else:
                fasta.append(line.strip().upper())
    file.close()
    return(fasta)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', action='store', dest='input',
                        required=True, help='path to HOMER matrix file')
    parser.add_argument('-o', '--output', action='store', dest='output',
                        required=True, help='path to PWM output file')
    parser.add_argument('-f', '--fasta', action='store', dest='fasta',
                        required=False, help='path to BED file, needed to calculate backgroun frequences for nucleotieds. \
                        Without this parametr background frequances = {A: 0.25, C: 0.25, G: 0.25, T: 0.25}')
    return(parser.parse_args())


def main():
    args = parse_args()
    input_homer = args.input
    output_pwm = args.output
    path_fasta = args.fasta

    pfm, inf = read_homer(input_homer)
    pcm = make_pcm_from_homer(pfm_homer=pfm)
    if not (path_fasta is None):
        fasta = read_fasta(path_fasta)
        background = background_freq(fasta, 'mono')
    else:
        background = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    pwm = make_pwm_from_pcm(pcm, background)

    with open(output_pwm, 'w') as file:
        file.write(
            '>{0} BEST_MATCH= {1}/{2} ALPHABET= ACGT\n'.format(inf['id_pfm'], inf['tf_name'], inf['tf_db']))
        for i in zip(pwm['A'], pwm['C'], pwm['G'], pwm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))
    file.close()


if __name__ == '__main__':
    main()
