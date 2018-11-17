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


if __name__ == '__main__':

    args = parse_args()
    input_homer = args.input
    output_pwm = args.output
    path_fasta = args.fasta

    pfm, inf = read_homer(input_homer)
    if not (path_fasta is None):
        fasta = read_fasta(path_fasta)
        background = background_freq(fasta, 'mono')
    else:
        background = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    pwm = make_pwm_from_homer(pfm, background)

    with open(output_pwm, 'w') as file:
        file.write('>{0}\t{1}/{2}\tALPHABET= ACGT\n'.format(inf['id_pfm'],
                                                            inf['tf_name'],
                                                            inf['tf_db']))
        for i in zip(pwm['A'], pwm['C'], pwm['G'], pwm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))
