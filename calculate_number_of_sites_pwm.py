import argparse
import sys
import functools
from multiprocessing import Pool


def read_fasta(path):
    fasta = list()
    append = fasta.append
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                line = line.strip().upper()
                append(line)
                append(complement(line))
    file.close()
    return(fasta)


def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(pwm)


def score(seq, pwm):
    position = 0
    score = 0
    for letter in seq:
        score += pwm[letter][position]
        position += 1
    return(score)


def complement(seq):
    seq = seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
    return(seq)


def scan_seq_by_pwm(seq, pwm, length_pwm, threshold):
    upper_threshold = 0
    counter = 0
    for i in range(len(seq) - length_pwm + 1):
        site = seq[i:length_pwm + i]
        if 'N' in site:
            continue
        s = score(site, pwm)
        counter += 1
        if s >= threshold:
            upper_threshold += 1
    return(upper_threshold, counter)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', action='store',
                        help='path to FASTA')
    parser.add_argument('pwm', action='store',
                        help='path to PWM file')
    parser.add_argument('threshold', action='store', type=float,
                        help='threshold for PWM')
    parser.add_argument('-t', '-threads', action='store', type=float, dest='threads',
                        default=2, help='Threads for paralell calculating')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    pwm_path = args.pwm
    fasta_path = args.fasta
    threshold = args.threshold
    threads = args.threads

    fasta = read_fasta(fasta_path)
    pwm = read_pwm(pwm_path)
    length_pwm = len(pwm['A'])
    total_upper_threshold = 0
    number_of_sites = 0
    # for record in fasta:
    #     upper_threshold, add_counter = scan_seq_by_pwm(record, pwm, threshold)
    #     total_upper_threshold += upper_threshold
    #     number_of_sites += add_counter
    with Pool(threads) as p:
        t = functools.partial(scan_seq_by_pwm, pwm=pwm, length_pwm=length_pwm, threshold=threshold)
        results = p.map(t, fasta)
    total_upper_threshold = sum([i[0] for i in results])
    number_of_sites = sum([i[1] for i in results])
    print(total_upper_threshold, number_of_sites, total_upper_threshold / number_of_sites)

if __name__ == '__main__':
    main()
