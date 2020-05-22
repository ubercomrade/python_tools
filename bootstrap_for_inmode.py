import subprocess
import csv
import itertools
import os
import sys
import argparse
import random
import math
import shutil


def read_sites(path):
    container = []
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                append(line.strip().upper())
            else:
                continue
    return(container)


def write_fasta(sites, tmp_dir, tag):
    with open('{0}/{1}'.format(tmp_dir, tag), 'w') as file:
        for index, site in enumerate(sites):
            file.write('>{0}\n{1}\n'.format(index, site))
    return(0)


def calculate_scores(path_to_inmode, path_to_java, motif_length, tmp_dir, tag):
    container = []
    args = [path_to_java, '-Xmx4096m', '-Xms1024m', 
            '--add-modules', 'java.xml.bind', '-jar',
            path_to_inmode, 'scan',
            'i={0}/Learned_DeNovo({1},2,2)_motif/XML_of_DeNovo({1},2,2)_motif.xml'.format(tmp_dir, motif_length),
            'id={}/{}.fa'.format(tmp_dir, tag), 'f=1.0', 'outdir={}'.format(tmp_dir), 'bs=false']
    r = subprocess.call(args)
    with open('{0}/{1}'.format(tmp_dir, "/Motif_hits_from_SequenceScan(1.0).BED")) as file:
        for line in file:
            container.append(math.log(float(line.split()[4]), 10))
    shutil.rmtree(tmp_dir)
    return(container)


def make_inmode(path_to_inmode, path_to_java, motif_length, order, tmp_dir)
    args = [path_to_java, '-Xmx4096m', '-Xms1024m', '--add-modules', 'java.xml.bind', '-jar', path_to_inmode,
    'denovo', 'i={}/train.fa'.format(tmp_dir), 'm={}'.format(motif_length), 'outdir={}'.format(tmp_dir),
    'mo={}'.format(order)]
    r = subprocess.call(args)
    return(0)


def bootstrap_inmode(sites, path_to_inmode, path_to_java, tmp_dir, order, size_of)
    true_scores = []
    false_scores = []
    number_of_sites = len(sites)
    len_of_site = len(sites[0])

    for i in range(10):
        train_sample = random.choices(sites, k=round(0.9 * number_of_sites))
        test_sample = [site for site in sites  if not site in train_sample]
        random_sample = [''.join(random.sample(list(random.choice(test_sample)), len_of_site)) for i in range(len(test_sample) * size_of_random_sample)]
        write_fasta(train_sample, tmp_dir, "train")
        write_fasta(test_sample, tmp_dir, "test")
        write_fasta(random_sample, tmp_dir, "shuffled")
        make_inmode(path_to_inmode, path_to_java, len_of_site, order, tmp_dir)
        for true_score in calculate_scores(path_to_inmode, path_to_java, len_of_site, tmp_dir, "test"):
            true_scores.append(true_score)
        for false_score in calculate_scores(path_to_inmode, path_to_java, len_of_site, tmp_dir, "shuffled"):
            false_scores.append(false_score)
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    table = []
    for tpr in [round(i * 0.01, 2) for i in range(0,105, 5)]:
        print(tpr)
        score = true_scores[round(true_length * tpr) - 1]
        actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
        fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
        table.append({'Scores': score, 'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
    return(table)


def write_table(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = data[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for line in data:
            writer.writerow(line)
    return(0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', action='store', help='path to write table with results')
    parser.add_argument('input', action='store', help='path to fasta file with sites')
    parser.add_argument('inmode', action='store', help='path to inmode program')
    parser.add_argument('len', action='store', type=int, help='len of TF site')
    parser.add_argument('out', action='store', help='path to write results')
    parser.add_argument('-j', '--java', action='store', dest='java',
                            required=False, default="java", help='path to java')
    parser.add_argument('-s', '--size', action='store', type=int, dest='size',
                            default= 1000, required=False, help='randome_sample times more than test_sample')
    parser.add_argument('-o', '--order', action='store', type=int, dest='order',
                            default= 2, required=False, help='order of BaMM model')
    parser.add_argument('-t', '--tmp', action='store', dest='tmp_dir',
                            default= './tmp', required=False, help='dir for writing tmp files')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    path = args.input
    out = args.output
    size_of_random_sample = args.size
    order = args.order    
    path_to_inmode = args.inmode
    path_to_java = args.java
    tmp_dir = args.tmp_dir

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    sites = read_sites(path)
    table = bootstrap_inmode(sites, path_to_inmode, path_to_java, tmp_dir, order, size_of_random_sample)
    write_table(out, table)
    shutil.rmtree(tmp_dir)
    return(0)


if __name__=='__main__':
    main()
