import random
import argparse
import sys

def read_fasta(path_in, path_out):
    fasta = list()
    append = fasta.append
    fasta_in = open(path_in, 'r')
    fasta_out = open(path_out, 'w')
    for index, line in enumerate(fasta_in):
        if not line.startswith('>'):
            line = line.strip().upper()
            line = clear_n(line)
            fasta_out.write(line + '\n')
        else:
            fasta_out.write('>{0}\n'.format(int(index / 2)))
                
    fasta_in.close()
    fasta_out.close()
    pass


def longest(ss):
    if len(ss[0]) > len(ss[1]):
        return(ss[0])
    else:
        return(ss[1])

    
def clear_n(string):
    while 1:
        position = string.find('N')
        if position == -1:
            break
        elif position == len(string) - 1:
            string = string[:position - 1]
            break
        elif string[position + 1] != 'N':
            string = string[:position] + random.choice('ACGT') + string[position + 1:]
        else:
            for index, n in enumerate(string[position:],position):
                if n != 'N':
                    string = longest([string[:position], string[index:]])
                    break
                elif index == len(string) - 1:
                    string = string[:position]
                    break
    return(string)
    
    
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('fasta_in', action='store', help='path to fasta in')
    parser.add_argument('fasta_out', action='store', help='path to fasta write results')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())
    

def main():
    args = parse_args()

    fasta_in = args.fasta_in
    fasta_out = args.fasta_out
    read_fasta(fasta_in, fasta_out)


if __name__ == '__main__':
    main()
            