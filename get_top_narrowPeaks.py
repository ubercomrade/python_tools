import argparse
import sys


def read_bed_cistrome(path):
    peaks = []
    with open(path, 'r') as file:
        for line in file:
            line = line.strip().split()
            line[-1] = float(line[-1])
            peaks.append(line)
    file.close()
    return(peaks)


def get_top_peaks(peaks, amount):
    sorted_peaks = sorted(peaks, key=lambda i: i[-1], reverse=True)
    return(sorted_peaks[:amount])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputBED', action='store', dest='input_bed',
                        required=True, help='path to BED file for reading')
    parser.add_argument('-o', '--outputBED', action='store', dest='output_bed',
                        required=True, help='path to BED file for writing')
    parser.add_argument('-a', '--amount', action='store', dest='amount',
                        default=100, required=False,
                        help='Amount of top peaks')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    path = args.input_bed
    output = args.output_bed
    amount = int(args.amount)

    peaks = read_bed_cistrome(path)
    peaks = get_top_peaks(peaks, amount)
    with open(output, 'w') as file:
        for i in peaks:
            i[-1] = str(i[-1])
            file.write('\t'.join(i) + '\n')
    file.close()
