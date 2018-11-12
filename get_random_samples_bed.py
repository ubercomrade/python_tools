import argparse
import random
import sys
import os


def read_bed(path):
    peaks = []
    with open(path, 'r') as file:
        for line in file:
            peaks.append(line)
    file.close()
    return(peaks)


def get_random_peaks(peaks, amount):
    random_peaks = random.sample(peaks, k=amount)
    other_peaks = [peak for peak in peaks if peak not in random_peaks]
    return(random_peaks, other_peaks)


def main(input, output, tag, amount_background, amount_train):
    peaks = read_bed(input)
    train_peaks, peaks = get_random_peaks(peaks, amount_train)
    background_peaks, other_peaks = get_random_peaks(peaks, amount_background)

    background_output = output + '/' + 'BACKGROUND_' + tag + '.bed'
    with open(background_output, 'w') as file:
        for i in background_peaks:
            file.write(i)
    file.close()

    train_output = output + '/' + 'TRAIN_' + tag + '.bed'
    with open(train_output, 'w') as file:
        for i in train_peaks:
            file.write(i)
    file.close()

    other_output = output + '/' + 'OTHER_' + tag + '.bed'
    with open(other_output, 'w') as file:
        for i in other_peaks:
            file.write(i)
    file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputBED', action='store', dest='input_bed',
                        required=True, help='path to BED file for reading')
    parser.add_argument('-o', '--output', action='store', dest='output_bed',
                        required=True, help='Folder for writing output BED files')
    parser.add_argument('--tag', action='store', dest='tag',
                        required=True, help='tag for outpot files')
    parser.add_argument('--train', action='store', dest='amount_train',
                        default=100, required=False,
                        help='Amount of train peaks')
    parser.add_argument('--background', action='store', dest='amount_background',
                        default=100, required=False,
                        help='Amount of background peaks')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    input = args.input_bed
    output = args.output_bed
    tag = args.tag
    amount_background = int(args.amount_background)
    amount_train = int(args.amount_train)

    if os.path.isdir(output):
        main(input, output, tag, amount_background, amount_train)
    else:
        os.mkdir(output)
        main(input, output, tag, amount_background, amount_train)
