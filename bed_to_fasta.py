import linecache
import argparse
import sys
import os
import pickle


class BioRecord:

    def __init__(self):
        self.recordName = None
        self.chromosome = None
        self.start = None
        self.end = None
        self.sequence = None
        self.strand = None
        self.length = None

    def setChromosome(self, chromosome):
        self.chromosome = chromosome

    def setStartPosition(self, start):
        try:
            start = int(start)
            if start >= 0:
                self.start = start
            else:
                raise Exception()
        except:
            print('start not numeric or start < 0!')
            self.start = None

    def setEndPosition(self, end):
        try:
            end = int(end)
            if end >= 0:
                self.end = end
            else:
                raise Exception()
        except:
            print('start not numeric or end < 0!')
            self.start = None

    def setName(self, name):
        self.recordName = name

    def setStrand(self, strand):
        try:
            strand = str(strand)
            if strand == '+' or strand == '-':
                self.strand = strand
            else:
                self.strand = None
                raise Exception('strand can be + or -')
        except:
            print('strand can be + or -')
            self.strand = None

    def setSequence(self, sequence):
        try:
            sequence = str(sequence)
            sequence = sequence.upper()
            control = {'A', 'C', 'G', 'T'}
            if len(set(list(sequence)) - control) != 0:
                raise Exception()
            else:
                self.sequence = sequence
        except:
            print('sequence contains wrong letter')
            self.sequence = None

    def setLength(self):
        if (not (self.start is None)) and (not (self.end is None)):
            self.length = self.end - self.start
        else:
            self.lenth = None

    def getChromosome(self):
        return(self.chromosome)

    def getSequence(self):
        return(self.sequence)

    def getStrand(self):
        return(self.strand)

    def getStart(self):
        return(self.start)

    def getEnd(self):
        return(self.end)

    def getPositions(self):
        try:
            if self.start != None and self.end != None:
                return([self.start, self.end])
            else:
                raise Exception()
        except:
            print('Start or End is None')
            return(None)

    def getName(self):
        return(self.recordName)

    def getLength(self):
        return(self.length)

    def __str__(self):
        return(self.recordName)

    def __repr__(self):
        return(self.recordName)


def read_bed3(path):
    output = list()
    with open(path, 'r') as file:
        counter = 1
        for line in file:
            line = line.strip().split()
            record = BioRecord()
            record.setName('Sequence_' + str(counter))
            record.setChromosome(line[0])
            record.setStartPosition(line[1])
            record.setEndPosition(line[2])
            record.setLength()
            output.append(record)
            counter += 1
    file.close()
    return(output)


def minimal_length(bio_records):
    '''
    Получить минимальную длину последовательности
    '''
    bio_copy = list(bio_records)
    bio_copy = sorted(bio_copy, key=lambda i: i.getLength(), reverse=False)
    length = bio_copy[0].getLength()
    return(length)


def maximal_length(bio_records):
    '''
    Получить максимальную длину последовательности
    '''
    bio_copy = list(bio_records)
    bio_copy = sorted(bio_copy, key=lambda i: i.getLength(), reverse=True)
    length = bio_copy[0].getLength()
    return(length)


def slice_to_size(bio_records, length):
    '''
    Привести все записи к единой длине
    '''
    bio_copy = list(bio_records)
    for i in bio_copy:
        if i.getLength() == length:
            continue
        else:
            start = i.getStart()
            end = i.getEnd()
            mid = (end - start) // 2 + start
            new_start = mid - length // 2
            new_end = mid + (length - length // 2)
            i.setStartPosition(new_start)
            i.setEndPosition(new_end)
    return(bio_copy)


def add_tail(bio_records, tail_length):
    bio_copy = list(bio_records)
    for i in bio_copy:
        start = i.getStart()
        end = i.getEnd()
        new_start = start - tail_length
        new_end = end + tail_length
        i.setStartPosition(new_start)
        i.setEndPosition(new_end)
    return(bio_copy)


def stat_of_fasta_file(path):
    '''
    start line == 1
    '''
    inf_data = path + '.pickle'
    flag = os.path.isfile(inf_data)
    if flag:
        with open(inf_data, 'rb') as file:
            output = pickle.load(file)
        file.close()
        return output
    else:
        stat = {}
        length = int()
        with open(path, 'r') as file:
            count = 0
            for line in file:
                line = line.strip().split()[0]
                if line[0] == '>':
                    stat[line[1:]] = count + 2
                else:
                    if count == 2:
                        length = len(line)
                count += 1
        file.close()
        output = (stat, length)
        with open(inf_data, 'wb') as file:
            pickle.dump(output, file)
        file.close()
        return output


def modify_bio_records(bio_records, to_min, to_max, to_size, tail):
    '''
    Модифицирование длины в зависимости от параметров
    '''
    bio_copy = list(bio_records)
    if to_min:
        min_length = minimal_length(bio_copy)
        bio_copy = slice_to_size(bio_copy, min_length)
        return(bio_copy)
    elif to_max:
        max_length = maximal_length(bio_copy)
        bio_copy = slice_to_size(bio_copy, max_length)
        return(bio_copy)
    elif not (to_size is None):
        to_size = int(to_size)
        bio_copy = slice_to_size(bio_copy, to_size)
        return(bio_copy)
    elif not (tail is None):
        tail = int(tail)
        bio_copy = add_tail(bio_copy, tail)
        return(bio_copy)
    else:
        return(bio_copy)


def bed_to_fasta(path_fasta, path_bed, to_min, to_max, to_size, tail):
    bed_peaks = read_bed3(path_bed)
    bed_peaks = modify_bio_records(bed_peaks, to_min, to_max, to_size, tail)
    chromosomes, length = stat_of_fasta_file(path_fasta)
    for peak in bed_peaks:
        seq = str()  # container for sequnce of peak
        chr_name = peak.getChromosome()
        try:
            chr_start_line = chromosomes[chr_name]
        except:
            print('{0} was not found in FASTA'.format(chr_name))
            continue
        peak_start_line = peak.getStart() // length + chr_start_line  # number of line where seq started
        peak_end_line = peak.getEnd() // length + chr_start_line  # number of line where seq ended
        position_start = peak.getStart() % length  # position number of nucleotide in first line
        # positionEnd = peaks.getEnd() % length  # position number of nucleotide in last line
        # positionStart = peakStartLine * length - peak.getStart()  # position number of nucleotide in first line
        # position number of nucleotide in last line
        position_end = peak.getEnd() - peak.getStart() + position_start
        for line_number in range(peak_start_line, peak_end_line + 1):
            seq += linecache.getline(path_fasta, line_number).strip()
        peak.setSequence(seq[position_start:position_end])
    linecache.clearcache()
    return(bed_peaks)


def write_fasta(bio_records, path_to_write):
    with open(path_to_write, 'w') as file:
        for record in bio_records:
            if record.getSequence() is not None:
                file.write('>' + record.getName() + '|' + record.getChromosome() + '|' +
                           str(record.getStart()) + '-' + str(record.getEnd()) + '\n')
                file.write(record.getSequence() + '\n')
            else:
                continue
    file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', '--inputFasta', action='store', dest='input_fasta',
                        required=True, help='path to FASTA file')
    parser.add_argument('-bed', '--inputBED', action='store', dest='input_bed',
                        required=True, help='path to BED file')
    parser.add_argument('-of', '--outputFasta', action='store', dest='output_fasta',
                        required=True, help='path to output Fasta')
    parser.add_argument('-m', '--minimize', action='store_true', dest='to_min',
                        default=False, required=False,
                        help='With this parametr all peaks will change their length to the length of the shortest sequence in set')
    parser.add_argument('-M', '--maximaze', action='store_true', dest='to_max',
                        default=False, required=False,
                        help='With this parametr all peaks will change their length to the length of the longest sequence in set')
    parser.add_argument('-l', '--length', action='store', dest='length',
                        default=None, required=False,
                        help='With this parametr all peaks will change their length to the length equal to the given parameter')
    parser.add_argument('-t', '--tail', action='store', dest='tail',
                        default=None, required=False,
                        help='With this parametr all peaks will get tail in start and end of sequence')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    input_fasta = args.input_fasta
    input_bed = args.input_bed
    output_fasta = args.output_fasta
    to_min = args.to_min
    to_max = args.to_max
    to_size = args.length
    tail = args.tail

    results = bed_to_fasta(input_fasta, input_bed, to_min, to_max, to_size, tail)
    write_fasta(results, output_fasta)
