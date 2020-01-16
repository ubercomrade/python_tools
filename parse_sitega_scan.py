'''
Copyright © 2018 Anton Tsukanov. Contacts: tsukanov@bionet.nsc.ru
License: http://www.gnu.org/licenses/gpl.txt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
'''

import argparse
import sys
import re
import pandas as pd


def parse_sitega(path):
    sitega = list()
    length = 30
    with open(path, 'r') as file:
        for line in file:
            #print(line)
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                name = line[0]
                chromosome = line[2]
                coordinates_strand = line[3]
                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                start = int(start)
                end = int(end)

            else:
                record = dict()
                line = line.strip().split()
                site = line[3].upper()
                strand = line[2]
                score = float(line[1])
                center = int(line[0])
                start_site = (center - (length // 2)) + int(start)
                end_site = int(start_site + length)
                    
                record['chr'] = chromosome
                record['start'] = start_site
                record['end'] = end_site
                record['name'] = name
                record['score'] = score
                record['strand'] = strand
                record['site'] = site
                sitega.append(record)
    file.close()
    return(sitega)


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('sitega', action='store',
                        help='path to sitega scan file')
    parser.add_argument('bed', action='store',
                        help='path to write file in bed format')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return(parser.parse_args())


def main():
    args = parse_args()
    path_in = args.sitega
    path_out = args.bed
    
    sitega = parse_sitega(path_in);
    sitega = pd.DataFrame(sitega)
    sitega = sitega[['chr', 'start', 'end', 'name', 'score', 'strand', 'site']]
    sitega.to_csv(path_out, sep="\t", index=False, header=False)
    

if __name__=="__main__":
    main()
