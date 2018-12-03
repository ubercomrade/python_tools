import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_bed_like_file(path):
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site'])
    return(df)


def read_peaks(path):
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
    if df['name'][0] == '.':
        name = 'peaks_'
        names = [name + str(i) for i in range(len(df))]
        df['name'] = names
    return(df)


def read_fasta_headers(path):
    '''
    Чтение фаста фаила и запись хедеров
    Шапка для FASTA: >uniq_id|chromosome|start-end|strand
    '''
    fasta_headers = list()
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                line = line[1:].strip().split('|')
                record = dict()
                record['name'] = line[0]
                record['chromosome'] = line[1]
                record['start'] = line[2].split('-')[0]
                record['end'] = line[2].split('-')[1]
                try:
                    record['strand'] = line[3]
                except:
                    #print('Record with out strand. Strand is +')
                    record['strand'] = '+'
                continue
            fasta_headers.append(record)
    file.close()
    df = pd.DataFrame(fasta_headers)
    df = df[['chromosome', 'start', 'end', 'name', 'strand']]
    return(df)


def overlap(site_1, site_2):
    '''
    Does the range (start1, end1) overlap with (start2, end2)?
    Based on De Morgan's laws
    '''

    return site_1['end'] >= site_2['start'] and site_2['end'] >= site_1['start']


def get_coord(site):
    return(np.array(site['start'], site['end']))


def calculate_fraction_of_sites(data, bamm_sites, pwm_sites):

    # Varable for write results
    results = {'peaks': int(), 'no_sites': int(), 'only_pwm_sites': int(), 'only_bamm_sites': int(),
               'only_overlap_sites': int(), 'added_sites_by_bamm': int(), 'added_sites_by_pwm': int(),
               'added_sites_by_both_methods': int()}
    results['peaks'] = len(data)
    results['no_sites'] = len(set(data['name']) - set(bamm_sites['name']) - set(pwm_sites['name']))
    results['only_bamm_sites'] = len(set(bamm_sites['name']) - set(pwm_sites['name']))
    results['only_pwm_sites'] = len(set(pwm_sites['name']) - set(bamm_sites['name']))

    # Caclculation fraction of sites in peaks
    both_method_peaks = set(bamm_sites['name']) & set(pwm_sites['name'])
    common_peaks_results = pd.DataFrame()

    # !!!Need to improve!!!
    # + strand
    for peak in both_method_peaks:
        subset_bamm = pd.DataFrame(bamm_sites[np.logical_and(
            bamm_sites['name'] == peak, bamm_sites['strand'] == '+')])
        subset_pwm = pd.DataFrame(pwm_sites[np.logical_and(
            pwm_sites['name'] == peak, pwm_sites['strand'] == '+')])
        for index_pwm, row_pwm in subset_pwm.iterrows():
            for index_bamm, row_bamm in subset_bamm.iterrows():
                intersection = overlap(row_pwm, row_bamm)
                subset_bamm.loc[index_bamm, 'overlaps'] = intersection
                subset_pwm.loc[index_pwm, 'overlaps'] = intersection
        common_peaks_results = pd.concat([common_peaks_results, subset_bamm, subset_pwm])

    # - strand
    for peak in both_method_peaks:
        subset_bamm = pd.DataFrame(bamm_sites[np.logical_and(
            bamm_sites['name'] == peak, bamm_sites['strand'] == '-')])
        subset_pwm = pd.DataFrame(pwm_sites[np.logical_and(
            pwm_sites['name'] == peak, pwm_sites['strand'] == '-')])
        for index_pwm, row_pwm in subset_pwm.iterrows():
            for index_bamm, row_bamm in subset_bamm.iterrows():
                intersection = overlap(row_pwm, row_bamm)
                subset_bamm.loc[index_bamm, 'overlaps'] = intersection
                subset_pwm.loc[index_pwm, 'overlaps'] = intersection
        common_peaks_results = pd.concat([common_peaks_results, subset_bamm, subset_pwm])

    # !!!Need to improve!!!
    # Ended calculation
    for peak in both_method_peaks:
        tmp = common_peaks_results[common_peaks_results['name'] == peak]
        tmp = set(tmp[tmp['overlaps'] == False]['type'])
        if 'pwm' in tmp and 'bamm' in tmp:
            results['added_sites_by_both_methods'] += 1
            continue
        elif 'pwm' in tmp:
            results['added_sites_by_pwm'] += 1
            continue
        elif 'bamm' in tmp:
            results['added_sites_by_bamm'] += 1
            continue
        else:
            results['only_overlap_sites'] += 1

    return(results)


def main():
    # MAIN
    bamm_path = '/home/anton/DATA/TF/CEPBA_mm10_GSM2845732/SCAN/CEBPA_all_BAMM.tsv'
    pwm_path = '/home/anton/DATA/TF/CEPBA_mm10_GSM2845732/SCAN/CEBPA_all_PWM.tsv'
    peaks_path = '/home/anton/DATA/TF/CEPBA_mm10_GSM2845732/BED/CEBPA_all.bed'

    # Read bamm sites
    bamm_sites = read_bed_like_file(bamm_path)
    bamm_sites['type'] = 'bamm'
    bamm_sites['overlaps'] = False

    # Read pwm sites
    pwm_sites = read_bed_like_file(pwm_path)
    pwm_sites['type'] = 'pwm'
    pwm_sites['overlaps'] = False

    # Read peaks and give names to every peask
    data = read_peaks(peaks_path)

    # Total results
    total_results = []

    # Get subdata by peak score
    bins = 7
    bins_width = (max(data['score']) - min(data['score'])) / bins
    fractions = []
    for step in range(bins):
        low = round(min(data['score']) + step * bins_width)
        high = round(min(data['score']) + (step + 1) * bins_width)
        fractions.append(str(low) + '-' + str(high))
        sub_data = data[np.logical_and(data['score'] >= low, data['score'] < high)]
        names = set(sub_data['name'])
        sub_pwm_sites = pwm_sites[[i in names for i in pwm_sites['name']]]
        sub_bamm_sites = bamm_sites[[i in names for i in bamm_sites['name']]]
        total_results.append(calculate_fraction_of_sites(sub_data, sub_bamm_sites, sub_pwm_sites))

    # Count data
    count = pd.DataFrame(total_results)
    count['fractions'] = fractions

    # Fraction data devided by total volume of peaks
    total_count = pd.DataFrame(total_results)
    total_count = total_count.drop(columns=["peaks"])
    total_count['fractions'] = fractions
    for column in total_count:
        if column == 'fractions':
            break
        total_count[column] = total_count[column] / sum(count['peaks'])

    # Fraction data devided by fractiomal volume of peaks
    fractional_count = pd.DataFrame(total_results)
    fractional_count['fractions'] = fractions
    for column in fractional_count:
        if column == 'fractions':
            break
        fractional_count[column] = fractional_count[column] / fractional_count['peaks']
    fractional_count = fractional_count.drop(columns="peaks")

    return(fractional_count, total_count, count)


fractional_count, total_count, count = main()


ax = fractional_count.plot.bar(stacked=True, x='fractions')
ax.legend(loc=[1.1, 0.45])
fig = ax.get_figure()
fig.savefig('/home/anton/DATA/TF/CEPBA_mm10_GSM2845732/FRACTIONS.JPG',
            dpi=150,  bbox_inches='tight')

ax = total_count.plot.bar(stacked=True, x='fractions')
ax.legend(loc=[1.1, 0.45])
fig = ax.get_figure()
fig.savefig('/home/anton/DATA/TF/CEPBA_mm10_GSM2845732/TOTAL_FRACTIONS.JPG',
            dpi=150,  bbox_inches='tight')

count.to_csv('/home/anton/DATA/TF/CEPBA_mm10_GSM2845732/count.tsv', sep='\t', index=False)
