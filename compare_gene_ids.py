import argparse
import sys
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_unweighted


def read_ids(path):
    container = list()
    with open(path, 'r') as file:
        for line in file:
            container.append(line.strip())
    file.close()
    return(set(container))


def write_ids(data, path):
    with open(path, 'w') as file:
        for ID in data:
            file.write('{0}\n'.format(ID))
    file.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-first', '--first_model', action='store', dest='first_model_path',
                        required=True, help='path to file that contain gene ids obtained by first_model in peaks file')
    parser.add_argument('-second', '--second_model', action='store', dest='second_model_path',
                        required=True, help='path to file that contain gene ids obtained by_second_model in peaks file')
    parser.add_argument('-third', '--third_model', action='store', dest='third_model_path',
                        required=True, help='path to file that contain gene ids obtained by_third_model in peaks file')
    parser.add_argument('-o', '--out', action='store', dest='out_path',
                        required=True, help='out_path')
    parser.add_argument('-fname', '--first_name', action='store', dest='fname',
                        required=False, default='First', help='First data name for plot')
    parser.add_argument('-sname', '--second_name', action='store', dest='sname',
                        required=False, default='Second', help='Second data name for plot')
    parser.add_argument('-tname', '--third_name', action='store', dest='tname',
                        required=False, default='Third', help='Third data name for plot')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():

    args = parse_args()
    first_model_path = args.first_model_path
    second_model_path = args.second_model_path
    third_model_path = args.third_model_path
    out_path = args.out_path
    fname = args.fname
    sname = args.sname
    tname = args.tname

    # if not os.path.isdir(out_path):
    #     os.mkdir(out_path)

    fm_ids = read_ids(first_model_path)
    sm_ids = read_ids(second_model_path)
    tm_ids = read_ids(third_model_path)

    plt.figure(figsize=(6, 4), dpi=150)
    venn3_unweighted([fm_ids, sm_ids, tm_ids], (fname, sname, tname))
    plt.savefig(out_path, dpi=150)

    only_fm_ids = fm_ids - sm_ids - tm_ids
    only_sm_ids = sm_ids - fm_ids - tm_ids
    only_tm_ids = tm_ids - fm_ids - sm_ids
    common_ids = fm_ids & sm_ids & tm_ids

    # write_ids(data=only_fm_ids,
    #           path=out_path + '/only.{0}.ids.txt'.format(fname))
    # write_ids(data=only_sm_ids,
    #           path=out_path + '/only.{0}.ids.txt'.format(sname))
    # write_ids(data=only_tm_ids,
    #           path=out_path + '/only.{0}.ids.txt'.format(tname))

if __name__=='__main__':
    main()
