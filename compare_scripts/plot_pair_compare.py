import argparse
import numpy as np
import csv
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import itertools

def read_counts(path):
    with open(path) as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter='\t')
        for row in csvreader:
            out = row
    
    for k in out.keys():
        out[k] = int(out[k])
    return(out)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('dir', action='store', help='dir with compare data')
    parser.add_argument('id', action='store', help='id for data')
    parser.add_argument('tools', action='store', choices=['PWM', 'BAMM', 'INMODE', 'SITEGA'],
                        nargs='+', help='list of used tools')
    #metavar='N'
    return(parser.parse_args())


def main():
    args = parse_args()
    main_path = args.dir
    ID = args.id
    tools = args.tools
    
    
    number_of_tools = len(tools)
    matrix = []
    coordinates = list(itertools.product(range(1, number_of_tools + 1), repeat=2))
    acc = 0
    for i in range(1, number_of_tools + 1):
        matrix.append(coordinates[acc:acc + number_of_tools])
        acc += number_of_tools

    coordinates = []
    for index, i in enumerate(matrix):
        coordinates += (i[:len(i) - index - 1])


    specs = []
    append = specs.append
    for i in range(number_of_tools):
        append([{'type':'domain'}] * number_of_tools)

    
    ftools = tools
    rtools = ftools[::-1]
    
    fig = make_subplots(rows=number_of_tools, cols=number_of_tools, specs=specs, column_titles=ftools[:-1], row_titles=rtools[:-1],
                   start_cell='top-left', horizontal_spacing = 0.4 / number_of_tools, vertical_spacing = 0.4 / number_of_tools)

    labels = ['Peaks with only first model sites',
              'Peaks with only second model sites',
              'Peaks with overlapped sites of models',
             'Peaks with both model sites but not overlapped',
             'Peaks with out sites']
    colors = ['#F66D44', '#FEAE65', '#E6F69D', '#AADEA7', '#64C2A6']

    for i, j in coordinates:
        tool1, tool2 = ftools[i - 1], rtools[j - 1]
        print(tool1, tool2, i, j)
        data = read_counts(main_path + '{0}_{1}.{2}_counts.tsv'.format(ID, tool1, tool2))
        print(data)
        #labels = list(data.keys())
        v = list(data.values())
        vals = np.array([v[0], v[1], v[2], v[5], v[6]])
        fig.add_trace(go.Pie(labels=labels, values=vals,
            name=str(index), sort=False,
            marker_colors=colors), j, i)


    fig.update_layout(height=700, width=700, legend_orientation="h",
                     font=dict(
                         family="Courier New, monospace",
                         size=14),
                      legend=dict(x=0.5, y=0.45))

    for i in fig['layout']['annotations']:
        if i['yanchor'] == 'middle':
            #i['x'] = i['x'] - 0.1
            i['x'] = -.12
        else:
            i['y'] += .07

    fig.write_image(main_path + '/compare_plot.pdf')
#    fig.show()

if __name__=="__main__":
    main()

