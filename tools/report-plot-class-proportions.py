#!/usr/bin/env python
from __future__ import division
import csv
import sys
import numpy as np
from collections import defaultdict

import matplotlib as mpl
mpl.use('Agg') # to generate plots without a X11 display
from matplotlib import pyplot as plt
from ecliptic.support.plotutils import iwork_colors, adjust_numbers_style

def load_original_stats(samples, class_mapping, default_class):
    reads = defaultdict(lambda: defaultdict(int))
    percentages = defaultdict(lambda: defaultdict(float))

    for smpspec in samples:
        sample = smpspec['name']

        for row in csv.DictReader(open(smpspec['file'])):
            for mapped_class in class_mapping[row['class']] or [default_class]:
                reads[mapped_class][sample] += int(row['reads'])
                percentages[mapped_class][sample] += float(row['proportion'])

    return dict(reads=reads, percentages=percentages)


def write_csv(options, stats, output, fmt):
    csvout = csv.writer(output)

    header = ['sample'] + options['table_classes']
    csvout.writerow(header)

    for sample in options['samples']:
        name = sample['name']
        datacolumns = [format(stats[cls][name], fmt) for cls in options['table_classes']]
        csvout.writerow([name] + datacolumns)


def plot_to_file(options, stats, outputs):
    figwidth = 2 + len(options['samples']) * len(options['plot_classes']) / 4

    fig = plt.figure(figsize=(figwidth, 3))

    ax = fig.add_subplot(1, 1, 1)

    class_names = options['plot_classes']
    xstep = 0.8 / len(options['samples'])
    xwidth = xstep * 0.8
    xcenterbase = np.arange(0, len(options['plot_classes']), 1)
    xleftbase = xcenterbase - 0.4

    for sampleno, sample in enumerate(options['samples']):
        name = sample['name']
        xlefts = xleftbase + xstep * sampleno
        data = [stats[cls][name] for cls in options['plot_classes']]
        plt.bar(xlefts, data, width=xwidth, facecolor=iwork_colors[sampleno],
                label=name, edgecolor='#666666', zorder=6)

    plt.xlim(-0.5, len(options['plot_classes']) - 0.5)
    plt.xticks(xcenterbase)
    ax.set_xticklabels(class_names)
    plt.setp(ax.get_xticklines(), visible=False)

    plt.legend(loc='best')
    plt.xlabel('Class')
    plt.ylabel('% Reads')

    if options['branding']:
        plt.annotate(options['branding'], (3, 3), fontsize=8, color='#888888',
            textcoords='figure points', horizontalalignment='left',
            verticalalignment='bottom')

    adjust_numbers_style(ax)
    plt.tight_layout()

    for output in outputs:
        plt.savefig(output)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Generate plots and tables of '
                    'assigned annotation proportions for the final report')
    parser.add_argument('samplestats', metavar='name:csvfile', type=str, nargs='+',
                        help='Input per-sample statistics.')
    parser.add_argument('--plot-classes', dest='plot_classes', action='store',
                        help='list of classes to plot (comma-separated)',
                        required=True)
    parser.add_argument('--table-classes', dest='table_classes', action='store',
                        help='list of classes to put into table (comma-separated)',
                        required=True)
    parser.add_argument('--write-reads-csv', dest='write_reads_csv', action='store',
                    default=None, help='file path to write a tabular output for read counts')
    parser.add_argument('--write-percentages-csv', dest='write_percentages_csv', action='store',
                        default=None, help='file path to write a tabular output for read '
                                           'percentages')
    parser.add_argument('--write-plot', dest='write_plot', action='append',
                        default=[], help='file path to write a bar chart')
    parser.add_argument('--class-mapping', dest='class_mappings', action='append',
                        default=[], help='rename or merge classes into a class')
    parser.add_argument('--default-class', dest='default_class', action='store',
                        default='Others', help='default class name for unspecified classes')
    parser.add_argument('--branding', dest='branding', action='store',
                        default=None, help='text for branding on the left bottom corner')
    args = parser.parse_args()

    samples = [dict(zip(('name', 'file'), token.split(':'))) for token in args.samplestats]
    class_mappings = defaultdict(list)
    for label, sources in [token.split('=') for token in args.class_mappings]:
        for src in sources.split(','):
            class_mappings[src].append(label)

    return {
        'samples': samples,
        'plot_classes': args.plot_classes.split(','),
        'table_classes': args.table_classes.split(','),
        'class_mappings': class_mappings,
        'default_class': args.default_class,
        'write_reads_csv': args.write_reads_csv,
        'write_percentages_csv': args.write_percentages_csv,
        'write_plot': args.write_plot,
        'branding': args.branding,
    }

if __name__ == '__main__':
    options = parse_arguments()

    stats = load_original_stats(options['samples'], options['class_mappings'],
                                options['default_class'])

    if options['write_reads_csv'] is not None:
        opt = options['write_reads_csv']
        outfile = sys.stdout if opt == '-' else open(opt, 'w')
        write_csv(options, stats['reads'], outfile, 'd')

    if options['write_percentages_csv'] is not None:
        opt = options['write_percentages_csv']
        outfile = sys.stdout if opt == '-' else open(opt, 'w')
        write_csv(options, stats['percentages'], outfile, '.3f')

    if options['write_plot']:
        plot_to_file(options, stats['percentages'], options['write_plot'])

