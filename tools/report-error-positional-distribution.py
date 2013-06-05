#!/usr/bin/env python
from __future__ import division
import csv
import sys
import os
import pickle
from ecliptic.support.plotutils import iwork_colors, adjust_numbers_style
import numpy as np

import matplotlib as mpl
mpl.use('Agg') # to enable plotting without X11
from matplotlib import pyplot as plt

PLOT_FORMATS = ['png', 'pdf']
BASES = 'ACGTX'
BASE_NONEXISTENT = 4
COLFORSTATS = 2 # index in pickled profile. 0: from 5' end  1: from 3' end  2: 12-division

# Make masking matrices for counting match, del, ins or substitutions.
EYE_MATRIX = np.identity(len(BASES), np.uint64)
MASK_MATCH = EYE_MATRIX.copy()
MASK_DEL = np.zeros([len(BASES)] * 2, np.uint64); MASK_DEL[:, BASE_NONEXISTENT] = 1
MASK_INS = np.zeros([len(BASES)] * 2, np.uint64); MASK_INS[BASE_NONEXISTENT, :] = 1
MASK_SUBST = 1 - MASK_MATCH - MASK_DEL - MASK_INS
for arr in (MASK_MATCH, MASK_DEL, MASK_INS, MASK_SUBST):
    arr[BASE_NONEXISTENT, BASE_NONEXISTENT] = 0

def write_tabular_output(options, data, fieldsize):
    w = csv.writer(open(options['csv'], 'w'))
    fieldheaders = [
        '%d %s-to-%s' % (fn + 1, basef, baset)
        for fn in range(fieldsize)
        for basef in BASES for baset in BASES
        if basef != 'X' or baset != 'X'
    ]
    w.writerow(['Sample'] + fieldheaders)

    for sample, _ in options['profiles']:
        profile = data[sample]

        fielddata = [
            profile[COLFORSTATS][fn, i, j]
            for fn in range(fieldsize)
            for i, basef in enumerate(BASES)
            for j, baset in enumerate(BASES)
            if basef != 'X' or baset != 'X'
        ]

        w.writerow([sample] + fielddata)


def plot_error_profiles_total(options, data, fieldsize):
    xpositions = np.arange(fieldsize)

    for sample, _ in options['profiles']:
        profile = data[sample][COLFORSTATS]

        fig = plt.figure(figsize=(3.5, 2.8))

        subst_freq = []
        ins_freq = []
        del_freq = []

        for fn in range(fieldsize):
            inserted = (profile[fn] * MASK_INS).sum()
            deleted = (profile[fn] * MASK_DEL).sum()
            substed = (profile[fn] * MASK_SUBST).sum()
            total = profile[fn].sum()

            subst_freq.append(substed / total * 100)
            ins_freq.append(inserted / total * 100)
            del_freq.append(deleted / total * 100)

        plt.plot(xpositions, subst_freq, label='subst', c=iwork_colors.yellow, linewidth=1.2)
        plt.plot(xpositions, ins_freq, label='ins', c=iwork_colors.red, linewidth=1.2)
        plt.plot(xpositions, del_freq, label='del', c=iwork_colors.blue, linewidth=1.2)

        ax = plt.gca()
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.xlim(0, fieldsize-1)
        plt.xlabel('Position in tag')
        plt.ylabel('Frequency (percent)')
        plt.title(sample)

        box = ax.get_position()
        ax.set_position([box.x0 + 0.1, box.y0 + box.height * 0.2,
                         box.width * 0.9, box.height * 0.77])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=3, columnspacing=0.5,
                  frameon=False, handletextpad=0.5)

        adjust_numbers_style(ax, xgrid=True)

        for format in PLOT_FORMATS:
            plt.savefig(options['plots'].format(format=format, type='total', sample=sample))

        plt.cla()
        plt.clf()


def plot_error_profiles_pertype(options, data, fieldsize, plottype):
    xpositions = np.arange(fieldsize)

    freqmask = {'del': MASK_DEL, 'subst': MASK_SUBST}[plottype]

    for sample, _ in options['profiles']:
        profile = data[sample][COLFORSTATS]

        fig = plt.figure(figsize=(3.5, 2.8))

        freqs = [[] for i in range(len(BASES) - 1)]

        for fn in range(fieldsize):
            profmasked = profile[fn] * freqmask

            for basei, freqtbl in enumerate(freqs):
                basereads = profile[fn, basei].sum()
                freqtbl.append(profmasked[basei].sum() / basereads)

        for base, freqtbl, color in zip(BASES, freqs, iwork_colors):
            plt.plot(xpositions, freqtbl, label=base, c=color, linewidth=1.2)

        ax = plt.gca()
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.xlim(0, fieldsize-1)
        plt.xlabel('Position in tag')
        plt.ylabel('Frequency (percent)')
        plt.title(sample)

        box = ax.get_position()
        ax.set_position([box.x0 + 0.1, box.y0 + box.height * 0.2,
                         box.width * 0.9, box.height * 0.77])
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, columnspacing=0.5,
                  frameon=False, handletextpad=0.5)

        adjust_numbers_style(ax, xgrid=True)

        for format in PLOT_FORMATS:
            plt.savefig(options['plots'].format(format=format, type=plottype, sample=sample))

        plt.cla()
        plt.clf()


def plot_error_profiles(options, data, fieldsize):
    plot_error_profiles_total(options, data, fieldsize)

    for plottype in ('subst', 'del'):
        plot_error_profiles_pertype(options, data, fieldsize, plottype)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Generate a table and plot of '
                                     'error frequency distribution for the report')
    parser.add_argument('profiles', metavar='name:pickle', type=str, nargs='+',
                        help='Input profiles with name.')
    parser.add_argument('--output-csv', dest='csv', metavar='PATH', type=str,
                        help='Tabular output file path')
    parser.add_argument('--output-plot', dest='plots', metavar='PATH', type=str,
                        help='Plot output file path (use {format}, {sample} and {type})')
    args = parser.parse_args()

    return {
        'profiles': [tok.split(':') for tok in args.profiles],
        'csv': args.csv,
        'plots': args.plots,
    }

if __name__ == '__main__':
    options = parse_arguments()

    data = dict((sample, pickle.load(open(filename)))
                for sample, filename in options['profiles'])
    fieldsize = data.itervalues().next()[COLFORSTATS].shape[0]
    write_tabular_output(options, data, fieldsize)
    plot_error_profiles(options, data, fieldsize)

