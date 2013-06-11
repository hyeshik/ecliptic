#!/usr/bin/env python
from __future__ import division
from StringIO import StringIO
from ecliptic.FormattedParsers import BindingSiteCatalogParser
from operator import attrgetter
from collections import Counter
from itertools import product, imap
import numpy as np
import math
import gzip
import csv

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from ecliptic.support.plotutils import iwork_colors, adjust_numbers_style

#<ParsedLine chrom='chr1' strand='+' pos=3532312 base='G' depth=10 del_score=0.7 mod_score=0.3 moddel_score=1.0 entropy_score=0.610864302055 t2c_score=0.0 del_fdr=0.0001 mod_fdr=1.0 moddel_fdr=0.0001 entropy_fdr=1.0 t2c_fdr=1.0 genome_seq='TGTTTGCCACCCTTGGTGTGGTAGTCACGTCCCTTTAACAAGTGTGGGTCGGAGCCAAGCAGCACGAACTGCAGGTGAACTCCATACAGCTACTTTACTAC' transcript_seq='' transcript_acc='' transcript_pos=None>


def get_center_seq(s, width):
    assert len(s) % 2 == 1
    return s[len(s) // 2 - width:len(s) // 2 + width + 1].replace('T', 'U')


def load_sequences(inp, scoretype, mindepth, minscore, seqwidth, is_transcriptome):
    getscore = attrgetter(scoretype + '_score')
    getfdr = attrgetter(scoretype + '_fdr')

    for i, row in enumerate(BindingSiteCatalogParser(inp)):
        if row.strand == '-':
            continue
        if row.depth >= mindepth and getscore(row) >= minscore:
        #if getfdr(row) < 0.0011:
            if not is_transcriptome:
                yield get_center_seq(row.genome_seq, seqwidth)
            elif row.transcript_seq:
                yield get_center_seq(row.transcript_seq, seqwidth)


class NMerCounter(object):

    def __init__(self, pattern_width, input_width):
        self.pattern_width = pattern_width
        self.input_width = input_width
        self.counter = Counter()
        self.input_left_indices = range(input_width - pattern_width + 1)

    def update(self, seq):
        assert len(seq) == self.input_width
        seq = seq.replace('U', 'T')
        subpatterns = [seq[i:i+self.pattern_width] for i in self.input_left_indices]
        self.counter.update(subpatterns)

    def summarized_table(self):
        total_combinations = 4 ** self.pattern_width
        total_incidents = max(sum(self.counter.itervalues()), 1)
        expected_even = math.log((total_incidents + 1) / (total_combinations + 1), 2)
        result = []

        for comb in imap(''.join, product(*['ACGT'] * self.pattern_width)):
            cnt = self.counter[comb]
            result.append({
                'seq': comb,
                'count': cnt,
                'freq': cnt / total_incidents * len(self.input_left_indices),
                'logenrichment': math.log(cnt + 1, 2) - expected_even,
            })

        enrichment_levels = [r['logenrichment'] for r in result]
        enrichment_avg = np.mean(enrichment_levels)
        enrichment_sd = np.std(enrichment_levels)

        for r in result:
            r['z'] = (r['logenrichment'] - enrichment_avg) / enrichment_sd

        result.sort(key=lambda x: x['z'], reverse=True)

        return result


def write_csv(options, result):
    w = csv.writer(open(options.csv_filename, 'w'))
    w.writerow(['Sequence', 'Count', 'Frequency (/site)', 'Log2 enrichment (z)'])

    for row in result:
        w.writerow([row['seq'], row['count'], row['freq'], row['z']])


def write_plot(result, pdf=None, png=None):
    fig = plt.figure(figsize=(3.8, 2))

    plt.hist([r['z'] for r in result], bins=20, facecolor='#444444', rwidth=0.8)

    plt.xlabel('Enrichment (z)')
    plt.ylabel('Number of seq.')
    plt.axvline(1.644854, color='red') # top 0.05 line

    adjust_numbers_style(plt.gca())

    plt.tight_layout()

    if pdf is not None:
        plt.savefig(pdf, format='pdf')

    if png is not None:
        plt.savefig(png, format='png')

    plt.clf()
    plt.cla()


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Generate list of enriched nmers in '
                                                 'flanking sequences of crosslinked sites')
    parser.add_argument('catalog', metavar='FILE', type=str,
                        help='Input binding site catalog file.')
    parser.add_argument('--nmer-width', metavar='WIDTH', type=int, required=True,
                        dest='nmer_width', help='Size of nmers to be tested')
    parser.add_argument('--score-type', metavar='TYPE', type=str, required=True,
                        dest='score_type', help='One of del, moddel, mod, entropy, or t2c')
    parser.add_argument('--minimum-depth', metavar='READS', type=int, required=True,
                        dest='minimum_depth', help='Minimum read count for reporting')
    parser.add_argument('--minimum-score', metavar='VALUE', type=float, required=True,
                        dest='minimum_score', help='Score cut-off for reporting')
    parser.add_argument('--flanking-window', metavar='WIDTH', type=int, required=True,
                        dest='flanking_window', help='Window width of weblogo')
    parser.add_argument('--seq-source', metavar='TYPE', type=str, required=True,
                        dest='source', help='One of genome or transcriptome')
    parser.add_argument('--output-csv', metavar='FILE', type=str,
                        dest='csv_filename', help='Path to CSV file for tabular output')
    parser.add_argument('--output-png', metavar='FILE', type=str,
                        dest='png_filename', help='Path to PNG file for plot output')
    parser.add_argument('--output-pdf', metavar='FILE', type=str,
                        dest='pdf_filename', help='Path to PDF file for plot output')

    return parser.parse_args()


if __name__ == '__main__':
    options = parse_arguments()

    seqiter = load_sequences(gzip.open(options.catalog),
                options.score_type, options.minimum_depth, options.minimum_score,
                options.flanking_window, options.source == 'transcriptome')

    counter = NMerCounter(options.nmer_width, options.flanking_window*2+1)

    for seq in seqiter:
        counter.update(seq)

    table_summary = counter.summarized_table()
    if options.csv_filename:
        write_csv(options, table_summary)

    if options.png_filename or options.pdf_filename:
        write_plot(table_summary, pdf=options.pdf_filename, png=options.png_filename)
