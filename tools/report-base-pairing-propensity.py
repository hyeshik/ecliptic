#!/usr/bin/env python
from __future__ import division
from ecliptic.FormattedParsers import BindingSiteCatalogParser
from operator import attrgetter
from collections import Counter
from itertools import product, imap, ifilter
import random
import futures
import numpy as np
import math
import gzip
import csv
import re
import tempfile
import zipfile

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from ecliptic.support.plotutils import iwork_colors, adjust_numbers_style
from matplotlib import cm


PSEUDOCOUNT_PROB = 0.001


def get_center_seq(s, width):
    assert len(s) % 2 == 1
    return s[len(s) // 2 - width:len(s) // 2 + width + 1].replace('T', 'U')


def load_sequences(inp, scoretype, mindepth, minscore, seqwidth, is_transcriptome):
    getscore = attrgetter(scoretype + '_score')
    getfdr = attrgetter(scoretype + '_fdr')

    for i, row in enumerate(BindingSiteCatalogParser(inp)):
        #if row.strand == '-':
        #    continue
        if row.depth >= mindepth and getscore(row) >= minscore:
        #if getfdr(row) < 0.0011:
            if not is_transcriptome:
                yield get_center_seq(row.genome_seq, seqwidth)
            elif row.transcript_seq:
                yield get_center_seq(row.transcript_seq, seqwidth)


def calculate_background_pairing(sequences):
    compositions = [Counter(bases) for bases in np.array(sequences, 'c').transpose()]

    pairing_cases = np.zeros([len(sequences[0])] * 2, np.uint64)
    for idx1, comp1 in enumerate(compositions):
        for idx2, comp2 in enumerate(compositions):
            if idx2 >= idx1:
                break

            propsum = sum((comp1[b1] * comp2[b2]) for b1, b2 in ('AU', 'UA', 'CG', 'GC'))
            pairing_cases[idx1, idx2] = propsum

    pairing_p = pairing_cases / (len(sequences) ** 2)

    return pairing_p


def calculate_actual_pairing(sequences):
    pairing_cases = np.zeros([len(sequences[0])] * 2, np.uint64)
    pairing = set(['AU', 'UA', 'CG', 'GC'])

    for seq in sequences:
        for idx1, base1 in enumerate(seq):
            for idx2, base2 in enumerate(seq):
                if idx2 >= idx1:
                    break

                if base1 + base2 in pairing:
                    pairing_cases[idx1, idx2] += 1
                    #pairing_cases[idx2, idx1] += 1

    pairing_p = pairing_cases / len(sequences)

    return pairing_p


def generate_basepairing_plot(options, allsequences, output, preview=False):
    # calculate basepairing propensity
    seqwidth = len(allsequences[0])
    background_pairing = calculate_background_pairing(allsequences)
    actual_pairing = calculate_actual_pairing(allsequences)
    pairing_enrichment = np.log2((actual_pairing + PSEUDOCOUNT_PROB)
                                 / (background_pairing + PSEUDOCOUNT_PROB))

    # draw it
    fig = plt.figure(figsize=((6, 6) if preview else (7, 9)))
    ax = fig.add_subplot(1, 1, 1, aspect=1)
    plt.pcolor(pairing_enrichment, cmap=cm.BrBG_r, vmin=-1, vmax=1)

    plt.xlim(0, seqwidth-1)
    plt.ylim(seqwidth, 1)

    # hide top and right axes
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()

    # format ticks
    tickinterval = 5
    tickmin = -options.flanking_window + (tickinterval -
                    (-options.flanking_window % tickinterval)) % tickinterval
    tickpositions = np.arange(tickmin, options.flanking_window, tickinterval)
    tickpositions_zeroi = tickpositions + options.flanking_window
    xticks = np.array([x for x in tickpositions_zeroi if 0 <= x <= seqwidth-1])
    yticks = np.array([y for y in tickpositions_zeroi if 1 <= y <= seqwidth])
    plt.xticks(xticks + 0.5)
    plt.yticks(yticks + 0.5)
    ax.xaxis.set_ticklabels(xticks - options.flanking_window)
    ax.yaxis.set_ticklabels(yticks - options.flanking_window)
    plt.setp(ax.xaxis.get_ticklabels(), rotation=-45)
    plt.setp(ax.yaxis.get_ticklabels(), rotation=-45)

    # format colorbar
    if not preview:
        cbar = plt.colorbar(orientation='horizontal', ticks=[-1, np.log2(0.75), 0, np.log2(1.5), 1],
                            shrink=0.5, aspect=10)
        cbar.set_label('Base-pairing propensity')
        cbar.ax.set_xticklabels(['50%', '75%', 'neutral', '150%', '200%'], position='top')
        cbar.ax.xaxis.set_ticks_position('top')

    plt.savefig(output, format=('png' if preview else 'pdf'))


def generate_base_weblogo(options, allsequences, output):
    import weblogolib as wl
    from weblogolib.colorscheme import ColorScheme, ColorGroup
    from corebio import seq_io
    from StringIO import StringIO

    ecliptic_color_scheme = ColorScheme([
        ColorGroup('A', iwork_colors.blue),
        ColorGroup('C', iwork_colors.green),
        ColorGroup('G', iwork_colors.yellow),
        ColorGroup('UT', iwork_colors.red),
    ])

    fastainput = StringIO(''.join('>%d\n%s\n' % s for s in enumerate(allsequences)))

    rna = wl.std_alphabets['rna']
    seqs = wl.read_seq_data(fastainput, alphabet=rna)
    logo = wl.LogoData.from_seqs(seqs)
    logooptions = wl.LogoOptions()
    logooptions.unit_name = 'probability'
    logooptions.color_scheme = ecliptic_color_scheme
    logooptions.show_fineprint = False
    logooptions.first_index = -options.flanking_window
    logooptions.title = "Ecliptic"

    format = wl.LogoFormat(logo, logooptions)
    wl.pdf_formatter(logo, format, open(output, 'w'))


def make_justified_plot(inputpng, outputpng, edgetrim=2, padding=10, tophalfpadding=25):
    # rotate 45 degree and crop tight.
    import Image

    im = Image.open(inputpng)
    rotated = im.rotate(45, expand=1)

    # put the image into large white background and rotate 45deg, then crop.
    white = Image.new('RGBA', (im.size[0]*2, im.size[1]*2), (255,)*4)
    newposition = ((white.size[0] - im.size[0]) // 2, (white.size[0] - im.size[1]) // 2)
    newposition += (newposition[0] + im.size[0], newposition[1] + im.size[1])
    white.paste(im, newposition)
    justified = white.rotate(45, Image.BICUBIC)
    centerbox = (justified.size[0] - rotated.size[0]) // 2, (justified.size[0] - rotated.size[1]) // 2
    centerbox += (centerbox[0] + rotated.size[0], centerbox[1] + rotated.size[1])
    center = justified.crop(centerbox)

    # find a minimal bounding box
    imagev = np.fromstring(center.tostring(), np.uint8).reshape(
                           (center.size[1], center.size[0], 4,))
    nonwhite = (imagev < 250).max(axis=2)[edgetrim:-edgetrim, edgetrim:-edgetrim]

    nonwhite_vert = np.where(nonwhite.max(axis=0))[0] + edgetrim
    nonwhite_horiz = np.where(nonwhite.max(axis=1))[0] + edgetrim
    top, bottom = nonwhite_vert.min(), nonwhite_vert.max()
    left, right = nonwhite_horiz.min(), nonwhite_horiz.max()

    # extend the bounding box for enough padding area
    top = max(0, top - padding)
    bottom = min(imagev.shape[1], bottom + padding)
    left = max(0, left - padding)
    right = min(imagev.shape[0], right + padding)
    top = int(bottom - (bottom - top) * 0.5 - tophalfpadding)

    # finalize the image and write output
    final = Image.fromstring('RGBA', (right - left, bottom - top),
                             imagev[top:bottom, left:right].tostring())
    final.save(outputpng, 'png')


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Produce a plot showing propensity to form '
                                                 'WC base pairing around crosslinked sites')
    parser.add_argument('catalog', metavar='FILE', type=str,
                        help='Input binding site catalog file.')
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
    parser.add_argument('--pattern', metavar='RE', type=str, default=None,
                        dest='pattern', help='Pattern that must be included in sequences in '
                             'regular expression (optional)')
    parser.add_argument('--max-sequences', metavar='COUNT', type=int, default=32768,
                        dest='max_seqs', help='Maximum number of sampled sequences')
    parser.add_argument('--output-preview-png', metavar='FILE', type=str,
                        dest='preview_png', default=None, help='Path for newly generated '
                        'preview image')
    parser.add_argument('--output-package', metavar='FILE', type=str,
                        dest='plot_package', default=None, help='Path for plot output '
                        'package (.zip)')

    return parser.parse_args()


INSTRUCTION = """\
This package contains two PDF files containing parts to be assembed
into a plot for base pairing propensity near cross-linked sites.

* basepairing-plot.pdf: the main heat map part and colormap
* sequence-logo.pdf: sequence logo for positional base probability

Place the parts into a single image using a vector image editor
such as Inkscape or Adobe Illustrator in this order:

1. Rotate the main heat map -45 degree and crop the upper half.
2. Place the color bar in appropriate place (lower right, maybe).
3. Open the logo.
4. Flip vertically x-axis and ticks on the logo, then move them
   to upper part of the logo.
5. Grab the whole logo objects to the heat map.
6. Align carefully.
7. Save the image.

Good luck!
""".replace('\n', '\r\n')

if __name__ == '__main__':
    options = parse_arguments()

    seqiter = load_sequences(gzip.open(options.catalog),
                options.score_type, options.minimum_depth, options.minimum_score,
                options.flanking_window, options.source == 'transcriptome')

    if options.pattern:
        seqiter = ifilter(re.compile(options.pattern).search, seqiter)

    allsequences = list(seqiter)
    if len(allsequences) >= options.max_seqs:
        random.shuffle(allsequences)
        del allsequences[options.max_seqs:]

    maketempfile = lambda suffix: tempfile.NamedTemporaryFile(suffix='.'+suffix)

    with maketempfile('pdf') as bppdf, maketempfile('pdf') as logopdf, \
            maketempfile('png') as previewpng:

        if options.plot_package is not None:
            generate_basepairing_plot(options, allsequences, bppdf.name)
            generate_base_weblogo(options, allsequences, logopdf.name)

            pkgzip = zipfile.ZipFile(options.plot_package, 'w', zipfile.ZIP_DEFLATED)
            pkgzip.write(bppdf.name, 'basepairing-plot.pdf')
            pkgzip.write(logopdf.name, 'sequence-logo.pdf')
            pkgzip.writestr('README.txt', INSTRUCTION)
            pkgzip.close()

        if options.preview_png is not None:
            generate_basepairing_plot(options, allsequences, previewpng.name, preview=True)
            make_justified_plot(previewpng.name, options.preview_png)

