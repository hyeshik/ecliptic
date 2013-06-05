#!/usr/bin/env python
from StringIO import StringIO
from ecliptic.FormattedParsers import BindingSiteCatalogParser
from ecliptic.support.plotutils import iwork_colors
from operator import attrgetter
import gzip
import csv

#<ParsedLine chrom='chr1' strand='+' pos=3532312 base='G' depth=10 del_score=0.7 mod_score=0.3 moddel_score=1.0 entropy_score=0.610864302055 t2c_score=0.0 del_fdr=0.0001 mod_fdr=1.0 moddel_fdr=0.0001 entropy_fdr=1.0 t2c_fdr=1.0 genome_seq='TGTTTGCCACCCTTGGTGTGGTAGTCACGTCCCTTTAACAAGTGTGGGTCGGAGCCAAGCAGCACGAACTGCAGGTGAACTCCATACAGCTACTTTACTAC' transcript_seq='' transcript_acc='' transcript_pos=None>


def get_center_seq(s, width):
    assert len(s) % 2 == 1
    return s[len(s) // 2 - width:len(s) // 2 + width + 1].replace('T', 'U')


def load_sequences(inp, scoretype, mindepth, minscore, seqwidth):
    getscore = attrgetter(scoretype + '_score')
    getfdr = attrgetter(scoretype + '_fdr')
    genome_seqs = StringIO()
    transcript_seqs = StringIO()
    genome_nseqs = transcript_nseqs = 0

    for i, row in enumerate(BindingSiteCatalogParser(inp)):
        if row.strand == '-':
            continue
        if row.depth >= mindepth and getscore(row) >= minscore:
        #if getfdr(row) < 0.0011:
            genome_seqs.write('>%d\n%s\n' % (i, get_center_seq(row.genome_seq, seqwidth)))
            genome_nseqs += 1
            if row.transcript_seq:
                transcript_seqs.write('>%d\n%s\n' % (i,
                                        get_center_seq(row.transcript_seq, seqwidth)))
                transcript_nseqs += 1

    genome_seqs.seek(0, 0)
    transcript_seqs.seek(0, 0)

    return genome_seqs, transcript_seqs, genome_nseqs, transcript_nseqs


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Generate weblogos for significantly '
                                                 'crosslinked sites for the report')
    parser.add_argument('catalog', metavar='FILE', type=str,
                        help='Input binding site catalog file.')
    parser.add_argument('--score-type', metavar='TYPE', type=str, required=True,
                        dest='score_type', help='One of del, moddel, mod, entropy, or t2c')
    parser.add_argument('--minimum-depth', metavar='READS', type=int, required=True,
                        dest='minimum_depth', help='Minimum read count for reporting')
    parser.add_argument('--minimum-score', metavar='VALUE', type=float, required=True,
                        dest='minimum_score', help='Score cut-off for reporting')
    parser.add_argument('--window', metavar='WIDTH', type=int, required=True,
                        dest='window', help='Window width of weblogo')
    parser.add_argument('--seq-source', metavar='TYPE', type=str, required=True,
                        dest='source', help='One of genome or transcriptome')
    parser.add_argument('--logo-type', metavar='TYPE', type=str, default='probability',
                        dest='logo_type', help='One of probability or bits')
    parser.add_argument('--output-pdf', metavar='FILE', type=str, default=None,
                        dest='output_pdf', help='Path to write a PDF output')
    parser.add_argument('--output-png', metavar='FILE', type=str, default=None,
                        dest='output_png', help='Path to write a PNG output')

    return parser.parse_args()


if __name__ == '__main__':
    from weblogolib import *
    from weblogolib.colorscheme import ColorScheme, ColorGroup
    from corebio import seq_io

    options = parse_arguments()

    ecliptic_color_scheme = ColorScheme([
        ColorGroup('A', iwork_colors.blue),
        ColorGroup('C', iwork_colors.green),
        ColorGroup('G', iwork_colors.yellow),
        ColorGroup('UT', iwork_colors.red),
    ])

    gseqs, trseqs, gseqn, trseqn = load_sequences(gzip.open(options.catalog),
                options.score_type, options.minimum_depth, options.minimum_score,
                options.window)
    if options.source == 'genome':
        oseqs, oseqn = gseqs, gseqn
    elif options.source == 'transcriptome':
        oseqs, oseqn = trseqs, trseqn
    else:
        raise ValueError('Unknown sequence source: %s' % repr(options.source))

    if oseqn <= 1:
        # generate empty graphic files to skip gracefully the toolchain
        if options.output_pdf is not None:
            open(options.output_pdf, 'w')
        if options.output_png is not None:
            open(options.output_png, 'w')
        raise SystemExit

    rna = std_alphabets['rna']
    seqs = read_seq_data(oseqs, alphabet=rna)
    logo = LogoData.from_seqs(seqs)
    #logo = LogoData.from_seqs(read_seq_data(gseqs, alphabet=rna))
    logooptions = LogoOptions()
    logooptions.unit_name = options.logo_type
    logooptions.color_scheme = ecliptic_color_scheme
    logooptions.show_fineprint = False
    logooptions.first_index = -options.window
    logooptions.title = "Ecliptic"

    format = LogoFormat(logo, logooptions)

    if options.output_pdf is not None:
        fout = open(options.output_pdf, 'w')
        pdf_formatter(logo, format, fout)

    if options.output_png is not None:
        fout = open(options.output_png, 'w')
        png_formatter(logo, format, fout)

