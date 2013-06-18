#!/usr/bin/env python
from ecliptic.support.bxwrap import MultiTrackSplitBinnedArray
import shelve
from itertools import groupby
from operator import itemgetter
from functools import partial
import futures
import pickle

BLOCKDEFS = {
    ('NM', '+'): [
        ('CDS', 'cdsBlocks'),
        ('5UTR', 'leftUtrBlocks'),
        ('3UTR', 'rightUtrBlocks'),
        ('Intron', 'intronBlocks'),
    ],
    ('NM', '-'): [
        ('CDS', 'cdsBlocks'),
        ('5UTR', 'rightUtrBlocks'),
        ('3UTR', 'leftUtrBlocks'),
        ('Intron', 'intronBlocks'),
    ],
    ('NR', '+'): [('Exon', 'exonBlocks'), ('Intron', 'intronBlocks')],
    ('NR', '-'): [('Exon', 'exonBlocks'), ('Intron', 'intronBlocks')],
}

def process(refdbfile, array_dir, nproc):
    refdb = shelve.open(refdbfile, 'r')
    orderedacc = sortbycoord(refdb)

    def iter_refseq():
        for acc in orderedacc:
            refinfo = refdb[acc]
            yield refinfo['chrom'], refinfo['strand'], acc

    counts_all = {}

    with futures.ProcessPoolExecutor(nproc) as executor:
        jobfun = partial(count_refseq_transcripts, array_dir, refdbfile)
        for jobret in executor.map(jobfun, (
                        list(grp) for k, grp in groupby(iter_refseq(), key=lambda x: x[:2]))):
            counts_all.update(jobret)

    return counts_all


def count_refseq_transcripts(array_dir, refdbfile, transcripts):
    # Open own database handler per process. Sharing database handler results unexpected
    # errors and various kinds of malfunctions.
    refdb = shelve.open(refdbfile, 'r')

    arr = MultiTrackSplitBinnedArray(array_dir)
    i_fivep = arr.tracks.index('5')
    i_threep = arr.tracks.index('3')
    i_insertion = arr.tracks.index('I')
    i_deletion = arr.tracks.index('D')
    r = {}

    for chrom, strand, acc in transcripts:
        refinfo = refdb[acc]

        blocks = BLOCKDEFS[acc[:2], strand]
        for settitle, label in blocks:
            if not refinfo[label]:
                continue

            cntarray = arr.get_blocks(chrom, refinfo[label], strand)
            cntsummary = (
                max(cntarray[i_fivep].sum(), cntarray[i_threep].sum()),
                cntarray[i_insertion].sum(),
                cntarray[i_deletion].sum())

            if sum(cntsummary) > 0:
                r[acc, settitle] = cntsummary

    return r

def sortbycoord(refdb):
    def sortkey(key):
        info = refdb[key]
        return info['chrom'], info['strand'], info['txStart'], info['txEnd']
    return sorted(refdb.keys(), key=sortkey)


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Counts reads mapped to each RefSeq '
                                                 'transcripts')
    parser.add_argument('--refflat-db', metavar='FILE', action='store', dest='refflat_db',
                        help='path to refflat database of ecliptic resources',
                        required=True)
    parser.add_argument('--array-dir', metavar='DIR', dest='array_dir', action='store',
                        help='path to directory containing BinnedArray track files',
                        required=True)
    parser.add_argument('--output', metavar='FILE', dest='output', action='store',
                        help='path to output pickle file to be created',
                        required=True)
    parser.add_argument('--parallel', metavar='N', dest='parallel', action='store',
                        type=int, help='path to output pickle file to be created',
                        default=4)
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    import sys

    options = parse_arguments()

    r = process(options.refflat_db, options.array_dir, options.parallel)
    pickle.dump(r, open(options.output, 'w'))

