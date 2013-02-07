#!/usr/bin/env python
from ecliptic.support.bxwrap import MultiTrackSplitBinnedArray
import shelve
import pickle

class ContinuousCache(object):

    def __init__(self, opener):
        self.curkey = None
        self.current = None
        self.opener = opener

    def get(self, key):
        if self.curkey != key:
            self.curkey = key
            self.current = self.opener(key)

        return self.current

BLOCKDEFS = {
    ('NM', '+'): [
        ('CDS', 'cdsBlocks'),
        ('5UTR', 'leftUtrBlocks'),
        ('3UTR', 'rightUtrBlocks'),
    ],
    ('NM', '-'): [
        ('CDS', 'cdsBlocks'),
        ('5UTR', 'rightUtrBlocks'),
        ('3UTR', 'leftUtrBlocks'),
    ],
    ('NR', '+'): [('Exon', 'exonBlocks')],
    ('NR', '-'): [('Exon', 'exonBlocks')],
}

def process(refdb, arr):
    i_fivep = arr.tracks.index('5')
    i_threep = arr.tracks.index('3')
    i_insertion = arr.tracks.index('I')
    i_deletion = arr.tracks.index('D')
    r = {}

    orderedacc = sortbycoord(rfdb)
    cachecleaner = ContinuousCache(lambda key: arr.clear_cache())

    for acc in orderedacc:
        refinfo = refdb[acc]
        cachecleaner.get((refinfo['chrom'], refinfo['strand']))

        blocks = BLOCKDEFS[acc[:2], refinfo['strand']]
        for settitle, label in blocks:
            if not refinfo[label]:
                continue

            cntarray = arr.get_blocks(refinfo['chrom'], refinfo[label],
                                      refinfo['strand'])
            cntsummary = (
                max(cntarray[i_fivep].sum(), cntarray[i_threep].sum()),
                cntarray[i_insertion].sum(),
                cntarray[i_deletion].sum())

            r[acc, settitle] = cntsummary

    return r

def sortbycoord(refdb):
    def sortkey(key):
        info = refdb[key]
        return info['chrom'], info['strand'], info['txStart'], info['txEnd']
    return sorted(refdb.keys(), key=sortkey)


if __name__ == '__main__':
    import sys

    output = sys.argv[1]
    refflatdb = sys.argv[2]
    arrprefix = sys.argv[3]

    rfdb = shelve.open(refflatdb, 'r')
    arr = MultiTrackSplitBinnedArray(arrprefix)
    r = process(rfdb, arr)
    pickle.dump(r, open(output, 'w'))

