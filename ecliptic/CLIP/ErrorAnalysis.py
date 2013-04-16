#!/usr/bin/env python
from __future__ import division
import pickle
import functools

__all__ = ['ErrorScoreCutoff']

def memoize(obj):
    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        if args not in cache:
            cache[args] = obj(*args, **kwargs)
        return cache[args]
    return memoizer


class ErrorScoreCutoff(object):

    EXTRAPOLATE_WIDTH = 5

    def __init__(self, cutoff_file):
        self.cutoffs = pickle.load(cutoff_file)

    def get_cutoff(self, fdr, base, depth):
        if base not in 'ACGT':
            return None
        elif depth in self.cutoffs[fdr][base]:
            return self.cutoffs[fdr][base][depth]
        else:
            # extrapolate cutoff using cutoffs for lower depths (conservative choose)
            return self.extrapolate(fdr, base, depth)

    @memoize
    def extrapolate(self, fdr, base, depth):
        pastpoints = []
        pastdepth = depth - 1
        cutoffpool = self.cutoffs[fdr][base]

        while pastdepth >= 1 and len(pastpoints) < self.EXTRAPOLATE_WIDTH:
            if pastdepth in cutoffpool:
                pastpoints.append(cutoffpool[pastdepth])
            pastdepth -= 1

        if len(pastpoints) < self.EXTRAPOLATE_WIDTH:
            return None # no enough past points to guess

        midpoints = sorted(pastpoints)[1:-1]
        return sum(midpoints) / len(midpoints)


if __name__ == '__main__':
    cutoffs = ErrorScoreCutoff(open('../../work/GSE37114-LIN28A-NKim/erroranalyses/CLIP-35L33G.entropy.fdr_cutoffs'))

    from matplotlib import pyplot as plt
    import numpy as np

    FDR = 0.0001
    BASE = 'G'
    Xrange = range(0, 200)

    xavail = [x for x in Xrange if x in cutoffs.cutoffs[FDR][BASE]]
    yavail = [cutoffs.cutoffs[FDR][BASE][x] for x in xavail]

    getcutoff = functools.partial(cutoffs.get_cutoff, FDR, BASE)
    xnotavail = [x for x in Xrange if x not in cutoffs.cutoffs[FDR][BASE]]
    ynotavail = [(getcutoff(x) or 0) for x in xnotavail]

    plt.figure(figsize=(30, 14))
    plt.scatter(xavail, yavail, s=10, c='blue', linewidth=0)
    plt.scatter(xnotavail, ynotavail, s=20, c='red', linewidth=0)
    plt.ylim(0, 1.6)
    plt.grid(True)
    plt.show()

