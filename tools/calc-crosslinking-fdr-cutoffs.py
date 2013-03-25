#!/usr/bin/env python
from __future__ import division
from itertools import groupby, dropwhile, takewhile
import numpy as np
import gzip

CUTOFF_MAXIMUM_DEPTH = 5000
CUTOFF_FDRS = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]

def load_dist_file(filename):
    for line in gzip.open(filename):
        depth, base, score, count = line.split()
        yield (int(depth), base, float(score), int(count))

def iterate_score_distribution(permfile, realfile):
    permiter = load_dist_file(permfile)
    realiter = load_dist_file(realfile)

    permnext = next(permiter)
    realnext = next(realiter)

    while permnext is not None or realnext is not None:
        if permnext is None:
            yield realnext[:3] + (0, realnext[3])
            realnext = next(realiter, None)
        elif realnext is None:
            yield permnext + (0,)
            permnext = next(permiter, None)
        else:
            pkey = permnext[:3]
            rkey = realnext[:3]
            if pkey == rkey:
                yield permnext + (realnext[3],)
                permnext = next(permiter, None)
                realnext = next(realiter, None)
            elif pkey < rkey:
                yield permnext + (0,)
                permnext = next(permiter, None)
            else:
                yield realnext[:3] + (0, realnext[3])
                realnext = next(realiter, None)

def calculate_fdr_curves(permfile, realfile, fdrlist):
    distiter = iterate_score_distribution(permfile, realfile)
    fdrcutoffs = dict((fdrlimit, dict((base, {}) for base in 'ACGT'))
                      for fdrlimit in fdrlist)

    for (depth, base), grp in groupby(distiter, key=lambda x: x[0:2]):
        if depth > CUTOFF_MAXIMUM_DEPTH:
            break
        if base not in 'ACGT':
            continue

        blockscores = list(grp)[::-1]

        # Unless drop first continuos zero count scores in the real count data,
        # they will make zero-division errors.
        try:
            firstnonzero = dropwhile(lambda x: x[1][4] < 1,
                                     enumerate(blockscores)).next()[0]
        except StopIteration:
            continue

        perm_counts = [score_count[3] for score_count in blockscores]
        real_counts = [score_count[4] for score_count in blockscores]
        counts = np.array([perm_counts, real_counts]).transpose()

        cumul_count = np.cumsum(counts, 0)
        # One of cumul_count (divisor) can be zero if no read is assigned to
        # a specific depth and base. So, clip it so that minimum value is 1,
        # this is not distorting the result because the cumulative fraction
        # is still zero for the datapoint.
        cumul_count_frac = (cumul_count /
            cumul_count[cumul_count.shape[0]-1].clip(1,
                np.max(cumul_count[cumul_count.shape[0]-1])+1))

        fdr = (cumul_count_frac[firstnonzero:, 0] /
               cumul_count_frac[firstnonzero:, 1]).clip(0, 1)

        cutoffs = [score_count[2] for score_count in blockscores][firstnonzero:]
        assert len(cutoffs) == len(fdr)

        for fdrlimit in fdrlist:
            try:
                validpositions = list(takewhile(lambda x: x < fdrlimit, fdr))
                if len(validpositions) < 1:
                    continue
                valid_cutoff = cutoffs[len(validpositions) - 1]
                fdrcutoffs[fdrlimit][base][depth] = valid_cutoff
            except StopIteration:
                pass

        #print depth, base, [fdrcutoffs[fdrlimit][base].get(depth, None) for fdrlimit in fdrlist]

    return fdrcutoffs

if __name__ == '__main__':
    import sys
    import cPickle as pickle

    prefix = sys.argv[1]
    output = sys.argv[2]

    #permfile = '../erroranalysis/CLIP-35L33G-entropy.perm.gz'
    #realfile = '../erroranalysis/CLIP-35L33G-entropy.real.gz'
    permfile = prefix + 'perm.gz'
    realfile = prefix + 'real.gz'

    cutoffs = calculate_fdr_curves(permfile, realfile, CUTOFF_FDRS)

    pickle.dump(cutoffs, open(output, 'w'))

