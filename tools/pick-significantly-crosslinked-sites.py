#!/usr/bin/env python

def process(inpfile, outfile, get_cutoff, scoretype):
    scorecol = {
        'del': 10, 'mod': 11, 'moddel': 12, 'entropy': 13, 't2c': 14,
    }[scoretype]

    for line in inpfile:
        fields = line[:-1].split()

        depth = int(fields[4])
        base = fields[3]
        cutoff = get_cutoff(base, depth)
        if cutoff is None:
            # insufficient reads or undecideable point without enough complexity
            continue

        score = float(fields[scorecol])
        if score >= cutoff > 0.0:
            print >> outfile, ' '.join(fields[:5] + [fields[scorecol]])

if __name__ == '__main__':
    import sys
    import gzip
    from functools import partial
    from ecliptic.CLIP.ErrorAnalysis import ErrorScoreCutoff

    nonzerosites = sys.argv[1]
    cutoffs = sys.argv[2]
    fdrlimit = sys.argv[3]
    scoretype = sys.argv[4]

    get_cutoff = partial(ErrorScoreCutoff(open(cutoffs)).get_cutoff, float(fdrlimit))
    process(gzip.open(nonzerosites), sys.stdout, get_cutoff, scoretype)

