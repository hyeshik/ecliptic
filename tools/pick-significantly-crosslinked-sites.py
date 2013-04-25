#!/usr/bin/env python
from ecliptic.Utils import TemporaryDirectory
from ecliptic.support.fileutils import ParallelMatchingReader, LineParser
from ecliptic.support.sequtils import GiantFASTAFile
from ecliptic.support.coordutils import LiftOverToTranscriptome
from itertools import izip
import os
import gzip

def generate_sites_bed(inpfile, get_cutoff, scoretype):
    scorecol = {
        'del': 10, 'mod': 11, 'moddel': 12, 'entropy': 13, 't2c': 14,
    }[scoretype]

    siteno = 0

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
            yield (fields[0], int(fields[2]), fields[1], base, depth, fields[scorecol])




class FlankingSequenceGetter(object):

    def __init__(self, reffile, flanking):
        self.reffile = GiantFASTAFile(reffile)
        self.flanking = flanking

    def get(self, chrom, pos, strand):
        chromsize = self.reffile.index[chrom][0]
        winleft = max(0, pos - self.flanking)
        winright = min(chromsize, pos + self.flanking + 1)
        leftfillsize = self.flanking - (pos - winleft)
        rightfillsize = self.flanking + 1 - (winright - pos)

        seq = self.reffile.get(chrom, winleft, winright, strand).upper()
        if leftfillsize != 0 or rightfillsize != 0:
            seq = leftfillsize * '-' + seq + rightfillsize * '-'

        return seq


def produce_output(outf, sites, liftover, gseqget, trseqget):
    converted_sites = liftover.liftover_points(sites)

    for gsite, trsite in izip(sites, converted_sites):
        # gsite: (chrom, pos, strand, base, depth, score)
        # trsite: (accession, pos) or None

        gseq = gseqget.get(gsite[0], int(gsite[1]), gsite[2])
        if trsite is None:
            trseq = traccession = trrefposition = ''
        else:
            trseq = trseqget.get(trsite[0], int(trsite[1]), '+')
            traccession = trsite[0]
            trrefposition = trsite[1]

        # chrom strand pos base depth score genome-seq transcriptome-seq tracc trpos
        print >> outf, '\t'.join([gsite[0], gsite[2], str(gsite[1]), # chrom strand pos
                                  gsite[3], str(gsite[4]), gsite[5], # base depth score
                                  gseq, trseq, traccession, str(trrefposition)])

bedparser = LineParser([
    ('chrom', None),
    ('start', int),
    ('end', int),
    ('name', None),
    ('score', None),
    ('strand', None)
])

if __name__ == '__main__':
    import sys
    from functools import partial
    from ecliptic.CLIP.ErrorAnalysis import ErrorScoreCutoff

    nonzerosites = sys.argv[1]
    cutoffs = sys.argv[2]
    fdrlimit = sys.argv[3]
    scoretype = sys.argv[4]
    alignmentdb = sys.argv[5]
    genomeref = sys.argv[6]
    nrrefseqref = sys.argv[7]
    flanking = int(sys.argv[8])
    nrrefseqbed = sys.argv[9]

    get_cutoff = partial(ErrorScoreCutoff(open(cutoffs)).get_cutoff, float(fdrlimit))
    liftover = LiftOverToTranscriptome(alignmentdb, nrrefseqbed)

    sites = list(generate_sites_bed(gzip.open(nonzerosites), get_cutoff, scoretype))

    gseqget = FlankingSequenceGetter(genomeref, flanking)
    trseqget = FlankingSequenceGetter(nrrefseqref, flanking)

    produce_output(sys.stdout, sites, liftover, gseqget, trseqget)

