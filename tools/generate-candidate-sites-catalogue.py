#!/usr/bin/env python
from ecliptic.Utils import TemporaryDirectory
from ecliptic.support.fileutils import ParallelMatchingReader, LineParser
from ecliptic.support.sequtils import GiantFASTAFile
from ecliptic.support.coordutils import LiftOverToTranscriptome
from itertools import izip
import pickle
import os
import gzip

SCORE_TYPES = ['del', 'mod', 'moddel', 'entropy', 't2c']
INPUT_SCORE_COLS = {'del': 10, 'mod': 11, 'moddel': 12, 'entropy': 13, 't2c': 14}

# nonzero.gz files (source)
#
# 0    1 2        3 4 5 6 7 8 9 10             11  12             13             14
# chr1 + 40672292 C 3 0 1 0 0 2 0.666666666667 0.0 0.666666666667 0.636514168295 0.0
# 0:chrom 1:strand 2:pos 3:base 4:depth 5:A 6:C 7:G 8:T 9:D
# 10:del 11:mod 12:moddel 13:entropy 14:t2c

# output (numbers are fictional)
#
# 0    1 2        3 4 5    6   7    8     9   10    11    12    13 14   15    16    17     18
# chr1 + 40672292 C 3 0.67 0.0 0.67 0.635 0.0 0.001 0.001 0.005 1 0.001 AGCDG AGCDF NM_123 142
# 0:chrom 1:strand 2:pos 3:base 4:depth 5:del 6:mod 7:moddel 8:entropy 9:t2c
# 10:fdr_del 11:fdr_mod 12:fdr_moddel 13:fdr_entorpy 14:fdr_t2c 15:flanking_genomic
# 16:flanking_transcriptomic 17:refseq_accession 18:refseq_pos
# 10-14 are the least FDR passing the site.

def prepare_original(inpfile, accfdrgetter):
    for line in inpfile:
        fields = line[:-1].split()

        depth = int(fields[4])
        base = fields[3]
        minfdrs = [accfdrgetter.min_accepted_fdr(scoretype, base, depth,
                                                 float(fields[INPUT_SCORE_COLS[scoretype]]))
                   for scoretype in SCORE_TYPES]
        yield [fields[0], int(fields[2]), fields[1]] + fields[:5] + fields[10:15] + map(str, minfdrs)
        #     |--> for liftover only <-------------|


class AcceptedFDRGetter(object):

    def __init__(self, fdrs, cutoffs_prefix):
        self.fdrs = sorted(fdrs)
        self.cutoffs = dict(
            (scoretype, ErrorScoreCutoff(open('{prefix}{type}.fdr_cutoffs'.format(
                                         prefix=cutoffs_prefix, type=scoretype))))
            for scoretype in SCORE_TYPES)

    def min_accepted_fdr(self, scoretype, base, depth, score):
        for fdr in self.fdrs:
            cutoff = self.cutoffs[scoretype].get_cutoff(fdr, base, depth)
            if cutoff is not None and score >= cutoff:
                return fdr
        else:
            return 1


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


def write_output(outf, gseqget, trseqget, gsite, trsite):
    # gsite: (chrom, pos, strand) + (chrom, strand, pos, base, depth, ....)
    # trsite: (accession, pos) or None

    gseq = gseqget.get(gsite[0], gsite[1], gsite[2])
    if trsite is None:
        trseq = traccession = trrefposition = ''
    else:
        trseq = trseqget.get(trsite[0], int(trsite[1]), '+')
        traccession = trsite[0]
        trrefposition = trsite[1]

    print >> outf, '\t'.join(gsite[3:] + [gseq, trseq, traccession, str(trrefposition)])


def process_by_blocks(inputgen, annotator, outputfunc, blocksize):
    while True:
        blk = [entry for entry, _ in izip(inputgen, xrange(blocksize))]
        if not blk:
            break

        annotations = annotator(blk)

        for orig, anno in izip(blk, annotations):
            outputfunc(orig, anno)


if __name__ == '__main__':
    import sys
    from functools import partial
    from ecliptic.CLIP.ErrorAnalysis import ErrorScoreCutoff

    nonzerosites = sys.argv[1]
    cutoffs_prefix = sys.argv[2]
    fdrs = sys.argv[3]
    alignmentdb = sys.argv[4]
    genomeref = sys.argv[5]
    nrrefseqref = sys.argv[6]
    nrrefseqbed = sys.argv[7]
    flanking = int(sys.argv[8])

    fdrs_getter = AcceptedFDRGetter(map(float, fdrs.split(',')), cutoffs_prefix)

    liftover = LiftOverToTranscriptome(alignmentdb, nrrefseqbed)
    gseqget = FlankingSequenceGetter(genomeref, flanking)
    trseqget = FlankingSequenceGetter(nrrefseqref, flanking)

    origstream = prepare_original(gzip.open(nonzerosites), fdrs_getter)
    writer = partial(write_output, sys.stdout, gseqget, trseqget)
    process_by_blocks(origstream, liftover.liftover_points, writer, 65536)

