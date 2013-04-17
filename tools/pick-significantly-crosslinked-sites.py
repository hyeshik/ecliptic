#!/usr/bin/env python
from ecliptic.Utils import TemporaryDirectory
from ecliptic.support.fileutils import ParallelMatchingReader, LineParser
from ecliptic.support.sequtils import GiantFASTAFile
import os
import gzip

def generate_sites_bed(inpfile, outfile, get_cutoff, scoretype):
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
            siteno += 1
            print >> outfile, '\t'.join((
                fields[0], fields[2], str(int(fields[2]) + 1),
                '%s.%08x:%s%s' % (scoretype, siteno, fields[3], fields[4]),
                fields[scorecol], fields[1]
            ))

def add_transcriptome_mapping(bedsource, chain, tmpdir, transcripts_by_strand):
    trmapped_file = os.path.join(tmpdir, 'mapped')
    trunmapped_file = os.path.join(tmpdir, 'unmapped') # unused; for debug only

    res = os.system('liftOver -multiple {bedsource} {chain} {trmapped_file} {trunmapped_file}'
                    .format(bedsource=bedsource, chain=chain, trmapped_file=trmapped_file,
                            trunmapped_file=trunmapped_file))
    if res != 0:
        raise OSError('liftOver failed with result code %d' % res)

    for gblk, trblk in ParallelMatchingReader(open(bedsource), open(trmapped_file), 3):
        assert gblk is not None
        gentry = next(gblk)

        # remove regions mapped to reverse to a transcript (e.g. crossmappings between
        # Tsix and Xist). The strand information in the liftOver output is unreliable
        # due to a bug in liftOver.
        trlist = transcripts_by_strand[gentry[5]]
        trblk = [trmap[:5] + ['+'] for trmap in (trblk or '') if trmap[0] in trlist]
#        if len(trblk) >= 2:
#            print '================='
#            print gentry
#            print trblk
#            raw_input('%s ===>' % tmpdir)
        assert len(trblk) <= 1

        yield (gentry, trblk[0] if trblk else None)


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

        seq = self.reffile.get(chrom, pos - self.flanking,
                               pos + self.flanking + 1, strand).upper()
        if leftfillsize != 0 or rightfillsize != 0:
            seq = leftfillsize * '-' + seq + rightfillsize * '-'

        return seq


def produce_output(outf, trmapping, gseqget, trseqget):
    for gmap, trmap in trmapping:
        gseq = gseqget.get(gmap[0], int(gmap[1]), gmap[5])
        trseq = '' if trmap is None else trseqget.get(trmap[0], int(trmap[1]), trmap[5])

        siteinfo = gmap[3].split(':')[1]
        base = siteinfo[0]
        depth = siteinfo[1:]
        # chrom strand pos base depth score genome-seq transcriptome-seq
        print >> outf, '\t'.join([gmap[0], gmap[5], gmap[1], base, depth, gmap[4], gseq, trseq])

bedparser = LineParser([
    ('chrom', None),
    ('start', int),
    ('end', int),
    ('name', None),
    ('score', None),
    ('strand', None)
])

def load_transcript_strand(bedfile):
    strands = {'+': set(), '-': set()}
    for row in bedparser.iter_parse(gzip.open(bedfile)):
        strands[row.strand].add(row.name.split('/')[0])
    return strands

if __name__ == '__main__':
    import sys
    from functools import partial
    from ecliptic.CLIP.ErrorAnalysis import ErrorScoreCutoff

    nonzerosites = sys.argv[1]
    cutoffs = sys.argv[2]
    fdrlimit = sys.argv[3]
    scoretype = sys.argv[4]
    transcriptome_chain = sys.argv[5]
    genomeref = sys.argv[6]
    nrrefseqref = sys.argv[7]
    flanking = int(sys.argv[8])
    nrrefseqbed = sys.argv[9]

    trstrands = load_transcript_strand(nrrefseqbed)
    get_cutoff = partial(ErrorScoreCutoff(open(cutoffs)).get_cutoff, float(fdrlimit))

    with TemporaryDirectory() as tmpdir:
        sitesbed = os.path.join(tmpdir, 'sites.bed')

        generate_sites_bed(gzip.open(nonzerosites), open(sitesbed, 'w'), get_cutoff, scoretype)

        trmappings = add_transcriptome_mapping(sitesbed, transcriptome_chain, tmpdir, trstrands)

        gseqget = FlankingSequenceGetter(genomeref, flanking)
        trseqget = FlankingSequenceGetter(nrrefseqref, flanking)

        produce_output(sys.stdout, trmappings, gseqget, trseqget)

