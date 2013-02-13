#!/usr/bin/env python
from ecliptic.support import sam
import os
import re
from collections import defaultdict
from functools import partial
import futures
import glob
import numpy as np

BASES = 'ACGT-'

class ReadMappingProfile(object):

    DIBASES = [n1+n2 for n1 in BASES for n2 in BASES]
    DIBASE2NUM = dict((bases, i) for i, bases in enumerate(DIBASES))

    def __init__(self, maxreadlength):
        self.cntarray = np.zeros([maxreadlength, len(self.DIBASES)], np.uint64)
        self.maxreadlength = maxreadlength

    def countup(self, pos, refbase, readbase, ndup):
        dibaseidx = self.DIBASE2NUM.get(refbase + readbase)
        if dibaseidx is not None:
            self.cntarray[pos, dibaseidx] += ndup

    def makematrix(self):
        # first: read position, second: from, third: to
        return self.cntarray.reshape([self.maxreadlength, len(BASES), len(BASES)])


def iter_cigar(cigar, pat=re.compile('(\d+)([A-Z])')):
    for cnt, cmd in pat.findall(cigar):
        yield (int(cnt), cmd)


def count_mutation_profile(cntarray, cigar, readseq, refseq, strand, ndup):
    readpos = refpos = 0

    readseq = readseq.upper()
    refseq = refseq.upper()

    icigar = iter_cigar(cigar)

    for cmdrepeat, command in (icigar if strand == '+' else reversed(list(icigar))):
        if command == 'M':
            read_fragment = readseq[readpos:readpos+cmdrepeat]
            ref_fragment = refseq[refpos:refpos+cmdrepeat]

            for i, (rdbase, rfbase) in enumerate(zip(read_fragment, ref_fragment)):
                cntarray.countup(readpos+i, rfbase, rdbase, ndup)

            readpos += cmdrepeat
            refpos += cmdrepeat
        elif command == 'N':
            refpos += cmdrepeat
        elif command == 'I':
            for i, rdbase in enumerate(readseq[readpos:readpos+cmdrepeat]):
                cntarray.countup(readpos+i, '-', rdbase, ndup)
            readpos += cmdrepeat
        elif command == 'D':
            for rfbase in refseq[refpos:refpos+cmdrepeat]:
                cntarray.countup(readpos, rfbase, '-', ndup)
            refpos += cmdrepeat
        else:
            raise ValueError('Unhandled cigar command %d%s' % (cmdrepeat, command))


def count_alignment_mutations(bamfile, seqfile, maxlength):
    samfile = os.popen('samtools view -h "%s"' % bamfile)
    parser = sam.SAMParser(samfile, seqfile=seqfile)

    cntarray = ReadMappingProfile(maxlength)

    for line in parser.iteralignments(withref=True, repl_softclip=True):
        assert len(line['mapped']) == 1

        aln = line['mapped'][0]
        ndup = int(line['qname'].split('-')[1])
        count_mutation_profile(cntarray, aln[5], line['seq'], aln[6], aln[3], ndup)

    return cntarray


def run_jobs(bamdir, reference, maxlength, maxworkers):
    with futures.ProcessPoolExecutor(maxworkers) as executor:
        bamfiles = glob.glob(os.path.join(bamdir, '*.bam'))

        runcount = partial(count_alignment_mutations, seqfile=reference, maxlength=maxlength)

        error_profile = reduce(lambda x, y: x + y,
                               (cnt.makematrix() for cnt in executor.map(runcount, bamfiles)))
        return error_profile


if __name__ == '__main__':
    import sys
    import pickle

    bamdir = sys.argv[1]
    reffasta = sys.argv[2]
    profile_out = sys.argv[3]
    maxlength = int(sys.argv[4])
    maxworkers = int(sys.argv[5])

    error_profile = run_jobs(bamdir, reffasta, maxlength, maxworkers)
    pickle.dump(error_profile, open(profile_out, 'w'))

