#!/usr/bin/env python
from __future__ import division
from ecliptic.support import sam
import os
import re
from collections import defaultdict
from functools import partial
import futures
import glob
import numpy as np
import struct

BASES = 'ACGT-'
PARTITIONS = 12 # divide position in a read as many as PARTITIONS

class ReadMappingProfile(object):

    DIBASES = [n1+n2 for n1 in BASES for n2 in BASES]
    DIBASE2NUM = dict((bases, i) for i, bases in enumerate(DIBASES))

    def __init__(self, maxreadlength):
        self.cntarray = np.zeros([maxreadlength, len(self.DIBASES)], np.uint64)
        self.cntarray_r = np.zeros([maxreadlength, len(self.DIBASES)], np.uint64)
        self.divcntarray = np.zeros([PARTITIONS, len(self.DIBASES)], np.uint64)
        self.maxreadlength = maxreadlength

    def countup(self, readlength, pos, refbase, readbase, ndup):
        dibaseidx = self.DIBASE2NUM.get(refbase + readbase)
        if dibaseidx is not None:
            self.cntarray[pos, dibaseidx] += ndup
            self.cntarray_r[readlength - 1 - pos, dibaseidx] += ndup
            partpos = min(int(pos / readlength * PARTITIONS), PARTITIONS-1)
            self.divcntarray[partpos, dibaseidx] += ndup

    def flatten(self):
        return np.hstack([self.cntarray.flatten(),
                          self.cntarray_r.flatten(),
                          self.divcntarray.flatten()])

    @classmethod
    def unpack(kls, maxreadlength, summarized):
        cntarraysize = maxreadlength * len(kls.DIBASES)
        divcntarraysize = PARTITIONS * len(kls.DIBASES)
        assert len(summarized) == cntarraysize * 2 + divcntarraysize

        # first: read position, second: from, third: to
        cntarray = summarized[:cntarraysize]
        cntarray = cntarray.reshape([maxreadlength, len(BASES), len(BASES)])
        cntarray_r = summarized[cntarraysize:cntarraysize*2]
        cntarray_r = cntarray_r.reshape([maxreadlength, len(BASES), len(BASES)])

        # first: partition position, second: from, third: to
        divcntarray = summarized[cntarraysize*2:len(summarized)]
        divcntarray = divcntarray.reshape([PARTITIONS, len(BASES), len(BASES)])

        return cntarray, cntarray_r, divcntarray


def iter_cigar(cigar, pat=re.compile('(\d+)([A-Z])')):
    for cnt, cmd in pat.findall(cigar):
        yield (int(cnt), cmd)


def count_mutation_profile(cntarray, cigar, readseq, refseq, strand, ndup):
    readpos = refpos = 0

    readseq = readseq.upper()
    refseq = refseq.upper()
    readseqlen = len(readseq)

    icigar = iter_cigar(cigar)

    for cmdrepeat, command in (icigar if strand == '+' else reversed(list(icigar))):
        if command == 'M':
            read_fragment = readseq[readpos:readpos+cmdrepeat]
            ref_fragment = refseq[refpos:refpos+cmdrepeat]

            for i, (rdbase, rfbase) in enumerate(zip(read_fragment, ref_fragment)):
                cntarray.countup(readseqlen, readpos+i, rfbase, rdbase, ndup)

            readpos += cmdrepeat
            refpos += cmdrepeat
        elif command == 'N':
            refpos += cmdrepeat
        elif command == 'I':
            for i, rdbase in enumerate(readseq[readpos:readpos+cmdrepeat]):
                cntarray.countup(readseqlen, readpos+i, '-', rdbase, ndup)
            readpos += cmdrepeat
        elif command == 'D':
            for rfbase in refseq[refpos:refpos+cmdrepeat]:
                cntarray.countup(readseqlen, readpos, rfbase, '-', ndup)
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

    return cntarray.flatten()


def run_jobs(bamdir, reference, maxlength, maxworkers):
    with futures.ProcessPoolExecutor(maxworkers) as executor:
        bamfiles = glob.glob(os.path.join(bamdir, '*.bam'))

        runcount = partial(count_alignment_mutations, seqfile=reference, maxlength=maxlength)

        error_profile = reduce(lambda x, y: x + y, executor.map(runcount, bamfiles))

        return ReadMappingProfile.unpack(maxlength, error_profile)


def write_binary_profile(profile, output):
    output.write(struct.pack('II', profile.shape[0], profile.shape[1]))
    profile.tofile(output)


if __name__ == '__main__':
    import sys
    import pickle

    bamdir = sys.argv[1]
    reffasta = sys.argv[2]
    profile_out_prefix = sys.argv[3]
    maxlength = int(sys.argv[4])
    maxworkers = int(sys.argv[5])

    error_profile = run_jobs(bamdir, reffasta, maxlength, maxworkers)
    pickle.dump(error_profile, open(profile_out_prefix + 'pickle', 'w'))
    write_binary_profile(error_profile[0], open(profile_out_prefix + 'arraydump', 'w'))

