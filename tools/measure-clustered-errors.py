#!/usr/bin/env python
from __future__ import division
from ecliptic.support import sam, sequtils
from ecliptic.Utils import TemporaryDirectory
import os
import sys
from struct import pack, unpack
import numpy as np
import glob
from itertools import groupby
import futures
import gzip
import re

BASES = 'ACGT-'
BASE2i = dict((b, i) for i, b in enumerate(BASES))
PATTERN_NONVALID_BASE = re.compile('[^ACGT]')
MAX_SPAN_BLOCKS = 16
MIN_CLUSTER_INTERVAL = 200


def shannon_entropy(counts):
    probs = counts[counts > 0] / counts.sum()
    return abs((np.log(probs) * probs).sum())

class SAMMeasureClusteredError(object):

    def __init__(self, samfile, seqfile, outnzname, outmname, outdname, outmdname, outename,
                       outt2cname):
        self.samfile = samfile
        self.seqs = sequtils.GiantFASTAFile(seqfile)

        # intermediate stats output for full non-zero read coordinates
        self.outnz = self.open_nonzero_raw_output(outnzname)

        # intermediate summarized outputs for deletion rate, modification rate or entropy
        self.outm = self.open_intermediate_summary_output(outmname)
        self.outd = self.open_intermediate_summary_output(outdname)
        self.outmd = self.open_intermediate_summary_output(outmdname)
        self.oute = self.open_intermediate_summary_output(outename)
        self.outt2c = self.open_intermediate_summary_output(outt2cname)

    def iter_sam(self, strand):
        filterarg = '-f 16' if strand == '-' else '-F 16'
        f = os.popen('samtools view {filter} {filename} | sort -k3,3 -k4,4n{strand}'.format(
                        filter=filterarg, filename=self.samfile,
                        strand={'+': '', '-': 'r'}[strand]))
        return sam.SAMParser(f).iteralignments(repl_softclip=True)

    def get_refseq(self, chrom, start, stop): # with automatic (-) strand detection
        if start >= 0:
            return self.seqs.get(chrom, start, stop).upper()

        start, stop = 1 - stop, 1 - start
        return sequtils.reverse_complement(self.seqs.get(chrom, start, stop).upper())

    def open_intermediate_summary_output(self, filename):
        return os.popen("sort -k1,1n -k2,2 -k3,3r | uniq -c >%s" % filename, 'w')

    def open_nonzero_raw_output(self, filename):
        return gzip.open(filename, 'w')

    def close(self):
        [f.close() for f in (self.outnz, self.outd, self.outm, self.outmd,
                             self.oute, self.outt2c)]

    def countup_seqreads(self, start, stop, refseq, readseq, cigar, nreads, posreads, strand,
                       cigar_pattern=sam.cigar_pattern):
        cigar_tokens = cigar_pattern.findall(cigar)
        if strand == '-':
            cigar_tokens = cigar_tokens[::-1]
        refleft = start
        readleft = 0

        for width, cmdtype in cigar_tokens:
            width = int(width)

            if cmdtype == 'M': # M spans the reference RNA
                for refi, readbase in zip(range(refleft, refleft+width),
                                          readseq[readleft:readleft+width]):
                    if readbase in BASE2i:
                        posreads[refi, BASE2i[readbase]] += nreads

                readleft += width
                refleft += width
            elif cmdtype == 'D': # deletions
                dbase = BASE2i['-']
                for refi in range(refleft, refleft+width):
                    posreads[refi, dbase] += nreads

                refleft += width
            elif cmdtype == 'I': # I doesn't exist in reference, ignore it
                readleft += width
            elif cmdtype == 'N':
                refleft += width
            else:
                raise ValueError("Unprocessed cigar command in %s." % cigar)

    def calculate_single(self, readlist, strand):
        leftend = min(read[1] for read in readlist)
        rightend = max(read[2] for read in readlist)
        refseq = self.get_refseq(readlist[0][0], leftend, rightend)
        posreads = np.zeros([len(refseq), 5], 'I')

        for chrom, start, stop, cigar, nreads, readseq in readlist:
            relstart, relstop = start - leftend, stop - leftend
            self.countup_seqreads(relstart, relstop, refseq, readseq, cigar,
                                  nreads, posreads, strand)

        chrom = readlist[0][0]
        nonzeropositions = np.where(posreads.sum(1) > 0)[0]

        for pos in nonzeropositions:
            refbase = refseq[pos]
            simcnt = posreads[pos]
            kcnt = simcnt.sum()
            del_ratio = simcnt[4] / kcnt
            if refbase in BASE2i:
                refcnt = simcnt[BASE2i[refbase]]
                moddel_ratio = (kcnt - refcnt) / kcnt
                mod_ratio = (kcnt - refcnt - simcnt[4]) / kcnt
            else:
                mod_ratio = moddel_ratio = 0.

            if refbase == 'T':
                t2c_ratio = simcnt[1] / kcnt
            else:
                t2c_ratio = 0.

            sentropy = shannon_entropy(simcnt)

            if 0 < del_ratio < 0.9 or 0 < mod_ratio < 0.9:
                abspos = leftend + pos if strand == '+' else 1 - leftend - pos
                print >> self.outnz, chrom, strand, abspos, refbase, kcnt, ' '.join(
                 map(str, simcnt)), del_ratio, mod_ratio, moddel_ratio, sentropy

            print >> self.outd, kcnt, refbase, '%.6f' % del_ratio
            print >> self.outm, kcnt, refbase, '%.6f' % mod_ratio
            print >> self.outmd, kcnt, refbase, '%.6f' % moddel_ratio
            print >> self.oute, kcnt, refbase, '%.6f' % sentropy
            print >> self.outt2c, kcnt, refbase, '%.6f' % t2c_ratio

    def process(self, strand):
        clustered = []
        clustered_limit = None
        clustered_chrom = None

        for l in self.iter_sam(strand):
            assert len(l['mapped']) == 1

            nreads = int(l['qname'].split('-')[1])
            chrom, start, stop, strand, mismatches, cigar = l['mapped'][0]
            if strand == '-':
                start, stop = 1 - stop, 1 - start
            markset = (chrom, start, stop, cigar, nreads, l['seq'])

            if not clustered: # at the beginning
                clustered.append(markset)
                clustered_chrom = chrom
                clustered_limit = stop + MIN_CLUSTER_INTERVAL
            elif chrom != clustered_chrom or clustered_limit < start: # start new cluster
                self.calculate_single(clustered, strand)
                clustered = [markset]
                clustered_chrom = chrom
                clustered_limit = stop + MIN_CLUSTER_INTERVAL
            else: # continue stacking to the current cluster
                clustered.append(markset)
                clustered_limit = max(clustered_limit, stop + MIN_CLUSTER_INTERVAL)

        if clustered:
            self.calculate_single(clustered, strand)


def run_job(samfile, refseqfile, outputprefix, strand):
    output_files = [outputprefix + '-' + suffix
                    for suffix in ['nonzero', 'mod', 'del', 'moddel', 'entropy', 't2c']]

    try:
        proc = SAMMeasureClusteredError(samfile, refseqfile, *output_files)
        proc.process(strand)
        proc.close()
    except:
        for fname in output_files:
            if os.path.exists(fname):
                os.unlink(fname)

        import traceback
        traceback.print_exc()
        return 1

    return 0


def merge_outputs(tmpdir, outputprefix, threads=8):
    # merge nonzero position dumps
    os.system("(cd '{tmpdir}' && zcat *-nonzero) | sort -k1,2 -k3,3n -k4,4 "
              "| pigz -p {threads} -c > '{output}'".format(
                tmpdir=tmpdir, threads=threads, output=outputprefix+'nonzero.real.gz'))

    # merge metric-specific summarized outputs
    for suffix in ['mod', 'del', 'moddel', 'entropy', 't2c']:
        cmdout = os.popen("cd '{tmpdir}' && sort -k2,2n -k3,3 -k4,4n *-{suffix}".format(
                          tmpdir=tmpdir, suffix=suffix))
        cmdoutiter = (line.split(None, 1) for line in cmdout)
        mergedoutput = gzip.open(outputprefix + suffix + '.real.gz', 'w')

        for data, grp in groupby(cmdoutiter, key=lambda v: v[1]):
            totalcnt = sum(int(cnt) for cnt, _ in grp)
            print >> mergedoutput, data[:-1], totalcnt


def run_parallel(bamsplitdir, refseqfile, finaloutputprefix, maxworkers=8):
    allfutures = []
    files = glob.glob(os.path.join(bamsplitdir, '*.bam'))

    with futures.ProcessPoolExecutor(maxworkers) as executor, TemporaryDirectory() as tmpdir:
        for fname in sorted(files, key=os.path.getsize, reverse=True):
            chrom = os.path.basename(fname).replace('.bam', '')
            for strand in '+-':
                outputprefix = os.path.join(tmpdir, chrom + strand)
                fobj = executor.submit(run_job, fname, refseqfile, outputprefix, strand)
                allfutures.append(fobj)

        nfailures = sum(f.result() for f in allfutures)

        if all(f.done() for f in allfutures) and nfailures < 1:
            merge_outputs(tmpdir, finaloutputprefix, max(2, maxworkers // 4))
        else:
            print >> sys.stderr, "Error on processing something."
            sys.exit(1)


if __name__ == '__main__':
    import sys

    bamsplitdir = sys.argv[1]   # directory with many bam files
    refseqfile = sys.argv[2]    # in fasta format
    outputprefix = sys.argv[3]  # output prefix
    threads = int(sys.argv[4])  # number of threads (or processes)

    run_parallel(bamsplitdir, refseqfile, outputprefix, threads)

