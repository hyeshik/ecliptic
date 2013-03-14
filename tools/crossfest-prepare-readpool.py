#!/usr/bin/env python
from ecliptic.support import sam, sequtils
from ecliptic.Utils import TemporaryDirectory
import os
import sys
from struct import pack, unpack
import glob
import futures
import gzip
import re


PATTERN_NONVALID_BASE = re.compile('[^ACGT]')
MAX_SPAN_BLOCKS = 16
MIN_CLUSTER_INTERVAL = 200

class SAMToReadPoolConverter(object):

    def __init__(self, samfile, seqfile, output):
        self.samfile = samfile
        self.seqs = sequtils.GiantFASTAFile(seqfile)
        self.output = output

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

    def write_serialized_cluster(self, readlist, strand):
        leftend = min(read[1] for read in readlist)
        rightend = max(read[2] for read in readlist)
        refseq = self.get_refseq(readlist[0][0], leftend, rightend)

        packed_readlist = []

        for chrom, start, stop, cigar, nreads in readlist:
            relstart, relstop = start - leftend, stop - leftend
            spanblocks = list(self.make_compact_spanning(relstart, cigar, strand))

            if len(spanblocks) > MAX_SPAN_BLOCKS:
                print >> sys.stderr, "Oops! Too many span blocks for %s:%d-%d %s - %s" % (
                        chrom, start, stop, cigar, spanblocks)

            # Skip sequences where reference include nonstandard bases.
            if any(PATTERN_NONVALID_BASE.search(refseq[blkstart:blkstop])
                    for blkstart, blkstop in spanblocks):
                continue

            formatted = pack('IB', nreads, len(spanblocks)) + ''.join(
                            pack('II', blkstart, blkstop) for blkstart, blkstop in spanblocks)

            #print chrom, start, stop, cigar, spanblocks
            #print '/'.join(self.get_refseq(chrom, leftend+blkstart, leftend+blkstop)
            #               for blkstart, blkstop in spanblocks)
            packed_readlist.append(formatted)

        self.output.write(readlist[0][0].ljust(32, '\0'))
        self.output.write(pack('iiII', leftend, rightend, len(packed_readlist), len(refseq)))
        self.output.write(refseq)
        for packed in packed_readlist:
            self.output.write(packed)

        return (leftend, rightend, rightend-leftend)

    def make_compact_spanning(self, leftend, cigar, strand, cigar_pattern=sam.cigar_pattern):
        cigar_tokens = cigar_pattern.findall(cigar)
        if strand == '-':
            cigar_tokens = cigar_tokens[::-1]
        frags = []
        refleft = leftend
        readleft = 0

        for width, cmdtype in cigar_tokens:
            width = int(width)

            if cmdtype in 'MD': # M or D spans the reference RNA
                if frags and frags[-1][1] == 'M':
                    frags[-1][0] += width # merge if the block is continued
                else:
                    frags.append([width, 'M'])
            elif cmdtype == 'I': # I doesn't exist in reference, ignore it
                pass
            elif cmdtype == 'N':
                if frags and frags[-1][1] == 'N':
                    frags[-1][0] += width
                else:
                    frags.append([width, 'N'])
            else:
                raise ValueError("Unprocessed cigar command in %s." % cigar)

        currleft = leftend
        for i, (width, cmdtype) in enumerate(frags):
            nextleft = currleft + width
            if i & 1 == 0:
                yield (currleft, nextleft) # yield matches only.
            currleft = nextleft

    def convert_to_readpool(self, strand):
        clustered = []
        clustered_limit = None
        clustered_chrom = None

        for l in self.iter_sam(strand):
            assert len(l['mapped']) == 1

            nreads = int(l['qname'].split('-')[1])
            chrom, start, stop, strand, mismatches, cigar = l['mapped'][0]
            if strand == '-':
                start, stop = 1 - stop, 1 - start
            markset = (chrom, start, stop, cigar, nreads)

            if not clustered: # at the beginning
                clustered.append(markset)
                clustered_chrom = chrom
                clustered_limit = stop + MIN_CLUSTER_INTERVAL
            elif chrom != clustered_chrom or clustered_limit < start: # start new cluster
                self.write_serialized_cluster(clustered, strand)
                clustered = [markset]
                clustered_chrom = chrom
                clustered_limit = stop + MIN_CLUSTER_INTERVAL
            else: # continue stacking to the current cluster
                clustered.append(markset)
                clustered_limit = max(clustered_limit, stop + MIN_CLUSTER_INTERVAL)

        if clustered:
            self.write_serialized_cluster(clustered, strand)


def run_job(samfile, refseqfile, outputfile, strand):
    output = open(outputfile, 'w')

    try:
        converter = SAMToReadPoolConverter(samfile, refseqfile, output)
        converter.convert_to_readpool(strand)
    except:
        os.unlink(outputfile)
        import traceback
        traceback.print_exc()
        return 1

    return 0


def merge_outputs(tmpdir, output, threads=8):
    inputfiles = ' '.join("'%s'" % fname for fname in sorted(os.listdir(tmpdir)))

    os.system("(cd '{tmpdir}' && cat {inputfiles}) | pigz -p {threads} -c > '{output}'".format(
                tmpdir=tmpdir, inputfiles=inputfiles, threads=threads, output=output))


def run_parallel(bamsplitdir, refseqfile, finaloutput, maxworkers=8):
    allfutures = []
    files = glob.glob(os.path.join(bamsplitdir, '*.bam'))

    with futures.ProcessPoolExecutor(maxworkers) as executor, TemporaryDirectory() as tmpdir:
        for fname in sorted(files, key=os.path.getsize, reverse=True):
            chrom = os.path.basename(fname).replace('.bam', '')
            for strand in '+-':
                outputpath = os.path.join(tmpdir, chrom + strand)
                fobj = executor.submit(run_job, fname, refseqfile, outputpath, strand)
                allfutures.append(fobj)

        nfailures = sum(f.result() for f in allfutures)

        if all(f.done() for f in allfutures) and nfailures < 1:
            merge_outputs(tmpdir, finaloutput, max(2, maxworkers // 4))
        else:
            print >> sys.stderr, "Error on processing something."
            sys.exit(1)


if __name__ == '__main__':
    import sys

    bamsplitdir = sys.argv[1]   # directory with many bam files
    refseqfile = sys.argv[2]    # in fasta format
    readpool = sys.argv[3]      # output file (new)
    threads = int(sys.argv[4])  # number of threads (or processes)

    run_parallel(bamsplitdir, refseqfile, readpool, threads)

