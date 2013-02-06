#!/usr/bin/env python
from bx.binned_array import BinnedArray
import os
import sys
import random
import futures
import gzip
import re
import os
from ecliptic.support import sam
from ecliptic.support.sequtils import reverse_complement


ENDTRIM = 2
cigar_pattern = re.compile('(\d+)([MIDNSHP])')
COMPTYPE = 'lzo'
TYPECODE = 'I' # maximum: 4,294,967,295 ~= 21x of maximum throughput from a HiSeq lane
newarray = lambda: BinnedArray(default=0, typecode=TYPECODE)
BASESPACE = ['_', '5', '3', 'A', 'C', 'G', 'T', 'I', 'D', 'N'] # '_' for depth (any mis/match)
BASESPACE2i = dict((n, i) for i, n in enumerate(BASESPACE))


def open_with_newdir(filename):
    dname = os.path.dirname(filename)
    if not os.path.isdir(dname):
        os.makedirs(dname)
    return open(filename, 'w')


def genrefblks(readseq, chrom, start, stop, strand, cigar, nreads):
    refpos = start
    readpos = 0
    if strand == '-':
        readseq = reverse_complement(readseq)

    tleftlim, trightlim = start + ENDTRIM, stop - ENDTRIM
    qleftlim, qrightlim = ENDTRIM, len(readseq) - ENDTRIM

    cigarcommands = cigar_pattern.findall(cigar)
    if cigarcommands[0][1] == 'S': # shift start site for the first soft clipping
        start -= int(cigarcommands[0][0])
    if len(cigarcommands) > 1 and cigarcommands[-1][1] == 'S': # last soft clipping
        stop += int(cigarcommands[-1][0])

    for num, cmd in cigarcommands:
        num = int(num)
        if cmd == 'M': # match
            mleft = max(qleftlim, readpos)
            mright = min(qrightlim, readpos + num)
            if mleft < mright:
                seq = readseq[mleft:mright]
                yield ('M', nreads, max(refpos, tleftlim), seq)
            refpos += num
            readpos += num
        elif cmd == 'S': # soft clip
            readpos += num
        elif cmd == 'N': # skip
            refpos += num
        elif cmd == 'D': # deletion
            if tleftlim <= refpos < trightlim:
                yield ('D', nreads, refpos, num)
            refpos += num
        elif cmd == 'I': # insertion
            ppos = (refpos if strand == '+' else (refpos-1))
            if tleftlim <= ppos < trightlim:
                yield ('I', nreads, ppos, num)
            readpos += num
        elif cmd == 'H': # hard clipping
            pass
        else:
            yield ('E', nreads, num, cmd, readseq)
            raise ValueError

    if strand == '+':
        fivep, threep = start, stop-1
    else:
        fivep, threep = stop-1, start

    yield ('5', nreads, fivep)
    yield ('3', nreads, threep)


def translate_sam(infile, strand):
    filteropt = '-f 16' if strand == '-' else '-F 16'
    inpf = os.popen('samtools view %s "%s"' % (filteropt, infile))
    insam = sam.SAMParser(inpf, zerobase=True)

    for aln in insam:
        nreads = int(aln['qname'].split('-')[1])
        for chrom, start, stop, strand, mm, cigar in aln['mapped']:
            for msg in genrefblks(aln['seq'], chrom, start, stop, strand, cigar, nreads):
                yield msg


def process(inputfile, prefix, strand):
    #print "%s starting." % os.path.basename(prefix)

    try:
        arr = [newarray() for i in BASESPACE]

        for fields in translate_sam(inputfile, strand):
            nread = int(fields[1])
            pos = int(fields[2])
            mtype = fields[0]

            if mtype == 'M':
                for relp, base in enumerate(fields[3]):
                    p = pos + relp
                    barr = arr[BASESPACE2i[base]]
                    arr[0].set(p, arr[0].get(p) + nread) # depth count
                    barr.set(p, barr.get(p) + nread) # per-base count
            elif mtype in '53ID':
                evarr = arr[BASESPACE2i[mtype]]
                evarr.set(pos, evarr.get(pos) + nread)
            else:
                raise ValueError("Unknown message type %s" % repr(mtype))

        for subtype, a in zip(BASESPACE, arr):
            a.to_file(open_with_newdir(prefix + subtype), COMPTYPE)
    except:
        import traceback
        traceback.print_exc()
        raise

    #print "%s finished." % os.path.basename(prefix)


def run_jobs(flatdir, dbdir, maxworkers):
    allfutures = []

    with futures.ProcessPoolExecutor(maxworkers) as executor:
        files = os.listdir(flatdir)
        random.shuffle(files) # randomize order to balance memory loads

        for f in files:
            inputfile = os.path.join(flatdir, f)

            for strand in '+-':
                outputprefix = os.path.join(dbdir, f.split('.')[1] + strand)
                fobj = executor.submit(process, inputfile, outputprefix, strand)
                allfutures.append(fobj)

        futures.wait(allfutures)

        if all(f.done() for f in allfutures):
            open('%s/done' % dbdir, 'w')
        else:
            print >> sys.stderr, "Error on processing something."


if __name__ == '__main__':
    import sys

    run_jobs(sys.argv[2], sys.argv[3], int(sys.argv[1]))
