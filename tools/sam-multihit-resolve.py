#!/usr/bin/env python
# emulates H. Guo's mapping strategy
from rnarry.sequtils import GiantFASTAFile
import re
import numpy as np

cigarpat = re.compile('(\d+)([A-Z])')

def read_sam_blks(inpf):
    stacked = []
    stackid = None
    for line in inpf:
        fields = line[:-1].split('\t')
        if not stacked or stackid != fields[0]:
            if stacked:
                yield stacked
            del stacked[:]
            stackid = fields[0]

        stacked.append(fields)

    if stacked:
        yield stacked

def get_alignment_mismatches(fields):
    clipping = sum(int(num)
             for num, tok in cigarpat.findall(fields[5]) if tok in 'HSID')
                # gsnap underestimates mismatches excluding all of those things
    nm = [int(tok.split(':')[2])
          for tok in fields[11:] if tok.startswith('NM:i:')]
    assert len(nm) == 1
    return clipping + nm[0]

def process(inpf, maxmm=1):
    for blk in read_sam_blks(inpf):
        if blk[0][0][0] == '@':
            for fields in blk:
                print '\t'.join(fields)
            continue

        if len(blk) == 1:
            if blk[0][2] != '*': # ignore unmapped reads
                if get_alignment_mismatches(blk[0]) <= 1:
                    print '\t'.join(blk[0])
            continue

        mismatches = np.array(map(get_alignment_mismatches, blk))

        for i in range(maxmm+1):
            mmaccepted = np.where(mismatches == i)[0]
            if len(mmaccepted) == 1:
                print '\t'.join(blk[int(mmaccepted[0])])
                break
            elif len(mmaccepted) > 1:
                break

if __name__ == '__main__':
    import sys

    maxmm = 1
    if len(sys.argv) > 1:
        maxmm = int(sys.argv[1])

    process(sys.stdin, maxmm)

