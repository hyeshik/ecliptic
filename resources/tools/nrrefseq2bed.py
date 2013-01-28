#!/usr/bin/env python
import shelve
import sys

def gen_bed_entries(flatentry):
    chrom = flatentry['chrom']
    strand = flatentry['strand']
    name = flatentry['name']

    if name.startswith('NR_'):
        for i, (blkstart, blkend) in enumerate(flatentry['exonBlocks']):
            yield (chrom, blkstart, blkend,
                   '%s/%d.e' % (name, i+1), '.', strand)
        return

    cdsstart, cdsend = flatentry['cdsStart'], flatentry['cdsEnd']

    if strand == '+':
        leftutr, rightutr = '5', '3'
    else:
        leftutr, rightutr = '3', '5'

    for i, (blkstart, blkend) in enumerate(flatentry['exonBlocks']):
        if cdsstart <= blkstart and blkend <= cdsend: # cds
            yield (chrom, blkstart, blkend,
                    '%s/%d.c' % (name, i+1), '.', strand)
            continue

        eff_cdsstart = min(blkend, max(blkstart, cdsstart))
        eff_cdsend = max(blkstart, min(blkend, cdsend))

        if blkstart < eff_cdsstart:
            yield (chrom, blkstart, eff_cdsstart,
                    '%s/%d.%s' % (name, i+1, leftutr), '.', strand)
        if eff_cdsstart < eff_cdsend:
            yield (chrom, eff_cdsstart, eff_cdsend,
                    '%s/%d.c' % (name, i+1), '.', strand)
        if eff_cdsend < blkend:
            yield (chrom, eff_cdsend, blkend,
                    '%s/%d.%s' % (name, i+1, rightutr), '.', strand)


db = shelve.open(sys.argv[1], 'r')
r = []
for k in db.itervalues():
    for fields in gen_bed_entries(k):
        r.append(fields)

r.sort()
for fields in r:
    print '\t'.join(map(str, fields))
