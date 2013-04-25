#!/usr/bin/env python
import shelve
import csv
import sys
import os
import tempfile
import bsddb
from itertools import groupby
from subprocess import PIPE, Popen
from collections import defaultdict

def read_nr_alignments(pslinput, nrlist, nrbedfile):
    inp = os.popen('sort -k10,10 {}'.format(pslinput))
    ambiggroups = {}
    for acc, line in groupby(inp, lambda L: L.split('\t')[9]):
        if acc not in nrlist:
            continue

        lines = list(line)
        if len(lines) == 1:
            yield (acc, lines[0][:-1])
            continue

        # multiple alignments available for a single RefSeq entry
        rows = list(csv.reader(lines, dialect='excel-tab'))
        ambiggroups[acc] = rows

    # overlap with nrRefSeq-genome.bed for the strict consistency
    overlapcnt = defaultdict(lambda: defaultdict(int))

    with tempfile.NamedTemporaryFile(suffix='.bed') as tmpbed:
        for acc, grpregions in ambiggroups.iteritems():
            for i, region in enumerate(grpregions):
                chrom, start, end, strand = region[13], region[15], region[16], region[8]
                pos = '{}\t{}\t{}\t{}/{}\t.\t{}'.format(chrom, start, end, acc, i, strand)
                print >> tmpbed, pos
        tmpbed.flush()

        with os.popen('bedtools intersect -wa -wb -a {} -b {}'.format(
                              tmpbed.name, nrbedfile)) as ch_stdout:
            for row in csv.reader(ch_stdout, dialect='excel-tab'):
                grpacc, grpidx = row[3].split('/')
                if grpacc == row[9].split('/')[0]:
                    overlapcnt[grpacc][int(grpidx)] += int(row[2]) - int(row[1])

    for acc, matchlengths in overlapcnt.iteritems():
        bestidx, _ = max(matchlengths.iteritems(), key=lambda (idx, overlap): overlap)
        yield (acc, '\t'.join(ambiggroups[acc][bestidx]))


if __name__ == '__main__':
    nrlistfile = sys.argv[1]
    nrbedfile = sys.argv[2]
    pslinput = sys.argv[3]
    dboutput = sys.argv[4]

    nrlist = set(open(nrlistfile).read().split())
    db = bsddb.hashopen(dboutput, 'c')
    for acc, pslentry in read_nr_alignments(pslinput, nrlist, nrbedfile):
        db[acc] = pslentry
    db.sync()

