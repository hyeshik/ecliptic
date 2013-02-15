#!/usr/bin/env python
#
# Copyright (c) 2011 Seoul National University RNA Biology Laboratory.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Hyeshik Chang <hyeshik@snu.ac.kr>
#

from collections import defaultdict
from ecliptic.support import sam
import gzip
import os
import csv

PRIORITY = """
miRNA rRNA tRNA Mt-tRNA snoRNA scRNA srpRNA snRNA lncRNA RNA ncRNA misc_RNA
Cis-reg ribozyme RC IRES frameshift_element
LINE SINE Simple_repeat Low_complexity Satellite DNA LTR CDS UTR3 UTR5
intron Other Unknown
""".split()
PRIORITYMAP = defaultdict(lambda: len(PRIORITYMAP))
#PRIORITYMAP = {}
PRIORITYMAP.update(dict((name, no) for no, name in enumerate(PRIORITY)))


class IntersectBedIterator(object):

    def __init__(self, inpfile):
        self.f = csv.reader(inpfile, dialect='excel-tab')
        self.queued = None
        self.queued_seqid = None

    def get(self, key):
        r = []

        seqid = int(key.split('-')[0])
        if self.queued_seqid is not None:
            if self.queued_seqid == seqid:
                r.append(self.queued)
            elif self.queued_seqid > seqid:
                return r

            self.queued = self.queued_seqid = None

        for fields in self.f:
            annoid = int(fields[3].split('-')[0])
            if annoid == seqid:
                r.append(fields)
            elif annoid > seqid:
                self.queued = fields
                self.queued_seqid = annoid
                break

        return [fields[15].split('|')[:2] + [int(fields[4])] for fields in r]


def removedup(L):
    r = []
    added = set()
    for v in L:
        if v not in added:
            r.append(v)
            added.add(v)
    return r


def open_bam(filepath):
    return os.popen('samtools view -h "{}"'.format(filepath))


def process(alnfile, annofile):
    iterintersection = IntersectBedIterator(gzip.open(annofile))
    alninp  = sam.SAMParser(open_bam(alnfile))

    for aln in alninp:
        seqname = aln['qname']
        annotations = iterintersection.get(seqname)
        annotations.sort(key=lambda (cls, name, score):
                            (PRIORITYMAP[cls], -score, name))

        if annotations:
            minscore = annotations[0][2] - 2
            annotations = filter(lambda (cls, name, score): score >= minscore,
                                 annotations)

        topanno = annotations[0][0] if annotations else ''

        print '%s\t%d\t%s\t%s' % (seqname, len(aln['mapped']), topanno,
                ','.join(removedup(['%s:%s' % (cls, name)
                                    for cls, name, score in annotations])))


if __name__ == '__main__':
    import sys

    alignments = sys.argv[1]
    annotations = sys.argv[2]

    process(alignments, annotations)

