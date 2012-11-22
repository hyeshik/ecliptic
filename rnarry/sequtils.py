# Copyright (c) 2011-2012 Hyeshik Chang
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

__all__ = [
    'iterseq_sequential',
    'GiantFASTAFile',
    'reverse_complement',
    'iter_windowed',
    'iter_windowed_str',
    'get_first_sequence_length',
]

from Bio import SeqIO
import string
import gzip
from collections import deque
import os
import re

def iterseq_sequential(fastapath):
    nextid = yield
    for seq in SeqIO.parse(open(fastapath, 'r'), format='fasta'):
        while seq.name == nextid:
            nextid = yield seq

import pysam
# provides an interface like 'samtools faidx'.
whitespace = re.compile('[ \t\r\n]')
class GiantFASTAFile(object):

    def __init__(self, filename):
        if not os.path.exists(filename + '.fai'):
            pysam.faidx(filename)

        self.fasta = open(filename)
        self.index = self.load_index(filename + '.fai')

    def load_index(self, filename):
        index = {}
        for line in open(filename):
            fields = line[:-1].split('\t')
            index[fields[0]] = tuple(map(int, fields[1:]))
        return index

    def get(self, seqid, start=None, stop=None): # zero-based, half-open
        length, filepos, colwidth, linesize = self.index[seqid]

        if start is None and stop is None:
            offset_st = filepos
            linenum_en = length // colwidth
            offset_en = filepos + length + linenum_en * (linesize - colwidth)
        else:
            linenum_st = start // colwidth
            offset_st = filepos + start + linenum_st * (linesize - colwidth)
            linenum_en = stop // colwidth
            offset_en = filepos + stop + linenum_en * (linesize - colwidth)

        self.fasta.seek(offset_st, 0)
        return whitespace.sub('', self.fasta.read(offset_en - offset_st))


revcmptrans = string.maketrans('ATUGCatugc', 'TAACGtaacg')
def reverse_complement(seq):
    return seq.translate(revcmptrans)[::-1]

def iter_windowed(it, width):
    i = iter(it)

    queue = deque()
    for _ in range(width):
        queue.append(i.next())

    yield tuple(queue)

    for next in i:
        queue.popleft()
        queue.append(next)
        yield tuple(queue)

def iter_windowed_str(it, width):
    for r in iter_windowed(it, width):
        yield ''.join(r)

def get_first_sequence_length(path, format='fastq-illumina', gzipped=True):
    if gzipped:
        opener = gzip.open
    else:
        opener = open
    return len(SeqIO.parse(opener(path), format=format).next().seq)

