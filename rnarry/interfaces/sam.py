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

"""
>>> calculate_cigar_length('11M')
(11, 11)
>>> calculate_cigar_length('3S11M2S')
(16, 16)
>>> calculate_cigar_length('5M2I7M')
(12, 14)
>>> calculate_cigar_length('11M2D5M3I5M')
(23, 24)

>>> generate_cigar_map('3M')
{0: 0, 1: 1, 2: 2}
>>> generate_cigar_map('2M1D3M')
{0: 0, 1: 1, 2: None, 3: 2, 4: 3, 5: 4}
>>> generate_cigar_map('1M1I3M')
{0: 0, 1: 2, 2: 3, 3: 4}
>>> CIGARMAP_REF2READ['1M1I3M']
(4, 5, {0: 0, 1: 2, 2: 3, 3: 4})

>>> get_read_refalign('AAAATTCTCCTAACAGTGGACGAGA', '13M1D12M')
'AAAATTCTCCTAADCAGTGGACGAGA'
"""

__all__ = [
    'calculate_cigar_length', 'SAMParser', 'trim_cigar',
    'generate_cigar_map', 'CIGARCachedRefMap', 'CIGARMAP_READ2REF',
    'CIGARMAP_REF2READ',
    'get_read_refalign',
    'F_PAIRED', 'F_PROPER_PAIR', 'F_UNMAPPED', 'F_MATE_UNMAPPED',
    'F_REVERSE_STRAND', 'F_REVERSE_MATE_STRAND', 'F_FIRST_READ',
    'F_SECOND_READ', 'F_NOT_PRIMARY', 'F_QC_FAILURE', 'F_TECH_DUPLICATE',
]

import re
from itertools import groupby
from collections import deque
from rnarry.sequtils import reverse_complement, GiantFASTAFile

F_PAIRED            = 0x0001
F_PROPER_PAIR       = 0x0002
F_UNMAPPED          = 0x0004
F_MATE_UNMAPPED     = 0x0008
F_REVERSE_STRAND    = 0x0010
F_REVERSE_MATE_STRAND = 0x0020
F_FIRST_READ        = 0x0040
F_SECOND_READ       = 0x0080
F_NOT_PRIMARY       = 0x0100
F_QC_FAILURE        = 0x0200
F_TECH_DUPLICATE    = 0x0400

cigar_pattern = re.compile('(\d+)([MIDNSHP])')
def calculate_cigar_length(cigar):
    lengths = [0, 0] # reference, read

    def found(match):
        size, action = match.groups()
        size = int(size)

        if action in 'MSPDN':
            lengths[0] += size
        if action in 'MSPI':
            lengths[1] += size

        return ''

    remaining = cigar_pattern.sub(found, cigar)
    if remaining != '':
        raise ValueError, 'Unknown pattern included: ' + remaining

    return tuple(lengths)

def trim_cigar(cigar, length, reverse=False):
    flat = ''.join(int(n) * op for n, op in cigar_pattern.findall(cigar))
    if reverse:
        flat = ''.join(reversed(flat))

    clen = 0
    aggregated = []
    for c in flat:
        if c in 'MSPI':
            clen += 1

        aggregated.append(c)
        if clen >= length:
            break

    r = []
    for k, it in groupby(aggregated):
        r.append('%d%s' % (len(list(it)), k))

    if not reverse:
        return ''.join(r)
    else:
        return ''.join(reversed(r))

def generate_cigar_map(cigar, read2ref=False):
    current = [0, 0] # reference, read
    trmap = {}

    if read2ref:
        readfwd, reffwd = 'DI'
    else:
        readfwd, reffwd = 'ID'

    def found(match):
        size, action = match.groups()
        size = int(size)

        if action in 'MS':
            for i in range(size):
                trmap[current[0] + i] = current[1] + i
            current[0] += size
            current[1] += size
        elif action == readfwd:
            current[1] += size
        elif action == reffwd:
            for i in range(size):
                trmap[current[0] + i] = None
            current[0] += size
        else:
            raise ValueError, 'Currently supported CIGAR actions are: MID only - %d%s' % (size, action)

        return ''

    remaining = cigar_pattern.sub(found, cigar)
    if remaining != '':
        raise ValueError, 'Unknown pattern included: ' + remaining

    return trmap


class CIGARCachedRefMap(object):

    def __init__(self, read2ref=False):
        self.read2ref = read2ref
        self.cache = {}

    def __getitem__(self, key):
        if key not in self.cache:
            self.cache[key] = calculate_cigar_length(key) + (generate_cigar_map(key, self.read2ref),)
        return self.cache[key]

CIGARMAP_REF2READ = CIGARCachedRefMap(read2ref=False)
CIGARMAP_READ2REF = CIGARCachedRefMap(read2ref=True)

def get_read_refalign(read, cigar):
    if cigar[:-1].isdigit() and cigar[-1] == 'M':
        return read # no insertion or deletion. map as it is.

    reflen, readlen, refmap = CIGARMAP_REF2READ[cigar]
    r = []
    for i in range(reflen):
        readpos = refmap[i]
        if readpos is not None:
            r.append(read[readpos])
        else:
            r.append('D')
    return ''.join(r)

class QueueableLineReader(object):
    def __init__(self, infile):
        self.infile = infile
        self.queue = deque()

    def __iter__(self):
        try:
            while True:
                if self.queue:
                    yield self.queue.popleft()
                else:
                    yield next(self.infile)
        except StopIteration:
            return

    def push(self, line):
        self.queue.append(line)

    def pushleft(self, line):
        self.queue.appendleft(line)


class SAMParser(object):
    def __init__(self, samfile, seqfile=None, zerobase=False, use_slider=False):
        if isinstance(samfile, basestring):
            self.samfile = open(samfile)
        else:
            self.samfile = samfile

        if use_slider:
            from rnarry import batchutils
            self.samfile = batchutils.slider_file(self.samfile)

        self.samfile = QueueableLineReader(self.samfile)
        self.seqfile = seqfile
        if seqfile is not None:
            self.seqidx = GiantFASTAFile(seqfile)
        self.zerobase = zerobase
        self.seqlen = {}

    def getsubseq(self, name, seqfrom, seqto): # [from, to) both 0-base
        if self.seqfile is None:
            raise ValueError, 'Sequence file is not given.'

        return self.seqidx.get(name, seqfrom, seqto)

    def __iter__(self):
        return self.iteralignments()

    # WARNING: 'mapped' coordiates are 1-based, both side inclusive by default
    def iteralignments(self, strands='+-', withref=False):
        geteditdist= lambda x: x[4]

        for line in self.samfile:
            fields = line[:-1].split('\t')

            if line[0] == '@':
                if line[:3] != '@SQ':
                    continue
                sqname = fields[1][3:]
                sqlen = [int(fl[3:]) for fl in fields[2:] if fl[:3] == 'LN:'][0]
                self.seqlen[sqname] = sqlen
                continue

            qname = fields[0]
            flags = int(fields[1])
            rname = fields[2]
            pos = int(fields[3]) # 1-based leftmost
            mapq = int(fields[4]) # phred-scaled
            cigar = fields[5]
            seq = fields[9]
            options = dict(v.split(':', 1) for v in fields[11:])

            if flags & F_REVERSE_STRAND:
                strand = '-'
                seq = reverse_complement(seq)
            else:
                strand = '+'

            editdist = int(options.get('NM', 'i:-1')[2:])

            if rname == '*' or strand not in strands:
                mapped = []
            else:
                reflen, _ = calculate_cigar_length(cigar)
                stop = pos + reflen - 1
                start = pos - 1 if self.zerobase else pos
                mapped = [(rname, start, stop, strand, editdist, cigar)]

            for altmatch in options.get('XA', 'Z:')[2:].split(';')[:-1]:
                altfields = altmatch.split(',')
                strand = altfields[1][0]
                pos = int(altfields[1][1:])
                rname = altfields[0]
                cigar = altfields[2]
                editdist = int(altfields[3])
                reflen, _ = calculate_cigar_length(cigar)
                stop = pos + reflen - 1
                start = pos - 1 if self.zerobase else pos

                if strand in strands:
                    mapped.append((rname, start, stop, strand, editdist,
                                   cigar))

            # search for alternative reads
            for altline in self.samfile:
                altfields = altline[:-1].split('\t')
                altqname = altfields[0]
                altflags = int(altfields[1])
                if altqname != qname:
                    self.samfile.push(altline)
                    break

                altrname = altfields[2]
                altpos = int(altfields[3]) # 1-based leftmost
                altmapq = int(altfields[4]) # phred-scaled
                altcigar = altfields[5]
                altseq = altfields[9]
                altoptions = dict(v.split(':', 1) for v in altfields[11:])

                if altflags & F_REVERSE_STRAND:
                    altstrand = '-'
                    altseq = reverse_complement(altseq)
                else:
                    altstrand = '+'

                alteditdist = int(altoptions.get('NM', 'i:-1')[2:])

                if altrname != '*' and altstrand in strands:
                    altreflen, _ = calculate_cigar_length(altcigar)
                    altstop = altpos + altreflen - 1
                    altstart = altpos - 1 if self.zerobase else altpos
                    mapped.append((altrname, altstart, altstop, altstrand,
                                   alteditdist, altcigar))

            mapped.sort(key=geteditdist)

            if withref:
                newmapped = []
                for m in mapped:
                    subseq = self.getsubseq(m[0], m[1]-1, m[2])
                    if m[3] == '-':
                        subseq = reverse_complement(subseq)
                    newmapped.append(m + (subseq,))
                mapped = newmapped

            yield {
                'qname': qname,
                'flags': flags,
                'mapq': mapq,
                'seq': seq,
                'options': options,
                'mapped': mapped, # positions are 1-based
            }


def _test():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _test()

