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

from __future__ import division
import csv
import gzip
from collections import defaultdict

def process(inpfile, outfile, additionals):
    counter = defaultdict(int)
    for title, reads in additionals:
        counter[title] += int(reads)

    for line in csv.reader(inpfile, dialect='excel-tab'):
        nreads = int(line[0].split('-')[1])
        nanno = int(line[1])

        annoclass = line[2]
        if not annoclass:
            if nanno == 0:
                annoclass = 'Unmapped'
            else:
                annoclass = 'Unannotated'

        counter[annoclass] += nreads

    orderedcount = sorted(counter.items(), key=lambda x: x[1], reverse=True)

    totalcnt = sum(counter.itervalues())

    out = csv.writer(outfile)
    out.writerow(('class', 'reads', 'proportion'))
    for clsname, cnt in orderedcount:
        out.writerow((clsname, cnt, cnt / totalcnt * 100))

if __name__ == '__main__':
    import sys

    additionals = [token.split(':') for token in sys.argv[2:]]
    process(gzip.GzipFile(sys.argv[1]), sys.stdout, additionals)
