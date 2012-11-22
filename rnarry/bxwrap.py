#!/usr/bin/env python
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

__all__ = ['SplitBinnedArray', 'MultiTrackSplitBinnedArray']

#from bx.binned_array import FileBinnedArray
from rnarry.bx_binnedarray_fixed import FileBinnedArray # bugfixed
import os
import numpy as np


class DummyArray(object):

    def get_range(self, start, end):
        return [0] * (end - start)


class SplitBinnedArray(object):

    def __init__(self, prefix):
        self.prefix = prefix
        self.opened = {}

    def get_array(self, chrom, strand):
        key = (chrom, strand)
        if key not in self.opened:
            fname = os.path.join(self.prefix, '%s.%s' % (chrom, strand))
            try:
                self.opened[key] = FileBinnedArray(open(fname))
            except IOError:
                self.opened[key] = DummyArray()
        return self.opened[key]

    def get(self, chrom, start, end, strand): # zero-based half-open
        localarray = self.get_array(chrom, strand)
        return localarray.get_range(start, end)

    def get_blocks(self, chrom, blocks, strand):
        pieces = []
        for start, end in blocks:
            pieces.append(self.get(chrom, start, end, strand))
        return np.hstack(pieces)

    def clear_cache(self):
        self.opened.clear()


class MultiTrackSplitBinnedArray(object):

    TRACKS = ('depth', '5', '3', 'A', 'C', 'G', 'T', 'I', 'D')

    def __init__(self, prefix, tracks=TRACKS):
        self.prefix = prefix
        self.tracks = tracks
        self.subarrays = [
            SplitBinnedArray(os.path.join(prefix, view))
            for view in tracks]

    def get(self, chrom, start, end, strand): # zero-based half-open
        return np.vstack([arr.get(chrom, start, end, strand)
                          for arr in self.subarrays])

    def get_blocks(self, chrom, blocks, strand):
        pieces = []
        for start, end in blocks:
            pieces.append(self.get(chrom, start, end, strand))
        return np.hstack(pieces)

    def clear_cache(self):
        for subarray in self.subarrays:
            subarray.clear_cache()


if __name__ == '__main__':
    chrom, start, stop = 'chr10', 127526013, 127526212
    chrom, start, stop = 'chr10', 127526103, 127526122

    a = MultiTrackSplitBinnedArray('PolyA-None.gspace')
    #print a.get('chr10', 127527070, 127527088, '+')
    nearpoint = a.get(chrom, start, stop, '+')
    print nearpoint

    print a.get_blocks('chr10', [(127526103, 127526112),
                        (127526112, 127526122)], '+')

#    import csv
#    w = csv.writer(open('polya-none.csv', 'w'))
#    w.writerow(['name'] + range(start, stop))
#    for view, data in zip(a.tracks, nearpoint):
#        w.writerow([view] + list(data))

