#
# Copyright (c) 2012 Hyeshik Chang
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

__all__ = [
    'TranscriptLevelReadsStats',
]

from rnarry.utils import LRU, lru_cached
import lzo as comp
import numpy as np
import bsddb

class TranscriptLevelReadsStats(object):

    NROW = 9
    ROWSPACE = ['depth', '5', '3', 'A', 'C', 'G', 'T', 'I', 'D', 'N']
    ROWSPACE2i = dict((tp, i) for i, tp in enumerate(ROWSPACE))

    def __init__(self, filepath, cachedepth=10):
        self.cache = LRU(cachedepth)
        self.db = bsddb.hashopen(filepath, 'r')
        self._get = lru_cached(self._get, self.cache)

    def _get(self, acc):
        flat = np.fromstring(comp.decompress(self.db[acc]), 'I')
        return flat.reshape([self.NROW, len(flat) // self.NROW])

    def get(self, acc, which=None):
        if which is None:
            return self._get(acc)
        else:
            return self._get(acc)[self.ROWSPACE2i[which]]

    def keys(self):
        return self.db.keys()

