#
# ecliptic.Paths
#  - Path configurations
#
#
# Copyright (C) 2013 Hyeshik Chang
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#

import os

__all__ = [
    'Paths', 'WORK_SUBDIRS', 'WORK_SYMLINKS',
]

WORK_SUBDIRS = ['sequences', 'alignments', 'annotations', 'bg_variations']
WORK_SYMLINKS = [('.', 'original'), ('PAIRS', 'PAIRS'), ('SAMPLES', 'SAMPLES')]

def pathgetter(name):
    def __getter(self, join=os.path.join):
        return join(self.basedir, name)
    return property(__getter)

# Singleton object for path configurations
class Paths:

    def __init__(self, basedir=None):
        # basedir may be reconfigured in later step of initialization
        self.basedir = basedir if basedir is not None else os.getcwd()

    @property
    def topdir(self):
        return self.basedir

    workdir         = pathgetter('work')
    datasourcedir   = pathgetter('originals')
    templatesdir    = pathgetter('templates')
    toolsdir        = pathgetter('tools')
    resourcesdir    = pathgetter('resources')


Paths = Paths()

WORK_SYMLINKS.extend([
    (Paths.resourcesdir, 'resources'),
])
