#
# ecliptic.DataSources
#  - Meta information for experimental data sources
#
#
# Copyright (C) 2012 Hyeshik Chang
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
import csv
from .Utils import TokensDictParser, iter_tab_separated
from .Paths import Paths

__all__ = [
    'DataSources', 'Project',
]


parse_samples_line = TokensDictParser(
    fields=['SRRno', 'GSMno', 'species', 'type', 'role', 'first_base', 'last_base',
            'quality_scale', 'strand', '3p_adapter', 'description'],
    intfields=['first_base', 'last_base'],
    nullstrfields=['3p_adapter']
)


class DataSources:

    def __init__(self, datadir):
        self.datadir = datadir
        self.projects = {}

        self.scan()

    def scan(self):
        for dirname in os.listdir(self.datadir):
            project_name = dirname
            self.projects[dirname] = Project(project_name, os.path.join(self.datadir, dirname))

    def __repr__(self):
        return '<DataSources at=%s nprojects=%d>' % (repr(self.datadir), len(self.projects))


class Project:

    def __init__(self, name, datadir):
        self.name = name
        self.datadir = datadir
        self.samples = {}
        self.pairs = {}

        self.scan()

    def scan(self):
        self.scan_samples()
        self.scan_pairs()

    def scan_samples(self):
        samples_list_fn = os.path.join(self.datadir, 'SAMPLES')

        for record in iter_tab_separated(open(samples_list_fn), parse_samples_line):
            self.samples[record['SRRno']] = record

    def scan_pairs(self):
        pass

    def __repr__(self):
        return '<Project at=%s nsamples=%d>' % (repr(self.datadir), len(self.samples))

    @property
    def workdir(self):
        return os.path.join(Paths.workdir, self.name)


if __name__ == '__main__':
    ds = DataSources('originals')
    print(ds)

