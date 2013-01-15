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
import yaml
from weakref import proxy
from .Utils import TokensDictParser, iter_tab_separated
from .Paths import Paths
from .Settings import settings

__all__ = [
    'DataSources', 'Project',
]

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

        self.samples = Sample.scan_samples(self, datadir)
        self.pairs = None #self.scan_pairs()

    def __repr__(self):
        return '<Project at=%s nsamples=%d>' % (repr(self.datadir), len(self.samples))

    @property
    def workdir(self):
        return os.path.join(Paths.workdir, self.name)


class Sample:

    attribute_names = [
        'label', 'SRAno', 'GEOno', 'species', 'routes', 'first_base', 'last_base',
        'minimum_length', 'quality_scale', 'strand', 'threep_adapter', 'description'
    ]

    def __init__(self, project, **attrs):
        self.project = proxy(project)

        for attrname, defaultvalue in settings['sample'].items():
            setattr(self, attrname, defaultvalue)

        for attrname in self.attribute_names:
            if attrname in attrs:
                setattr(self, attrname, attrs[attrname])

        # TODO: check if any mandatory attribute is missing here.

    @classmethod
    def scan_samples(cls, project, datadir):
        samples_list_fn = os.path.join(datadir, 'SAMPLES')

        samples = {}
        for label, record in yaml.load(open(samples_list_fn)).items():
            samples[label] = cls(project, label=label, **record)

        return samples


if __name__ == '__main__':
    ds = DataSources('originals')
    print(ds)

