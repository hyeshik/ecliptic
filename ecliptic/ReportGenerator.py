#
# ecliptic.ReportGenerator
#  - Generates analysis report from templates
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
import pwd
import csv
from datetime import datetime
from collections import OrderedDict
from .Paths import Paths
from .DataSources import *
from .Utils import readlink_recursive

__all__ = [
    'ReportTemplateVariablesProvider',
    'ProjectWithReporting',
    'AnalysisResultReporting',
]


class ReportTemplateVariablesProvider:

    def __init__(self, work_dir):
        self.work_dir = os.path.abspath(work_dir)
        project_name = os.path.basename(self.work_dir)
        self.project = ProjectWithReporting(project_name, self.work_dir)

    def toplevel(self):
        return {'project': self.project, 'table': self.table,
                'join_tables': self.join_tables, 'samples': self.samples,
                'filesize': self.filesize}

    def table(self, name):
        filename = os.path.join(self.project.datadir, 'report', 'tables', name)
        return list(csv.DictReader(open(filename)))

    def join_tables(self, join_name, *names):
        # Join multiple tables `names' by key `join_name' by presented orders in
        # the earliest existing table of each key.

        filenames = [os.path.join(self.project.datadir, 'report', 'tables', name)
                     for name in names]
        input_contents = [
            (i, OrderedDict((row[join_name], row) for row in csv.DictReader(open(name))))
            for i, name in enumerate(filenames)]

        while input_contents:
            pivot_i, pivot = input_contents.pop(0)
            for key, pivotv in pivot.items():
                merged_value = dict(((pivot_i, k), v) for k, v in pivotv.items())
                for o_i, other in input_contents:
                    if key in other:
                        merged_value.update(dict(
                            ((o_i, k), v) for k, v in other[key].items()))
                        del other[key]

                merged_value['key'] = key
                yield merged_value

    def samples(self, workflow=None):
        for sample in self.project.samples.values():
            if workflow is None or (
                    workflow in sample.workflows if isinstance(workflow, str)
                    else any(w in sample.workflows for w in workflow)):
                yield sample

    def filesize(self, filename):
        if os.path.isfile(filename):
            return os.path.getsize(filename)
        elif os.path.isdir(filename):
            return sum(
                sum(os.path.getsize(os.path.join(dirname, f)) for f in files)
                for dirname, _, files in os.walk(filename))
        else:
            return 0


class ProjectWithReporting(Project):

    @property
    def sample_description_time(self):
        sampledeffn = readlink_recursive(os.path.join(self.workdir, 'SAMPLES'))
        return datetime.fromtimestamp(os.path.getmtime(sampledeffn))

    @property
    def sample_description_owner(self):
        sampledeffn = readlink_recursive(os.path.join(self.workdir, 'SAMPLES'))
        return pwd.getpwuid(os.stat(sampledeffn).st_uid).pw_name

    @property
    def pair_description_time(self):
        pairdeffn = readlink_recursive(os.path.join(self.workdir, 'PAIRS'))
        return datetime.fromtimestamp(os.path.getmtime(pairdeffn))

    @property
    def pair_description_owner(self):
        pairdeffn = readlink_recursive(os.path.join(self.workdir, 'PAIRS'))
        return pwd.getpwuid(os.stat(pairdeffn).st_uid).pw_name

    @property
    def analysis_begin_time(self):
        tm_begin = float(open(os.path.join(self.workdir, '.time.begin')).read())
        return datetime.fromtimestamp(tm_begin)

    @property
    def analysis_finish_time(self):
        tm_finish = float(open(os.path.join(self.workdir, '.time.finish')).read())
        return datetime.fromtimestamp(tm_finish)

