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
from datetime import datetime
from .Paths import Paths
from .DataSources import *
from .Utils import readlink_recursive

__all__ = [
    'ReportTemplateVariablesProvider',
]


class ReportTemplateVariablesProvider:

    def __init__(self, work_dir):
        self.work_dir = os.path.abspath(work_dir)
        project_name = os.path.basename(self.work_dir)
        self.project = ProjectWithReporting(project_name, self.work_dir)

    def toplevel(self):
        return {'project': self.project}


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

