#
# ecliptic.ScriptGenerator
#  - Generates scripts from templates to process data
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

from jinja2 import Template, FileSystemLoader, Environment
import os
import re
from .Paths import Paths

__all__ = [
    'ScriptGenerator',
]

class ScriptGenerator:

    tagline = re.compile(r'ecliptic:([A-Za-z0-9 \t]*)')

    def __init__(self, templatedir):
        self.templatedir = templatedir
        self.templates = self.scan_templates()

        template_global_variables = {
            'Paths': Paths,
        }
        self.env = Environment(loader=FileSystemLoader(templatedir))
        self.env.globals.update(template_global_variables)

    def scan_templates(self):
        templates = []

        for fname in os.listdir(self.templatedir):
            if fname.startswith('.'):
                continue

            fullpath = os.path.join(self.templatedir, fname)
            tagspec = self.tagline.findall(open(fullpath).read(200)) + ['']
            tags = tagspec[0].split() if tagspec else []

            if 'work' not in tags:
                continue

            templates.append(fname)

        return templates

    def generate_project_scripts(self, project):
        if not os.path.isdir(project.workdir):
            os.makedirs(project.workdir)

        for name in self.templates:
            outputpath = os.path.join(project.workdir, name)
            tmpl = self.env.get_template(name)
            open(outputpath, 'w').write(tmpl.render())

