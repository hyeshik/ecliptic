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
import time
from itertools import chain
from .Paths import Paths, WORK_SUBDIRS, WORK_SYMLINKS
from .DataSources import *
from .TemplateFilters import register_filters

__all__ = [
    'ScriptGenerator',
]

class ScriptGenerator:

    tagline = re.compile(r'ecliptic:([A-Za-z0-9 \t]*)')

    def __init__(self, chainsdir):
        self.templatedir = chainsdir
        self.templates = {
            'work': self.scan_templates('work'),
            'setting': self.scan_templates('setting'),
        }

        template_global_variables = {
            'Paths': Paths,
        }
        # change variable block markers to enable escaped braces "{{ }}" in str.format strings
        self.env = Environment(loader=FileSystemLoader(self.templatedir),
                               variable_start_string='{(',
                               variable_end_string=')}')
        self.env.globals.update(template_global_variables)
        register_filters(self.env)

    def scan_templates(self, typetag):
        templates = []

        for fname in os.listdir(self.templatedir):
            if fname.startswith('.'):
                continue

            fullpath = os.path.join(self.templatedir, fname)
            maximum_tag_head = 200
            tagspec = self.tagline.findall(open(fullpath).read(maximum_tag_head)) + ['']
            tags = tagspec[0].split() if tagspec else []

            if typetag not in tags:
                continue

            templates.append(fname)

        return templates

    def generate_project_scripts(self, project):
        if isinstance(project, DataSources):
            ds = project
            for project in ds.projects.values():
                self.generate_project_scripts(project)
            return

        self.prepare_workdirs(project)

        variables = self.prepare_template_variables(project)

        for name in self.templates['work']:
            outputpath = os.path.join(project.workdir, name)
            tmpl = self.env.get_template(name)
            outputfile = open(outputpath, 'w')
            print('# DO NOT EDIT THIS FILE MANUALLY. Generated by ecliptic on {}'.format(
                    time.asctime()), file=outputfile)
            outputfile.write(tmpl.render(**variables))

        for name in self.templates['setting']:
            outputpath = os.path.join(project.workdir, name)
            # Skip generation when the setting file already exists.
            if os.path.exists(outputpath):
                continue

            tmpl = self.env.get_template(name)

            outputfile = open(outputpath, 'w')
            print('# Uncomment and change the following options here as you need.\n',
                  file=outputfile)

            # Comment out non-comment non-blank lines.
            output = '\n'.join(('#'+line if line and line[:1] != '#' else line)
                               for line in tmpl.render(**variables).splitlines())

            outputfile.write(output + '\n')

    def prepare_workdirs(self, project):
        for subdir in WORK_SUBDIRS:
            subdir_fullpath = os.path.join(project.workdir, subdir)
            if not os.path.isdir(subdir_fullpath):
                os.makedirs(subdir_fullpath)

        for orig, name in WORK_SYMLINKS:
            linktarget = os.path.join(project.workdir, name)
            if callable(orig):
                orig = orig()
            linksource = os.path.abspath(os.path.join(project.datadir, orig))

            try:
                os.unlink(linktarget)
            except OSError:
                pass
            os.symlink(linksource, linktarget)

    unallowed_letters_in_name = re.compile('[^A-Za-z0-9_]')
    unallowed_letter_replacement = '_'
    def prepare_template_variables(self, project):
        variables = {}

        # add shortcuts to available tools
        toolsdefs = []
        for fname in os.listdir(Paths.toolsdir):
            if fname.startswith('.'):
                continue

            fullname = os.path.join(Paths.toolsdir, fname)
            if os.path.isfile(fullname) and os.access(fullname, os.X_OK):
                namebase = os.path.basename(fullname).rsplit('.', 1)[0]
                safenamebase = self.unallowed_letters_in_name.sub(
                        self.unallowed_letter_replacement, namebase).upper() + '_CMD'
                toolsdefs.append('{} = _cmd({})'.format(safenamebase, repr(fullname)))

        variables['INTERNAL_TOOLS'] = '\n'.join(toolsdefs)

        variables['SAMPLES'] = project.samples
        variables['PAIRS'] = project.pairs
        variables['WORKFLOWS'] = set(chain(*[sample.workflows
                                     for sample in project.samples.values()]))
        variables['NAME'] = project.name
        variables['REFERENCES'] = project.references
        variables['DATADIR'] = project.datadir
        variables['WORKDIR'] = project.workdir
        variables['THREADS'] = 32 # XXX: change this to detect machine configuration

        from ecliptic import VERSION_STRING
        variables['ECLIPTIC_VERSION'] = VERSION_STRING

        return variables

