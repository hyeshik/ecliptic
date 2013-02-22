#
# ecliptic.ScriptRunner
#  - Runs sequentially scripts to run jobs
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
import sys
import time
import glob
try:
    from colorama import Fore, Style
except ImportError:
    Fore = Style = None
#from .Paths import Paths
#from .DataSources import *

__all__ = [
    'ScriptRunner',
]

class ScriptRunner:

    def __init__(self, datasources):
        self.datasources = datasources

    def scan_scripts(self):
        scripts = []

        for project in self.datasources.projects.values():
            scriptpattern = os.path.join(project.workdir, 'Snakefile*')
            for filename in glob.glob(scriptpattern):
                tokens = os.path.basename(filename).split('.', 2)
                if len(tokens) >= 2 and tokens[1].isdigit():
                    priority = int(tokens.pop(1))
                else:
                    priority = float('inf') # lowest priority for unprioritized tasks

                scripts.append((priority, project.name, tuple(tokens), project, filename))

        return sorted(scripts)

    def run(self, jobs=None):
        for prio, _, _, proj, snakefile in self.scan_scripts():
            if jobs is not None and proj.name not in jobs:
                continue

            if Fore is not None:
                tokens = snakefile.rsplit(os.path.sep, 2)
                tokens[-2] = (Style.BRIGHT + Fore.RED + tokens[-2] +
                              Style.RESET_ALL + Fore.RESET)
                snakefile_printing = os.path.sep.join(tokens)
            else:
                snakefile_printing = snakefile

            print("==> {}".format(snakefile_printing))

            ret = os.system('cd "{}" && snakemake -s "{}" {} -j'.format(
                            proj.workdir, snakefile, 'all'))
            if ret != 0:
                print('Error code {}. Terminating.'.format(ret), file=sys.stderr)
                sys.exit(1)

        print('\ndone.')

    # XXX: add runner for submitting jobs to doligo or openpbs

