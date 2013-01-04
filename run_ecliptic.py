#!/usr/bin/env python3
from ecliptic.DataSources import *
from ecliptic.ScriptGenerator import *

ds = DataSources('originals')
scriptgen = ScriptGenerator('templates')
scriptgen.generate_project_scripts(ds.projects['GSE37685-HNRNPU'])

