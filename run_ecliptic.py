#!/usr/bin/env python3
from ecliptic.Paths import *
from ecliptic.DataSources import *
from ecliptic.ScriptGenerator import *
from ecliptic.ScriptRunner import *

ds = DataSources(Paths.datasourcedir)
ScriptGenerator(Paths.templatesdir).generate_project_scripts(ds)
ScriptRunner(ds).run()

