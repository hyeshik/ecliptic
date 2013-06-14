#!/usr/bin/env python3
from __future__ import division
import shutil
from jinja2 import Template, FileSystemLoader, Environment
from ecliptic.ReportGenerator import ReportTemplateVariablesProvider
from ecliptic.Paths import Paths
from ecliptic.TemplateFilters import register_filters
import os


REPORT_HTML_PATH = 'report.html'
REPORT_CSS_PATH = 'report/style.css'


def process(options):
    templatepath = lambda filename: os.path.join(Paths.reportingdir, filename)
    outputpath = lambda filename: os.path.join(options.work_dir, filename)

    variables = ReportTemplateVariablesProvider(options.work_dir)

    # Copy resource files
    shutil.copy(templatepath('style.css'), outputpath(REPORT_CSS_PATH))

    # Generate HTML file from template
    tmplenv = Environment(loader=FileSystemLoader(Paths.reportingdir))
    register_filters(tmplenv)
    tmpl = tmplenv.get_template('report.html')
    open(outputpath(REPORT_HTML_PATH), 'w').write(tmpl.render(**variables.toplevel()))


def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Generates a HTML file that compiles all '
                                                 'plots and tables for reporting')
    parser.add_argument('--ecliptic-dir', metavar='DIR', type=str, required=True,
                        dest='ecliptic_dir', help='Path to ecliptic data and templates.')
    parser.add_argument('--work-dir', metavar='DIR', type=str, required=True,
                        dest='work_dir', help='Directory containing project-specific files')
    parser.add_argument('--meta-info', metavar='FILE', type=str, required=True,
                        dest='meta_info', help='Pickle file that is passed from pipeline')

    return parser.parse_args()


if __name__ == '__main__':
    options = parse_arguments()
    Paths.basedir = os.path.abspath(options.ecliptic_dir)
    process(options)

