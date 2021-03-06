#
# ecliptic.TemplateFilters
#  - Support filters for script generators
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
import re

__all_filters__ = [
    'rule_name', 'thousand_sep', 'human_readable_bytes',
]

__all__ = [
    'register_filters',
] + __all_filters__


def rule_name(value, illegal=re.compile('[^A-Za-z0-9_]')):
    return illegal.sub('_', value)

def thousand_sep(value):
    return '{:,d}'.format(int(value))

def human_readable_bytes(value):
    num = value
    for x in ('bytes', 'KiB', 'MiB', 'GiB', 'TiB', 'PiB', 'EiB', 'ZiB', 'YiB'):
        if num < 1024.:
            return '{:.1f} {}'.format(num, x)
        num /= 1024.

    raise ValueError('Bytes too big: {}'.format(value))

def register_filters(env):
    for name in __all_filters__:
        env.filters[name] = eval(name)

