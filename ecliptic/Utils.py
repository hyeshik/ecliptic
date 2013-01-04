#
# ecliptic.Utils
#  - Support functions
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

import csv

__all__ = [
    'TokensDictParser', 'iter_tab_separated',
]


class TokensDictParser:

    def __init__(self, fields, intfields=(), nullstrfields=()):
        self.fields = fields
        self.intfields = intfields
        self.nullstrfields = nullstrfields

    def __call__(self, tokens):
        record = dict(zip(self.fields, tokens))

        for fn in self.intfields:
            record[fn] = int(record[fn])

        for fn in self.nullstrfields:
            record[fn] = record[fn] if record[fn] != 'none' else None

        return record


def iter_tab_separated(file, parser):
    for line in csv.reader(file, dialect='excel-tab'):
        yield parser(line)

