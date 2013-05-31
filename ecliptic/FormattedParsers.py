#
# ecliptic.FormattedParsers
#  - parsers for internal file formats in ecliptic pipeline
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

from .support.fileutils import LineParser

__all__ = [
    'BindingSiteCatalogParser',
]

# $JOBDIR/erroranalyses/*.catalog.gz
#
BindingSiteCatalogParser = LineParser([
    ('chrom', None),
    ('strand', None),
    ('pos', int),
    ('base', None),
    ('depth', int),
    ('del_score', float),
    ('mod_score', float),
    ('moddel_score', float),
    ('entropy_score', float),
    ('t2c_score', float),
    ('del_fdr', float),
    ('mod_fdr', float),
    ('moddel_fdr', float),
    ('entropy_fdr', float),
    ('t2c_fdr', float),
    ('genome_seq', None),
    ('transcript_seq', None),
    ('transcript_acc', None),
    ('transcript_pos', int),
])

