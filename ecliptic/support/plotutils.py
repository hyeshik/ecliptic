# Copyright (c) 2011 Seoul National University RNA Biology Laboratory.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Written by: Hyeshik Chang
# Written on: March 6th, 2012

from __future__ import division

__all__ = [
    'adjust_numbers_style',
    'LinearToHeatRGB',
    'iwork_colors',
    'prepare_cumulative',
    'ntile_color',
]

import matplotlib.lines as mpllines
from matplotlib import pyplot as plt
import numpy as np
import colorsys


def adjust_numbers_style(ax, spines=('left', 'bottom'), xgrid=False, smart_bounds=False):
    for loc, spine in ax.spines.iteritems():
        if loc in spines:
            spine.set_smart_bounds(smart_bounds)
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
        ax.yaxis.set_tick_params(pad=9)
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_tick_params(pad=7)
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])
    for line in ax.get_xticklines():
        line.set_marker(mpllines.TICKDOWN)

    for line in ax.get_yticklines():
        line.set_marker(mpllines.TICKLEFT)

    ax.yaxis.grid(True, linestyle='-', alpha=0.3)
    if xgrid:
        ax.xaxis.grid(True, linestyle='-', alpha=0.3)
    else:
        ax.xaxis.grid(False)

    leg = ax.get_legend()
    if leg is not None:
        ltext  = leg.get_texts()
        llines = leg.get_lines()
        frame  = leg.get_frame()

        #frame.set_facecolor('0.80')
        plt.setp(ltext, fontsize='12')
        plt.setp(llines, linewidth=1.5)

class ColorSet(object):
    def __init__(self, colors):
        self.bykey = dict(colors)
        self.byiter = [color for k, color in colors]

    def __getitem__(self, key):
        return self.byiter[key % len(self.byiter)]

    def __iter__(self):
        return iter(self.byiter)

    def __getattr__(self, name):
        if name in self.bykey:
            return self.bykey[name]
        else:
            return object.__getattr__(self, name)

iwork_colors = ColorSet([
    ('blue', '#2e578c'),
    ('green', '#5d9648'),
    ('yellow', '#e7a13d'),
    ('red', '#bc2d30'),
    ('violet', '#6f3d79'),
    ('mint', '#01ae94'),
    ('gray', '#7d807f'),
    ('orange', '#ff8000'),
])

class LinearToHeatRGB(object):
    def __init__(self, vmin, vmax):
        self.vmin = vmin
        self.vmax = vmax
        self.vinterval = vmax - vmin

    def __call__(self, v):
        raise NotImplementedError


def prepare_cumulative(grp, width=100, reverse=False):
    grp = sorted(grp)
    xvalues = [grp[int(i / width * (len(grp)-1))] for i in range(width+1)]
    yvalues = np.arange(width+1) / width
    if reverse:
        yvalues = yvalues[::-1]
    return xvalues, yvalues


def ntile_color(n):
    colors = []
    for i in range(n):
        r, g, b = colorsys.hsv_to_rgb(i / (n - 1) * 0.95, 0.92, 0.75)
        colors.append('#%02x%02x%02x' % (int(r * 255), int(g * 255),
                                         int(b * 255)))
    return colors


if __name__ == '__main__':
    f = LinearToHeatRGB(0, 1)
    print map(f, [0, 0.3, 0.5, 0.8, 1])

