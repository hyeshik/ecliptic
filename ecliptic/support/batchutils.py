# ----------------------------------------------------------------------------
# "THE BEER-WARE LICENSE" (Revision 42, (c) Poul-Henning Kamp):
# Hye-Shik Chang <perky@FreeBSD.org> wrote this file.  As long as you retain
# this notice you can do whatever you want with this stuff. If we meet some
# day, and you think this stuff is worth it, you can buy me a beer in return.
# ----------------------------------------------------------------------------

from __future__ import division
import os, sys, re
import time

__all__ = ['TerminalSlider', 'slider', 'slider_file']

def dotint(v):
    dotstr = str(v)
    left = dotstr[:len(dotstr)%3]
    right = re.findall('...', dotstr[len(dotstr)%3:])
    if left:
        right.insert(0, left)
    return ','.join(right)

class TerminalSlider:
    def __init__(self, maxpos, initpos=0):
        self.barfill = -1
        self.maxpos = maxpos
        self.initpos = initpos
        self.sessmax = maxpos - initpos
        try:
            import fcntl, termios, struct
            s = struct.pack("HHHH", 0, 0, 0, 0)
            self.cols = struct.unpack("HHHH", fcntl.ioctl(sys.stderr.fileno(),
                                termios.TIOCGWINSZ, s))[1]
        except:
            self.cols = 78

        # 4(percent) + 2([]) + 4(spaces) + 2(endmark)
        self.barsize = self.cols - 12
        self.barsize -= len(dotint(maxpos)) # done value
        self.barfmt = '\r%%-4s [%%s>%%s] %%-%ds ' % len(dotint(maxpos))
        self.barsize -= 23
        self.st = time.time()
        self.update(0)

    def update(self, value, *ext):
        if self.maxpos > 0:
            newfill = (value + self.initpos) * self.barsize // self.maxpos
            percent = str(value * 100 // self.maxpos)+'%'
        else:
            newfill = self.barsize
            percent = '0%'
        self.barfill = newfill
        sys.stderr.write(self.barfmt %
            (percent,
             '='*(newfill - 1), ' '*(self.barsize-newfill), dotint(value)))

        if value * 30 >= self.sessmax: # over 3%
            elapsed = time.time() - self.st
            if value > 0:
                estimated = int(self.sessmax / value * elapsed - elapsed)
            else:
                estimated = self.sessmax
            sys.stderr.write("   ETA %02d:%02d" %
                             (estimated // 60, estimated % 60))
        sys.stderr.flush()

    def end(self):
        print >> sys.stderr

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.end()

def slider(it, length=None):
    if length is None:
        try:
            length = len(it)
        except ValueError:
            it = list(it)
            length = len(it)

    s = TerminalSlider(length, 0)
    try:
        for i, el in enumerate(it):
            s.update(i)
            yield el
        else:
            s.update(length)
    finally:
        s.end()

def slider_file(f, updateevery=1000):
    f.seek(0, os.SEEK_END)
    filesize = f.tell()
    f.seek(0, os.SEEK_SET)

    s = TerminalSlider(filesize, 0)
    try:
        for i, it in enumerate(f):
            if i % updateevery == 0:
                s.update(f.tell())
            yield it
        else:
            s.update(filesize)
    finally:
        s.end()

if __name__ == '__main__':
    f = TerminalSlider(100, 0)
    for i in range(100):
        time.sleep(0.2)
        f.update(i)
