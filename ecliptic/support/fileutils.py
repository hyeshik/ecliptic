from __future__ import division
from itertools import groupby

__all__ = [
    'LineParser', 'ParsedLine', 'ParsedLineComment',
]

class LineParser(object):

    def __init__(self, fields_spec, separator='\t', linefeed='\n', comment=None):
        self.field2index = dict((field[0], i) for i, field in enumerate(fields_spec))
        self.fields_spec = fields_spec
        self.separator = separator
        self.linefeed = linefeed
        self.comment = comment

    def parse(self, line):
        if self.comment is not None and line.startswith(self.comment):
            return ParsedLineComment(line.rstrip(self.linefeed), line)

        fields = line.rstrip(self.linefeed).split(self.separator)
        return ParsedLine(self.field2index, [
            (data if adapter is None else (adapter(data) if data else None))
            for data, (name, adapter) in zip(fields, self.fields_spec)], line)

        return ParsedLine(self, line)

    def iter_parse(self, it):
        for line in it:
            yield self.parse(line)

    def __call__(self, it):
        return self.iter_parse(it)


class ParsedLineComment(object):

    def __init__(self, data, line):
        self.data = data
        self.line = line

    def __str__(self):
        return self.line


class ParsedLine(object):

    def __init__(self, field2index, data, line):
        self.field2index = field2index
        self.data = data
        self.line = line

    def __getitem__(self, index):
        return self.data[index]

    def __len__(self):
        return len(self.data)

    def __getattr__(self, name):
        if name.startswith('_'):
            return object.__getattr__(self, name)
        else:
            return self.data[self.field2index[name]]

    def __getitem__(self, name):
        return self.data.__getitem__(name)

    def __repr__(self):
        return '<ParsedLine %s>' % ' '.join(
            '%s=%s' % (name, repr(self.data[idx]))
            for name, idx in sorted(self.field2index.iteritems(), key=lambda x: x[1]))

    def __str__(self):
        return self.line


class ParallelMatchingReader(object):

    def __init__(self, src1, src2, src1key, src2key=None, separator='\t', linefeed='\n'):
        if src2key is None:
            src2key = src1key

        self.src1 = groupby(self.splititer(src1, separator, linefeed),
                            key=lambda x: x[src1key])
        self.src2 = groupby(self.splititer(src2, separator, linefeed),
                            key=lambda x: x[src2key])
        self.src1key = src1key
        self.src2key = src2key

    @staticmethod
    def splititer(src, sep, lf):
        for line in src:
            yield line.rstrip(lf).split(sep)

    def __iter__(self):
        src1ontap = src2ontap = None

        while True:
            if src1ontap is None:
                src1ontap = next(self.src1, None)
            if src2ontap is None:
                src2ontap = next(self.src2, None)

            if src1ontap is None:
                if src2ontap is None:
                    # both inputs are null
                    break

                yield (None, src2ontap[1])
                for _, v in self.src2:
                    yield (None, v)

                break
            elif src2ontap is None:
                yield (src1ontap[1], None)
                for _, v in self.src1:
                    yield (v, None)

                break

            # inputs from both sources are available
            if src1ontap[0] < src2ontap[0]:
                yield (src1ontap[1], None)
                src1ontap = None
            elif src1ontap[0] > src2ontap[0]:
                yield (None, src2ontap[1])
                src2ontap = None
            else:
                yield (src1ontap[1], src2ontap[1])
                src1ontap = src2ontap = None


if __name__ == '__main__':
    for a, b in ParallelMatchingReader(open('t2ctemp.bed'), open('new-t2c-genome.bed'), 3, 3):
        print list(a or ''), list(b or '')


