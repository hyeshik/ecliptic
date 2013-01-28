#!/usr/bin/env python
import sys, os
import shelve
import gzip
import pickle
from numpy import array, cumsum
from collections import defaultdict

REFGENE_HEADER = [
    'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd',
    'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds',
    'id', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames',
]

def parse_refgene(annofile):
    for line in annofile:
        if line.startswith('#'):
            continue

        fields = dict(zip(REFGENE_HEADER, line[:-1].split('\t')))

        yield {
            'txStart': int(fields['txStart']),
            'txEnd': int(fields['txEnd']),
            'exonStarts': map(int, fields['exonStarts'][:-1].split(',')),
            'exonEnds': map(int, fields['exonEnds'][:-1].split(',')),
            'cdsStart': int(fields['cdsStart']),
            'cdsEnd': int(fields['cdsEnd']),
            'name': fields['name'],
            'chrom': fields['chrom'],
            'strand': fields['strand'],
            'geneName': fields['name2'],
        }

REFSEQ_CLSMAP = {'5': 'UTR5', '3': 'UTR3', 'i': 'intron', 'c': 'CDS'}
def load_refgene(bedout, refgene):
    def add(e, start, stop, ei, regtype):
        instancename = '%s.%d' % (e['name'], serial)
        annoid = '%s#%d%s' % (instancename, ei, regtype)

        if instancename.startswith('NR_') or instancename.startswith('XR_'):
            cls = 'ncRNA'
        else:
            cls = REFSEQ_CLSMAP[regtype]

        commonname = e['geneName'] or e['name']
        if not commonname.upper().startswith('MIR'):
            print >> bedout, '\t'.join((
                e['chrom'], str(start), str(stop),
                '%s|%s|%s' % (cls, commonname, annoid),
                '.', e['strand']
            ))

        lengthstat[instancename][regtype] += stop - start
        intervals[instancename][regtype].append((start, stop))
        intervals[instancename]['strand'] = e['strand']
        intervals[instancename]['chrom'] = e['chrom']

    lengthstat = defaultdict(lambda: {'i': 0, '5': 0, '3': 0, 'c': 0})
    intervals = defaultdict(lambda: {'i': [], '5': [], '3': [], 'c': []})
    idserial = defaultdict(int)
    instances = defaultdict(list)

    for i, e in enumerate(parse_refgene(refgene)):
        if i and i % 5000 == 0:
            print '%d done.' % i

        if e['strand'] == '+':
            firstutr, lastutr = '53'
        else:
            firstutr, lastutr = '35'

        idserial[e['name']] += 1
        serial = idserial[e['name']]
        instances[e['name']].append('%s.%d' % (e['name'], serial))

        for ei, (estart, estop) in enumerate(zip(e['exonStarts'], e['exonEnds'])):
            if e['strand'] == '-':
                ei = len(e['exonStarts']) - 1 - ei

            have_cdsstart = (estart <= e['cdsStart'] < estop)
            have_cdsend = (estart <= e['cdsEnd'] < estop)

            if have_cdsstart and have_cdsend:
                if estart != e['cdsStart']:
                    add(e, estart, e['cdsStart'], ei, firstutr)
                if e['cdsStart'] != e['cdsEnd']:
                    add(e, e['cdsStart'], e['cdsEnd'], ei, 'c')
                if estop != e['cdsEnd']:
                    add(e, e['cdsEnd'], estop, ei, lastutr)
            elif have_cdsstart:
                if estart != e['cdsStart']:
                    add(e, estart, e['cdsStart'], ei, firstutr)
                add(e, e['cdsStart'], estop, ei, 'c')
            elif have_cdsend:
                add(e, estart, e['cdsEnd'], ei, 'c')
                if estop != e['cdsEnd']:
                    add(e, e['cdsEnd'], estop, ei, lastutr)
            else:
                if estop <= e['cdsStart']:
                    regtype = firstutr
                elif estart >= e['cdsEnd']:
                    regtype = lastutr
                else:
                    regtype = 'c'

                add(e, estart, estop, ei, regtype)

        for ei, (estop, estart) in enumerate(zip(e['exonEnds'][:-1],
                                                 e['exonStarts'][1:])):
            if e['strand'] == '-':
                ei = len(e['exonStarts']) - 2 - ei
            add(e, estop, estart, ei, 'i')

    return lengthstat, intervals, instances


def main(args):
    if len(args) < 2:
        print 'Usage: %s dbprefix refgenefile' % args[0]
        sys.exit(1)

    dbprefix = args[1]
    refgenefile = args[2]

    refgene = gzip.open(refgenefile) # XXX
    bedout = gzip.open(os.path.join(dbprefix, 'cat.refseq.bed.gz'), 'w')
    lengthstat, intervals, instances = load_refgene(bedout, refgene)

    pickle.dump(dict(lengthstat),
                open(os.path.join(dbprefix, 'refseq-lengths.pickle'), 'w'))
    pickle.dump(dict(intervals),
                open(os.path.join(dbprefix, 'refseq-intervals.pickle'), 'w'))
    pickle.dump(dict(instances),
                open(os.path.join(dbprefix, 'refseq-instances.pickle'), 'w'))

if __name__ == '__main__':
    main(sys.argv)

