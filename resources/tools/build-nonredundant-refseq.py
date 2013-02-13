#!/usr/bin/env python
import sys

LABELS = [
    'geneName', 'name', 'chrom', 'strand', 'txStart', 'txEnd',
    'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds'
]
INTFIELDS = ['txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount']

def parse_refflat_line(inpf):
    for line in inpf:
        d = dict(zip(LABELS, line[:-1].split('\t')))

        for f in INTFIELDS:
            d[f] = int(d[f])

        d['exonBlocks'] = zip(map(int, d['exonStarts'].split(',')[:-1]),
                              map(int, d['exonEnds'].split(',')[:-1]))
        del d['exonStarts']
        del d['exonEnds']

        yield d

def get_exon_len(entry):
    return sum(blkend-blkstart for blkstart, blkend in entry['exonBlocks'])

def does_overlap(A, B, C, D):
    return (A <= D <= B) or (C <= B <= D)

def eliminate_redunant_isoforms(cluster):
    survived = set(range(len(cluster)))

    for i, entry in enumerate(cluster):
        if i not in survived:
            continue

        for j in list(survived):
            if i == j:
                continue

            other = cluster[j]
            if any(does_overlap(astart, aend-1, bstart, bend-1)
                   for astart, aend in entry['exonBlocks']
                   for bstart, bend in other['exonBlocks']):
                survived.remove(j)

    return [cluster[i] for i in sorted(survived)]

def process(inpf, db, outtxt, anno):
    cluster = []
    clusterstart = -1
    clusterend = -1
    clusterchrom = None
    clusterstrand = None

    def flush():
        if clusterchrom != row['chrom']:
            print >> sys.stderr, 'Processing', row['chrom']

        if not cluster:
            return

        if len(cluster) > 1:
            cluster.sort(key=lambda e: (e['name'][:2], -get_exon_len(e),
                                        int(e['name'][3:]))) # NM_ is preferred
            for entry in eliminate_redunant_isoforms(cluster):
                entry['description'] = anno.get(entry['name'], '')
                add_mRNA_block_annotation(entry)
                db[entry['name']] = entry
        else:
            entry = cluster[0]
            entry['description'] = anno.get(entry['name'], '')
            add_mRNA_block_annotation(entry)
            db[entry['name']] = entry

        del cluster[:]

    for row in parse_refflat_line(inpf):
        if (clusterchrom != row['chrom']
                or clusterstrand != row['strand']
                or clusterend <= row['txStart']):
            flush()

            clusterstart = row['txStart']
            clusterend = row['txEnd']
            clusterchrom = row['chrom']
            clusterstrand = row['strand']
        elif clusterend < row['txEnd']:
            clusterend = row['txEnd']

        cluster.append(row)

    flush()

    print >> outtxt, '\n'.join(sorted(db.keys()))

def split_blocks(entry):
    cdsstart, cdsend = entry['cdsStart'], entry['cdsEnd']
    leftutrblocks, rightutrblocks, cdsblocks = [], [], []

    for exstart, exend in entry['exonBlocks']:
        if exend <= cdsstart:
            leftutrblocks.append((exstart, exend))
        elif exstart >= cdsend:
            rightutrblocks.append((exstart, exend))
        else:
            eff_cdsstart = min(exend, max(exstart, cdsstart))
            eff_cdsend = max(exstart, min(exend, cdsend))

            if eff_cdsstart != exstart:
                assert eff_cdsstart > exstart
                leftutrblocks.append((exstart, eff_cdsstart))
            if eff_cdsstart != eff_cdsend:
                assert eff_cdsstart < eff_cdsend
                cdsblocks.append((eff_cdsstart, eff_cdsend))
            if eff_cdsend != exend:
                assert eff_cdsend < exend
                rightutrblocks.append((eff_cdsend, exend))

    entry['leftUtrBlocks'] = leftutrblocks
    entry['rightUtrBlocks'] = rightutrblocks
    entry['cdsBlocks'] = cdsblocks

def block_length(blk):
    return sum(en - st for st, en in blk)

def add_mRNA_block_annotation(entry):
    name = entry['name']

    if name.startswith('NM_'):
        split_blocks(entry)
        lenleft = block_length(entry['leftUtrBlocks'])
        lenright = block_length(entry['rightUtrBlocks'])
        lencds = block_length(entry['cdsBlocks'])

        if entry['strand'] == '+':
            entry['partLengths'] = (lenleft, lencds, lenright)
        else:
            entry['partLengths'] = (lenright, lencds, lenleft)

    entry['totalLength'] = block_length(entry['exonBlocks'])

    entry['intronBlocks'] = [
        (leftend, rightstart) for (_, leftend), (rightstart, _)
        in zip(entry['exonBlocks'][:-1], entry['exonBlocks'][1:])]
    entry['intronLength'] = block_length(entry['intronBlocks'])


if __name__ == '__main__':
    import os, shelve
    import sys

    reflinkfile = sys.argv[1]
    output = sys.argv[2]
    listoutput = sys.argv[3]

    outdb = shelve.open(output, 'c')
    outtxt = open(listoutput, 'w')

    import gzip, csv
    anno = {}
    for rl in csv.reader(gzip.open(reflinkfile), dialect='excel-tab'):
        if rl[1]:
            anno[rl[2]] = rl[1]

    process(sys.stdin, outdb, outtxt, anno)

