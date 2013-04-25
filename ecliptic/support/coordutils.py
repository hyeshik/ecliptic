#!/usr/bin/env python
import bsddb
import os
import csv
import tempfile
from itertools import groupby

__all__ = ['LiftOverToTranscriptome']

class LiftOverToTranscriptome(object):
    """
    Alternative implementation of UCSC GenomeBrowser's liftOver.
    This implementation is better (less buggy) at handling many short alignments
    on "negative" strand.
    """

    def __init__(self, database, reference_bed):
        self.db = bsddb.hashopen(database, 'r')
        self.reference_bed = reference_bed

    def liftover_points(self, points):
        with tempfile.NamedTemporaryFile(suffix='bed') as bedtmp:
            for i, point in enumerate(points): # point: (chrom, pos, strand, ....)
                print >> bedtmp, '{chrom}\t{start}\t{end}\t{name}\t.\t{strand}'.format(
                    chrom=point[0], start=point[1], end=point[1]+1, name=i, strand=point[2])

            bedtmp.flush()
            result = os.popen('bedtools intersect -loj -s -wa -wb -a {tmp} -b {reference}'.format(
                        tmp=bedtmp.name, reference=self.reference_bed))
            for i, grp in groupby(csv.reader(result, dialect='excel-tab'), lambda x: x[3]):
                allmatches = list(grp)
                assert len(allmatches) == 1
                match = allmatches[0]
                if match[9] == '.': # no matched transcript
                    yield None
                else:
                    pos = int(match[1])
                    strand = match[5]
                    transcript = match[9].split('/')[0]
                    trpos = self.convert_pos(pos, transcript)
                    yield None if trpos is None else (transcript, trpos)

    def convert_pos(self, pos, transcript):
        try:
            pslentry = self.db[transcript].split('\t')
        except KeyError:
            return None

        qSize = int(pslentry[10])
        blockSizes = map(int, pslentry[18][:-1].split(','))
        qStarts = map(int, pslentry[19][:-1].split(','))
        tStarts = map(int, pslentry[20][:-1].split(','))
        strand = pslentry[8]

        for qst, tst, blksize in zip(qStarts, tStarts, blockSizes):
            relpos = pos - tst
            if relpos < 0: # indel
                break
            elif relpos < blksize:
                fwdpos = qst + relpos
                return fwdpos if strand == '+' else qSize - 1 - fwdpos

#0       1    2      3    4    5    6    7       8    9               10      11   12      13      14           15         16         17   18                       19                        20
['3286', '0', '915', '0', '0', '0', '4', '8326', '-', 'NM_001177658', '4216', '0', '4201', 'chr1', '195471971', '4773199', '4785726', '5', '3602,124,166,155,154,', '15,3617,3741,3907,4062,', '4773199,4777524,4782567,4783950,4785572,']


if __name__ == '__main__':
    lotr = LiftOverToTranscriptome('../../resources/mm10/nrRefSeq-alignment.db',
                                   '../../resources/mm10/nrRefSeq-genome.bed.gz')
    print list(lotr.liftover_points([('chr1', 4775204, '-'),
                                    ('chr4', 134008024, '-'),
                                    ('chrX', 1512312, '+')
                                    ]))

