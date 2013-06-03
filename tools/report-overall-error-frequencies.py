#!/usr/bin/env python
from __future__ import division
import csv
import sys
import futures
import os
from ecliptic.support import sam

def scan_alignment(bamfile):
    inpsam = os.popen('samtools view "{}"'.format(bamfile))

    totaltags = totalbases = 0
    cnt_subst = cnt_del = cnt_ins = 0

    for record in sam.parse_sam_simple(inpsam):
        # <ParsedLine qname='6-3799' flag=0 rname='chr9' pos=120960330 mapq=40 cigar='36M' rnext='*' pnext='0' tlen=0 seq='AAAAACTTTATCGGGGATACATGCGGTAGGGTAAAT' qual='*'>
        nreads = int(record.qname.split('-')[1])
        totaltags += nreads

        reflen, readlen = sam.calculate_cigar_length(record.cigar)
        totalbases += nreads * reflen

        nmtag = [int(tok.split(':')[2]) for tok in str(record).split('\t')[11:]
                 if tok.startswith('NM:i:')][0]
        rec_inserts = rec_deletions = 0
        for size, cigarcmd in sam.cigar_pattern.findall(record.cigar):
            if cigarcmd == 'D':
                rec_deletions += int(size)
            elif cigarcmd == 'I':
                rec_inserts += int(size)

        rec_substs = nmtag - rec_inserts - rec_deletions
        assert rec_substs >= 0

        cnt_subst += rec_substs * nreads
        cnt_del += rec_deletions * nreads
        cnt_ins += rec_inserts * nreads

    return {
        'tags': totaltags,
        'bases': totalbases,
        'substitutions': cnt_subst,
        'deletions': cnt_del,
        'insertions': cnt_ins,
    }

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Generate error frequency table '
                                     'for per-tag per-nt statistics in the report')
    parser.add_argument('alignment', metavar='name:bam', type=str, nargs='+',
                        help='Input alignment file with name.')
    parser.add_argument('--parallel', metavar='N', type=int, default=12,
                        help='Number of parallel processes.')
    args = parser.parse_args()

    return {
        'alignments': [tok.split(':') for tok in args.alignment],
        'parallel': args.parallel,
    }

if __name__ == '__main__':
    options = parse_arguments()

    final_output = sys.stdout
    csvout = csv.writer(final_output)
    csvout.writerow(['Sample', 'Subst/tag', 'Del/tag', 'Ins/tag',
                     'Subst/nt', 'Del/nt', 'Ins/nt'])

    with futures.ProcessPoolExecutor(options['parallel']) as executor:
        results = {}

        for name, bamfile in options['alignments']:
            results[name] = executor.submit(scan_alignment, bamfile)

        for name, bamfile in options['alignments']:
            result = results[name].result()
            csvout.writerow([
                name,
                result['substitutions'] / result['tags'],
                result['deletions'] / result['tags'],
                result['insertions'] / result['tags'],
                result['substitutions'] / result['bases'],
                result['deletions'] / result['bases'],
                result['insertions'] / result['bases'],
            ])

