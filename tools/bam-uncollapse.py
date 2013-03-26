#!/usr/bin/env python
import sys
import os

def process(inputfile, outputfile):
    with os.popen('samtools view -bS - > {}'.format(outputfile), 'w') as output:
        header = os.popen('samtools view -H {}'.format(inputfile)).read()
        output.write(header)

        body = os.popen('samtools view {}'.format(inputfile))
        for line in body:
            name, remaining = line.split('\t', 1)
            nreads = int(name.split('-')[1])
            for i in xrange(nreads):
                output.write(name + '.{}\t'.format(i) + remaining)

if __name__ == '__main__':
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
    process(inputfile, outputfile)
