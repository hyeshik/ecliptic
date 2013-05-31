#!/usr/bin/env python
from ecliptic.support.sam import (
    parse_sam_simple, F_REVERSE_STRAND, calculate_cigar_length,
    cigar_replace_softclips, cigar_pattern)
from ecliptic.support.sequtils import GiantFASTAFile

def process(inpfile, output, refseqfile, delpos='5'):
    if delpos == '5':
        is_leftbind = lambda line: not line.flag & F_REVERSE_STRAND
    elif delpos == '3':
        is_leftbind = lambda line: line.flag & F_REVERSE_STRAND
    else:
        raise ValueError('Unknown binding oriention %s' % repr(delpos))

    for line in parse_sam_simple(inpfile):
        if str(line)[0] == '@' or 'D' not in line.cigar:
            output.write(str(line))
            continue

        leftpos = line.pos - 1 # to zero-based coord
        reflen, seqlen = calculate_cigar_length(line.cigar)
        refseq = refseqfile.get(line.rname, leftpos, leftpos + reflen).upper()
        assert len(refseq) == reflen

        cigar_tokens = map(list, cigar_pattern.findall(line.cigar))
        leftclip = int(cigar_tokens[0][0]) if cigar_tokens[0][1] == 'S' else 0
        rightclip = int(cigar_tokens[-1][0]) if cigar_tokens[-1][1] == 'S' else 0

        refleft = 0
        tainted = False
        for toki, (width, cmd) in enumerate(cigar_tokens):
            width = int(width)
            if cmd != 'D':
                if cmd in 'MPN':
                    refleft += width
                continue

            deleted = refseq[refleft:refleft+width]

            if len(set(deleted)) > 1:
                # do not break deletions with two or more heteronucleotide, e.g.:
                # AGGGGAGTGGTGCAGTGTGATTCCCGCCGGT
                # AGGG--GTGGTGCAGTGTGATTCCCGCCGGT
                #     ^^ 2D
                refleft += width
                continue

            lefttok = cigar_tokens[toki-1]
            righttok = cigar_tokens[toki+1]

            if is_leftbind(line): # push to the left as possible
                spanleft = refleft
                if lefttok[1] == 'M':
                    lefttokwidth = int(lefttok[0])

                    for i in xrange(spanleft - 1, spanleft - lefttokwidth - 1, -1):
                        if refseq[i] == deleted[0]:
                            spanleft = i
                        else:
                            break

                    leftexpand = refleft - spanleft

                    if leftexpand > 0:
                        lefttok[0] = str(int(lefttok[0]) - leftexpand)
                        righttok[0] = str(int(righttok[0]) + leftexpand)
                        tainted = True

            else: # push to the right as possible
                spanright = refleft + width - 1
                if righttok[1] == 'M':
                    righttokwidth = int(righttok[0])

                    for i in xrange(refleft + width, refleft + width + righttokwidth):
                        if refseq[i] == deleted[0]:
                            spanright = i
                        else:
                            break

                    rightexpand = spanright - refleft - width + 1

                    if rightexpand > 0:
                        lefttok[0] = str(int(lefttok[0]) + rightexpand)
                        righttok[0] = str(int(righttok[0]) - rightexpand)
                        tainted = True

            refleft += width

        if not tainted:
            output.write(str(line))
        else:
            fields = str(line).split('\t')
            fields[5] = ''.join(w+cmd for w, cmd in cigar_tokens)
            output.write('\t'.join(fields))


if __name__ == '__main__':
    import os, sys

    #inputfile = 'work/GSE37114-LIN28A-NKim/alignments/CLIP-35L33G.bam'
    #refgenome = 'resources/mm10/genome.fa'
    #saminput = os.popen('samtools view "{}"'.format(inputfile))

    refgenome = sys.argv[1]
    delbinding = sys.argv[2]

    process(sys.stdin, sys.stdout, GiantFASTAFile(refgenome), delbinding)

