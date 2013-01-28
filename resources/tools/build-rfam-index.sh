#!/bin/sh
TMPDIR="$1/rfam-index"
RFAMREF="$2"
RFAMFULL="$3"
SPECIES="$4"
GENOME2BIT="$5"
OUTPUT="$6"

mkdir -p $TMPDIR

zgrep -E '^#=GF (AC|ID|DE|TP)' $RFAMFULL \
 | sed -e 's/#=GF \([A-Z][A-Z]\)   \(.*\)$/\1	\2/' \
 | awk -F'	' '{
  if ($1 == "AC" || $1 == "ID" || $1 == "DE") {
    printf "%s	", $2;
  }
  else if ($1 == "TP") {
    print $2;
  }
}' > $TMPDIR/Rfam.catalogue

zcat $RFAMREF > $TMPDIR/Rfam.fasta

grep -i "$SPECIES" $TMPDIR/Rfam.fasta | awk '{print $1}' | sed -e 's/^>//' \
  > $TMPDIR/sequenceids

faSomeRecords $TMPDIR/Rfam.fasta $TMPDIR/sequenceids $TMPDIR/Rfam.excerpt.fasta

blat -q=rna -minIdentity=85 $GENOME2BIT \
  $TMPDIR/Rfam.excerpt.fasta $TMPDIR/Rfam.psl

pslCDnaFilter -globalNearBest=0.01 -minCover=0.85 $TMPDIR/Rfam.psl \
  $TMPDIR/Rfam.bestHits.psl

pslToBed $TMPDIR/Rfam.bestHits.psl $TMPDIR/rfam.bed

cat $TMPDIR/rfam.bed | python tools/rfam-catalog-classify.py $TMPDIR/Rfam.catalogue \
  | gzip - > $OUTPUT

rm -rf $TMPDIR
