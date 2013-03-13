#!/bin/sh  
executable() {
  if [ ! -x "`which $1`" ]; then
    if [ -z "$2" ]; then
      PROGNAME="$1"
    else
      PROGNAME="$2"
    fi

    echo "An executable \`$1' is not found. Please install $PROGNAME."
    return 2
  fi
}

python2mod() {
  if [ -z "$2" ]; then
    MODNAME="$1"
  else
    MODNAME="$2"
  fi

  if ! python -c "import $1" >/dev/null 2>&1; then
    echo "A Python 2.x module \`$1' is not found. Please install $MODNAME."
  fi
}

python3mod() {
  if [ -z "$2" ]; then
    MODNAME="$1"
  else
    MODNAME="$2"
  fi

  if ! python3 -c "import $1" >/dev/null 2>&1; then
    echo "A Python 3.x module \`$1' is not found. Please install $MODNAME."
  fi
}

#========================
# Check for executables
#========================

executable "pigz"
executable "wget"
executable "snakemake"
executable "echidna"
executable "bamtools"
executable "bam2fastx" "TopHat"

# UCSC Genome Browser tools
for name in pslCDnaFilter faSomeRecords twoBitToFa pslCDnaFilter pslToBed blat; do
  executable "$name" "Jim Kent's tools"
done

# GMAP/GSNAP
for name in gmap_build psl_splicesites iit_store gsnap snpindex; do
  executable "$name" "gmap/gsnap"
done

# A. Gordon's FASTX_Toolkit
for name in fastx_clipper fastq_quality_trimmer fastq_quality_filter fastx_collapser \
            fastx_trimmer fastx_artifacts_filter; do
  executable "$name" "FASTX_Toolkit"
done

# bedtools
for name in bedtools intersectBed; do
  executable "$name" "bedtools"
done

# samtools
for name in samtools bcftools vcfutils.pl; do
  executable "$name" "samtools"
done


#==============================
# Check for Python 2.x modules
#==============================

if executable "python" "Python 2.x"; then
  python2mod lzo
  python2mod Bio BioPython
  python2mod bx "James Taylor's bx-python"
  python2mod numpy
  python2mod matplotlib
  python2mod pysam
  python2mod rpy
fi


#==============================
# Check for Python 3.x modules
#==============================

if executable "python3" "Python 3.x"; then
  python3mod colorama
  python3mod jinja2
  python3mod numpy
  python3mod yaml
fi
