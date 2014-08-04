#!/usr/bin/env bash

#PBS -N bowtie2
#PBS -l nodes=1:ppn=8
#PBS -j oe

# Check if we're in a Torque environment.
if [[ -z $PBS_ENVIRONMENT ]]; then
  # Not in a Torque environment - submit to queue and exit.
  [ -z "$1" ] && echo "At least one file must be provided." && exit 1
  [ ! -z "$2" ] && qsub -v F1=$1,F2=$2 $0 || qsub -v F1=$1 $0
  exit
fi

# Global settings.
GENEDB=/glusterfs/users/mark/data/hg19/bowtie/Homo_sapiens.GRCh37.73.fil.gtf
BOWTIE=/glusterfs/users/mark/src/bowtie2/bowtie2
STOOLS=/glusterfs/users/mark/src/samtools/samtools

# Bowtie2 arguments.
B2ARGS=()
B2ARGS+=(-q) # FastQ inputs.
B2ARGS+=(--phred33) # Utilizing Phred33 quality scores.
B2ARGS+=(--threads $PBS_NUM_PPN) # Number of threads to utilize.
B2ARGS+=(--non-deterministic) # Seed random number generator with current time.

cd $PBS_O_WORKDIR
if [ -z "$F2" ]; then
  $BOWTIE $B2ARGS -x $GENEDB -U <(zcat ${F1})
else
  $BOWTIE $B2ARGS -x $GENEDB -1 <(zcat ${F1}) -2 <(${F2})
fi
