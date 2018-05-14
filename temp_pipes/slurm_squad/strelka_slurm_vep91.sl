#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
export PERL5LIB=$PERL5LIB:/cephfs/users/mbrown/PIPELINES/TOOLS/cpanm/lib/perl5
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cephfs/software/samtools-1.6/htslib-1.6/htslib/lib
source /cephfs/users/mbrown/PIPELINES/TOOLS/miniconda3/bin/activate strelka
source /cephfs/users/mbrown/.virtualenvs/DNAseq/bin/activate
$strelka_pipe --tumor $tumor --normal $normal -j $j