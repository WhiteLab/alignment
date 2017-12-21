#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
export PERL5LIB=$PERL5LIB:/cephfs/users/mbrown/PIPELINES/TOOLS/cpanm/lib/perl5
source /cephfs/users/mbrown/.virtualenvs/DNAseq/bin/activate
$pipeline --tumor $tumor --normal $normal --json $j