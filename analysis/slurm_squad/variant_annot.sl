#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
source /cephfs/users/mbrown/.virtualenvs/DNAseq/bin/activate
$pipeline --tumor $tumor --normal $normal --json $j