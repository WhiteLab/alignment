#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu

source /cephfs/users/mbrown/PIPELINES/TOOLS/miniconda3/bin/activate strelka
source /cephfs/users/mbrown/.virtualenvs/DNAseq/bin/activate
$strelka_pipe --tumor $tumor --normal $normal -j $j