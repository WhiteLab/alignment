#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
source /cephfs/users/mbrown/.virtualenvs/bin/activate
$pipeline --file1 $f1 --file2 $f2 --seqtype $t --json $j
