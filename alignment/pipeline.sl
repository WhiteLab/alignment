#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
#SBATCH -c $cores
$SBATCH --mem $mem
#SBATCH -o $log
$pipeline $f1 $f2 $t $j
