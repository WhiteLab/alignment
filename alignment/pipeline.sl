#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
$pipeline -f1 $f1 -f2 $f2 -t $t -j $j
