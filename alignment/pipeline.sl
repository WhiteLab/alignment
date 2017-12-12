#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
echo $pipeline $f1 $f2 $t $j
sh -c "$pipeline --file1 $f1 --file2 $f2 --seqtype $t --json $j"
