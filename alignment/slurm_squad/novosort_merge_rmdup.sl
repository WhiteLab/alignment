#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
source /cephfs/users/mbrown/.virtualenvs/DNAseq/bin/activate
$novosort -c $threads -m $ram --rd -o $sample.merged.final.bam --index --tmpdir ./TMP $bam_string 2>> $loc
