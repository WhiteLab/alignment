#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
source /cephfs/users/mbrown/.virtualenvs/DNAseq/bin/activate
mkdir TMP
$novosort -c $threads -m $ram --rd -o $out_bam --index --tmpdir TMP $bam_string 2>> $loc
rmdir TMP
