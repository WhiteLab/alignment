#!/usr/bin/env bash

#SBATCH --mail-user=miguelb@uchicago.edu
source /cephfs/users/mbrown/.virtualenvs/DNAseq/bin/activate
mkdir TMP $tmp
$novosort -c $threads -m $ram -o $in_bam --assumesorted --index --tmpdir ./TMP $bam_string 2>> $loc
$java_tool -Xmx$ram -jar $picard_tool MarkDuplicates CREATE_INDEX=true TMP_DIR=$tmp REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_RECORDS_IN_RAM=$recs +  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 INPUT=$in_bam OUTPUT=$out_bam METRICS_FILE=$mets VALIDATION_STRINGENCY=LENIENT
rmdir TMP $tmp
rm $in_bam $in_bai