#!/usr/bin/env bash

cd $PBS_O_WORKDIR

/usr/lib/jvm/java-6-openjdk-amd64/bin/java -jar /home/mark/src/mutect/muTect-1.1.4.jar      \
  --analysis_type MuTect                              \
  --reference_sequence hg19.fa                        \
  --input_file:normal $1                              \
  --input_file:tumor $2                               \
  --out ${3}.calls                                    \
  --coverage_file ${3}.coverage                       \
  --intervals nimblegen.merged.intervals
