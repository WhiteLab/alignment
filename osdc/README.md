DNAseq Paired End pipeline
===========================
Adapted from Jason Grundstad's pipeline to run on PDC by Miguel Brown, 2015 Februrary

## MAIN:

#### basic_node_setup.py
usage: basic_node_setup.py [-h] [-id BID] [-j CONFIG_FILE] [-w WAIT]

VM spawner for pipeline processes. Sets up vm for sample analysis and attaches
storage references

optional arguments:
  -h, --help            show this help message and exit
  -id BID, --BID BID    Project Bionimbus ID
  -j CONFIG_FILE, --json CONFIG_FILE
                        JSON config file with snapshot ids and set up params
  -w WAIT, --wait WAIT  Wait time before giving up on spawning an image.
                        Reommended value 300 (in seconds)

#### parse_qc.pl - run at end of pipeline to gather qc stats

#### pipeline_wrapper.py 
usage: pipeline_wrapper.py [-h] [-f FN] [-j CONFIG]

Pipeline wrapper script to process multiple paired end set serially.

optional arguments:
  -h, --help            show this help message and exit
  -f FN, --file FN      File with bionimbus ID, seqtype and sample lane list
  -j CONFIG, --json CONFIG
                        JSON config file with tool locations, reference
                        locations, and data staging informationlane list

## MODULES:
#### date_time.py
Simple helper module that prints the current timestamp

#### pipeline.py
usage: pipeline.py [-h] [-f1 END1] [-f2 END2] [-t SEQTYPE] [-j CONFIG_FILE]

DNA alignment paired-end QC pipeline

optional arguments:
  -h, --help            show this help message and exit
  -f1 END1, --file1 END1
                        First fastq file
  -f2 END2, --file2 END2
                        Second fastq file
  -t SEQTYPE, --seqtype SEQTYPE
                        Type of sequencing peformed. Likely choices are
                        genome, exome, and capture
  -j CONFIG_FILE, --json CONFIG_FILE
                        JSON config file containing tool and reference
                        locations

##### Runs the following submodules in order:
1. fastx
2. bwa_mem_pe
3. novosort_sort_pe
(picard_sort_pe - deprecated in favor of novosort, but still avaialable)
4. picard_rmdup
5. flagstats
6. picard_insert_size 
7. coverage
8. upload_to_swift

## Pipeline submodule descriptions:
#### fastx.py
usage: fastx.py [-h] [-f FASTX_TOOL] [-sa SAMPLE] [-f1 END1] [-f2 END2]

FASTX quality stats module. Provides quality stats for fastq file and is
independent of alignment.

optional arguments:
  -h, --help            show this help message and exit
  -f FASTX_TOOL, --fastx FASTX_TOOL
                        Location of fastx_quality_stats tool. Version 0.0.13.2
                        preferred.
  -sa SAMPLE, --sample SAMPLE
                        Sample/location name prefix
  -f1 END1, --file1 END1
                        First of paired-end fastq file
  -f2 END2, --file2 END2
                        Second of paired-end fastq file

#### bwa_mem_pe.py
usage: bwa_mem_pe.py [-h] [-b BWA_TOOL] [-rg RGRP] [-br BWA_REF] [-f1 END1]
                     [-f2 END2] [-s SAMTOOLS_TOOL] [-sr SAMTOOLS_REF]
                     [-sa SAMPLE] [-l LOG_DIR]

BWA paired-end alignment module. Typically run first in pipeline.

optional arguments:
  -h, --help            show this help message and exit
  -b BWA_TOOL, --bwa BWA_TOOL
                        Location of bwa alignment tool. Version 0.7.8
                        preferred.
  -rg RGRP, --RGRP RGRP
                        SAM header read group string
  -br BWA_REF, --bwa_reference BWA_REF
                        Location of bwa reference file
  -f1 END1, --file1 END1
                        First of paired-end fastq file
  -f2 END2, --file2 END2
                        Second of paired-end fastq file
  -s SAMTOOLS_TOOL, --samtools SAMTOOLS_TOOL
                        Location of samtools tool. Version 1.19 preferred.
  -sr SAMTOOLS_REF, --samtools_reference SAMTOOLS_REF
                        Location of samtools reference
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location
#### novosort_sort_pe.py 
usage: novosort_sort_pe.py [-h] [-n NOVOSORT] [-sa SAMPLE] [-l LOG_DIR]

novosort tool to sort BAM module.

optional arguments:
  -h, --help            show this help message and exit
  -n NOVOSORT, --novosort NOVOSORT
                        novosort binary location
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location

##### picard_sort_pe.py - available, but deprecated in favor of novosort
usage: picard_sort_pe.py [-h] [-j JAVA_TOOL] [-p PICARD_TOOL] [-pt PICARD_TMP]
                         [-sa SAMPLE] [-l LOG_DIR]

Picard tools to sort BAM module.

optional arguments:
  -h, --help            show this help message and exit
  -j JAVA_TOOL, --java JAVA_TOOL
                        Java location directory, version jdk1.7.0_45 preferred
  -p PICARD_TOOL, --picard PICARD_TOOL
                        Picard jar file location
  -pt PICARD_TMP, --picard_temp PICARD_TMP
                        Picard temp folder location to create
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location

#### picard_rmdup.py
usage: picard_rmdup.py [-h] [-j JAVA_TOOL] [-p PICARD_TOOL] [-pt PICARD_TMP]
                       [-sa SAMPLE] [-l LOG_DIR]

Picard tools remove duplicate module. Removes duplicates from BAM file, run
after sorting BAM.

optional arguments:
  -h, --help            show this help message and exit
  -j JAVA_TOOL, --java JAVA_TOOL
                        Java location directory, version jdk1.7.0_45 preferred
  -p PICARD_TOOL, --picard PICARD_TOOL
                        Picard jar file location
  -pt PICARD_TMP, --picard_temp PICARD_TMP
                        Picard temp folder location to create
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location

#### flagstats.py
usage: flagstats.py [-h] [-s SAMTOOLS_TOOL] [-sa SAMPLE]

Flag stats from samtools module. Assumes bwa alignent and picard tools have
been run to sort bam and remove duplicates.

optional arguments:
  -h, --help            show this help message and exit
  -s SAMTOOLS_TOOL, --samtools SAMTOOLS_TOOL
                        Location of samtools tool. Version 1.19 preferred.
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix

#### picard_insert_size.py
usage: picard_insert_size.py [-h] [-j JAVA_TOOL] [-p PICARD_TOOL] [-sa SAMPLE]
                             [-l LOG_DIR]

Picard collect insert size metrics module. Gathers insert size metrics, run
after removing BAM duplicates.

optional arguments:
  -h, --help            show this help message and exit
  -j JAVA_TOOL, --java JAVA_TOOL
                        Java location directory, version jdk1.7.0_45 preferred
  -p PICARD_TOOL, --picard PICARD_TOOL
                        Picard jar file location
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -l LOG_DIR, --log LOG_DIR
                        LOG directory location

#### coverage.py 
usage: coverage.py [-h] [-bt BEDTOOLS2_TOOL] [-sa SAMPLE] [-c COVERAGE]
                   [-bf BED_FILE]

Bedtools coverage calculation module. Typically run last in pipeline. See
coverage parameter.

optional arguments:
  -h, --help            show this help message and exit
  -bt BEDTOOLS2_TOOL, --bedtools BEDTOOLS2_TOOL
                        Location of bedtools2 tool.
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -c COVERAGE, --coverage COVERAGE
                        Name of submodule to run. Choose from genome, exome,
                        capture or all.
  -bf BED_FILE, --bed_file BED_FILE
                        Bedfile list. If running all, list as string in order
                        format 'exome,genome,capture'. Else, just list the one
                        bed file

#### upload_to_swift.py
usage: upload_to_swift.py [-h] [-o OBJ] [-c CONT]

Uploads current directory contents to specified object and container

optional arguments:
  -h, --help            show this help message and exit
  -o OBJ, --object OBJ  Swift object name to uplaod to. i.e. PANCAN
  -c CONT, --container CONT
                        Swfit container name to upload to. i.e.
                        ALIGN/2015-1234

#### merge_qc_stats.py 
usage: merge_qc_stats.py [-h] [-o OBJ] [-c CONT] [-l LANE_LIST]

Uses pipeline lane list to create a summary table of qc stats

optional arguments:
  -h, --help            show this help message and exit
  -o OBJ, --object OBJ  Swift object name, i.e. PANCAN
  -c CONT, --container CONT
                        Swift container prefix, i.e. RAW/2015-1234
  -l LANE_LIST, --lane_list LANE_LIST
                        Original lane list used to run pipeline

novosort_merge_pe.py                                                                                                          
usage: novosort_merge_pe.py [-h] [-sa SAMPLE] [-o OBJ] [-j CONFIG_FILE]
                            [-w WAIT]

novosort tool to merge BAM files module.

optional arguments:
  -h, --help            show this help message and exit
  -sa SAMPLE, --sample SAMPLE
                        Sample/project name prefix
  -o OBJ, --obj OBJ     Swift object name
  -j CONFIG_FILE, --json CONFIG_FILE
                        JSON config file with tool and ref locations
  -w WAIT, --wait WAIT  Wait time to download bam files. 900 (seconds)
                        recommended
## UTILITY:
#### attach_cinder.py
usage: attach_cinder.py [-h] [-sid SID] [-vid VID] [-id BID] [-s SIZE]
                        [-ip IP] [-w WAIT]

Attaches cinder volume with references to existing vm

optional arguments:
  -h, --help            show this help message and exit
  -sid SID, --snapshot-id SID
                        ID of snapshot. Use cinder to find
  -vid VID, --virt-mach VID
                        Virtual machine id. Use Nova to find
  -id BID, --BID BID    Bionimbpus project id
  -s SIZE, --size SIZE  Cinder reference size. Recommended value 200 (in GB)
  -ip IP, --ip_add IP   VM IP address
  -w WAIT, --wait WAIT  Wait time before giving up on spawning an image.
                        Recommended value 300 (in seconds)

#### cleanup.py                                                                                                                                               
usage: cleanup.py [-h] [-cid CID] [-vid VID] [-id BID] [-ip VIP]

Breaks down a vm built with the standard of having the bionimbus project ID as
part of it's name and attached to a reference

optional arguments:
  -h, --help            show this help message and exit
  -cid CID, --cinder-id CID
                        ID of attached cinder volume. Use cinder to find
  -vid VID, --virt-mach VID
                        Virtual machine id. Use Nova to find
  -id BID, --BID BID    Bionimbpus project id
  -ip VIP, --ip_add VIP
                        VM IP address

##### hg19_pe_config.json
JSON config file with standard references and tools locations

#### mount.sh
Command basic_node_setup.py uses to mount a reference to a vm.

#### setup_vm.py
usage: setup_vm.py [-h] [-id BID] [-im IMAGE] [-w WAIT] [-f FLAVOR] [-k KEY]

VM spawner for pipeline processes

optional arguments:
  -h, --help            show this help message and exit
  -id BID, --BID BID    Project Bionimbus ID
  -im IMAGE, --image IMAGE
                        Image id to spawn
  -w WAIT, --wait WAIT  Wait time before giving up on spawning an image.
                        Reommended value 300 (in seconds)
  -f FLAVOR, --flavor FLAVOR
                        Image "flavor" to spawn
  -k KEY, --key KEY     Image key-pair to use

##### std_vm_config.json
JSON configuration parameters for creating a pipeline vm and attaching reference storage to it

#### unmount.sh
Command cleanup.py uses to unmount a reference from a vm.
