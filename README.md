DNAseq analysis pipeline in two parts, alignment QC and variant calling.  Tools featured are cutadapt, bwa, mutect, scalpel, platypus, variant effect predictor, and, picard.
Designed to run on an hpc with slurm job controller.
Disclaimer:  There are other modes and tools in progress outside of what's outlined here.  Use at your own risk!

## Quick Start - (ain't nobody got time for that!)
## CRASH COURSE RUN:
You don't have time to read through what each script does, and you're a BAMF.  However, being familiar with the Linux environment is highly recommended.

##### 1) Create job run files
###### a) Get list of fastq files to process from swift, one file per line, using a new-line separated list of bionimbus ids:
```
cat <bnids_list> | xargs -IFN find <fastq file home directory>/BN -name '*.gz' > fastq_list 
```
###### b) Use this list to create a run file, in this example for custom capture:
```
<DNA pipe home dir>/utility/fq2lane.py -f fastq_list -s <seq type> > lane_list
```
Typical fastq list:

2018-254_180309_NB551089_0045_AHJGNFBGX5_1_1_sequence.txt.gz
2018-254_180309_NB551089_0045_AHJGNFBGX5_1_2_sequence.txt.gz
2016-2199_180302_K00216R_0042_BHNVTKBBXX_1_1_sequence.txt.gz
2016-2199_180302_K00216R_0042_BHNVTKBBXX_1_2_sequence.txt.gz
2016-2209_180302_K00216R_0042_BHNVTKBBXX_1_1_sequence.txt.gz
2016-2209_180302_K00216R_0042_BHNVTKBBXX_1_2_sequence.txt.gz

Resultant lane_list:

2018-254	capture	180309_NB551089_0045_AHJGNFBGX5_1
2016-2199	capture	180302_K00216R_0042_BHNVTKBBXX_1
2016-2209	capture	180302_K00216R_0042_BHNVTKBBXX_1

##### c) Check config file - this file is typically in <RNA pipe hoe dir>/utility/config_files/v2.5_config.json, and can be copied and modified.  Fields that are likely to be adjusted:

    "refs":{
    "project_dir": "/cephfs/PROJECTS/",
    "project": "PANCAN", # typically as a subdir if the project dir
    "align_dir": "ALIGN_RNASEQ", # typically as a subdir of project
    "analysis": "ANALYSIS",
    "annotation": "ANNOTATION",
    "config": "/cephfs/users/mbrown/RNAseq/utility/config_files/slurm_config.json" # full path to config file being used

    },
    "params":{
	"threads":"8",
	"ram":"30",
	"user": "mbrown", # this and next param to set acls
    "group": "CPCI",
    }
See example config in utility/config_files/utility/slurm_DNAseq_vep91_config.json
##### 2) Pipeline run - QC:

```
<DNAseq pipe dir>/alignment/pipeline_wrapper.py -f lane_list.txt -j modified_config.json 2> run.log
```

The pipeline will iterate throught the list, Create run directories with the align_dir specified above, trim fastqs, align with bwa, and qc the fastq and bam files.  Logs track most of the steps.  Multiple qc tables can ba concatenated for convenience after run using:
```
<DNAseq pipe dir>/alignment/merge_qc_stats.py -h
usage: merge_qc_stats.py [-h] [-p P_DIR] [-d F_DIR] [-l LANE_LIST]

Uses pipeline lane list to create a summary table of qc stats

optional arguments:
  -h, --help            show this help message and exit
  -p P_DIR, --project_dir P_DIR
                        Project directory, i.e. /cephfs/PROJECTS/PANCAN
  -d F_DIR, --file_dir F_DIR
                        file directory prefix with qc files
  -l LANE_LIST, --lane_list LANE_LIST
                        Original lane list used to run pipeline
```
Resultant files are organized into subdirectories: BAM, LOGS, and QC

##### 3) Pipeline run - Variant calling:
This will take the same pipeline run files from alignment along with a parameter indicating where to start the process/what part of the process to execute and run mutect, scalpel, and platypus for variant calling, and variant effect predictor for annotation

```
<DNAseq pipe dir>/analysis/variant_annot_pipeline_wrapper.py -h
usage: variant_annot_pipeline_wrapper.py [-h] [-sp SAMPLE_PAIRS]
                                         [-j CONFIG_FILE] [-e ESTEP]

Pipeline for variant calls and annotation using mutect and snpEff

optional arguments:
  -h, --help            show this help message and exit
  -sp SAMPLE_PAIRS, --sample-pairs SAMPLE_PAIRS
                        Tumor/normal sample pair list
  -j CONFIG_FILE, --json CONFIG_FILE
                        JSON config file with tool and ref locations, USE FULL
                        PATH!
  -e ESTEP, --execute ESTEP
                        Steps to start at, valid entries are start, indel,
                        snv-annot, snv-indel, germ-call, germ-annot
```
Output will be in the ANALYSIS and ANNOTATION subdirectories named as specificed in the config file. Notably, variant call summary files will be generated in .xls format.