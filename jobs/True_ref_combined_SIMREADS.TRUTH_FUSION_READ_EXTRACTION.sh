#!/bin/bash

#$ -N FUSION_READ_EXTRACTION_JOBNAME
#$ -cwd
#$ -t 1::1
#$ -pe threaded 8
#$ -l h_vmem=5G
#$ -l h_rt=0:10:0
#$ -e .//logs/True_ref_combined_SIMREADS.TRUTH_FUSION_READ_EXTRACTION.err
#$ -o .//logs/True_ref_combined_SIMREADS.TRUTH_FUSION_READ_EXTRACTION.log

echo `date`: Running on `uname -n`

if [ -e .//logs/True_ref_combined_SIMREADS.TRUTH_SPLIT_VCF.done ]; then
  if [ ! -e True_ref_combined_SIMREADS.TRUTH_FUSION_READ_EXTRACTION.done];then
    bash /hpc/cog_bioinf/kloosterman/users/mroosmalen/NanoFG/github/NanoFG/scripts/fusion_gene_read_extraction.py     -b True_ref_combined_SIMREADS.sorted.bam     -v .//split_vcf/$SGE_TASK_ID.vcf     -o .//split_vcf

