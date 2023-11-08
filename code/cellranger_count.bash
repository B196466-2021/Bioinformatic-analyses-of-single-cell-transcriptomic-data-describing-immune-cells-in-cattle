#!/bin/bash
# Grid Engine options
#$ -N CellRanger_count
#$ -cwd
#$ -o /exports/eddie/scratch/s2170612/CellRanger_count.out
#$ -e /exports/eddie/scratch/s2170612/CellRanger_count.err
#$ -m bea
#$ -pe sharedmem 8
#$ -l h_vmem=20G
#$ -l h_rt=48:00:0
#$ -R y
# Initialise the modules framework
. /etc/profile.d/modules.sh
module load igmm/apps/cellranger/5.0.0

#### -pe sharedmem 6 --- Ani set it to 6 core, 20G. Setting a smaller amout for now to make sure it runs

cellranger count --id=11858HJPool01-N__11858HJ0001_cellranger \
--transcriptome=/exports/eddie/scratch/s2170612/ref \
--fastqs=/exports/eddie/scratch/s2170612/fastq/11858HJPool01-N__11858HJ0001 \
--sample=11858HJPool01-N__11858HJ0001 \