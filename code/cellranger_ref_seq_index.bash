#!/bin/bash
# Grid Engine options
#$ -N CellRanger_reference
#$ -cwd
#$ -o /exports/cmvm/eddie/eb/groups/Hope_Group/Tiancheng/data/CellRanger_reference.out
#$ -e /exports/cmvm/eddie/eb/groups/Hope_Group/Tiancheng/data/CellRanger_reference.err
#$ -pe sharedmem 4
#$ -l h_vmem=20G
#$ -l h_rt=1:00:0
#$ -R y
# Initialise the modules framework
. /etc/profile.d/modules.sh
module load igmm/apps/cellranger/5.0.0

WORKING_DIR=$(pwd)
DATA_DIR=${data_dir}
IN_GFP=${in_gfp}
IN_FASTA=${in_fasta}

CELLRANGER_DIR=$WORKING_DIR/$DATA_DIR/ref/cellranger
GTF_BASENAME=$(basename $IN_GFP)
GTF_BASENAME_FILTERED=${GTF_BASENAME/.gtf/_filtered.gtf}
FASTA_BASENAME=$(basename $IN_FASTA)
FASTA_BASENAME_short=${FASTA_BASENAME/.fa/}


mkdir -p $(dirname $CELLRANGER_DIR) $CELLRANGER_DIR


# Making own ref file
# GTF file
cd $CELLRANGER_DIR

# Filter GTF
cellranger mkgtf $GTF_BASENAME $GTF_BASENAME_FILTERED \
                   --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_LV_gene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_C_gene

# Make reference
cellranger mkref --genome=$FASTA_BASENAME_short --fasta=$FASTA_BASENAME --genes=$GTF_BASENAME_FILTERED --nthreads=3