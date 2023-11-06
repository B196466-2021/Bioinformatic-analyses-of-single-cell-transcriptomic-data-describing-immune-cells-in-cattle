#!/bin/bash
#!/bin/bash
# Grid Engine options
#$ -N CellRanger_count
#$ -cwd
#$ -o /exports/cmvm/eddie/eb/groups/Hope_Group/Barbara_Shih/2021-04-28-_9649_-Jayne_Hope_scRNA-Seq_bovine_lymph/log/CellRanger_count.out
#$ -e /exports/cmvm/eddie/eb/groups/Hope_Group/Barbara_Shih/2021-04-28-_9649_-Jayne_Hope_scRNA-Seq_bovine_lymph/log/CellRanger_count.err
#$ -m bea
#$ -pe sharedmem 6
#$ -l h_vmem=20G
#$ -l h_rt=120:00:0
#$ -R y
# Initialise the modules framework
. /etc/profile.d/modules.sh
module load igmm/apps/cellranger/5.0.0

cellranger count --id=sample1 --transcriptome=/exports/cmvm/eddie/eb/groups/Hope_Group/Tiancheng/data/ref --fastqs=/exports/cmvm/eddie/eb/groups/Hope_Group/Tiancheng/data/11858HJPool01-N__11858HJ0001 --sample=11858HJPool01-N__11858HJ0001
