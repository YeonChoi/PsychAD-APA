#!/bin/bash

#BSUB -J snakemake
#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -n 1
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=3000]
#BSUB -W 140:00
#BSUB -eo ./logs/snakemake.err
#BSUB -oo ./logs/snakemake.out
#BSUB -L /bin/bash

ml purge
ml samtools/1.13

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate base_APA

# For many node
snakemake --profile lsf --restart-times 1 --latency-wait 10 --use-conda

# for single node
#snakemake -j --restart-times 1 --latency-wait 10 #--batch all=1/3

# Below are from PNM
#snakemake --profile prachu_lsf --restart-times 2 --latency-wait 10 #--batch all=1/3
#snakemake -j 4 --restart-times 2 --latency-wait 10
