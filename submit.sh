#!/bin/bash
#SBATCH -J prepare_Rdata
#SBATCH -o %j_prepare.out
#SBATCH -N 1
#SBATCH --cpus-per-task=5
#SBATCH --mem=50G
#SBATCH -p Cnode_all


			 
set -e
module unload R

source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3

echo "[`date`] Using R at: $(which Rscript)"
echo "[`date`] R version: $(Rscript -e 'cat(R.version.string)')"

Rscript 1.data_prepare.R 
Rscript 2.analysis.R

echo "calculate done!"
