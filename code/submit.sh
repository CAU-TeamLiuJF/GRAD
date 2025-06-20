#!/bin/bash
#SBATCH -J GRAD
#SBATCH -o %j_GRAD.out
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH -p Cnode_all


			 
set -e
module unload R

source /public/home/liujf/software/program/Miniconda3-py311_23.5.0-3/bin/activate /public/home/liujf/workspace/xueyh/software/r4.3

echo "[`date`] Using R at: $(which Rscript)"
echo "[`date`] R version: $(Rscript -e 'cat(R.version.string)')"

Rscript data_prepare.R 
Rscript main.R

echo "calculate done!"
