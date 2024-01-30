#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=0:20:00

#SBATCH --account=def-mhbrice

module load gcc/9.3.0 r/4.3.1

Rscript R-scripts/03-Analysis/INLA-test.R > output_inla_test.txt
