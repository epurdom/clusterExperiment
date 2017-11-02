#!/bin/bash
#SBATCH --job-name=StandardSmall
#SBATCH --cpus-per-task 32

R CMD BATCH --vanilla testStandardSmall.R testStandardSmall_Oct31.Rout