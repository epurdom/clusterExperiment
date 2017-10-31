#!/bin/bash
#SBATCH --job-name=StandardSmall_single
#SBATCH --cpus-per-task 10

R CMD BATCH --vanilla testStandardSmall_singleCore.R testStandardSmall_singleCore_Oct31.Rout