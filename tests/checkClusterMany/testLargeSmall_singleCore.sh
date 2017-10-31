#!/bin/bash
#SBATCH --job-name=LargeSmall_single
#SBATCH --cpus-per-task 10
R CMD BATCH --vanilla testLargeSmall_singleCore.R testLargeSmall_singleCore_Oct31.Rout