#!/bin/bash
#SBATCH --job-name=LargeProfile
#SBATCH --cpus-per-task 32

R CMD BATCH --vanilla testLargeProfile_singleCore.R testLargeProfile_singleCore_Oct31.Rout