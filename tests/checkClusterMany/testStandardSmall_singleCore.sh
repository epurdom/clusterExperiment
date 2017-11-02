#!/bin/bash
#SBATCH --job-name=StandardSmall_single
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_LargeProfile_singleCore.txt; sleep 15; done &

R CMD BATCH --vanilla testStandardSmall_singleCore.R testStandardSmall_singleCore_Oct31.Rout