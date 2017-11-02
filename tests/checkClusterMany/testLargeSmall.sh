#!/bin/bash
#SBATCH --job-name=LargeSmall
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_LargeSmall.txt; sleep 15; done &
R CMD BATCH --vanilla testLargeSmall_singleCore.R testLargeSmall_singleCore_Oct31.Rout