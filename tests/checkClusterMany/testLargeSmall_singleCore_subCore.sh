#!/bin/bash
#SBATCH --job-name=LargeSub
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_LargeSmall_singleCore_subCore.txt; sleep 15; done &

R CMD BATCH --vanilla testLargeSmall_singleCore_subCore.R memTestDirectory/testLargeLargeSmall_singleCore_subCore_Oct31.Rout