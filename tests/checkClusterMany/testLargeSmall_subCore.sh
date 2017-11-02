#!/bin/bash
#SBATCH --job-name=LargeSmall
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_LargeSmall_subCore.txt; sleep 15; done &
R CMD BATCH --vanilla testLargeSmall_subCore.R  memTestDirectory/testLargeSmall_subCore_Oct31.Rout