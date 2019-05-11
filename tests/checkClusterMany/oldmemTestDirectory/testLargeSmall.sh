#!/bin/bash
#SBATCH --job-name=LargeSmall
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_LargeSmall.txt; sleep 15; done &
R CMD BATCH --vanilla testLargeSmall.R  memTestDirectory/testLargeSmall_Oct31.Rout