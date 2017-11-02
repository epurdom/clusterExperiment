#!/bin/bash
#SBATCH --job-name=LargeProfile
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_LargeProfile_singleCore.txt; sleep 15; done &

R CMD BATCH --vanilla testLargeProfile_singleCore.R memTestDirectory/testLargeProfile_singleCore_Oct31.Rout