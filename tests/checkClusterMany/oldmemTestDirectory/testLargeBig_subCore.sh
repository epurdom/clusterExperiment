#!/bin/bash
#SBATCH --job-name=LgBigSub
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_LargeBig_subCore.txt; sleep 15; done &
R CMD BATCH --vanilla testLargeBig_subCore.R memTestDirectory/testLargeBig_subCore_Oct31.Rout