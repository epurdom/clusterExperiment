#!/bin/bash
#SBATCH --job-name=LargeBig
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_LargeBig.txt; sleep 15; done &
R CMD BATCH --vanilla testLargeBig.R memTestDirectory/testLargeBig_Oct31.Rout