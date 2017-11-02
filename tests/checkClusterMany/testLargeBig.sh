#!/bin/bash
#SBATCH --job-name=LargeBig
#SBATCH --cpus-per-task 32
while true; do free -h >> memoryLogger_LargeBig.txt; sleep 2; done &
R CMD BATCH --vanilla testLargeBig.R testLargeBig_Oct31.Rout