#!/bin/bash
#SBATCH --job-name=StandardSmall
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_StandardSmall.txt; sleep 15; done &

R CMD BATCH --vanilla testStandardSmall.R testStandardSmall_Oct31.Rout