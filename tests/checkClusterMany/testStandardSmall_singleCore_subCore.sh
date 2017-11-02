#!/bin/bash
#SBATCH --job-name=StandSubcore
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_StandardSmall_singleCore_subCore.txt; sleep 15; done &

R CMD BATCH --vanilla testStandardSmall_singleCore_subCore.R testStandardSmall_singleCore_subCore_Oct31.Rout