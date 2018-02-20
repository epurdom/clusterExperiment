#!/bin/bash
#SBATCH --job-name=StBigSub
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_StandardBig_subCore.txt; sleep 15; done &
R CMD BATCH --vanilla testStandardBig_subCore.R  memTestDirectory/testStandardBig_subCore_Oct31.Rout