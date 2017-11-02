#!/bin/bash
#SBATCH --job-name=StandBig
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_StandardBig.txt; sleep 15; done &
R CMD BATCH --vanilla testStandardBig.R  memTestDirectory/testStandardBig_Oct31.Rout