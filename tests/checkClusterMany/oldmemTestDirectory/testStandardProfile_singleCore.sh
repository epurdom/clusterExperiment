#!/bin/bash
#SBATCH --job-name=StandardProfile
#SBATCH --cpus-per-task 32
while true; do free -h >> memTestDirectory/memoryLogger_StandardProfile_singleCore.txt; sleep 15; done &

R CMD BATCH --vanilla testStandardProfile_singleCore.R memTestDirectory/testStandardProfile_singleCore_Oct31.Rout