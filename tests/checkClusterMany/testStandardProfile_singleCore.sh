#!/bin/bash
#SBATCH --job-name=StandardProfile
#SBATCH --cpus-per-task 10

R CMD BATCH --vanilla testStandardProfile_singleCore.R testStandardProfile_singleCore_Oct31.Rout