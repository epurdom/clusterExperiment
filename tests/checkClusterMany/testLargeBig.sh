#!/bin/bash
#SBATCH --job-name=LargeBig
#SBATCH --cpus-per-task 32
R CMD BATCH --vanilla testLargeBig.R testLargeBig_Oct31.Rout