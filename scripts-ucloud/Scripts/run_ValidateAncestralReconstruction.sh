#!/bin/bash
#
#SBATCH --account sduvarcall_slim	# account
#SBATCH --nodes 1			# number of nodes
#SBATCH --time 24:00:00			# max time (HH:MM:SS)
#SBATCH --output /work/sduvarcall/haplotype-project/logs/slurm-%j.out

Rscript /work/sduvarcall/haplotype-project/Scripts/run_ValidateAncestralReconstruction.R
