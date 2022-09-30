#! /bin/bash
#
#SBATCH --account sduvarcall_slim	# account
#SBATCH --nodes 1			# number of nodes
#SBATCH --time 24:00:00			# max time (HH:MM:SS)

Rscript classify-mutationAge.R
