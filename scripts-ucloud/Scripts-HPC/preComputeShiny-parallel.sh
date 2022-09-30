#!/bin/bash
#
#SBATCH --account sduvarcall_slim		# account
#SBATCH --nodes 1						# number of nodes
#SBATCH --time 24:00:00					# max time (HH:MM:SS)

cores=24

#for i in {1..1}
#for i in {1..4}
for i in $(seq 1 1 $cores)
do
	echo $i
	Rscript --vanilla preComputeShinyData.R &
	sleep 5
done
sleep 86400 # Keep alive 24 hours - perhaps make more advanced at some point