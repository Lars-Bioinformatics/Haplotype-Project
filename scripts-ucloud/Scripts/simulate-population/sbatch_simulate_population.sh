#!/bin/bash
#SBATCH --account sduvarcall_slim
#SBATCH --nodes 1
#SBATCH --time 24:00:00

python2 $1 --start $2 --stop $3 --step $4 --samples $5 --sims $6 --seed $7 --start_folder $8 --stop_folder $9 --genealogy ${10} --ancestral_method ${11} --simulation_type ${12} --simulate_population ${13} --cs_correction ${14}
