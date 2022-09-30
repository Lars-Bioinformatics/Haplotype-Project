#!/bin/bash
#SBATCH --account sduvarcall_slim		# account
#SBATCH --nodes 1						# number of nodes
#SBATCH --time 24:00:00					# max time (HH:MM:SS)

# Script dir
script_dir="/work/sduvarcall/haplotype-project/Scripts/simulate-population/"

# cmd="python2 ${script_dir}simulate_population.py --stop 2010 --sims 50 --step 50 --sims_type famHaplotypes --ancestral_method simulatedFounder --genealogy star"

cmd="python2 ${script_dir}simulate_population.py --stop 510 --sims 50 --step 10 --sims_type famHaplotypes --ancestral_method simulatedFounder --genealogy star --samples 20 --extra_option none"

echo $cmd

$cmd