import os, platform

# Script dir
script_dir = "/work/sduvarcall/haplotype-project/Scripts/simulate-population/"

# simulate population?
#simulate_pop = False
simulate_pop = True

# genealogy
genealogy = "star"
#genealogy = "correlated"

# Ancestral method - should only be run in parallel if simulate_pop = False
#ancestral = "mostFreqBase"
#ancestral = "branchBound"
#ancestral = "branchBoundIndep"
ancestral = "simulatedFounder"
#ancestral = "knownBreaks"

# Simulation type
# sim_type = "knownBreaks"
# sim_type = "popFreqs"
sim_type = "alleleFreqs"
#sim_type = "haplotypes"
#sim_type = "famHaplotypes"

# Chance sharing correction
cs_correction = "None"
#cs_correction = "adjustLengths"
#cs_correction = "gandolfo"

# Generation to test
first_gen = 10
start_gen = first_gen
seed = 42
samples = 100

# jump = 100
# stop_gen = start_gen+jump
# last_gen = 510
# step_gen = 20
# sims = 10

jump = 100
stop_gen = start_gen+jump
last_gen = 2510
step_gen = 50
sims = 10

# jump = 500
# stop_gen = start_gen+jump
# last_gen = 10010
# step_gen = 50
# sims = 10

# jump = 100
# stop_gen = start_gen+jump
# last_gen = 3010
# step_gen = 20
# sims = 20
while start_gen < last_gen:

	# Make last interval smaller if needed
	if stop_gen > last_gen:
		stop_gen = last_gen

	cmd = ("sbatch " + script_dir + "sbatch_simulate_population.sh " + script_dir + "simulate_population.py " + 
			str(start_gen) + " " + str(stop_gen) + " " + str(step_gen) + " " + str(samples) + " " + str(sims) + " " + str(seed) + " " + str(first_gen) + " " + str(last_gen) + " " + str(genealogy) + " " + str(ancestral) + " " + str(sim_type) + " " + str(simulate_pop) + " " + str(cs_correction))

	# print cmd
	# Print arguments
	print str(start_gen) + " " + str(stop_gen) + " " + str(step_gen) + " " + str(samples) + " " + str(sims) + " " + str(seed) + " " + str(first_gen) + " " + str(last_gen) + " " + str(genealogy) + " " + str(ancestral) + " " + str(sim_type) + " " + str(simulate_pop) + " " + str(cs_correction)

	if platform.system() == "Linux":
		os.system(cmd)
	
	print start_gen, stop_gen
	
	start_gen += jump + step_gen
	stop_gen += jump + step_gen
	
	# if samples < 100:
	# 	samples += 10
	# else: # samples >= 100
	# 	samples += 50
	
	# Not sure if needs to change
	#seed += 1

print "Jobs submitted to slurm"