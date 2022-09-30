import os, platform

# Script dir
script_dir = "/work/sduvarcall/haplotype-project/Scripts/simulate-population/"
#script_dir = "/Users/lars/Documents/PhD/haplotype-project/Scripts/simulate-population/"

# simulate population?
# simulate_pop = False
simulate_pop = True

# genealogy
genealogy = "star"
# genealogy = "correlated"

# Ancestral method
#ancestral = "mostFreqBase"
#ancestral = "branchBound"
# ancestral = "branchBoundIndep"
ancestral = "simulatedFounder"
#ancestral = "knownBreaks"

# Simulation type
# sim_type = "knownBreaks"
# sim_type = "popFreqs"
sim_type = "famHaplotypes"
# sim_type = "haplotypes"

# Generation to test
test_gen = 50
seed = 42
sims = 50

# Testing sample size
samples = 10
while samples < 1001:
	# cmd = ("sbatch --account sduvarcall_slim --nodes 1 --time 1:0 python2 " + script_dir + "simulate_population.py " +
	# 		  "--start " + str(test_gen) + " --stop " + str(test_gen) + " --step " + str(test_gen) +
	# 		  " --samples " + str(samples) + " --seed " + str(seed) + " --sims " + str(sims))
	# cmd = ("sbatch " + script_dir + "sbatch_simulate_population.sh " + script_dir + "simulate_population.py " +
	# 		str(test_gen) + " " + str(test_gen) + " " + str(test_gen) + " " + str(samples) + " " + str(sims) + " " + str(seed) + " " + str(-1) + " " + str(-1) + " " + str(genealogy) + " " + str(ancestral) + " " + str(sim_type) + " " + str(simulate_pop))
	
	cmd="python2 "+script_dir+"simulate_population.py --start " + str(test_gen) + " --stop " + str(test_gen) + " --step " + str(test_gen) + " --sims " + str(sims) + " --sims_type " + sim_type + " --ancestral_method " + ancestral + " --genealogy " + genealogy + " --samples " + str(samples) + " --extra_option none &"
	
	print cmd
	# Print arguments
	print str(test_gen) + " " + str(test_gen) + " " + str(test_gen) + " " + str(samples) + " " + str(sims) + " " + str(seed) + " " + str(genealogy) + " " + str(ancestral) + " " + str(sim_type) + " " + str(simulate_pop)

	if platform.system() == "Linux":
		os.system(cmd)
	
	if samples < 100:
		samples += 10
	else: # samples >= 100
		samples += 50
	
	# Not sure if needs to change
	#seed += 1

# print "Jobs submitted to slurm"
print "Done"