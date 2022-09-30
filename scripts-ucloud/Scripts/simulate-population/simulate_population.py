#!python2

import random, sys, copy, os, platform, argparse, time
import numpy as np
from itertools import izip
from multiprocessing import Pool

# FIELDS
#star_genealogy = True
#star_genealogy = False
#genealogy = "starGenealogy" if star_genealogy else "correlatedGenealogy"
#gene = "BRCA1"
#custom_filename_string = "test" # "training" # "test"
#custom_filename_string = "training"
# if True, walk from mutation to outmost mutation. Will stop when break is found.
walk_outwards = True
# if True, allow breaks on both sides of a mutation in the same generation
allowBreaksBothSidesSameGen = True
# start_generations = 10
# stop_generations = 500
# step_generations = 10
# sim_num = 20
# max_samples = 100
#seed = 42 # training
#seed = 24 # test
#seed = 1024 # test2
#seed = 9012

# Arguments parser
def str2bool(v):
	if v.lower() in ('yes', 'true', 't', 'y', '1'):
		return True
	elif v.lower() in ('no', 'false', 'f', 'n', '0'):
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected e.g. true, false, yes, no, t, f, 1, 0.')

parser = argparse.ArgumentParser(description='Simulate population with BRCA1/2 mutation')
parser.add_argument("--gene", type=str, default="BRCA1", help="Select gene to simulate: BRCA1/2 (default: BRCA1)")
parser.add_argument("--genealogy", type=str, default="star", help="Select genealogy: star or correlated (default: star)")
parser.add_argument("--start", type=int, default=10, help="First generation to simulate (default: 10)")
parser.add_argument("--stop", type=int, default=510, help="Last generation to simulate (default: 100)")
parser.add_argument("--step", type=int, default=10, help="steps between simulated generation (default: 10)")
parser.add_argument("--start_folder", type=int, default=-1, help="Name of first generation in location name, only needed if splitting simulations on multiple runs")
parser.add_argument("--stop_folder", type=int, default=-1, help="Name of last generation in location name, only needed if splitting simulations on multiple runs")
parser.add_argument("--samples", type=int, default=100, help="Number of samples (default: 100)")
parser.add_argument("--sims", type=int, default=20, help="Number of simulations per generation (default: 20)")
parser.add_argument("--seed", type=int, default=42, help="Set seed for random number generator (default: 42)")
parser.add_argument("--sims_type", type=str, default="famHaplotypes", help="Select simulation type: knownBreaks, popFreqs, alleleFreqs, haplotypes, famHaplotypes (default: famHaplotypes)")
parser.add_argument("--ancestral_method", type=str, default="branchBoundIndep", help="Set ancestral method strategy: mostFreqBase, branchBound, branchBoundIndep, simulatedFounder (default: branchBoundIndep)")
parser.add_argument("--cs_correction", type=str, default="None", help="Select chance sharing correction approach: None, gandolfo, adjustLengths (default: None)")
parser.add_argument("--simulate_population", type=str2bool, default=True, help="Simulate population (default: True)", metavar="True/False")
parser.add_argument("--extra_option", type=str, default="compare-age-methods", help="Select option: compare-age-methods, dist-mat, haplo-length, none (default: compare-age-methods)")
args = parser.parse_args()
print args
print args.start

gene = args.gene
genealogy = args.genealogy #temp value
start_generations = args.start
stop_generations = args.stop
step_generations = args.step
sim_num = args.sims
max_samples = args.samples
seed = args.seed
simulatePop = args.simulate_population
option = args.extra_option
sim_type = args.sims_type
cs_correction = args.cs_correction

print max_samples

# If simulate known breaks, we do not need to construct ancestral haplotype
if sim_type == "knownBreaks":
	ancestral_method = "knownBreaks"
else:
	ancestral_method = args.ancestral_method

if args.start_folder == -1 or args.stop_folder == -1:
	start_folder = start_generations
	stop_folder = stop_generations
else:
	start_folder = args.start_folder
	stop_folder = args.stop_folder

# Correct values for genealogy
star_genealogy = True if genealogy == "star" else False
genealogy = "starGenealogy" if star_genealogy else "correlatedGenealogy"

random.seed(seed)
np.random.seed(seed)

if not simulatePop:
	print "Not simulating population, only running extra option"

# sys.exit()

# Count stats
counts = dict()
count_left = 0
count_right = 0
count_randvals = dict()

#gene = "test"
#stop_generations = 10
#sim_num = 1
#max_samples = 20

# WORKDIR
if platform.system() == "Linux":
	workdir = "/work/sduvarcall/haplotype-project/classify-simulated-population-age/"
	rscript_dir = "/work/sduvarcall/haplotype-project/Scripts/"
else: # Mac
	workdir = "/Users/lars/Documents/PhD/haplotype-project/classify-simulated-population-age/"
	rscript_dir = "/Users/lars/Documents/PhD/haplotype-project/Scripts/"

print start_folder
print stop_folder
print step_generations
if start_folder == stop_folder and start_folder == step_generations:
	data_dir = gene+"_simulated-test_sample_size-"+ \
    		   sim_type+"-"+genealogy+"-"+str(sim_num)+ \
			   "_simulations-generations_"+str(start_folder)+ \
			   "_"+str(stop_folder)+"_step"+str(step_generations)+ \
			   "-"+"seed_"+str(seed)+"/"
else:
	data_dir = gene+"_simulated-"+sim_type+"-"+genealogy+ \
    		   "-"+str(max_samples)+"_samples-"+str(sim_num)+ \
			   "_simulations-generations_"+str(start_folder)+ \
			   "_"+str(stop_folder)+"_step"+str(step_generations)+ \
			   "-"+"seed_"+str(seed)+"/"

## Best case ancestral
# 	data_dir = gene+"_simulated-"+genealogy+ \
# 	           "-"+str(max_samples)+"_samples-"+str(sim_num)+ \
# 	           "_simulations-generations_"+str(start_folder)+ \
# 	           "_"+str(stop_folder)+"_step"+str(step_generations)+ \
# 	           "-"+"seed_"+str(seed)+"/"

## When running in intervals on multiple computers - fix data_dir
# data_dir = gene+"_simulated-"+genealogy+ \
#            "-"+str(max_samples)+"_samples-"+str(sim_num)+ \
#            "_simulations-generations_10"+ \
#            "_10010_step"+str(step_generations)+ \
#            "-"+"seed_"+str(seed)+"/"


# Testing multicore processing
# def work(gen):
# 	print gen, gene, genealogy
# 	time.sleep(2)
# 	return(gen)
#
# pool = Pool()
# res = pool.map(work, [10,20,30])
# pool.close()
# pool.join()
# print res
# sys.exit()

def get_recomb_frac_accum(izip_iter):
	recomb_fract_accum = []
	snp_index = 0
	for snp1, snp2 in izip_iter:
		if len(recomb_fract_accum) == 0:
			recomb_fract_accum.append((snp2-snp1)/100)
		else:
			recomb_fract_accum.append((snp2-snp1)/100 + recomb_fract_accum[snp_index-1])
		snp_index += 1
	return(recomb_fract_accum)

def get_recomb_frac(izip_iter):
	recomb_fract_accum = []
	for snp1, snp2 in izip_iter:
		recomb_fract_accum.append((snp2-snp1)/100)
	return(recomb_fract_accum)

if gene == "BRCA1":
	brca_snps = 6834
	brca_snps_left = 2269
	
	# test 15 snps on either side
	# brca_snps = 32
	# brca_snps_left = 16

	snp_coords = [float(x.strip().split('\t')[3]) for x in open(workdir+"../Scripts/simulate-population/brca1_snp_coords.txt").readlines()]
	snp_coords_mb = [float(x.strip().split('\t')[2]) for x in open(workdir+"../Scripts/simulate-population/brca1_snp_coords.txt").readlines()]

	### Get accumulated freqs
	# if (walk_outwards):
	# 	# Go from mutation towars outer area
	# 	recomb_fract_accum_left = get_recomb_frac_accum(izip(reversed(snp_coords[:brca_snps_left-1]), reversed(snp_coords[1:brca_snps_left])))
	# 	recomb_fract_accum_right = get_recomb_frac_accum(izip(snp_coords[brca_snps_left+1:brca_snps-1], snp_coords[brca_snps_left+2:brca_snps]))
	# else:
	# 	# Go from outer area towards mutation
	# 	recomb_fract_accum_left = get_recomb_frac_accum(izip(snp_coords[:brca_snps_left-1], snp_coords[1:brca_snps_left]))
	# 	recomb_fract_accum_right = get_recomb_frac_accum(izip(reversed(snp_coords[brca_snps_left:brca_snps-1]), reversed(snp_coords[brca_snps_left+1:brca_snps])))

	### Get total interval prob
	#brca_prob_left = 0.1695771
	#brca_prob_right = 0.1477736
	# brca_prob_left = recomb_fract_accum_left[-1]
	# brca_prob_right = recomb_fract_accum_right[-1]
	# print brca_prob_left
	# print brca_prob_right

	### Print for inspection
	# print len(snp_coords)
	# print len(recomb_fract_accum_right)
	# print len(recomb_fract_accum_left), len(recomb_fract_accum_right), len(recomb_fract_accum_left)+len(recomb_fract_accum_right)
	# print snp_coords[brca_snps-10:brca_snps]
	# print [i for i in reversed(snp_coords)]
	# print recomb_fract_accum_left
	# print recomb_fract_accum_right
	# exit()

	# Get interval freqs
	if (walk_outwards):
		# Go from mutation towars outer area
		recomb_fract_left = get_recomb_frac(izip(reversed(snp_coords[:brca_snps_left-1]), reversed(snp_coords[1:brca_snps_left])))
		recomb_fract_right = get_recomb_frac(izip(snp_coords[brca_snps_left+1:brca_snps-1], snp_coords[brca_snps_left+2:brca_snps]))
	else:
		# Go from outer area towards mutation
		recomb_fract_left = get_recomb_frac(izip(snp_coords[:brca_snps_left-1], snp_coords[1:brca_snps_left]))
		recomb_fract_right = get_recomb_frac(izip(reversed(snp_coords[brca_snps_left:brca_snps-1]), reversed(snp_coords[brca_snps_left+1:brca_snps])))
	
	# test 15 snps on either side
	# recomb_fract_left = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]
	# recomb_fract_right = recomb_fract_left

	#print len(snp_coords)
	#print len(recomb_fract_accum_right)
	print len(recomb_fract_left), len(recomb_fract_right), len(recomb_fract_left)+len(recomb_fract_right)
	#print snp_coords[brca_snps-10:brca_snps]
	#print [i for i in reversed(snp_coords)]
	# print recomb_fract_left
	# print recomb_fract_right
	# sys.exit()
	
	# Get allele frequencies for population
	with open(workdir+"pop_freqs_minorAllele_brca1_ordered.txt") as f:
		allele_headers = f.readline().split("\t")
		allele_freqs = [float(i) for i in f.readline().split("\t")]
		# print allele_headers[:10]
		# print snp_freqs[:10]
	
	# Get genotype frequencies for population
	freqs_list = []
	with open(workdir+"pop_freqs_genotypes_brca1_ordered.txt") as f:
		header = f.readline()
		for line in f:
			freqs = line.strip().split("\t")
			#print freqs
			#freqs_list.append({0:float(freqs[0]), 1:float(freqs[1]), 2:float(freqs[2])})
			freqs_list.append([float(i) for i in freqs])
		#print freqs_list

	if sim_type == "haplotypes":
		# Get haplotypes
		haplotypes = []
		with open(workdir+"../Scripts/simulate-population/haplotypes_brca1_geno_plink.txt") as f:
			header_tmp = f.readline()
			for line in f:
				haplotype = line.strip().split("\t")
				haplotypes.append(haplotype[1:])
				# haplotypes.append(map(int, line.strip().split("\t")[1:]))
		# Get Mut1HGVS
		with open(workdir+"../Scripts/simulate-population/haplotypes_brca1_pheno_plink.txt") as f:
			header_tmp = f.readline().strip().split("\t")
			index_hgvs = header_tmp.index("Mut1HGVS")
			mut1hgvs_list = [line.strip().split("\t")[index_hgvs] for line in f]
				
	if sim_type == "famHaplotypes":
		# Get famHaplotypes
		haplotypes = []
		with open(workdir+"../Scripts/simulate-population/" +
				  "famHaplotypes_simulationBasis_brca1_geno.txt") as f:
			header_tmp = f.readline()
			for line in f:
				haplotype = line.strip().split("\t")
				haplotypes.append(haplotype[1:]) # OBS!: Beware of columns before genotype data
		# Get Mut1HGVS
		#with open(workdir+"../Scripts/simulate-population/famHaplotypes_brca1_pheno_plink.txt") as f:
		with open(workdir+"../Scripts/simulate-population/" + 
				  "famHaplotypes_simulationBasis_brca1_pheno.txt") as f:
			header_tmp = f.readline().strip().split("\t")
			index_hgvs = header_tmp.index("Mut1HGVS")
			mut1hgvs_list = [line.strip().split("\t")[index_hgvs] for line in f]
	
	if sim_type == "famHaplotypesNoHet":
		# Get famHaplotypes
		haplotypes = []
		with open(workdir+"../Scripts/simulate-population/" +
				  "famHaplotypes_simulationBasis_brca1_geno_no-het.txt") as f:
			header_tmp = f.readline()
			for line in f:
				haplotype = line.strip().split("\t")
				haplotypes.append(haplotype[1:]) # OBS!: Beware of columns before genotype data
		# Get Mut1HGVS
		#with open(workdir+"../Scripts/simulate-population/famHaplotypes_brca1_pheno_plink.txt") as f:
		with open(workdir+"../Scripts/simulate-population/" + 
				  "famHaplotypes_simulationBasis_brca1_pheno.txt") as f:
			header_tmp = f.readline().strip().split("\t")
			index_hgvs = header_tmp.index("Mut1HGVS")
			mut1hgvs_list = [line.strip().split("\t")[index_hgvs] for line in f]


	# Get snp names i.e. first line of file
	with open(workdir+"../Scripts/simulate-population/brca1_snps_names.txt") as f:
		snp_header = f.readline()
		# print len(snp_header.split("\t"))

elif gene == "BRCA2":
	brca_snps = 4597
	brca_snps_left = 2693
	
	snp_coords = [float(x.strip().split('\t')[3]) for x in open(workdir+"brca2_snp_coords.txt").readlines()]
	
	if (walk_outwards):
		# Go from mutation towars outer area
		recomb_fract_accum_left = get_recomb_frac(izip(reversed(snp_coords[:brca_snps_left-1]), reversed(snp_coords[1:brca_snps_left])))
		recomb_fract_accum_right = get_recomb_frac(izip(snp_coords[brca_snps_left:brca_snps-1], snp_coords[brca_snps_left+1:brca_snps]))
	else:
		# Go from outer area towards mutation
		recomb_fract_accum_left = get_recomb_frac(izip(snp_coords[:brca_snps_left-1], snp_coords[1:brca_snps_left]))
		recomb_fract_accum_right = get_recomb_frac(izip(reversed(snp_coords[brca_snps_left:brca_snps-1]), reversed(snp_coords[brca_snps_left+1:brca_snps])))
	
	# brca_prob_left = 0.2903979
	# brca_prob_right = 0.1597097
	brca_prob_left = recomb_fract_accum_left[-1]
	brca_prob_right = recomb_fract_accum_right[-1]
	
	# Get genotype frequencies for population
	freqs_list = []
	with open(workdir+"pop_freqs_genotypes_brca2_ordered.txt") as f:
		header = f.readline()
		for line in f:
			freqs = line.strip().split("\t")
			freqs_list.append([float(i) for i in freqs])
		print freqs_list
	
	with open(workdir+"brca2_snps_names.txt") as f:
		snp_header = f.readline()
	
elif gene == "test":
	# generations = 100
	# sim_num = 1
	# max_samples = 1
	
	brca_snps = 10
	brca_snps_left = 5
	
	recomb_fract_accum_left = [0.1,0.2,0.3,0.4]
	recomb_fract_accum_right = [0.1,0.2,0.3,0.4]
	
	brca_prob_left = 0.4
	brca_prob_right = 0.4
	
	with open(workdir+"brca1_snps_names.txt") as f:
		for line in f:
			snp_header = line
			break
	
else:
	print("Gene unknown! Use either BRCA1 or BRCA2 (or test)")
	sys.exit()


def computeGenerations_old(ancestors):

	descendants = []
	print "Generation:", 0
	print len(ancestors)

	for g in xrange(generations):
		print "Generation:", g+1
		count = 0
		for ancestor in ancestors:
			if star_genealogy:
				num_descendants = 1
			else:
				num_descendants = 1 if random.random()>0.3 else 2 #random.randrange(1,3)
			for i in xrange(num_descendants):
				descendant = ancestor
				# OBS: Might bug in index range!
				if random.random() < brca_prob_left:
					bp = random.randrange(0,brca_snps_left+1)
					#q = bp-1; p = bp+2
					#descendant = descendant[:q] + [-1 for i in descendant[q:p]] + descendant[p:]
					descendant = [0 if random.random()>0.3 else 2 for i in descendant[:bp]] + descendant[bp:]
				if random.random() < brca_prob_right:
					bp = random.randrange(brca_snps_left+1, len(descendant)+1)
					# q = bp-1; p = bp+2
					# descendant = descendant[:q] + [-1 for i in descendant[q:p]] + descendant[p:]
					descendant = descendant[:bp] + [0 if random.random()>0.3 else 2 for i in descendant[bp:]]

				count += 1
				descendants.append(descendant)
				
		print count

		# Sample subset of data
		if len(descendants) > max_samples:
			samples = random.sample(xrange(len(descendants)), max_samples)
			descendants = [descendants[int(i)] for i in samples]
	
		ancestors = descendants
		descendants = []
		print len(ancestors)
	
	return ancestors

def computeGenerations(ancestors, generations):

	descendants = []
	#print "Generation:", 0
	#print len(ancestors)

	for g in xrange(1,generations+1):
		#print "Generation:", g
		count = 0
		for ancestor in ancestors:
			#print ancestor[brca_snps_left-1000:brca_snps_left+1000]
			if star_genealogy:
				num_descendants = 1
			else:
				num_descendants = 1 if random.random()>0.3 else 2 #random.randrange(1,3)
			for i in xrange(num_descendants):
				#descendant = copy.deepcopy(ancestor)
				descendant = ancestor
				
				### Left side ###
				if (walk_outwards):
					bp = brca_snps_left  # breakpoint index (going from mutations towards outmost snp)
				else:
					bp = 0  # breakpoint index (reversed i.e. going from outmost snp towards mutation)
				rand_val = random.random()
				#print brca_prob_left, brca_prob_right, (brca_prob_right+brca_prob_left)
				#print rand_val
				if rand_val < brca_prob_left:
					for val in recomb_fract_accum_left:
						if (walk_outwards):
							bp -= 1
						else:
							bp += 1
						if rand_val < val:
							#print "Left:", rand_val, val, bp
							global count_left
							count_left += 1
							#descendant = [0 if random.random()>0.3 else 2 for i in descendant[:bp]] + descendant[bp:]
							#descendant = [g for i in descendant[:bp]] + descendant[bp:]
							#descendant[bp] = 2
							descendant = descendant[:bp-1] + [2] + descendant[bp:]
							counts[bp] = counts.get(bp, 0) + 1
							count_randvals[-rand_val] = count_randvals.get(-rand_val, 0) + 1
							break
				
				# If we allow breaks to happen on both sides of mutation in the same generation, then draw a new random number
				if allowBreaksBothSidesSameGen:
					rand_val = random.random()
				
				
				### Right side ###
				if (walk_outwards):
					bp = brca_snps_left+1  # breakpoint index (going from mutations towards outmost snp)
				else:
					bp = brca_snps # breakpoint index (reversed i.e. going from outmost snp towards mutation)
				rand_val = rand_val - brca_prob_left
				#print rand_val
				if rand_val < brca_prob_right and rand_val > 0:
					for val in recomb_fract_accum_right:
						# print snp1, snp2, snp2-snp1
						# print len(snp_coords), brca_snps_left
						# print len(snp_coords[(brca_snps_left+1):(brca_snps)]), len(snp_coords[(brca_snps_left):(brca_snps-1)])
						# print len(snp_coords[:brca_snps_left-1]), len(snp_coords[1:brca_snps_left])
						if (walk_outwards):
							bp += 1
						else:
							bp -= 1
						if rand_val < val:
							#print "Right:", rand_val, val, bp
							#sys.exit()
							global count_right
							count_right += 1
							#descendant = descendant[:bp] + [0 if random.random()>0.3 else 2 for i in descendant[bp:]]
							#descendant = descendant[:bp] + [g for i in descendant[bp:]]
							#descendant[bp] = 2
							descendant = descendant[:bp] + [2] + descendant[bp+1:]
							counts[bp] = counts.get(bp, 0) + 1
							count_randvals[rand_val] = count_randvals.get(rand_val, 0) + 1
							break

				count += 1
				descendants.append(descendant)
				
		#print count

		# Sample subset of data
		if len(descendants) > max_samples:
			samples = random.sample(xrange(len(descendants)), max_samples)
			descendants = [descendants[int(i)] for i in samples]
	
		ancestors = descendants
		descendants = []
		#print len(ancestors)
	
	return ancestors

'''
Best case scenario i.e. all breaks are known and identifiable
'''
def computeGenerationsBinomial(ancestors, generations):

	descendants = []
	#print "Generation:", 0
	#print len(ancestors)

	for g in xrange(1,generations+1):
		#print "Generation:", g
		count = 0
		for ancestor in ancestors:
			#print ancestor[brca_snps_left-1000:brca_snps_left+1000]
			if star_genealogy:
				num_descendants = 1
			else:
				num_descendants = 1 if random.random()>0.3 else 2 #random.randrange(1,3)
			for i in xrange(num_descendants):
				#descendant = copy.deepcopy(ancestor)
				descendant = ancestor
				
				### Left side ###
				if (walk_outwards):
					bp = brca_snps_left  # breakpoint index (going from mutations towards outmost snp)
				else:
					bp = 0  # breakpoint index (reversed i.e. going from outmost snp towards mutation)

				val_accum = 0
				for val in recomb_fract_left:
					val_accum += val
					if (walk_outwards):
						bp -= 1
					else:
						bp += 1
					rand_val = random.random()
					if rand_val < val:
						#print "Left:", rand_val, val, bp
						global count_left
						count_left += 1
						descendant = descendant[:bp-1] + [2] + descendant[bp:]
						counts[bp] = counts.get(bp, 0) + 1
						count_randvals[-val_accum] = count_randvals.get(-val_accum, 0) + 1
						break
				
				### Right side ###
				if (walk_outwards):
					bp = brca_snps_left+1 # breakpoint index (going from mutations towards outmost snp)
				else:
					bp = brca_snps # breakpoint index (reversed i.e. going from outmost snp towards mutation)

				val_accum = 0
				for val in recomb_fract_right:
					val_accum += val
					if (walk_outwards):
						bp += 1
					else:
						bp -= 1
					rand_val = random.random()
					if rand_val < val:
						#print "Right:", rand_val, val, bp
						global count_right
						count_right += 1
						descendant = descendant[:bp-1] + [2] + descendant[bp:]
						counts[bp] = counts.get(bp, 0) + 1
						count_randvals[val_accum] = count_randvals.get(val_accum, 0) + 1
						#print bp-1, val_accum
						break

				count += 1
				descendants.append(descendant)
				
		#print count

		# Sample subset of data
		if len(descendants) > max_samples:
			samples = random.sample(xrange(len(descendants)), max_samples)
			descendants = [descendants[int(i)] for i in samples]
	
		ancestors = descendants
		descendants = []
		#print len(ancestors)
	
	return ancestors

def combinedGenoAndPheno(start_generations, stop_generations, step_generations):

    os.system("mkdir -p "+workdir+data_dir)
    os.system("mkdir -p "+workdir+data_dir+"simulated_population/")

    for generations in xrange(start_generations, stop_generations+1, step_generations):
		print "Generation:", generations
		
		file_name = workdir+data_dir+"simulated_population/"+gene+ \
                    "_simulated-"+genealogy+"-"+str(max_samples)+"_samples-"+ \
                    str(sim_num)+"_simulations-generations_"+str(generations)+ \
                    "-seed_"+str(seed)
        
		# Genotype file
		out_geno = file_name+"-geno.txt"
		# Phenotype file
		out_pheno = file_name+"-pheno.txt"
		# Founder file
		out_founder = file_name+"-founder.txt"
		
		if os.path.isfile(out_geno):
			sys.exit("File already exists: " + out_geno)

		#founder = [0 if random.random()>0.3 else 2 for i in xrange(brca_snps)]
		founder = [0 for i in xrange(brca_snps)]

		with open(out_geno, "w") as g:
			with open(out_pheno, "w") as p:
				with open(out_founder, "w") as f:
					# Write file headers
					g.write(snp_header)
					f.write(snp_header)
					p.write("Onc_ID\tFamCode\tCountry\tMut1HGVS\n")
				
					# Write founder genotypes
					f.write("founder\t"+"\t".join(str(i) for i in founder) + "\n")

					sid = 1
					for pop_num in xrange(1,sim_num+1):
						# Being overwrited in each run, needs to be re-defined
						founder = [0 for i in xrange(brca_snps)]
						if star_genealogy:
							ancestors = [founder] * max_samples
						else:
							ancestors = [founder]
                        
						#ancestors = computeGenerations(ancestors, generations)
						ancestors = computeGenerationsBinomial(ancestors, generations)
						
						# number of breaks in generation
						#print(sum(sum(ancestor)/2 for ancestor in ancestors))
						
						if gene == "test":
							print ancestors
							break
						mut_name = "mutSim_"+gene+"_"+str(generations)+"_"+str(pop_num)
						for ancestor in ancestors:
							g.write(str(sid)+"\t"+"\t".join(str(a) for a in ancestor) + "\n")
							p.write(str(sid)+"\tFam"+str(pop_num)+"\tsimuland"+str(pop_num)+"\t"+mut_name+"\n")
							sid += 1

# if simulatePop:
# 	combinedGenoAndPheno(start_generations, stop_generations, step_generations)

'''
Simulated using population frequencies from CIMBA data
'''

def get_haplotypes_subset(haplotypes, hap_num):
	ref = mut1hgvs_list[hap_num]
	#print ref
	ref_indexes = [i for i, e in enumerate(mut1hgvs_list) if e != ref]
	# print len(ref_indexes)
	# print len(haplotypes)
	# print len(mut1hgvs_list)
	haplotypes_subset = [haplotypes[i] for i in ref_indexes]
	# print len(haplotypes_subset)
	return(haplotypes_subset)

def getHaplotype(start, stop, hap_subset=None, hap_num=None):
	haplotype = []

	if sim_type == "haplotypes" or sim_type == "famHaplotypes" or sim_type == "famHaplotypesNoHet":
		if hap_num == None:
			hap_num = random.randint(0, len(hap_subset)-1)
		#haplotype = map(int, haplotypes[hap_num][start:stop])
		haplotype = hap_subset[hap_num][start:stop]

	elif sim_type == "popFreqs": # or genotype freqs
		for freq in freqs_list[start:stop]:
			rand_val = random.random()
			if rand_val < freq[0]:
				haplotype.append(0)
			elif rand_val < (freq[1]+freq[0]):
				haplotype.append(1)
			else:
				haplotype.append(2)

	else: #elif sim_type == "alleleFreqs":
		for freq in allele_freqs[start:stop]:
		# for freq in allele_freqs[2212:2215]:
			rand_val = random.random()
			# print rand_val
			# print freq
			if rand_val < freq:
				# print "minor"
				haplotype.append(2)
			else:
				# print "major"
				haplotype.append(0)
		#sys.exit()
	

	return(haplotype)

def computeGenerationsBinomial_popFreqs(ancestors, generations, hap_subset=None):
	descendants = []
	#print "Generation:", 0
	#print len(ancestors)

	for g in xrange(1,generations+1):
		#print "Generation:", g
		count = 0
		for ancestor in ancestors:
			#print ancestor[brca_snps_left-1000:brca_snps_left+1000]
			if star_genealogy:
				num_descendants = 1
			else:
				#num_descendants = 1 if random.random()>0.3 else 2 #random.randrange(1,3)
				# if len(ancestors) > 2:
				if len(ancestors) > 4:
					num_descendants = np.random.poisson(1.25,1)[0]
				else:
					num_descendants = 0
					while num_descendants == 0:
						num_descendants = np.random.poisson(1.25,1)[0]
					#num_descendants = 1 if random.random()>0.25 else 2 #random.randrange(1,3)
					# num_descendants = 2
					# for i in xrange(num_descendants):
					# 	descendants.append(ancestor)
					# continue
			for i in xrange(num_descendants):
				#descendant = copy.deepcopy(ancestor)
				descendant = ancestor
				# if (num_descendants == 2):
				# 	print descendant
				
				### Left side ###
				if (walk_outwards):
					bp = brca_snps_left  # breakpoint index (going from mutations towards outmost snp)
				else:
					bp = 0  # breakpoint index (reversed i.e. going from outmost snp towards mutation)

				val_accum = 0
				for val in recomb_fract_left:
					val_accum += val
					if (walk_outwards):
						bp -= 1
					else:
						bp += 1
					rand_val = random.random()
					if rand_val < val:
						#print "Left:", rand_val, val, bp
						global count_left
						count_left += 1
						if (sim_type == "knownBreaks"):
							descendant = descendant[:bp-1] + [2] + descendant[bp:]
						else: # popFreqs
							descendant = getHaplotype(0,bp,hap_subset) + descendant[bp:]
							#descendant = getHaplotype(0,bp-1) + [5] + descendant[bp:]
						counts[bp] = counts.get(bp, 0) + 1
						count_randvals[-val_accum] = count_randvals.get(-val_accum, 0) + 1
						break
				
				### Right side ###
				if (walk_outwards):
					bp = brca_snps_left+1 # breakpoint index (going from mutations towards outmost snp)
				else:
					bp = brca_snps # breakpoint index (reversed i.e. going from outmost snp towards mutation)

				val_accum = 0
				for val in recomb_fract_right:
					val_accum += val
					if (walk_outwards):
						bp += 1
					else:
						bp -= 1
					rand_val = random.random()
					if rand_val < val:
						#print "Right:", rand_val, val, bp
						global count_right
						count_right += 1
						if (sim_type == "knownBreaks"):
							descendant = descendant[:bp-1] + [2] + descendant[bp:]
						else: # popFreqs
							descendant = descendant[:bp-1] + getHaplotype(bp-1,len(descendant),hap_subset)
							#descendant = descendant[:bp-1] + [5] + getHaplotype(bp,len(descendant))
						counts[bp] = counts.get(bp, 0) + 1
						count_randvals[val_accum] = count_randvals.get(val_accum, 0) + 1
						#print bp-1, val_accum
						break

				count += 1
				#print len(descendant)
				descendants.append(descendant)
				
		#print count

		# Sample subset of data
		#if not star_genealogy and len(descendants) > (max_samples*0.2):
		if not star_genealogy and len(descendants) > max_samples:
			# sample descendants: random.sample(choose_from, n)
			#samples = sorted(random.sample(xrange(len(descendants)), int(round(min(len(descendants)*0.8,max_samples),0))))
			samples = sorted(random.sample(xrange(len(descendants)), max_samples))
			descendants = [descendants[int(i)] for i in samples]
	
		ancestors = descendants
		descendants = []
		# for a in ancestors:
		# 	print a[0:10]
		# print len(ancestors)
	
	return ancestors

def combinedGenoAndPheno_popFreqs(start_generations, stop_generations, step_generations):

    os.system("mkdir -p "+workdir+data_dir)
    os.system("mkdir -p "+workdir+data_dir+"simulated_population/")

    for generations in xrange(start_generations, stop_generations+1, step_generations):
		print "Generation:", generations
		
		file_name = workdir+data_dir+"simulated_population/"+gene+ \
                    "_simulated-"+genealogy+"-"+str(max_samples)+"_samples-"+ \
                    str(sim_num)+"_simulations-generations_"+str(generations)+ \
                    "-seed_"+str(seed)
        
		# Genotype file
		out_geno = file_name+"-geno.txt"
		# Phenotype file
		out_pheno = file_name+"-pheno.txt"
		# Founder file
		out_founder = file_name+"-founder.txt"
		
		if os.path.isfile(out_geno):
			#sys.exit("File already exists: " + out_geno)
			print "File already exists: " + out_geno
			#continue

		# org_founder = getHaplotype(0, brca_snps)
		# print org_founder[1:10]
		# sys.exit()

		with open(out_geno, "w") as g:
			with open(out_pheno, "w") as p:
				with open(out_founder, "w") as f:
					# Write file headers
					g.write(snp_header)
					f.write(snp_header)
					p.write("Onc_ID\tFamCode\tCountry\tMut1HGVS\n")
				
					# Write founder genotypes
					# if (sim_type == "knownBreaks"):
					# 	founder = [0 for i in xrange(brca_snps)]
					# 	f.write("founder\t"+"\t".join(str(i) for i in founder) + "\n")

					sid = 1 # sample id
					for pop_num in xrange(1,sim_num+1):
						haplotypes_subset = None
						# Being overwrited in each run, needs to be re-defined
						if (sim_type == "knownBreaks"):
							founder = [0 for i in xrange(brca_snps)]
							f.write("founder\t"+"\t".join(str(i) for i in founder) + "\n")
						elif (sim_type == "haplotypes" or 
							  sim_type == "famHaplotypes" or 
							  sim_type == "famHaplotypesNoHet"):
							hap_num = random.randint(0, len(haplotypes)-1)
							founder = getHaplotype(0, brca_snps, haplotypes, hap_num)
							f.write("mutSim_"+gene+"_"+str(generations)+"_"+str(pop_num)+"\t"+
									"\t".join(str(i) for i in founder) + "\n")
							haplotypes_subset = get_haplotypes_subset(haplotypes, hap_num)
						else: # popFreqs
							founder = getHaplotype(0, brca_snps)
							f.write("mutSim_"+gene+"_"+str(generations)+"_"+str(pop_num)+"\t"+
									"\t".join(str(i) for i in founder) + "\n")
						
						if star_genealogy:
							ancestors = [founder] * max_samples
						else:
							ancestors = [founder]
							
						# print len(ancestors)
						ancestors = computeGenerationsBinomial_popFreqs(ancestors, generations,
							 											haplotypes_subset)
						# print len(ancestors)
						
						# number of breaks in generation
						#print(sum(sum(ancestor)/2 for ancestor in ancestors))

						mut_name = "mutSim_"+gene+"_"+str(generations)+"_"+str(pop_num)
						for ancestor in ancestors:
							g.write(str(sid)+"\t"+"\t".join(str(a) for a in ancestor) + "\n")
							p.write(str(sid)+"\tFam"+str(pop_num)+"\tsimuland"+str(pop_num)+"\t"+mut_name+"\n")
							sid += 1


def simulate_generations(generations):
	print "Generation:", generations
	
	file_name = workdir+data_dir+"simulated_population/"+gene+ \
                "_simulated-"+genealogy+"-"+str(max_samples)+"_samples-"+ \
                str(sim_num)+"_simulations-generations_"+str(generations)+ \
                "-seed_"+str(seed)
    
	# Genotype file
	out_geno = file_name+"-geno.txt"
	# Phenotype file
	out_pheno = file_name+"-pheno.txt"
	# Founder file
	out_founder = file_name+"-founder.txt"
	
	if os.path.isfile(out_geno):
		#sys.exit("File already exists: " + out_geno)
		print "File already exists: " + out_geno
		#return()

	# org_founder = getHaplotype(0, brca_snps)
	# print org_founder[1:10]
	# sys.exit()

	with open(out_geno, "w") as g:
		with open(out_pheno, "w") as p:
			with open(out_founder, "w") as f:
				# Write file headers
				g.write(snp_header)
				f.write(snp_header)
				p.write("Onc_ID\tFamCode\tCountry\tMut1HGVS\n")
			
				# Write founder genotypes
				# if (sim_type == "knownBreaks"):
				# 	founder = [0 for i in xrange(brca_snps)]
				# 	f.write("founder\t"+"\t".join(str(i) for i in founder) + "\n")

				sid = 1 # sample id
				for pop_num in xrange(1,sim_num+1):
					haplotypes_subset = None
					# Being overwrited in each run, needs to be re-defined
					if (sim_type == "knownBreaks"):
						founder = [0 for i in xrange(brca_snps)]
						f.write("founder\t"+"\t".join(str(i) for i in founder) + "\n")
					elif (sim_type == "haplotypes" or 
						  sim_type == "famHaplotypes" or 
						  sim_type == "famHaplotypesNoHet"):
						hap_num = random.randint(0, len(haplotypes)-1)
						founder = getHaplotype(0, brca_snps, haplotypes, hap_num)
						f.write("mutSim_"+gene+"_"+str(generations)+"_"+str(pop_num)+"\t"+
								"\t".join(str(i) for i in founder) + "\n")
						haplotypes_subset = get_haplotypes_subset(haplotypes, hap_num)
					else: # popFreqs
						founder = getHaplotype(0, brca_snps)
						f.write("mutSim_"+gene+"_"+str(generations)+"_"+str(pop_num)+"\t"+
								"\t".join(str(i) for i in founder) + "\n")
					
					if star_genealogy:
						ancestors = [founder] * max_samples
					else:
						ancestors = [founder]
						
					# print len(ancestors)
					ancestors = computeGenerationsBinomial_popFreqs(ancestors, generations,
						 											haplotypes_subset)
					# print len(ancestors)
					
					# number of breaks in generation
					#print(sum(sum(ancestor)/2 for ancestor in ancestors))

					mut_name = "mutSim_"+gene+"_"+str(generations)+"_"+str(pop_num)
					for ancestor in ancestors:
						g.write(str(sid)+"\t"+"\t".join(str(a) for a in ancestor) + "\n")
						p.write(str(sid)+"\tFam"+str(pop_num)+"\tsimuland"+str(pop_num)+"\t"+mut_name+"\n")
						sid += 1

def combinedGenoAndPheno_popFreqs_parallel(start_generations, stop_generations, step_generations):

	os.system("mkdir -p "+workdir+data_dir)
	os.system("mkdir -p "+workdir+data_dir+"simulated_population/")

	gens = [generations for generations in xrange(start_generations, stop_generations+1, step_generations)]
	#gens = [generations for generations in xrange(stop_generations, start_generations-1, -step_generations)]
	print gens
	
	if len(gens) == 1:
		simulate_generations(gens[0])
	else:
		pool = Pool()
		res = pool.map(simulate_generations, gens)
		pool.close()
		pool.join()
		#print res

if simulatePop:
	#combinedGenoAndPheno_popFreqs(start_generations, stop_generations, step_generations)
	combinedGenoAndPheno_popFreqs_parallel(start_generations, stop_generations, step_generations)
	


if option != "none":
	os.system("Rscript --vanilla "+rscript_dir+"prepare_simulated_data.R "+ \
          #str(custom_filename_string)+" "+ \
          str(max_samples)+" "+ \
          str(sim_num)+" "+ \
          str(gene)+" "+ \
          str(start_generations)+" "+ \
          str(stop_generations)+" "+ \
          str(step_generations)+" "+ \
          str(seed)+" "+ \
          str(genealogy)+" "+ \
          str(workdir+data_dir)+" "+ \
		  str(option)+" "+ \
		  str(ancestral_method)+" "+ \
		  str(cs_correction))


run_stats = False
if platform.system() == "Darwin" and run_stats:
	import matplotlib.pyplot as plt
	# Statistics
	print counts
	print len(counts)
	#x = [i for i in counts]

	x = snp_coords
	# x = snp_coords_mb   # physical coords
	print len(x)
	# Inject positions with no breaks
	for i in xrange(len(snp_coords)):
		if counts.get(i, 0) == 0:
			counts[i] = 0

	y = [counts[i] for i in counts]
	print len(y)
	print sum(y)
	# plt.plot(x,y, "bo")
	# plt.show()

	print count_left
	print count_right

	#print count_randvals
	print len(count_randvals)

	#plt.plot([i for i in xrange(len(recomb_fract_accum_left))],recomb_fract_accum_left, "bo")
	#plt.show()

	#plt.plot([i for i in xrange(len(recomb_fract_accum_right)+brca_snps_left, brca_snps_left,-1)], recomb_fract_accum_right, "bo")
	#plt.show()

	# brca1: snp index in intervals from -0.17 to 0.15 with steps 0.01. 
	snp_index = [146,270,417,522,594,648,733,843,929,969,1050,1369,1593,1688,1909,2044,2266,2483,3047,4085,5129,5571,5621,5711,5831,5966,6095,6295,6590,6711,6805,6830]
	# test 15 snps on either side
	#snp_index = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]

	# brca1
	x = [-0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.10, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15]
	# test 15 snps on either side
	#x = [-0.14, -0.13, -0.12, -0.11, -0.10, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15]
	y = [0 for i in x]
	print len(snp_index), len(x), len(y)
	for k in counts:
		for i, c in enumerate(snp_index):
			if k <= c:
				y[i] = y[i] + counts[k]
				break
	plt.plot(x,y,'bo')
	# plt.show()

	y = [0 for i in x]
	for k in count_randvals:
		for i, v in enumerate(x):
			if v <= 0:
				#if round(k,8) < v:
				if k < v:
					y[i] = y[i] + count_randvals[k]
					break
			else:
				#if round(k,2) <= v:
				if k <= v:
					y[i] = y[i] + count_randvals[k]
					break

	plt.plot(x,y,'ro')
	plt.show()