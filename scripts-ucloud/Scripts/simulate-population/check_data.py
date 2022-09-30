#infile = "/Users/lars/Documents/PhD/haplotype-project/Scripts/simulate-population/BRCA2_simulated_population_200_samples.txt"
#infile = "/Users/lars/Documents/PhD/haplotype-project/Scripts/simulate-population/BRCA2_simulated_200-samples_10-muts_20-generations_geno.txt"
infile = "/Users/lars/Documents/PhD/haplotype-project/Scripts/simulate-population/BRCA2_simulated_200-samples_10-muts_20-generations_pheno.txt"

with open(infile) as f:
	for line in f:
		l = line.split("\t")
		print(len(l))