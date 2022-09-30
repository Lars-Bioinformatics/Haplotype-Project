library(readxl)
library(dplyr)
library(forcats)
library(fastmatch)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(gtable)
library(data.table)

removeSNPsInGene = FALSE

# setwd("~/Documents/PhD/haplotype-project")

# Test cases
### BRCA1 ###
#mut = "c.5266dupC"; gene = "BRCA1"; pop_ref = brca2_geno_merged; pop_factor = 1.2
#mut = "c.1016dupA"; gene = "BRCA1"
#mut = "c.1016delA"; gene = "BRCA1"
#mut = "c.-200-?_80+?del"; gene = "BRCA1"
#mut = "c.1687C>T"; gene = "BRCA1"
#mut = "c.3481_3491del11"; gene = "BRCA1"


### BRCA2 ###
### mut = "c.7617+1G>A"; gene = "BRCA2"


# Test if pattern occurs due to country haplotype or if real
# validateNoCountryPattern()

# mut = "c.1687C>T"; gene = "BRCA1"
# haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "GERMANY")
# haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "AUSTRIA")
# haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "ITALY")
# haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "USA")
# haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "SWEDEN")
#
# mut = "c.181T>G"; gene = "BRCA1"
# haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "GERMANY")
# haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "ITALY")
#
# mut = "c.7617+1G>A"; gene = "BRCA2"
# haplotype_mutation(mut, pop_factor = pop_factor, country = "DENMARK")
# haplotype_mutation(mut, pop_factor = pop_factor, country = "SPAIN")

# Compute haplotype for specific mutation - no country consensus
### haplotype_mutation(mut, pop_ref, pop_factor = pop_factor)

# Compute haplotypes for all country consensuses
### haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, allCountryConsensus = T)

# Run single mutation, on one or more countries.
runMutation <- function(mut, gene, country, allCountryConsensus) {
    gene <<- gene
    # Compute haplotypes
    haplotype_mutation(mut, country = country, allCountryConsensus = allCountryConsensus)
}

# Run single mutation, on one or more countries.
runMutations <- function(muts, gene, allCountryConsensus) {
    gene <<- gene
    # Compute haplotypes
    for (mut in muts){
        haplotype_mutation(mut, allCountryConsensus = allCountryConsensus)
    }
}


# Run BRCA1 mutations with more than 100 patients
# Runs with all consensus countries over 0.1 of total samples for that mutation
runBRCA1_100samples <- function(){

    gene = "BRCA1"
    # muts = c("c.181T>G", "c.190T>C", "c.5266dupC", "c.1016dupA", "c.1687C>T", "c.1961delA",
    #          "c.3481_3491del11", "c.4065_4068delTCAA", "c.68_69delAG")
    muts = c("c.-200-?_80+?del", "c.1687C>T", "c.181T>G", "c.211A>G", "c.2475delC",
             "c.2681_2682delAA", "c.3319G>T", "c.3331_3334delCAAG", "c.3481_3491del11",
             "c.3700_3704del5", "c.3756_3759delGTCT", "c.4035delA", "c.4065_4068delTCAA",
             "c.4186-?_4357+?dup", "c.427G>T", "c.4327C>T", "c.5123C>A", "c.5266dupC",
             "c.5333-36_5406+400del510", "c.5503C>T", "c.68_69delAG")

    # Compute haplotype breaks for each mutation
    for (mut in muts){
        haplotype_mutation(mut, allCountryConsensus = T)
    }
}

# Run BRCA2 mutations with more than 100 patients
# Runs with all consensus countries over 0.1 of total samples for that mutation
runBRCA2_100samples <- function(){

    gene = "BRCA2"
    muts = c("c.1310_1313delAAGA", "c.1813dupA", "c.2808_2811delACAA", "c.3847_3848delGT",
             "c.4478_4481delAAAG", "c.5645C>A", "c.5682C>G", "c.5722_5723delCT", "c.5946delT",
             "c.6275_6276delTT", "c.7069_7070delCT", "c.7617+1G>A", "c.771_775del5",
             "c.7934delG", "c.8537_8538delAG")

    for (mut in muts){
        haplotype_mutation(mut, allCountryConsensus = T)
    }
}

## ---- haplotype_mutation
# Main function to run everyting
haplotype_mutation <- function(mut, pop_ref = NULL, pop_factor = 1, country = NULL, 
                               distinct = F, allCountryConsensus = F, prepare = T, 
                               group = NULL, cutoff=0.5, ancestral_method = "branchBoundIndep",
                               min_samples=2){
    "
	Haplotype samples with a given brca mutation
	"

    start.time = Sys.time()

    if (prepare){
        #prepare_dataframe(mut, gene, distinct)
        prepare_dataframe(mut, gene)
    }

    if (!is.null(pop_ref)){
        matched2 <- select_rare_snps(matched, pop_ref, pop_factor)
    }

    if (ancestral_method == "mostFreqBase"){
        if (!is.null(country)){
            # If haplotype should be based on samples from specific country
            #haplotypes <<- findConsensus(filter(matched, SNP %in% subset(brca_data, Country == country)$Onc_ID))
            print("LOL")
            country <<- country
            haplotypes <<- findConsensus_plink(filter(matched_plink, SNP %in% subset(brca_data, Country == country)$Onc_ID), cutoff)
        } else if (!is.null(group)) {
            print("Kage")
            # brca_data <<- brca_data.filtered
            # haplotypes <<- findConsensus(filter(matched, SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID))
            haplotypes <<- findConsensus_plink(filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID), cutoff)
        } else {
            print("Kage2")
            # If the haplotype should be defined across all samples
            #haplotypes <<- findConsensus(matched)
            haplotypes <<- findConsensus_plink(matched_plink, cutoff)
        }
        #haplotypes <<- haplotypes
    
        # Find SNPs breaking the haplotype
        # breaks <<- findHaplotypeBreaks(matched, haplotypes)
        breaks <<- findHaplotypeBreaks_plink(matched_plink, haplotypes)
    } else if (ancestral_method == "simulatedFounder"){
    #} else if (grepl("simulated", gene)){
        print("FOUNDER")
        haplotypes = read.table(paste0(filename_sim, "-founder.txt"), header = T)
        haplotypes <- as.integer(haplotypes[haplotypes$SNP == mut,-1])
        breaks <<- findHaplotypeBreaks_plink(matched_plink, haplotypes)
        #haplotypes = rep(0, ncol(matched_plink)-1)
        #haplotypes <<- findConsensus_plink(filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID), cutoff)
    } else {
        if (!is.null(country)){
            matched_plink_subset = matched_plink[matched_plink$SNP %in% subset(brca_data, Country == country)$Onc_ID,]
        } else {
            matched_plink_subset = matched_plink[matched_plink$SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID,]
        }
        if (ancestral_method == "branchBoundIndep"){
            print("BranchBoundIndep")
            res = branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = cutoff, indep_sides = TRUE, min_samples = min_samples)
        } else if (ancestral_method == "branchBound"){ #branchBound
            print("BranchBound")
            res = branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = cutoff, indep_sides = FALSE, min_samples = min_samples)
        } else if (ancestral_method == "branchBoundIndep_alleleFreq"){
            res = branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = cutoff, indep_sides = T, 
                                             min_samples = min_samples, alleleFreqs = alleleFreqs)
        } else if (ancestral_method == "branchBound_alleleFreq"){
            res = branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = cutoff, indep_sides = F, 
                                             min_samples = min_samples, alleleFreqs = alleleFreqs)
        } else {
            stop("Ancestral method does not exist")
        }
        haplotypes <<- res[[1]]
        breaks <<- res[[2]]
    }

    # If no samples breaking haplotype, then stop analysis of current mutation
    if (length(breaks) == 0){
        print(paste("No haplotype breaks found in families for mut:", mut))
        return()
    }

    # Map found breaks (SNPs) to their genomic position
    # chr_pos <<- mapSNPs(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    chr_pos <<- mapSNPs_plink(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    #chr_pos <<- mapSNPs_plink_new(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)

    # Save information required to plot and find nearest breaks
    # save(chr_pos, file = paste0("cache/chr_pos-", brca_name, "-", mut_name, ".RData"))
    # Load information, if already computed (SHOULD BE MOVED BEFORE COMPUTATIONS)
    # load("cache/chr_pos-BRCA1-c.1016delA.RData")

    if (!prepare) {
        print(Sys.time()-start.time)
        return(chr_pos)
    }

    # Find the nearest break on each side of the genome
    # and compute the mean length of the breaks for each country
    # haplotypeDistance <- nearestBreaksStatistics(chr_pos, mut, country)

    # Visualize the breaks around the BRCA gene
    if (!is.null(group)){
        p <- plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
        p2 <- plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
    } else {
        p <- plot_breakpoints_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut)
        p2 <- plot_nearest_brca_break_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut)
    }

    # ggplot(haplotypeDistance, aes(x=sample_id, y=positive_pos)) + geom_bar()

    return(list(p,p2))

    # Computing all haplotypes for different ancestral references (countries), higher than 10% of total data size
    if (allCountryConsensus == T){
        quantile = 0.1
        countries <- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        for (country in countries){
            haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = country)
        }
    }


}

prepare_dataframe <- function(mut, gene, distinct=F, gene_info_obj=NULL){
    # Read pheno and geno data
    ## readData()

    # Get BRCA genes coordinate information
    getBRCAinfo(gene, mut, gene_info_obj)

    # fixed mut name for saving files on windows
    mut_name <- gsub(">", "_", mut)
    mut_name <<- gsub("\\?", "_", mut_name)

    # Computes general statistics of dataset
    # computeStatistics()

    # Prints for checking output
    print(dim(brca_data)) # 8x67
    #print(dim(haplotypes_brca)) # 15679x12468
    print(dim(chr_coords_all)) # 12467x3
    print(mut) # c.1016delA

    # Filter out the coordinate of SNPs needed i.e. on current chromosome (13 or 17)
    chr_coords <<- filter(chr_coords_all, Chr_numeric == brca_chr) %>% 
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(brca_geno_plink[,-1])) %>% 
        distinct(position_b37, .keep_all = T) %>% 
        # sort snps according to location
        arrange(position_b37)
    dim(chr_coords)
    dim(chr_coords_all)

    # Extract samples with given mutation and remove SNPs not related to actual BRCA gene
    #matched <<- extractSamples(brca_data, haplotypes_brca, chr_coords, distinct) # Should be commented out
    matched_plink <<- extractSamples(brca_data, brca_geno_plink, chr_coords, distinct)
    matched <<- matched_plink # Here for now to ensure working code. should be fixed!
    #matched <<- extractSamples(brca_data, haplotypes_brca, chr_coords, F)
    # matched <- matched[, c(1, na.omit(fmatch(chr_coords$SNP, colnames(matched))))]
    
    # Save reference for easy compatibility with plotting
    matched_single <<- matched
    matched_plink_single <<- matched_plink
    brca_data.single <<- brca_data
    
    return(matched_plink)
}

###############################################################################
### Read Data
###############################################################################
## Read input data into global variables
readData <- function(){
    # Read phenotype data
    brca1_pheno <- read.csv("input/139_ouh_june_2017/B1_Onco_phenotype_distribution_311215.csv")
    brca2_pheno <- read.csv("input/139_ouh_june_2017/B2_Onco_phenotype_distribution_180816.csv")
    brca1_pheno_hisp <- read.csv("input/139_ouh_june_2017/B1_Hispanic_OncoArray_phenotypes_230817.csv")
    brca2_pheno_hisp <- read.csv("input/139_ouh_june_2017/B2_Hispanic_OncoArray_phenotypes_230817.csv")

    brca1_pheno_hisp <<- brca1_pheno_hisp %>% mutate(Country="Hispanic")
    brca2_pheno_hisp <<- brca2_pheno_hisp %>% mutate(Country="Hispanic")

    # Note there's 20 more columns in brca1_pheno compared to brca1_pheno_hisp.
    dim(brca2_pheno)
    dim(brca2_pheno_hisp)
    colnames(brca2_pheno[,1:47])
    colnames(brca2_pheno_hisp)
    # I think these are not important and therefore removed in the merge step.
    brca1_pheno_merged <<- rbind(brca1_pheno[, 1:47], brca1_pheno_hisp) %>% arrange(FamCode)

    # Suppress warnings of NAs introduced - non-important data
    brca2_pheno_merged <<- suppressWarnings(rbind(brca2_pheno[, 1:47], brca2_pheno_hisp) %>% arrange(FamCode))

    #brca2_pheno_merged <<- rbind(brca2_pheno[, 32:32], brca2_pheno_hisp[, 32:32])

    # Hisp is broken, thus they are excluded in the analyses
    # brca1_pheno_merged <<- brca1_pheno
    # brca2_pheno_merged <<- brca2_pheno

    # Read genotype data
    # if(!file.exists("cache/brca1_geno_merged.RData") || !file.exists("cache/brca2_geno_merged.RData")){
    #     # Read and merge genotype data - then save R object
    #     brca1_geno <- read.table("input/139_ouh_june_2017/139_ouh_brca1_onco_geno.txt", header = T, stringsAsFactors = T)
    #     brca2_geno <- read.table("input/139_ouh_june_2017/139_ouh_brca2_onco_geno.txt", header = T, stringsAsFactors = T)
    #     brca1_geno_hisp <- read.table("input/139_ouh_june_2017/139_ouh_female_hisp_brca1_onco_geno.txt", header = T, stringsAsFactors = T)
    #     brca2_geno_hisp <- read.table("input/139_ouh_june_2017/139_ouh_female_hisp_brca2_onco_geno.txt", header = T, stringsAsFactors = T)
    #
    #     brca1_geno_merged <- rbind(brca1_geno, brca1_geno_hisp) %>% arrange(SNP)
    #     brca2_geno_merged <- rbind(brca2_geno, brca2_geno_hisp) %>% arrange(SNP)
    #
    #     # brca1_geno_merged <- brca1_geno
    #     # brca2_geno_merged <- brca2_geno
    #
    #     save(brca1_geno_merged, file = "cache/brca1_geno_merged.RData")
    #     save(brca2_geno_merged, file = "cache/brca2_geno_merged.RData")
    # }
    # load("cache/brca1_geno_merged.RData", .GlobalEnv)
    # load("cache/brca2_geno_merged.RData", .GlobalEnv)
    # load("cache/brca1_geno_merged_noFactors.RData", .GlobalEnv)
    # load("cache/brca2_geno_merged_noFactors.RData", .GlobalEnv)

    # Read genotype data in plink format
    if (!file.exists("cache/brca1_geno_plink.RData") || !file.exists("cache/brca2_geno_plink.RData")){
        brca1_geno_eu_plink <- read.table2("input/139_ouh_june_2017/139_ouh_brca1_onco_geno_plink_format.txt", header = T)
        brca2_geno_eu_plink <- read.table2("input/139_ouh_june_2017/139_ouh_brca2_onco_geno_plink_format.txt", header = T)
        brca1_geno_hisp_plink <- read.table2("input/139_ouh_june_2017/139_ouh_female_hisp_brca1_onco_geno_plink_format.txt", header = T)
        brca2_geno_hisp_plink <- read.table2("input/139_ouh_june_2017/139_ouh_female_hisp_brca2_onco_geno_plink_format.txt", header = T)

        brca1_geno_plink <- rbind(brca1_geno_eu_plink, brca1_geno_hisp_plink)
        brca2_geno_plink <- rbind(brca2_geno_eu_plink, brca2_geno_hisp_plink)

        save(brca1_geno_plink, file = "cache/brca1_geno_plink.RData")
        save(brca2_geno_plink, file = "cache/brca2_geno_plink.RData")
    }
    load("cache/brca1_geno_plink.RData", envir = .GlobalEnv)
    load("cache/brca2_geno_plink.RData", envir = .GlobalEnv)

    # Read SNP genomic coordinates
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")

    # Read Genetic coordinates (centimorgan)
    morgan.brca1 <<- read.table2("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    morgan.brca2 <<- read.table2("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map

    # BRCA mutation coordinates
    brca1_mutation_coords <<- read.table("input/139_ouh_june_2017/BRCA1_mut_position.txt", header = T, sep = "\t")
    brca2_mutation_coords <<- read.table("input/139_ouh_june_2017/BRCA2_mut_position.txt", header = T, sep = "\t")
    
    # genotype frequencies in population
    popFreqs_brca1 <<- read.table2("cache/pop_freqs_genotypes_brca1_ordered.txt", header = T)
    popFreqs_brca2 <<- read.table2("cache/pop_freqs_genotypes_brca2_ordered.txt", header = T)
    
    # allele frequencies in population
    alleleFreqs_brca1 <<- as.numeric(read.table2("cache/pop_freqs_minorAllele_brca1_ordered.txt", header = T))
    alleleFreqs_brca2 <<- as.numeric(read.table2("cache/pop_freqs_minorAllele_brca2_ordered.txt", header = T))
}

readSimulatedData <- function(){
    
    sim_simulations <<- 20; sim_generations <<- 100; sim_samples <<- 100
    brca1_simulated_pheno_merged <<- read.table(file = paste0("classify-simulated-population-age/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_42/simulated_population/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_100-seed_42-pheno.txt"), header = T)
    brca1_simulated_geno_plink <<- read.table(file = paste0("classify-simulated-population-age/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_42/simulated_population/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_100-seed_42-geno.txt"), header = T)
    
    brca1_simulated_pheno_merged <<- read.table(file = paste0("classify-simulated-population-age/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_42/simulated_population/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_150-seed_42-pheno.txt"), header = T)
    brca1_simulated_geno_plink <<- read.table(file = paste0("classify-simulated-population-age/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_42/simulated_population/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_150-seed_42-geno.txt"), header = T)
    
    brca1_simulated_pheno_merged <<- read.table(file = paste0("classify-simulated-population-age/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_100_100_step10-seed_42/simulated_population/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_100-seed_42-pheno.txt"), header = T)
    brca1_simulated_geno_plink <<- read.table(file = paste0("classify-simulated-population-age/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_100_100_step10-seed_42/simulated_population/training-BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_100-seed_42-geno.txt"), header = T)
    
    # brca2_pheno_merged <<- read.table(file = "Scripts/simulate-population/BRCA2_simulated_population_200-samples_population-1_pheno.txt", header = T)
    # brca2_geno_plink <<- read.table(file = "Scripts/simulate-population/BRCA2_simulated_population_200-samples_population-1_geno.txt", header = T)
    # brca2_pheno_merged <<- read.table(file = "Scripts/simulate-population/BRCA2_simulated_population_200-samples_population-2_pheno.txt", header = T)
    # brca2_geno_plink <<- read.table(file = "Scripts/simulate-population/BRCA2_simulated_population_200-samples_population-2_geno.txt", header = T)
    #brca2_simulated_pheno_merged <<- read.table(file = "Scripts/simulate-population/BRCA2_simulated_starGenealogy_200-samples_2-muts_20-generations_pheno.txt", header = T)
    #brca2_simulated_geno_plink <<- read.table(file = "Scripts/simulate-population/BRCA2_simulated_starGenealogy_200-samples_2-muts_20-generations_geno.txt", header = T)

    # simulations = 2
    # simulations = 5; generations = 100
    #sim_simulations <<- 5; sim_generations <<- 1000; sim_samples <<- 200
    #sim_simulations <<- 10; sim_generations <<- 20; sim_samples <<- 200
    #simulations = 10; generations = 10
    #simulations = 10; generations = 20
    # simulations = 20
    #brca1_simulated_pheno_merged <<- read.table(file = paste0("Scripts/simulate-population/training3_BRCA1_simulated_starGenealogy_",sim_samples,"-samples_",sim_simulations,"-muts_",sim_generations,"-generations_pheno.txt"), header = T)
    #brca1_simulated_geno_plink <<- read.table(file = paste0("Scripts/simulate-population/training3_BRCA1_simulated_starGenealogy_",sim_samples,"-samples_",sim_simulations,"-muts_",sim_generations,"-generations_geno.txt"), header = T)
}

read.table2 <- function(file, header){
    return(as.data.frame(fread(file = file, header = header)))
}

getBRCAinfo <- function(brca_name, mut, gene_info_obj=NULL){
    # Gene information for custom gene in shiny app
    if (!is.null(gene_info_obj)){
        brca_name <<- brca_name
        row_with_info <- which(gene_info_obj[1,] == mut)
        brca_start <<- gene_info_obj[row_with_info, 4]
        brca_stop <<- gene_info_obj[row_with_info, 5]
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_chr <<- gene_info_obj[row_with_info, 3]
        
        morgan.coords <<- read.table2(paste0("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr",brca_chr,"_combined_b37.txt"), header = T) # phase3 map
        brca_cM_middle <<- morgan.coords[which.min(abs(morgan.coords$position-brca_middle)),3]
        
        haplotypes_brca <<- open_projects[[brca_name]]$project_geno
        brca_geno_plink <<- open_projects[[brca_name]]$project_geno
        brca_data <<- filter(open_projects[[brca_name]]$project_pheno, Mut1HGVS %in% mut)
        
        chr_coords_all <<- project_coords
    }
    
    # BRCA gene information
    if (brca_name == "BRCA1"){
        brca_name <<- "BRCA1"
        #brca_start <<- 41197695; brca_stop <<- 41276113
        brca_start <<- brca1_mutation_coords[which(brca1_mutation_coords$Mut1HGVS == mut), 2]
        brca_stop <<- brca1_mutation_coords[which(brca1_mutation_coords$Mut1HGVS == mut), 3]
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_cM_middle <<- morgan.brca1[which.min(abs(morgan.brca1$position-brca_middle)),3]
        brca_chr <<- 17
        
        # Assume c.427G>T and c.3319G>T has same mutation - to verify clustering
        # brca1_pheno_merged[brca1_pheno_merged$Mut1HGVS == "c.427G>T", "Mut1HGVS"] = "c.3319G>T"
        alleleFreqs <<- alleleFreqs_brca1
        
        # Legacy names
        haplotypes_brca <<- brca1_geno_plink #brca1_geno_merged
        brca_geno_plink <<- brca1_geno_plink
        if (is.null(mut)){
            brca_data <<- brca1_pheno_merged
        } else {
            brca_data <<- filter(brca1_pheno_merged, Mut1HGVS %in% mut)
        }
    } else if (brca_name == "BRCA2") { # BRCA2
        brca_name <<- "BRCA2"
        #brca_start <<- 32890598; brca_stop <<- 32972907
        brca_start <<- brca2_mutation_coords[which(brca2_mutation_coords$Mut1HGVS == mut), 2]
        brca_stop <<- brca2_mutation_coords[which(brca2_mutation_coords$Mut1HGVS == mut), 3]
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_cM_middle <<- morgan.brca2[which.min(abs(morgan.brca2$position-brca_middle)),3]
        brca_chr <<- 13
        
        alleleFreqs <<- alleleFreqs_brca2
        
        # Legacy names
        haplotypes_brca <<- brca2_geno_plink #brca2_geno_merged
        brca_geno_plink <<- brca2_geno_plink
        if (is.null(mut)){
            brca_data <<- brca2_pheno_merged
        } else {
            brca_data <<- filter(brca2_pheno_merged, Mut1HGVS %in% mut)
        }
    } else if (brca_name == "BRCA1_simulated") {
        brca_name <<- "BRCA1"
        brca_start <<- 41197695; brca_stop <<- 41276113
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_cM_middle <<- morgan.brca1[which.min(abs(morgan.brca1$position-brca_middle)),3]
        brca_chr <<- 17
        #brca_start <<- brca_middle; brca_stop <<- brca_middle; 
        
        alleleFreqs <<- alleleFreqs_brca1
        
        # Legacy names
        haplotypes_brca <<- brca1_simulated_geno_plink #brca1_geno_merged
        brca_geno_plink <<- brca1_simulated_geno_plink
        if (is.null(mut)){
            brca_data <<- brca1_simulated_pheno_merged
        } else {
            brca_data <<- filter(brca1_simulated_pheno_merged, Mut1HGVS %in% mut)
        }
        #brca_data <<- filter(brca1_simulated_pheno_merged, Mut1HGVS %in% mut)
    } else if (brca_name == "BRCA2_simulated") { # BRCA2
        brca_name <<- "BRCA2"
        brca_start <<- 32890598; brca_stop <<- 32972907
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_cM_middle <<- morgan.brca2[which.min(abs(morgan.brca2$position-brca_middle)),3]
        brca_chr <<- 13

        alleleFreqs <<- alleleFreqs_brca2
        
        # Legacy names
        haplotypes_brca <<- brca2_simulated_geno_plink #brca2_geno_merged
        brca_geno_plink <<- brca2_simulated_geno_plink
        if (is.null(mut)){
            brca_data <<- brca2_simulated_pheno_merged
        } else {
            brca_data <<- filter(brca2_simulated_pheno_merged, Mut1HGVS %in% mut)
        }
    } else if (brca_name == "BRCA2_toySample") { # BRCA2
        brca_name <<- "BRCA2"
        #brca_start <<- 32890598; brca_stop <<- 32972907
        brca_start <<- brca2_mutation_coords[which(brca2_mutation_coords$Mut1HGVS == mut), 2]
        brca_stop <<- brca2_mutation_coords[which(brca2_mutation_coords$Mut1HGVS == mut), 3]
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_cM_middle <<- morgan.brca2[which.min(abs(morgan.brca2$position-brca_middle)),3]
        brca_chr <<- 13
        
        alleleFreqs <<- alleleFreqs_brca2
        
        # Legacy names
        haplotypes_brca <<- brca2_toy_geno
        brca_geno_plink <<- brca2_toy_geno
        if (is.null(mut)){
            brca_data <<- brca2_toy_pheno
        } else {
            brca_data <<- filter(brca2_toy_pheno, Mut1HGVS %in% mut)
        }
    }
    brca.mut.cM <<- brca_cM_middle
}

computeStatistics <- function(){
    ## Compute general info about dataset ##
    # Number of families with specific mutation
    #MutFamCount <- brca1_pheno_merged %>% distinct(Mut1HGVS, FamCode) %>% group_by(Mut1HGVS) %>% summarise(n=n())
    MutFamCount <<- brca_data %>% distinct(Mut1HGVS, FamCode) %>% group_by(Mut1HGVS) %>% summarise(n=n())

    # Number of countries with specific mutation
    #MutCountryCount <- brca1_pheno_merged %>% distinct(Mut1HGVS, Onc_ID, Country) %>% group_by(Country) %>% summarise(n=n())
    MutCountryCount <<- brca_data %>% distinct(Mut1HGVS, Onc_ID, Country) %>% count(Country)

    # Number of patients with specific mutation
    #MutPatientCount <- count(brca1_pheno_merged, Mut1HGVS) %>% filter(n>1) %>% summarise(sum(n))
    #MutPatientCount <- count(brca2_pheno_merged, Mut1HGVS) %>% filter(n>1) %>% summarise(sum(n))
    MutPatientCount <<- count(brca_data, Mut1HGVS)

    # Count the number of family members
    #FamMemberCount <- brca1_pheno_merged %>% distinct(Onc_ID, FamCode) %>% count(FamCode) %>% count(FamMembers=n)
    #FamMemberCount <- brca2_pheno_merged %>% distinct(Onc_ID, FamCode) %>% count(FamCode) %>% count(FamMembers=n)
    FamMemberCount <<- brca_data %>% distinct(Onc_ID, FamCode) %>% count(FamCode) %>% count(FamMembers=n)

    # Count number of families with at least two members
    MultiFamCount <<- brca_data %>% count(FamCode) %>% filter(n>1)

    # Number of different mutations with at least 2 samples
    MutCount_brca1 <- brca1_pheno_merged %>% count(Mut1HGVS) %>% filter(n>2) %>% nrow()
    MutCount_brca2 <- brca2_pheno_merged %>% count(Mut1HGVS) %>% filter(n>2) %>% nrow()

    hispStatistics <- function(gene = "BRCA1"){
        if (gene=="BRCA1"){
            brca_pheno = brca1_pheno
            brca_pheno_hisp = brca1_pheno_hisp
        } else {
            brca_pheno = brca2_pheno
            brca_pheno_hisp = brca2_pheno_hisp
        }

        brca_shared <- filter(brca_pheno, Mut1HGVS %in% brca_pheno_hisp$Mut1HGVS)
        total = brca_shared %>% distinct(Mut1HGVS, Onc_ID, Country) %>% count(Mut1HGVS) %>% rename(samples_total_eu=n)

        total$spanish = 0
        total$non_spanish = 0

        spanish = brca_shared %>% distinct(Mut1HGVS, Onc_ID, Country) %>% filter(Country=="SPAIN") %>% count(Mut1HGVS)
        index = match(spanish$Mut1HGVS, total$Mut1HGVS)
        total[index, "spanish"] = spanish$n

        non_spanish = brca_shared %>% distinct(Mut1HGVS, Onc_ID, Country) %>% filter(Country != "SPAIN") %>% count(Mut1HGVS)
        index = match(non_spanish$Mut1HGVS, total$Mut1HGVS)
        total[index, "non_spanish"] = non_spanish$n

        total$samples_total_hispanic <- filter(brca_pheno_hisp, Mut1HGVS %in% brca_pheno$Mut1HGVS) %>% distinct(Mut1HGVS, Onc_ID, Country) %>% count(Mut1HGVS) %>% .$n

        total = total %>% select("Mut1HGVS", "spanish", "non_spanish", "samples_total_eu", "samples_total_hispanic")

        write.table(total, file = paste0("Hispanic-statistics-", gene, ".txt"), quote = F, row.names = F, sep = "\t")
    }
    #hispStatistics(gene = "BRCA1")
    #hispStatistics(gene = "BRCA2")
}

computeDensityPlots <- function(option=1, centiMorgan=F, one_fam_member=F, plotBreaks=F, fam=F){
    # Compute common ancestral haplotype for all samples - Nul fordeling
    # Paramteter: Option: select value between 1-6
    # plus:
    # one_fam_member: keep only one individual from each family
    # centiMorgan: Use centiMorgan (genetic) distance instead of physical distance

    # Pick SNPs
    if (option %in% c(1,3,5)){
        SNPs_gene = "BRCA1"
    } else {
        SNPs_gene = "BRCA2"
    }
    getBRCAinfo(SNPs_gene, NULL)
    prepare_dataframe("all", SNPs_gene)

    brca_data <<- rbind(brca1_pheno_merged, brca2_pheno_merged)
    if (one_fam_member) brca_data <<- distinct(brca_data, FamCode, .keep_all = T)

    if (fam){
        brca1_geno_w_brca1_snps_plink <- famHaplotypes_brca1_geno_plink[na.omit(fmatch(brca_data$Onc_ID, famHaplotypes_brca1_geno_plink$SNP)),c(2, na.omit(fmatch(chr_coords$SNP, colnames(famHaplotypes_brca1_geno_plink))))]
        brca2_geno_w_brca1_snps_plink <- famHaplotypes_brca2_geno_plink[na.omit(fmatch(brca_data$Onc_ID, famHaplotypes_brca2_geno_plink$SNP)),c(2, na.omit(fmatch(chr_coords$SNP, colnames(famHaplotypes_brca2_geno_plink))))]
    } else {
        brca1_geno_w_brca1_snps_plink <- brca1_geno_plink[na.omit(fmatch(brca_data$Onc_ID, brca1_geno_plink$SNP)),c(1, na.omit(fmatch(chr_coords$SNP, colnames(brca1_geno_plink))))]
        brca2_geno_w_brca1_snps_plink <- brca2_geno_plink[na.omit(fmatch(brca_data$Onc_ID, brca2_geno_plink$SNP)),c(1, na.omit(fmatch(chr_coords$SNP, colnames(brca2_geno_plink))))]
    }

    ## Pick one group of samples
    if (option %in% c(1,2)){
        matched_plink <- brca1_geno_w_brca1_snps_plink # BRCA1
    } else if (option %in% c(3,4)){
        matched_plink <- brca2_geno_w_brca1_snps_plink # BRCA2
    } else if (option %in% c(5,6)){
        matched_plink <- rbind(brca1_geno_w_brca1_snps_plink, brca2_geno_w_brca1_snps_plink) # BRCA1 and BRCA2
    }

    haplotypes <<- findConsensus_plink(matched_plink)
    breaks <<- findHaplotypeBreaks_plink(matched_plink, haplotypes)
    chr_pos <<- mapNearestSNPs_plink(matched_plink, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    if (centiMorgan){
        d <- density(chr_pos$haplo_length_cM)
        cm.text = " - centiMorgan"
    } else {
        d <- density(chr_pos$haplo_length)
        cm.text = ""
    }
    if (option == 1) caption = "BRCA1_SNPs-BRCA1_samples"
    if (option == 3) caption = "BRCA1_SNPs-BRCA2_samples"
    if (option == 5) caption = "BRCA1_SNPs-BRCA1-2_samples"
    if (option == 2) caption = "BRCA2_SNPs-BRCA1_samples"
    if (option == 4) caption = "BRCA2_SNPs-BRCA2_samples"
    if (option == 6) caption = "BRCA2_SNPs-BRCA1-2_samples"
    if (one_fam_member){
        png(paste0("Density plot - ", caption, " - one fam member", cm.text, ".png"))
    }else if (fam){
        png(paste0("Density plot - ", caption, " - families", cm.text, ".png"))
    } else {
        png(paste0("Density plot - ", caption, cm.text, ".png"))
    }
    plot(d, main = paste("Haplotype length -", caption))
    polygon(d, col="red")
    dev.off()
    print("Done with density plot")

    # Number of samples breaking in brca gene
    print(nrow(filter(chr_pos, positive_pos < brca_stop, negative_pos > brca_start)))

    if (plotBreaks){
        chr_pos_minmax = chr_pos %>% mutate(position_pos_adjusted = positive_pos-brca_middle, position_neg_adjusted = negative_pos-brca_middle)
        plots = plot_nearest_brca_break_country(chr_pos_minmax, brca_chr, brca_name, brca_start, brca_stop, mut, fixObject=F)
        num_samples = length(unique(chr_pos_minmax$sample_id))

        ofm.text = if (one_fam_member) "-one_fam_member" else ""
        fm.text = if (fam) "-families" else ""
        ggsave(plot = plots[[1]], filename = paste0(caption, "-", "nearest_break-country-", num_samples, "_samples", omf.text, fm.text, ".png"), width = 16, height = 9)
        ggsave(plot = plots[[2]], filename = paste0(caption, "-", "nearest_break-positive_order-", num_samples, "_samples", omf.text, fm.text, ".png"), width = 16, height = 9)
        ggsave(plot = plots[[3]], filename = paste0(caption, "-", "nearest_break-negative_order-", num_samples, "_samples", omf.text, fm.text, ".png"), width = 16, height = 9)
        ggsave(plot = plots[[4]], filename = paste0(caption, "-", "nearest_break-haplot_length_order-", num_samples, "_samples", omf.text, fm.text, ".png"), width = 16, height = 9)
    }
}

select_rare_snps <- function(matched, brca_geno_reference, pop_factor = 1){
    ### generate population frequencies ###

    pop_freqs = c()
    for (pos in 2:ncol(brca_geno_reference)){
        # Combine a SNP column to one string and split string into one-character vector
        seq <- strsplit(paste0(brca_geno_reference[, pos], collapse = ""), split = "")[[1]]
        # Counts the number of occurences of each element
        counts <- table(seq)
        if (names(counts[1]) == "-"){
            counts <- counts[2:length(counts)]
        }
        freq = c()
        for (i in 1:length(counts)){
            freq = c(freq, counts[i]/sum(counts))
        }
        #print(paste("pos:", pos))
        #print(counts)
        #print(freq)
        pop_freqs = c(pop_freqs, list(freq))
    }
    head(pop_freqs)

    # List of all snps
    snp_names <- colnames(brca_geno_reference)[2:ncol(brca_geno_reference)]

    ### Get snps
    rare_snp_names = snp_names[na.omit(fmatch(chr_coords$SNP, snp_names))]
    pop_freqs_rare_snps = pop_freqs[na.omit(fmatch(chr_coords$SNP, snp_names))]
    length(rare_snp_names)
    length(pop_freqs)
    length(pop_freqs_rare_snps)
    length(chr_coords$SNP)
    length(matched)

    # BRCA rare snps
    rare_snps <<- c()
    #pos = 11
    for (pos in 2:ncol(matched)){
        # Combine a SNP column to one string and split string into one-character vector
        seq <- strsplit(paste0(matched[, pos], collapse = ""), split = "")[[1]]
        # Counts the number of occurences of each element
        counts <- table(seq)
        if (names(counts[1]) == "-"){
            counts <- counts[2:length(counts)]
        }
        fam_freq = c()
        for (i in 1:length(counts)){
            fam_freq = c(fam_freq, counts[i]/sum(counts))
        }
        fam_freq

        fam_max_name <<- names(which.max(fam_freq))
        fam_max <<- fam_freq[fam_max_name]
        pop_freq <<- pop_freqs_rare_snps[[pos-1]][fam_max_name]

        # Keep snp if no population frequency or if population
        # frequency * given factor is smaller than family frequency
        if (is.na(pop_freq) || fam_max > pop_freq*pop_factor){
        #if (!is.na(pop_freq) && fam_max > (pop_freq)){
            rare_snps = c(rare_snps, rare_snp_names[pos-1])
        }

    }
    print("Done")
    length(rare_snps)

    matched = matched %>% select(SNP, rare_snps)

    return(matched)
}

extractSamples <- function(brca_data, brca_geno_plink, chr_coords, distinct = F){
    "
    Extract relevant samples from table of all genotypes
    "

    # Keep only one member of a family and remove the others, if chosen
    if (distinct){
        brca_data <- brca_data[!duplicated(brca_data$FamCode), ]
        # Equivalent to above
        # brca_data <- subset(brca_data, !duplicated(FamCode))
        # brca_data <- distinct(brca_data, FamCode, .keep_all = T)
    }

    # Extract samples (with given brca mutation) from the big haplotypes dataset
    # Faster to use fmatch than match
    index <- fmatch(brca_data$Onc_ID, brca_geno_plink$SNP)
    matched <- brca_geno_plink[index, ]
    # as above, but slower
    # matched <- filter(brca_geno, SNP %in% brca_data$Onc_ID)
    #print(dim(matched))

	# Only keep SNPs (columns) that are related to current chromosome
    snp_col <- which(names(matched)=="SNP")
    matched <- matched[, c(snp_col, na.omit(fmatch(chr_coords$SNP, colnames(matched))))]
	#print(dim(matched))

    return(matched)
}

###############################################################################
### Find consensus (ancestor) haplotype
###############################################################################
findConsensus_plink <- function(matched_plink, cutoff=0.5){
    #col = matched_plink[,100]
    ancestral_locus <- function(col){
        #print(col)
        #print(cutoff)
        #if (length(table(col)) == 3 && table(col)[3] > table(col)[1]) print(table(col))
        #col=col2
        #col = table(col)/length(col)
        col = table(col)
        if (length(col) == 1){
            return(as.integer(names(col[1])))
        }

        zeros = max(0, col[names(col)==0])
        ones = max(0, col[names(col)==1])
        twos = max(0, col[names(col)==2])
        freqs <- c("0"=(2*zeros+ones)/(2*zeros+2*ones+2*twos), "2"=(2*twos+ones)/(2*zeros+2*ones+2*twos))
        max_allele = which.max(freqs)
        if (freqs[max_allele] > cutoff){
            return(as.integer(names(max_allele)))
        } else(
            return(1)
        )

        # col = col[names(col) != 1]
        # col_max = which(col == max(col)) # important not to use which.max
        # if (length(col_max) == 1 && col_max){
        #     # if (col[col_max] < 0.3){
        #     #     return(1)
        #     # }
        #     return(as.integer(names(col_max[1])))
        # } else { # length == 2
        #     return(1)
        # }
    }

    #m2 <- matched_plink[1:10,2:10]
    # col2=filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID) %>% .$chr13_32975487_C_T
    # ancestral_locus(col2)
    cutoff <<- max(cutoff, 0.5) # should be at least 0.5 (and at most 1)
    #consensus <- apply(matched_plink[,2:ncol(matched_plink)], 2, ancestral_locus)
    consensus <- apply(matched_plink[, -1], 2, ancestral_locus)
    return(consensus)
}

findConsensus <- function(matched, cutoff = 0.5){
    "
    Find haplotype consensus for each SNP
    Haplotype consensus is defined as the most
    frequent letter (A,C,G,T,I,D) in all samples
    "

    start.time = Sys.time()
    cutoff = max(cutoff, 0.5) # should be at least 0.5 (and at most 1)
    print("Find haplotype consensus for each SNP")
    haplotypes <<- c("SNP", rep("", ncol(matched)-1))
    for (pos in 2:ncol(matched)){
        if (pos %% 1000 == 0) print(paste("Process:", pos, "of", ncol(matched), "SNPs"))

        # Combine a SNP column to one string
        #seq <- paste0(matched[, pos], collapse = "")
        # Split string into one-character vector
        #seq2 <- strsplit(seq, split = "")[[1]]
        # Counts the number of occurences of each element
        #counts <- table(seq2)
        # print(counts)

        # Works for famHaplotypes - might not work for phased data!
        counts <- table(unlist(strsplit(as.character(matched[, pos]), split = "")))

        #counts <- table(unlist(lapply(matched[, pos], as.character)))
        #counts <- table(unlist(strsplit(matched[, pos], split = "")))

        # counts <- table(c(substr(matched[, pos], 1, 1),substr(matched[, pos], 2, 2)))

        # if only '-' character in string, return no consensus
        if (length(counts) == 1 && names(counts[1]) == "-"){
            haplotypes[pos] <- "-"
            next
        }

        # Skip the '-' character - as it's unknown value
        if (names(counts[1]) == "-") {
            counts = counts[2:length(counts)]
        }

        # If only one element (except -), return that as consensus
        if (length(counts) == 1){
            haplotypes[pos] <- names(counts[1])
            next
        }

        # Check if more than one element has the highest frequency
        freqs = counts/sum(counts)
        max_allele = which.max(freqs)
        print(freqs)
        print(max_allele)
        print(cutoff)
        if (freqs[max_allele] > cutoff){
            haplotypes[pos] <- names(freqs[max_allele])
        } else {
            haplotypes[pos] <- paste0(names(counts), collapse = "")
        }

        # # List of candidates for most frequent haplotype
        # candidates = names(counts[start_index])
        # # Highest number of occurences
        # best = counts[start_index]
        #
        # for (i in 1:length(counts)){
        #     if (counts[i] == best){
        #         candidates = c(candidates, names(counts[i]))
        #     }
        #     if (counts[i] > best){
        #         candidates = names(counts[i])
        #         best = counts[i]
        #     }
        #     # print(candidates)
        # }
        #
        # # If only one candidate, return as consensus
        # if (length(candidates) == 1){
        #     new = candidates
        #     #names(new) <- names(matched[pos])
        #     haplotypes[pos] <- new
        # # When multiple possible haplotypes, combine and return
        # } else {
        #     new = paste0(names(counts[1]), names(counts[2]))
        #     #print(paste0(names(counts[1]), names(counts[2])))
        #     #new = "F"
        #     #names(new) <- names(matched[pos])
        #     haplotypes[pos] <- new
        # }
    }
    names(haplotypes) <- colnames(matched)
    print(paste0("Time used to find consensus: ", Sys.time() - start.time))
    # Return list representing consensus haplotype
    return(haplotypes)
}

ancestral_locus <- function(col, cutoff = 0.5, alt_allele_freq = NULL){
    col = table(col)
    if (length(col) == 1){
        return(as.integer(names(col[1])))
    }
    
    zeros = max(0, col[names(col)==0])
    ones = max(0, col[names(col)==1])
    twos = max(0, col[names(col)==2])
    total = (2*zeros+2*ones+2*twos)
    
    if (is.null(alt_allele_freq) || alt_allele_freq == 0){
        freqs <- c("0"=(2*zeros+ones)/total, "2"=(2*twos+ones)/total)
        max_allele = which.max(freqs)
        if (freqs[max_allele] > cutoff){
            return(as.integer(names(max_allele)))
        } else(
            return(1)
        )
    } else {
        p.val = 0.4
        total = (zeros+2*ones+twos)
        # Check if observed is higher than expected
        if (alt_allele_freq*total<(twos+ones)){ # && alt_allele_freq>0.05){# && twos > 0){
            p.val = binom.test(x = (twos+ones), n = total, p = alt_allele_freq, conf.level = 0.95)$p.value
            #p.val = binom.test(x = (2*twos+ones), n = total, p = as.numeric(alt_allele_freq), conf.level = 0.95)$p.value
            # p.val = binom.test(x = (2*zeros+ones), n = total, p = 1-alt_allele_freq, conf.level = 0.95)$p.value
            # binom.test(x = 139, n = total, p = 1-alt_allele_freq, conf.level = 0.95)$conf.int[1]*total
            # binom.test(x = 140, n = total, p = 1-alt_allele_freq, conf.level = 0.95)$conf.int[2]*total
        }
        if (twos>=zeros || p.val < 0.05){
            return(2)
        # If close to expected, then set 1 since high uncertainty of correct ancestral locus
        } else if (p.val>0.5){ # Might need tuning!
            return(1)
        } else {
            return(0)
        }
    }
}

branch_and_bound_ancestral <- function(matched_plink, matched_plink_subset, cutoff=0.5, indep_sides = T, min_samples = 5, alleleFreqs = NULL, verbose = F){
    snps.left = sum(chr_coords$position_b37 < brca_middle)
    
    matched_all = matched_plink[,-1]
    breaks = matrix(FALSE, nrow=nrow(matched_all), ncol=ncol(matched_all))
    rownames(breaks) = rownames(matched_all)
    colnames(breaks) = colnames(matched_all)
    
    rows_kept = rep(F, nrow(matched_all))
    rows_kept[fmatch(matched_plink_subset$SNP, matched_plink$SNP)] = T
    
    haplotypes = rep(1, ncol(matched_all))
    names(haplotypes) = names(matched_all)
    
    # Least required number of samples needed to define ancestral haplotype
    # if (nrow(matched_plink_subset) < 3){
    #     min_samples = 1
    # } 
    # else {
    #     min_samples = 2
    # }
    
    if (min_samples < 1){
        min_samples = max(nrow(matched_plink_subset) * min_samples, 5)
    }
    
    if (nrow(matched_plink_subset) < min_samples){
        message("HELLO")
        #return(-1)
    }
    
    if (indep_sides) {
        # Right side
        rows_right = rows_kept
        t=Sys.time()
        r = sum(rows_right)
        for (i in (snps.left+1):length(chr_coords$position_b37)){
            break_right = i
            #print(matched_all[rows_right,i])
            
            #col = c(1, 0, 2, 1, 2, 1)
            #col = matched_kept_right[,i]
            if (is.null(alleleFreqs)){
                haplotypes[i] = ancestral_locus(matched_all[rows_right,i], cutoff)
            } else {
                haplotypes[i] = ancestral_locus(matched_all[rows_right,i], cutoff, alt_allele_freq = alleleFreqs[i])
            }
            # print(paste("index:", i))
            # print(matched_all[rows_left,i])
            # print(paste("hap2:", haplotypes[i]))
            # print(paste("founder:", sim_founder[i]))
            
            # adj = 1L-min(matched_kept_right[,i])
            # haplotypes[i] = as.integer(names(which.max(setNames(tabulate(matched_kept_right[,i]+adj), min(matched_kept_right[,i]):max(matched_kept_right[,i])))))
            # 
            # adj = 1L-min(col)
            # freqs = setNames(tabulate(col+adj), min(col):max(col)) / (nrow(matched_kept_right)*2)
            # freqs
            # max_allele = which.max(freqs)
            # if (freqs[max_allele] > cutoff){
            #     haplotypes[i] = as.integer(names(max_allele))
            # } else {
            #     haplotypes[i] = 1
            # }
            breaks[,i] = abs(matched_all[,i]-haplotypes[i])==2
            #print(breaks[rows_right,i])
            #matched_kept_right = matched_kept_right[which(breaks[row.names(breaks) %in% row.names(matched_kept_right),i]==F),]
            bak = rows_right
            rows_right[breaks[,i]] = F
            #dim(matched_kept_right)
            #i = i + 1
            
            if (verbose == "right" && sum(rows_right) != r){
                print(paste("index:", i))
                #print(matched_all[rows_left,i])
                print(matched_all[bak,i])
                print(paste("hap2:", haplotypes[i]))
                print(paste("founder:", sim_founder[i]))
                r=sum(rows_right)
            }
            
            # Special case if two samples left
            if (sum(rows_right) == 2 && sum(matched_all[rows_right,i]) == 2 && matched_all[rows_right,i][1] != 1){
                breaks[rows_right, i] = T
                #print(i-snps.left)
                print(i)
                #break_right = i
                break
            }
            
            #if (sum(rows_right) < 3){
            if (sum(rows_right) < min_samples){
            #if (sum(rows_right) < 5 || sum(rows_right) < nrow(matched_plink_subset)*0.25){
                #breaks[row.names(breaks) %in% row.names(matched_kept_right),i]=T
                breaks[rows_right, i] = T
                #print(i-snps.left)
                print(i)
                #break_right = i
                break
            }
        }
        Sys.time()-t
        
        # Left side
        #matched_kept_left = matched_kept
        rows_left = rows_kept
        t=Sys.time()
        l = sum(rows_left)
        for (i in snps.left:1){
        #for (i in snps.left:2578){
            break_left = i
            #haplotypes[i] = ancestral_locus(matched_kept_left[,i], cutoff)
            if (is.null(alleleFreqs)){
                haplotypes[i] = ancestral_locus(matched_all[rows_left,i], cutoff)
            } else {
                haplotypes[i] = ancestral_locus(matched_all[rows_left,i], cutoff, alt_allele_freq = alleleFreqs[i])
            }
            
            #haplotypes[i] = as.integer(names(which.max(setNames(tabulate(matched_kept_left[,i]+1), sort(unique(matched_kept_left[,i]))))))
            # adj = 1L-min(matched_kept_left[,i])
            # freqs = setNames(tabulate(matched_kept_left[,i]+adj), min(matched_kept_left[,i]):max(matched_kept_left[,i])) / nrow(matched_kept_left)
            # max_allele = which.max(freqs)
            # if (freqs[max_allele] > cutoff){
            #     haplotypes[i] = as.integer(names(max_allele))
            # } else {
            #     haplotypes[i] = 1
            # }
            breaks[,i] = abs(matched_all[,i]-haplotypes[i])==2
            #matched_kept_left = matched_kept_left[which(breaks[row.names(breaks) %in% row.names(matched_kept_left),i]==F),]
            bak = rows_left
            rows_left[breaks[,i]] = F
            
            if (verbose == "left" && sum(rows_left) != l){
                print(paste("index:", i))
                #print(matched_all[rows_left,i])
                print(matched_all[bak,i])
                print(paste("hap2:", haplotypes[i]))
                print(paste("founder:", sim_founder[i]))
                l=sum(rows_left)
            }
            
            # Special case if two samples left
            if (sum(rows_left) == 2 && sum(matched_all[rows_left,i]) == 2 && matched_all[rows_left,i][1] != 1){
                breaks[rows_right, i] = T
                print(i)
                #break_left = i
                break
            }
            
            #if (sum(rows_left) < 3){
            if (sum(rows_left) < min_samples){
            #if (sum(rows_left) < 5 || sum(rows_left) < nrow(matched_plink_subset)*0.25){
                breaks[rows_left,i]=T
                print(i)
                #break_left = i
                break
            }
        }
        Sys.time()-t
    } else {
        for (i in 0:max(snps.left,length(chr_coords$position_b37)-snps.left)){
            if (i < snps.left){
                index = snps.left-i
                # Left side
                if (is.null(alleleFreqs)){
                    haplotypes[index] = ancestral_locus(matched_all[rows_kept,index], cutoff)
                } else {
                    haplotypes[index] = ancestral_locus(matched_all[rows_kept,index], cutoff, alt_allele_freq = alleleFreqs[index])
                }
                # Breaks
                breaks[,index] = abs(matched_all[,index]-haplotypes[index])==2
                # Eliminate samples with sub-optimal ancestral haplotype
                rows_kept[breaks[,index]] = F
            }
            
            if (i < length(chr_coords$position_b37)-snps.left){
                index = snps.left+1+i
                # Right side
                if (is.null(alleleFreqs)){
                    haplotypes[index] = ancestral_locus(matched_all[rows_kept,index], cutoff)
                } else {
                    haplotypes[index] = ancestral_locus(matched_all[rows_kept,index], cutoff, alt_allele_freq = alleleFreqs[index])
                }
                # Breaks
                breaks[,index] = abs(matched_all[,index]-haplotypes[index])==2
                # Eliminate samples with sub-optimal ancestral haplotype
                rows_kept[breaks[,index]] = F
            }
            if (sum(rows_kept) < min_samples){
                #if sum(rows_kept) < 5 || sum(rows_kept) < nrow(matched_subset)*0.25){
                breaks[rows_kept,i]=T
                print(i)
                break_left = i
                break_right = i
                break
            }
        }
    }
    
    #breaks <<- breaks
    #haplotypes <<- haplotypes
    
    return(list(haplotypes, breaks, break_left, break_right))
}

###############################################################################
### Find breaks (or recombinations) 
###############################################################################
findHaplotypeBreaks_plink <- function(matched_plink, consensus){
    #c <- c(1:9)
    #breaks <- t(apply(matched_plink[,-1], 1, function(row) abs(row-consensus)==2 | row == -1))
    breaks <- t(apply(matched_plink[,-1], 1, function(row) abs(row-consensus)==2))
    return(breaks)
}

findHaplotypeBreaks_single_snp <- function(matched, consensus){
    get_breaks <- function(row){
        breaks = !(row == consensus[2:length(consensus)])
        breaks[nchar(haplotypes[2:length(haplotypes)]) > 1] = F
        return(breaks)
    }
    #c <- c(1:9)
    #breaks <- t(apply(matched[,2:ncol(matched)], 1, get_breaks))
    breaks <- t(apply(matched[, -1], 1, get_breaks))
    return(breaks)
}

## ---- findHaplotypeBreaks
findHaplotypeBreaks <- function(matched, haplotypes){
    start.time = Sys.time()
    breaks <- list()
    for (pos in 2:ncol(matched)){

        if (pos %% 500 == 0) print(paste("Process:", pos, "of", ncol(matched), "SNPs"))

        haplotype = haplotypes[pos]
        # print(haplotype)
        if (nchar(haplotype) > 1){
            next
        }

        split_rows <- strsplit(as.character(matched[, pos]), split = "")
        #split_rows <- strsplit(matched[, pos], split = "")

        for (i in 1:nrow(matched)){
            # i = 1
            #row <- 1
            #row <- strsplit(as.character(matched[i,pos]), split = "")[[1]]
            #row <- c("A", "A")
            # print(row)
            row <- split_rows[[i]]
            if (row[1] == "-" || row[2] == "-"){
                next
            }
            if (row[1] != haplotype && row[2] != haplotype){
                #print(row)
                #print(haplotype)
                index = as.character(matched[i,1])
                breaks[[index]] = c(breaks[[index]], colnames(matched[pos]))
            }
        }
    }
    print(paste("Num samples:", length(breaks)))
    print(paste("Num breaks:", length(breaks[[1]])))
    print(paste0("Time used to find breaks: ", Sys.time() - start.time))
    return(breaks)
}

###############################################################################
### Map (or link) the found breaks (snps) to their genomic position
###############################################################################
### Map found SNPs to their genomic position ###
mapSNPs_plink <- function(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name){
    #mapSNPs <- function(){
    chr_pos = NULL
    # For each sample (patient) with given mutation
    for (i in 1:nrow(breaks)){
        # If no breaks for sample, then skip and continue
        if (length(which(breaks[i,]==T)) == 0) next
        # Extract rows with SNPs in brca1 chromosome (17) and where SNP is homozygote for that sample
        #chr_pos2 = filter(chr_coords, chr_coords$Chr_numeric == brca_chr, chr_coords$SNP %in% breaks[[i]])
        chr_pos2 = filter(chr_coords, SNP %in% names(which(breaks[i,]==T)))
        # print(dim(chr_pos2))
        # Removes SNPs in the brca gene
        if (removeSNPsInGene == TRUE){
            chr_pos2 = filter(chr_pos2, as.numeric(as.character(position_b37)) < brca_start | as.numeric(as.character(position_b37)) > brca_stop)
        }
        # Add sample_id column
        chr_pos2["sample_id"] = matched_plink[i,1]
        # Add position_b37_adjusted column, where brca1 = 0
        chr_pos2 = mutate(chr_pos2, position_b37_adjusted = as.numeric(as.character(position_b37)) - brca_middle)
        # Add cM_adjusted column, where brca1 = 0
        chr_pos2 = mutate(chr_pos2, cM_adjusted = cM - brca_cM_middle)
        # Combine all samples into big data frame for plotting (likely not most efficient way)
        chr_pos = bind_rows(chr_pos, chr_pos2)
    }
    # Add columns Country and FamCode to data frame
    index <- fmatch(chr_pos$sample_id, brca_data$Onc_ID)
    # chr_pos[c("Country", "FamCode")] <- brca_data[index, ] %>% select(Country, FamCode)
    chr_pos[c("Country", "FamCode", "cluster_groups")] <- brca_data[index, ] %>% select(Country, FamCode, cluster_groups)

    # Add fake-break at max and min dist for samples with no breaks, so they show in plots
    for (sample in unique(matched_plink$SNP)[!(unique(matched_plink$SNP) %in% unique(filter(chr_pos, position_b37_adjusted < 0)$sample_id))]){
        if (is.na(sample)) next
        #df_max = data.frame("SNP"="None", Chr_numeric=brca_chr, position_b37=max(chr_pos$position_b37), cM=max(chr_pos$cM), sample_id=sample, position_b37_adjusted=max(chr_pos$position_b37_adjusted), cM_adjusted=max(chr_pos$cM_adjusted))
        df_min = data.frame("SNP"="None", Chr_numeric=brca_chr, position_b37=min(chr_pos$position_b37), cM=min(chr_pos$cM), sample_id=sample, position_b37_adjusted=min(chr_pos$position_b37_adjusted), cM_adjusted=min(chr_pos$cM_adjusted))
        #df <- rbind(df_max,df_min)
        df = cbind(df_min, brca_data %>% filter(Onc_ID == sample) %>% select(Country, FamCode, cluster_groups))
        chr_pos = rbind(chr_pos, df)
    }

    for (sample in unique(matched_plink$SNP)[!(unique(matched_plink$SNP) %in% unique(filter(chr_pos, position_b37_adjusted > 0)$sample_id))]){
        if (is.na(sample)) next
        df_max = data.frame("SNP"="None", Chr_numeric=brca_chr, position_b37=max(chr_pos$position_b37), cM=max(chr_pos$cM), sample_id=sample, position_b37_adjusted=max(chr_pos$position_b37_adjusted), cM_adjusted=max(chr_pos$cM_adjusted))
        #df_min = data.frame("SNP"="None", Chr_numeric=brca_chr, position_b37=min(chr_pos$position_b37), cM=min(chr_pos$cM), sample_id=sample, position_b37_adjusted=min(chr_pos$position_b37_adjusted), cM_adjusted=min(chr_pos$cM_adjusted))
        #df <- rbind(df_max,df_min)
        df = cbind(df_max, brca_data %>% filter(Onc_ID == sample) %>% select(Country, FamCode, cluster_groups))
        chr_pos = rbind(chr_pos, df)
    }

    # Print some data information - here counting samples by Country
    print(chr_pos %>% distinct(sample_id, FamCode, Country) %>% count(Country))
    print(dim(chr_pos))
    # unique(chr_pos$sample_id)
    # print(chr_pos)

    return(chr_pos)
}

### Map found SNPs to their genomic position ###
mapSNPs_plink_new <- function(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name){
    #mapSNPs <- function(){
    # chr_pos = NULL
    # # For each sample (patient) with given mutation
    # for (i in 1:nrow(breaks)){
    #     # If no breaks for sample, then skip and continue
    #     if (length(names(which(breaks[i,]==T))) == 0) next
    #     # Extract rows with SNPs in brca1 chromosome (17) and where SNP is homozygote for that sample
    #     #chr_pos2 = filter(chr_coords, chr_coords$Chr_numeric == brca_chr, chr_coords$SNP %in% breaks[[i]])
    #     chr_pos2 = filter(chr_coords, SNP %in% names(which(breaks[i,]==T)))
    #     # print(dim(chr_pos2))
    #     # Removes SNPs in the brca1 gene
    #     # chr_pos2 = filter(chr_pos2, as.numeric(as.character(position_b37)) < brca_start | as.numeric(as.character(position_b37)) > brca_stop)
    #     # Add sample_id column
    #     chr_pos2["sample_id"] = matched_plink[i,1]
    #     # Add position_b37_adjusted column, where brca1 = 0
    #     chr_pos2 = mutate(chr_pos2, position_b37_adjusted = as.numeric(as.character(position_b37)) - brca_middle)
    #     # Combine all samples into big data frame for plotting (likely not most efficient way)
    #     chr_pos = bind_rows(chr_pos, chr_pos2)
    # }

    breaks2 = mutate(as.data.frame(breaks), SNP=matched_plink$SNP)
    chr_pos <- apply(breaks2, 1, function(row){
        if (length(which(row==T)) == 0) return(row)
        chr_pos2 = filter(chr_coords, SNP %in% names(row[which(row==T)])) %>% mutate(sample_id = row[length(row)])
        #chr_pos2["sample_id"] = row[length(row)]
    })
    chr_pos = bind_rows(chr_pos)
    chr_pos = mutate(as.data.frame(chr_pos), position_b37_adjusted = position_b37 - brca_middle)

    # Add columns Country and FamCode to data frame
    index <- fmatch(chr_pos$sample_id, brca_data$Onc_ID)
    # chr_pos[c("Country", "FamCode")] <- brca_data[index, ] %>% select(Country, FamCode)
    chr_pos[c("Country", "FamCode", "cluster_groups")] <- brca_data[index, c("Country", "FamCode", "cluster_groups")]

    # Add fake-break at max and min dist for samples with no breaks, so they show in plots
    for (sample in unique(matched_plink$SNP)[!(unique(matched_plink$SNP) %in% unique(chr_pos$sample_id))]){
        if (is.na(sample)) next
        df_max = data.frame("SNP"="None", Chr_numeric=brca_chr, position_b37=max(chr_pos$position_b37), sample_id=sample, position_b37_adjusted=max(chr_pos$position_b37_adjusted))
        df_min = data.frame("SNP"="None", Chr_numeric=brca_chr, position_b37=min(chr_pos$position_b37), sample_id=sample, position_b37_adjusted=min(chr_pos$position_b37_adjusted))
        df <- rbind(df_max,df_min)
        df = cbind(df, brca_data %>% filter(Onc_ID == sample) %>% select(Country, FamCode, cluster_groups))
        chr_pos = rbind(chr_pos, df)
    }

    # Print some data information - here counting samples by Country
    print(chr_pos %>% distinct(sample_id, FamCode, Country) %>% count(Country))
    print(dim(chr_pos))
    # unique(chr_pos$sample_id)
    # print(chr_pos)

    return(chr_pos)
}

### Map nearest found breaks to their genomic position ###
mapNearestSNPs_plink <- function(matched_plink, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name){
    chr_pos = NULL

    # remove snps in brca gene
    if (removeSNPsInGene == TRUE){
        chr_coords2 <- filter(chr_coords, position_b37 < brca_start | position_b37 > brca_stop)
    } else {
        chr_coords2 <- chr_coords
    }

    # For each sample (patient) with given mutation
    for (i in 1:nrow(breaks)){
        # If no breaks for sample, then skip and continue
        if (length(which(breaks[i,]==T)) == 0) next
        # Extract rows with homozygote SNPs
        x = filter(chr_coords2, SNP %in% names(which(breaks[i,]==T)))
        # x = filter(chr_coords2, SNP %in% names(breaks[i,which(breaks[i,]==T)]))
        no_pos_break=F
        if (length(which(x$position_b37 > brca_middle)) == 0){
            #pos = which(chr_coords2$position_b37==max(chr_coords2$position_b37))
            no_pos_break=T
        } else {
            pos = which(x$position_b37 == min(x$position_b37[x$position_b37 > brca_middle]))
        }
        no_neg_break=F
        if (length(which(x$position_b37 < brca_middle)) == 0){
            #neg = which(chr_coords2$position_b37==min(chr_coords2$position_b37))
            no_neg_break=T
        } else {
            neg = which(x$position_b37 == max(x$position_b37[x$position_b37 < brca_middle]))
        }
        #positive_pos = min(x$position_b37[which(x$position_b37 > brca_middle)])
        #negative_pos = max(x$position_b37[which(x$position_b37 < brca_middle)])

        # If breaks on both sides
        if (!no_neg_break && !no_pos_break){
            chr_pos2 = data.frame(sample_id=matched_plink[i,1], positive_pos=x[pos,3], negative_pos=x[neg,3], positive_cM=x[pos,4], negative_cM=x[neg,4])
        # if no positive break
        } else if (no_pos_break) {
            chr_pos2 = data.frame(sample_id=matched_plink[i,1], positive_pos=NA, negative_pos=x[neg,3], positive_cM=NA, negative_cM=x[neg,4])
        # if no negative break
        } else if (no_neg_break) {
            chr_pos2 = data.frame(sample_id=matched_plink[i,1], positive_pos=x[pos,3], negative_pos=NA, positive_cM=x[pos,4], negative_cM=NA)
        }
        
        # Remove SNPs in brca gene
        # chr_pos2 = filter(chr_pos2, as.numeric(as.character(negative_pos)) < brca_start | as.numeric(as.character(positive_pos)) > brca_stop)

        # Combine all samples into big data frame for plotting (likely not most efficient way)
        chr_pos = bind_rows(chr_pos, chr_pos2)
    }
    
    if (is.null(chr_pos)){
        return(NULL)
    }
    
    # Add missing values
    # chr_pos$positive_pos[is.na(chr_pos$positive_pos)] = max(chr_pos$positive_pos, na.rm = T)
    # chr_pos$negative_pos[is.na(chr_pos$negative_pos)] = min(chr_pos$negative_pos, na.rm = T)
    # chr_pos$positive_cM[is.na(chr_pos$positive_cM)] = max(chr_pos$positive_cM, na.rm = T)
    # chr_pos$negative_cM[is.na(chr_pos$negative_cM)] = min(chr_pos$negative_cM, na.rm = T)
    
    # Add haplo_length fields
    chr_pos <- mutate(chr_pos, haplo_length = positive_pos - negative_pos, haplo_length_cM = positive_cM - negative_cM)
    
    # Add columns Country and FamCode to data frame
    index <- fmatch(chr_pos$sample_id, brca_data$Onc_ID)
    chr_pos[c("Country", "FamCode")] <- brca_data[index, ] %>% select(Country, FamCode)

    # Add fake-break at max and min dist for samples with no breaks, so they show in plots
    for (sample in unique(matched_plink$SNP)[!(unique(matched_plink$SNP) %in% unique(chr_pos$sample_id))]){
        if (is.na(sample)) next
        #df_max = data.frame("SNP"="None", Chr_numeric=brca_chr, position_b37=max(chr_coords$position_b37), cM=max(chr_coords$cM), sample_id=sample, position_b37_adjusted=max(chr_coords$position_b37)-brca_middle)
        #df_min = data.frame("SNP"="None", Chr_numeric=brca_chr, position_b37=min(chr_coords$position_b37), cM=min(chr_coords$cM), sample_id=sample, position_b37_adjusted=min(chr_coords$position_b37)-brca_middle)
        #df <- rbind(df_max,df_min)
        df <- data.frame(sample_id=sample, positive_pos=max(chr_coords2$position_b37), negative_pos=min(chr_coords2$position_b37),
                         positive_cM=max(chr_coords2$cM), negative_cM=min(chr_coords2$cM),
                         haplo_length=max(chr_coords2$position_b37)-min(chr_coords2$position_b37),
                         haplo_length_cM=max(chr_coords2$cM)-min(chr_coords2$cM))
        df = cbind(df, brca_data %>% filter(Onc_ID == sample) %>% select(Country, FamCode))#, cluster_groups))
        chr_pos = rbind(chr_pos, df)
    }

    # Add fake-break if only max dist has no break
    for (sample in which(is.na(chr_pos$positive_pos))){
        ## Max value as seen in other
        chr_pos[sample,][is.na(chr_pos[sample,])] = c(max(chr_coords2$position_b37),
                                              max(chr_coords2$cM),
                                              max(chr_coords2$position_b37)-chr_pos[sample,"negative_pos"],
                                              max(chr_coords2$cM)-chr_pos[sample,"negative_cM"])
    }

    # Add fake-break if only min dist has no break
    for (sample in which(is.na(chr_pos$negative_pos))){
        chr_pos[sample,][is.na(chr_pos[sample,])] = c(min(chr_coords2$position_b37),
                                                    min(chr_coords2$cM),
                                                    chr_pos[sample,"positive_pos"]-min(chr_coords2$position_b37),
                                                    chr_pos[sample,"positive_cM"]-min(chr_coords2$cM))
    }

    # Print some data information - here counting samples by Country
    #print(chr_pos %>% distinct(sample_id, FamCode, Country) %>% count(Country))
    print(dim(chr_pos))
    # unique(chr_pos$sample_id)
    # print(chr_pos)

    return(chr_pos)
}

## ---- mapSNPs
### Map found SNPs to their genomic position ###
mapSNPs <- function(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name){
#mapSNPs <- function(){
    chr_pos = NULL
    # For each sample (patient) with given mutation
    for (i in 1:length(breaks)){
        # Extract rows with SNPs in brca1 chromosome (17) and where SNP is homozygote for that sample
        #chr_pos2 = filter(chr_coords, chr_coords$Chr_numeric == brca_chr, chr_coords$SNP %in% breaks[[i]])
		chr_pos2 = filter(chr_coords, chr_coords$SNP %in% breaks[[i]])
		# print(dim(chr_pos2))
        # Removes SNPs in the brca gene
		if (removeSNPsInGene == TRUE){
            chr_pos2 = filter(chr_pos2, as.numeric(as.character(position_b37)) < brca_start | as.numeric(as.character(position_b37)) > brca_stop)
		}
        # Add sample_id column
        chr_pos2["sample_id"] = names(breaks[i])
        # Add position_b37_adjusted column, where brca1 = 0
        chr_pos2 = mutate(chr_pos2, position_b37_adjusted = as.numeric(as.character(position_b37)) - brca_middle)
        # Combine all samples into big data frame for plotting (likely not most efficient way)
        chr_pos = bind_rows(chr_pos, chr_pos2)
    }
    # Add columns Country and FamCode to data frame
    index <- fmatch(chr_pos$sample_id, brca_data$Onc_ID)
    # chr_pos[c("Country", "FamCode")] <- brca_data[index, ] %>% select(Country, FamCode)
    chr_pos[c("Country", "FamCode", "cluster_groups")] <- brca_data[index, ] %>% select(Country, FamCode, cluster_groups)

    # Print some data information - here counting samples by Country
    print(chr_pos %>% distinct(sample_id, FamCode, Country) %>% count(Country))
    print(dim(chr_pos))
    # unique(chr_pos$sample_id)
    # print(chr_pos)

    return(chr_pos)
}

###############################################################################
### Plotting
###############################################################################
### Define colors for plot
colors_plot <- function(pheno_data = brca_data, palette = "Set1", color_names = NULL, otherCountry = F){
    # Create colors for use in plot
    # palette defines the colorpalette to use e.g. Set1, Dark2, Set2, Set3 etc.
    # if color_names is given, then produce colors equal to the length of that vector
    # Otherwise, produce colors for each country in pheno_data
    # if otherCountry set, then combine the country with frequencies below 0.05 of total dataset

    # Get number of groups
    if (!is.null(color_names)){
        country_freq = color_names
    } else if (otherCountry){
        country_freq <- table(fct_lump(pheno_data$Country, prop = 0.05))
    } else {
        country_freq <- sort(table(as.character(pheno_data$Country)), decreasing = T)
    }


    if (palette %in% c("kelly", "sasha") && length(country_freq) > 20){
        print(gene)
        print(mut)
        print(pheno_data)
        print(country_freq)
        print("Too many countries for kelly or sasha colors! Switching to Dark2")
        palette = "Dark2"
    }

    if (palette == "sasha"){
        # Sasha's 20 colors - https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
        sasha_colors = c('#e6194b', '#3cb44b', '#ffe119', '#0082c8', '#f58231',
                         '#911eb4', '#46f0f0', '#f032e6', '#d2f53c', '#fabebe',
                         '#008080', '#e6beff', '#aa6e28', '#fffac8', '#800000',
                         '#aaffc3', '#808000', '#ffd8b1', '#000080', '#808080',
                         '#000000','#FFFFFF') # Black and white
        cols <- sasha_colors[1:length(country_freq)]
    } else if (palette == "sasha2"){
        # Sasha's 20 colors rearranged for heatmap annotation
        sasha_colors = c('#e6194b', '#3cb44b', '#0082c8', '#f58231', '#911eb4', 
                         '#46f0f0', '#f032e6', '#ffe119', '#d2f53c', '#fabebe',
                         '#008080', '#e6beff', '#aa6e28', '#fffac8', '#800000',
                         '#aaffc3', '#808000', '#ffd8b1', '#000080', '#808080',
                         '#000000','#FFFFFF') # Black and white
        cols <- sasha_colors[1:length(country_freq)]
    } else if (palette == "kelly"){
        # Kelly's 22 colors of maximum contrast
        # https://gist.github.com/ollieglass/f6ddd781eeae1d24e391265432297538
        kelly_colors = c('#F3C300', '#875692', '#F38400',
                         '#A1CAF1', '#BE0032', '#C2B280', '#848482', '#008856',
                         '#E68FAC', '#0067A5', '#F99379', '#604E97', '#F6A600',
                         '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522',
                         '#E25822', '#2B3D26', '#222222')
        cols <- kelly_colors[1:length(country_freq)]
    } else if (palette == "ggplot"){
        n = length(country_freq)
        hues = seq(15, 375, length = n + 1)
        cols = hcl(h = hues, l = 65, c = 100)[1:n]
    } else {
        # Get number of colors according to number of groups
        n_palette_col = if (palette == "Set1") 9 else 8 # Dark2=8, Set1=9
        if (length(country_freq) <= n_palette_col){
            cols <- brewer.pal(n = max(length(country_freq), 3), name = palette)[1:(length(country_freq))]
        } else {
            cols <- colorRampPalette(brewer.pal(n = n_palette_col, name = palette))(length(country_freq))
        }
    }

    #cols <- rainbow_hcl(length(country_freq))

    # Assign color to country or group
    if (is.null(color_names)){
        names(cols) <- names(country_freq)
    } else {
        names(cols) <- color_names
    }

    # Make it global (for convenience) and return
    #cols <<- cols
    return(cols)
}

### Visual output of haplotypes and breaks ###
plot_breakpoints_country <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, otherCountry = F){ #, country = F){

    # Define y axis step distance
    tick_steps = 3000000
    # Set the adjusted beginning and ending position of brca gene
    brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
    # Number of samples
    num_samples <<- length(unique(chr_pos$sample_id))

    if (otherCountry){
        quantile = 0.05
        countries <- chr_pos %>% distinct(sample_id, Country) %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        chr_pos2 = chr_pos %>% mutate(Country = fct_other(Country, keep = countries))
        # chr_pos2 = chr_pos %>% mutate(Country = as.character(Country))
        # index = !(chr_pos2$Country %in% countries)
        # chr_pos2[index,"Country"] = "Other"
        # chr_pos2$Country = as.factor(chr_pos2$Country)
    } else {
        chr_pos2 = chr_pos
    }
    # print(unique(chr_pos2$Country))
    # print(unique(chr_pos$Country))

    cols = colors_plot(brca_data, palette = "sasha", otherCountry = otherCountry)
    print(cols)

    #lab_cols = cols[as.character(distinct(chr_pos, sample_id, Country)$Country)]
    #lab_cols = cols[lab_cols[lab_cols != "Other"]]
    # Plot countries together that has more than 5% samples of total datasets
    p = chr_pos2 %>%
    #p = chr_pos %>% mutate(Country = fct_lump(Country, prop = 0.05)) %>% # group_by(Country) %>%

        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        ggplot(aes(x=interaction(sample_id, Country), y=position_b37_adjusted)) +
        #ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_b37_adjusted)) +

        # Plots the breakpoints around the brca gene
        # geom_point(aes(color = Country), size = 0.8) +
        geom_point(aes(color = Country), size = 1.2) +
        #geom_point(aes(color = cols[Country]), size = 0.8) +
        #geom_point(aes(color = as.character(cluster_groups)), size = 0.8) +
        #geom_point(size = 0.8, color = "darkblue") +

        # Plots the area of the brca gene
        # Plot before points to keep behind
        #geom_hline(yintercept = 0, color = "red", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.4) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.4) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +
        
        theme_minimal() + # theme
        
        # Defines the color of the dots
        #scale_color_brewer(palette="Set1") +
        scale_color_manual(values = cols) +

        #geom_hline(aes(yintercept = 5*tick_steps), color = "white", size = 0.2) +
        #geom_hline(aes(yintercept = -5*tick_steps), color = "white", size = 0.2) +

        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", mut)) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title

        # Bold face for legend title
        theme(legend.title = element_text(face="bold")) +
        
        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols)) +  # vertical x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[as.character(distinct(chr_pos2, sample_id, Country)$Country)])) +  # vertical x-labels
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = distinct(chr_pos2, sample_id, Country)$Country)) +  # vertical x-labels
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = lab_cols)) +  # vertical x-labels

        # Define x-ticks and their label names
        scale_x_discrete(breaks=interaction(chr_pos2$sample_id, chr_pos2$Country), labels = chr_pos2$FamCode) +
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country)) +

        # Define y-ticks and their label names (sets BRCA name in y axis)
        scale_y_continuous(breaks = c(seq(from = tick_steps, to = 5*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -5*tick_steps, by = -tick_steps), 0),
                           labels = c(seq(from = tick_steps, to = 5*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -5*tick_steps, by = -tick_steps), brca_name)) + #, limits = c(-9000000,9000000)) +


        # Define visible y-axis
        # ylim(-10000000,10000000) +
        #coord_cartesian(ylim=c(-2000000,2000000)) + # preferred

        # Set x and y axis label names
        ylab(paste("Physical distance (bp) from", gene)) + xlab("Sample")

    print(p)
    #p = p #+ theme(axis.text.x = element_text(colour = cols[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]))
    #print(p)
    # print(p + ylim(-1000000,1000000))

    # Define 'consensus' for filename
    #consensus <- if (is.null(country)) "All_consensus" else paste0(country, "_consensus")

    # Save plots to SVG and PNG files
    # ggsave(paste0("plots/", brca_name, "-", mut_name, "-", length(num_samples), "_samples.svg"))
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", num_samples, "_samples.png"), scale = 1.5)
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", num_samples, "_samples.png"), scale = 1, width = 16, height = 10)

    return(p)
}

### Produces four plots, ordered by country, by haplotype length, by positve genomic points and by negative points - all breakpoints coloured by country.
plot_nearest_brca_break_country <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, otherCountry=F, fixObject=T){ #, country){

    if (fixObject){
        # Extract (positive position) breakpoint closest to brca gene
        positive <- chr_pos %>%
            group_by(sample_id) %>%
            filter(position_b37_adjusted >= 0) %>%
            slice(which.min(position_b37_adjusted)) %>%
            rename(position_pos = position_b37, position_pos_adjusted = position_b37_adjusted, snp_pos = SNP, cM_pos = cM, cM_pos_adjusted = cM_adjusted)
        # Extract (negative position) breakpoint closest to brca gene
        negative <- chr_pos %>%
            group_by(sample_id) %>%
            filter(position_b37_adjusted <= 0) %>%
            slice(which.max(position_b37_adjusted)) %>%
            rename(position_neg = position_b37, position_neg_adjusted = position_b37_adjusted, snp_neg = SNP, cM_neg = cM, cM_neg_adjusted = cM_adjusted)

        # chr_pos_minmax <- rbind(positive, negative)
        chr_pos_minmax <- merge(positive, negative, by = intersect(names(positive), names(negative))) %>% mutate(haplo_length = position_pos-position_neg)
    } else {
        chr_pos_minmax <- chr_pos
    }

    # Define y axis step distance
    tick_steps = 1000000
    # Set the adjusted beginning and ending position of brca gene
    brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
    # Number of samples
    num_samples = length(unique(chr_pos_minmax$sample_id))
    # Point size
    if (num_samples < 20){
        point_size <- 2.5
    } else if (num_samples < 100){
        point_size <- 2
    } else if (num_samples < 200){
        point_size <- 1.5
    } else {
        point_size <- 0.8
    }

    #otherCountry = T
    if (otherCountry){
        cols = colors_plot(chr_pos_minmax, palette = "sasha", otherCountry = T)
        chr_pos_minmax <- chr_pos_minmax %>% mutate(Country = fct_lump(Country, prop = 0.05))
    } else {
        cols = colors_plot(chr_pos_minmax, palette = "sasha", otherCountry = F)
    }

    # Plot countries together that has less than 5% samples of total datasets
    #p = chr_pos_minmax %>% #mutate(Country = fct_lump(Country, prop = 0.1)) %>% # group_by(Country) %>%
    p = chr_pos_minmax %>% #group_by(Country) %>%

        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        ggplot(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_pos_adjusted)) +

        # Plots the area of the brca gene
        #geom_hline(yintercept = 0, color = "darkblue", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.1) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.1) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +

        # Plots the breakpoints around the brca gene
        geom_point(aes(color = Country), size = point_size) +

        # OBS: Negative points need special care, depending on plot needed. See printing of plots below.
        #geom_point(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_neg_adjusted, color = Country), size = 0.8) +

        # Defines the color of the dots
        #scale_color_brewer(palette="Set1") +
        scale_color_manual(values = cols) +

        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", mut_name)) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title

        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)] )) +  # vertical x-labels
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols)) +  # vertical x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[chr_pos_minmax$Country])) +  # vertical x-labels

        # Define x-ticks and their label names
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country), labels = chr_pos$FamCode) +
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country)) +

        # Define y-ticks and their label names (sets BRCA name in y axis)
        scale_y_continuous(breaks = c(seq(from = tick_steps, to = 15*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -15*tick_steps, by = -tick_steps), 0),
                           labels = c(seq(from = tick_steps, to = 15*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -15*tick_steps, by = -tick_steps), brca_name)) +


        # Define visible y-axis
        #ylim(-1000000,1000000) +

        # Set x and y axis label names
        ylab(paste("Physical distance (bp) from", gene)) + xlab("Sample")

    # Decide consensus for filename
    #consensus <- if (is.null(country)) "All_consensus" else paste0(country, "_consensus")

    # Negative points added when sorting by country
    negative_positions <- geom_point(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_neg_adjusted,
                                         color = Country), size = point_size)

    # Order points by country
    p1 = p + negative_positions + #geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +
        scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$Country), labels = chr_pos_minmax$FamCode) #+
        #theme(axis.text.x = element_text(colour = cols[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]))
    p1

    # Save plots to PNG file
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-country-", num_samples, "_samples.png"))


    # Negative points added when sorting by genomic position
    negative_positions <- geom_point(aes(x=reorder(sample_id, position_pos_adjusted), y=position_neg_adjusted,
                                         color = Country), size = point_size)
    # negative_positions <- geom_point(aes(x=sample_id, y=position_neg_adjusted,
    #                                      color = Country), size = point_size)

    # Order points by positive genomics positions
    p2 = p + aes(x=reorder(sample_id, position_pos_adjusted)) + negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$position_pos_adjusted), labels = chr_pos_minmax$FamCode) #+
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]))
    p2

    # Order points by positive genomics positions
    #print(p + aes(x=reorder(sample_id, position_pos_adjusted)) + negative_positions + xlab("Sample"))
    # Save plots to PNG file
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-positive_order-", num_samples, "_samples.png"))


    # Order points by negative genomics positions
    p3 = p + aes(x=reorder(sample_id, -position_neg_adjusted)) + negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$position_neg_adjusted), labels = chr_pos_minmax$FamCode) #+
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]))
    p3

    p4 = p + aes(x=reorder(sample_id, haplo_length)) + negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$haplo_length), labels = chr_pos_minmax$FamCode) #+
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]))
    p4

    # Order points by negative genomics positions
    #print(p + aes(x=reorder(sample_id, -position_neg_adjusted)) + negative_positions + xlab("Sample"))
    # Save plots to PNG file
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-negative_order-", num_samples, "_samples.png"))

    return(list(p1,p2,p3,p4))
}


### Visual output of haplotypes and breaks ###
plot_breakpoints_group <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group, col = "position_b37_adjusted", zoom = F){

    # Number of samples
    num_samples <<- length(unique(chr_pos$sample_id))

    # Set up colors for plot
    label_cols <- colors_plot(palette = "sasha")[as.character(distinct(chr_pos, sample_id, Country)$Country)]
    cols_groups <- colors_plot(palette = "Set1", color_names = names(sort(table(cluster_groups), decreasing = T)) )
    #cols_groups <- colors_plot(palette = "Set1", color_names = 1:max(cluster_groups))
    # lol_cols <- c("#7570B3", "#1B9E77", "#E6AB02", "#66A61E", "#D95F02","#E7298A")
    # names(lol_cols) <- c("USA", "DENMARK", "GERMANY", "FRANCE", "SPAIN", "CANADA")

    #
    if (col %in% c("position_b37_adjusted", "physical_distance_from_brca", "physical_dist")) {
        # Define y axis step distance
        tick_steps = 3000000; tick_steps_from = tick_steps; tick_steps_to = 5*tick_steps
        tick_steps2 = -3000000; tick_steps_from2 = tick_steps2; tick_steps_to2 = 5*tick_steps2
        tick_steps3 = 3000000; tick_steps_from3 = tick_steps3; tick_steps_to3 = 5*tick_steps3
        tick_steps4 = -3000000; tick_steps_from4 = tick_steps4; tick_steps_to4 = 5*tick_steps4
        brca_pos = 0
        # Set the adjusted beginning and ending position of brca gene
        brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max

        p = chr_pos %>%
            # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
            # The genomic position defines the y axis
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_b37_adjusted)) +
            ylab(paste("Physical distance (bp) from", gene))

    } else if (col %in% c("position_b37", "physical_distance", "absolute_physical_dist")) {
        # Define y axis step distance
        brca_middle2 = round(brca_middle/10^6)*10^6
        tick_steps = 3000000; tick_steps_from = tick_steps; tick_steps_to = 5*tick_steps
        tick_steps2 = -3000000; tick_steps_from2 = tick_steps2; tick_steps_to2 = 5*tick_steps2
        tick_steps3 = 3000000; tick_steps_from3 = brca_middle2+tick_steps3; tick_steps_to3 = brca_middle2+5*tick_steps3
        tick_steps4 = -3000000; tick_steps_from4 = brca_middle2+tick_steps4; tick_steps_to4 = brca_middle2+5*tick_steps4

        brca_pos = 0
        # Set the adjusted beginning and ending position of brca gene
        brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
        #brca_max = brca_stop; brca_min = brca_start

        p = chr_pos %>%
            # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
            # The genomic position defines the y axis
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_b37_adjusted)) +
            ylab("Physical distance (bp)")

    } else if (col %in% c("centimorgan", "genetic_distance", "absolute_genetic_dist")) {

        brca_pos2 = round(brca_cM_middle)
        tick_steps = 3; tick_steps_from = tick_steps; tick_steps_to = 15*tick_steps
        tick_steps2 = -3; tick_steps_from2 = tick_steps2; tick_steps_to2 = 15*tick_steps2
        tick_steps3 = 3; tick_steps_from3 = brca_pos2+tick_steps3; tick_steps_to3 = brca_pos2+15*tick_steps3
        tick_steps4 = -3; tick_steps_from4 = brca_pos2+tick_steps4; tick_steps_to4 = brca_pos2+15*tick_steps4

        brca_pos = 0
        brca_max = 0.1; brca_min = -0.1

        p = chr_pos %>%
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=cM_adjusted)) +
            ylab("Genetic distance (cM)")
    } else if (col %in% c("genetic_distance_from_gene", "genetic_dist", "cm", "cM")) {
        brca_pos = 0
        brca_max = 0.1; brca_min = -0.1

        tick_steps = 3; tick_steps_from = tick_steps; tick_steps_to = 15*tick_steps
        tick_steps2 = -3; tick_steps_from2 = tick_steps2; tick_steps_to2 = 15*tick_steps2
        tick_steps3 = 3; tick_steps_from3 = tick_steps3; tick_steps_to3 = 15*tick_steps3
        tick_steps4 = -3; tick_steps_from4 = tick_steps4; tick_steps_to4 = 15*tick_steps4

        p = chr_pos %>%
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=cM_adjusted)) +
            ylab(paste("Genetic distance (cM) from", gene))
    }

    if (zoom && col %in% c("genetic_distance_from_gene", "genetic_dist")){
        tick_steps = 1; tick_steps_from = tick_steps
        tick_steps2 = -1; tick_steps_from2 = tick_steps2
        tick_steps3 = 1; tick_steps_from3 = tick_steps3
        tick_steps4 = -1; tick_steps_from4 = tick_steps4
    } else if (zoom && col %in% c("position_b37_adjusted", "physical_distance_from_brca", "physical_dist")){
        tick_steps = 10^6; tick_steps_from = tick_steps
        tick_steps2 = -10^6; tick_steps_from2 = tick_steps2
        tick_steps3 = 10^6; tick_steps_from3 = tick_steps3
        tick_steps4 = -10^6; tick_steps_from4 = tick_steps4
    }

    # Plot countries together that has less than 5% samples of total datasets
    #p = chr_pos %>% #mutate(Country = fct_lump(Country, prop = 0.05)) %>% # group_by(Country) %>%


    p = p +
        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        #ggplot(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_b37_adjusted)) +
        #ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_b37_adjusted)) +
        #ggplot(aes(x=interaction(FamCode, cluster_groups), y=position_b37_adjusted)) +

        # Plots the area of the brca gene - plot before points to keep behind
        #geom_hline(yintercept = 0, color = "red", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.4) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.4) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +

        # Plots the breakpoints around the brca gene
        #geom_point(aes(color = Country), size = 0.8) +
        geom_point(aes(color = factor(cluster_groups)), size = 0.8) +
        #geom_point(size = 0.8, color = "darkblue") +

        theme_minimal() + # theme
        
        # Defines the color of the dots
        #scale_color_brewer(palette="Set1") +
        scale_color_manual(values = cols_groups) +
        #scale_color_hue() +

        #scale_fill_manual(values = c("#7570B3", "#1B9E77", "#E6AB02", "#66A61E", "#D95F02","#E7298A"), breaks = c("USA", "DENMARK", "GERMANY", "FRANCE", "SPAIN", "CANADA")) +

        #geom_hline(aes(yintercept = 5*tick_steps), color = "white", size = 0.2) +
        #geom_hline(aes(yintercept = -5*tick_steps), color = "white", size = 0.2) +

        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", mut, "- Consensus: Group", group)) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title
        
        # Bold face for legend title
        theme(legend.title = element_text(face="bold")) +

        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = label_cols)) +  # vertical x-labels

        # Define x-ticks and their label names
        scale_x_discrete(breaks = interaction(chr_pos$sample_id, chr_pos$cluster_groups),
                         labels = chr_pos$FamCode) +
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country)) +

        scale_y_continuous(breaks = c(seq(from = tick_steps_from, to = tick_steps_to, by = tick_steps),
                                      seq(from = tick_steps_from2, to = tick_steps_to2, by = tick_steps2), brca_pos),
                           labels = c(seq(from = tick_steps_from3, to = tick_steps_to3, by = tick_steps3),
                                      seq(from = tick_steps_from4, to = tick_steps_to4, by = tick_steps4), brca_name)) +

        # Define visible y-axis
        #ylim(-1000000,1000000) +
        #coord_cartesian(ylim = c(-10^6,10^6)) +

        # Legend title
        labs(colour="Groups") +

        # Set x and y axis label names
        #ylab("Genomic position") +
        xlab("Sample")

    if (zoom){
        p = p + coord_cartesian(ylim = c(-zoom,zoom))
    } else {
        # Define y-ticks and their label names (sets BRCA name in y axis)
        p = p
    }

    print(p)
    # print(p + ylim(-1000000,1000000))

    # Define 'consensus' for filename
    #consensus <- if (is.null(country)) "All_consensus" else paste0(country, "_consensus")

    # Save plots to SVG and PNG files
    # ggsave(paste0("plots/", brca_name, "-", mut_name, "-", length(num_samples), "_samples.svg"))
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", num_samples, "_samples.png"), scale = 1.5)
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", num_samples, "_samples.png"), scale = 1, width = 16, height = 10)

    return(p)
}

plot_nearest_brca_break_group <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group, col = "position_b37_adjusted"){

    # Extract (positive position) breakpoint closest to brca gene
    pos <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted >= 0) %>%
        slice(which.min(position_b37_adjusted)) %>% #slice(1) %>%
        rename(position_pos = position_b37, position_pos_adjusted = position_b37_adjusted,
               cM_pos_adjusted = cM_adjusted, snp_pos = SNP, cM_pos = cM)
    # Extract (negative position) breakpoint closest to brca gene
    neg <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted <= 0) %>%
        arrange(sample_id, desc(position_b37_adjusted)) %>%
        slice(which.max(position_b37_adjusted)) %>% # slice(1) %>%
        rename(position_neg = position_b37, position_neg_adjusted = position_b37_adjusted,
               cM_neg_adjusted = cM_adjusted, snp_neg = SNP, cM_neg = cM)

    chr_pos_minmax <- merge(pos, neg, by = intersect(names(pos), names(neg))) %>%
        mutate(haplo_length = position_pos-position_neg)

    # Define y axis step distance
    #tick_steps = 1000000
    # Set the adjusted beginning and ending position of brca gene
    #brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
    # Number of samples
    num_samples = length(unique(chr_pos_minmax$sample_id))
    # Point size
    if (num_samples < 20){
        point_size <- 2.5
    } else if (num_samples < 100){
        point_size <- 2
    } else if (num_samples < 200){
        point_size <- 1.5
    } else {
        point_size <- 0.8
    }

    line_size = 0.1


    # Set up genomic distance measure
    if (col %in% c("position_b37_adjusted", "physical_distance_from_brca", "physical_dist")) {
        # Define y axis step distance
        tick_steps = 1000000; tick_steps_from = tick_steps; tick_steps_to = 15*tick_steps
        tick_steps2 = -1000000; tick_steps_from2 = tick_steps2; tick_steps_to2 = 15*tick_steps2
        tick_steps3 = 1000000; tick_steps_from3 = tick_steps3; tick_steps_to3 = 15*tick_steps3
        tick_steps4 = -1000000; tick_steps_from4 = tick_steps4; tick_steps_to4 = 15*tick_steps4

        brca_pos = 0
        # Set the adjusted beginning and ending position of brca gene
        brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max

        p = chr_pos_minmax %>%
            # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
            # The physical position defines the y axis
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
            ylab(paste("Physical distance (bp) from", gene))

    }
    else if (col %in% c("position_b37", "physical_distance", "absolute_physical_dist")) {
        # Define y axis step distance
        brca_middle2 = round(brca_middle/10^6)*10^6 # Rounding brca_middle for nice output
        tick_steps = 1000000; tick_steps_from = tick_steps; tick_steps_to = 15*tick_steps
        tick_steps2 = -1000000; tick_steps_from2 = tick_steps2; tick_steps_to2 = 15*tick_steps2
        tick_steps3 = 1000000; tick_steps_from3 = brca_middle2+tick_steps3; tick_steps_to3 = brca_middle2+15*tick_steps3
        tick_steps4 = -1000000; tick_steps_from4 = brca_middle2+tick_steps4; tick_steps_to4 = brca_middle2+15*tick_steps4

        brca_pos = 0
        # Set the adjusted beginning and ending position of brca gene
        brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
        #brca_max = brca_stop; brca_min = brca_start

        # Dirty fix to incorporate different genomic measures
        # chr_pos_minmax <<- chr_pos_minmax %>% rename(pos_old=position_pos_adjusted, neg_old=position_neg_adjusted,
        #                           position_pos_adjusted=position_pos, position_neg_adjusted=position_neg)

        p = chr_pos_minmax %>%
            # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
            # The physical position defines the y axis
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
            ylab("Physical distance (bp)")

        # negative_positions <- geom_col(aes(x=interaction(sample_id, cluster_groups),
        #                                    y=position_neg, fill = as.character(cluster_groups)))
        #
        # negative_pos_reorder <- geom_col(aes(x=reorder(sample_id, position_pos), y=position_neg,
        #                                      fill = as.character(cluster_groups)))

    }
    else if (col %in% c("centimorgan", "genetic_distance", "absolute_genetic_dist")) {
        brca_pos = 0
        brca_max = 0.1; brca_min = -0.1

        brca_pos2 = round(brca_cM_middle)
        tick_steps = 1; tick_steps_from = tick_steps; tick_steps_to = 15*tick_steps
        tick_steps2 = -1; tick_steps_from2 = tick_steps2; tick_steps_to2 = 15*tick_steps2
        tick_steps3 = 1; tick_steps_from3 = brca_pos2+tick_steps3; tick_steps_to3 = brca_pos2+15*tick_steps3
        tick_steps4 = -1; tick_steps_from4 = brca_pos2+tick_steps4; tick_steps_to4 = brca_pos2+15*tick_steps4

        # Dirty fix to incorporate different genomic measures
        chr_pos_minmax <- chr_pos_minmax %>% rename(pos_old=position_pos_adjusted, neg_old=position_neg_adjusted,
                                                    position_pos_adjusted=cM_pos_adjusted, position_neg_adjusted=cM_neg_adjusted)

        p = chr_pos_minmax %>%
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
            ylab("Genetic distance (cM)")
    }
    else if (col %in% c("genetic_distance_from_gene", "genetic_dist")) {
        brca_pos = 0
        brca_max = 0.1; brca_min = -0.1

        # Dirty fix to incorporate different genomic measures
        chr_pos_minmax <- chr_pos_minmax %>% rename(pos_old=position_pos_adjusted, neg_old=position_neg_adjusted,
                                                    position_pos_adjusted=cM_pos_adjusted, position_neg_adjusted=cM_neg_adjusted)

        tick_steps = 1; tick_steps_from = tick_steps; tick_steps_to = 15*tick_steps
        tick_steps2 = -1; tick_steps_from2 = tick_steps2; tick_steps_to2 = 15*tick_steps2
        tick_steps3 = 1; tick_steps_from3 = tick_steps3; tick_steps_to3 = 15*tick_steps3
        tick_steps4 = -1; tick_steps_from4 = tick_steps4; tick_steps_to4 = 15*tick_steps4

        p = chr_pos_minmax %>%
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
            ylab(paste("Genetic distance (cM) from", gene))
    }

    # Set up colors for plot
    # label_cols <- colors_plot(palette = "sasha")[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]
    # cols_groups <- colors_plot(palette = "Set1", color_names = names(sort(table(brca_data$cluster_groups), decreasing = T)) )
    label_cols <- colors_plot(palette = "sasha")[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]
    cols_groups <- colors_plot(palette = "Set1", color_names = names(sort(table(cluster_groups), decreasing = T)) )

    # Negative points added when sorting by country
    negative_positions_sampleOrder <- geom_point(aes(x=interaction(sample_id, cluster_groups),
                                                     y=position_neg_adjusted, color = as.character(cluster_groups)),
                                                     size = point_size)

    # Negative points added when sorting by genomic position
    negative_positions <- geom_point(aes(x=reorder(sample_id, position_pos_adjusted), y=position_neg_adjusted,
                                         color = as.character(cluster_groups)), size = point_size)

    # Plot countries together that has less than 5% samples of total datasets
    #p = chr_pos_minmax %>% #mutate(Country = fct_lump(Country, prop = 0.1)) %>% # group_by(Country) %>%
    #p = chr_pos_minmax %>% # mutate(Country = fct_lump(Country, prop = 0.05)) %>% #group_by(Country) %>%

    p = p +
        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        #ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +


        # Plots the area of the brca gene - plot before points to keep behind
        #geom_hline(yintercept = 0, color = "darkblue", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = line_size) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = line_size) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +


        # Plots the breakpoints around the brca gene
        geom_point(aes(color = as.character(cluster_groups)), size = point_size) +
        #geom_point(aes(color = Country), size = point_size) +

        # OBS: Negative points need special care, depending on plot needed. See printing of plots below.
        #geom_point(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_neg_adjusted, color = Country), size = 0.8) +


        # Defines the color of the dots
        scale_color_brewer(palette="Set1") +
        #scale_color_manual(values = cols_groups) +
        #scale_fill_manual(values = cols_groups) +

        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", mut, "- Group consensus:", group)) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title

        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = label_cols)) +  # vertical x-labels

        # Define x-ticks and their label names
        #scale_x_discrete(breaks=interaction(chr_pos_minmax$sample_id, chr_pos_minmax$Country), labels = chr_pos_minmax$sample_id) +
        #scale_x_discrete(breaks=interaction(chr_pos_minmax$sample_id, chr_pos_minmax$Country)) +
        # scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$cluster_groups),
        #                 labels = chr_pos_minmax$FamCode) +

        # Define y-ticks and their label names (sets BRCA name in y axis)
        # scale_y_continuous(breaks = c(seq(from = tick_steps, to = 15*tick_steps, by = tick_steps),
        #                               seq(from = -tick_steps, to = -15*tick_steps, by = -tick_steps), 0),
        #                    labels = c(seq(from = tick_steps, to = 15*tick_steps, by = tick_steps),
        #                               seq(from = -tick_steps, to = -15*tick_steps, by = -tick_steps), brca_name)) +
       scale_y_continuous(breaks = c(seq(from = tick_steps_from, to = tick_steps_to, by = tick_steps),
                                     seq(from = tick_steps_from2, to = tick_steps_to2, by = tick_steps2), brca_pos),
                          labels = c(seq(from = tick_steps_from3, to = tick_steps_to3, by = tick_steps3),
                                     seq(from = tick_steps_from4, to = tick_steps_to4, by = tick_steps4), brca_name)) +



        # Define visible y-axis
        #ylim(-1000000,1000000) +

        # Legend title
        labs(colour="Groups") +

        # Set x and y axis label names
        #ylab("Genomic position") +
        xlab("Sample")


    # Negative points added when sorting by position
    # negative_positions <- geom_point(aes(x=interaction(sample_id, cluster_groups),
    #                                     y=position_neg_adjusted, color = as.character(chr_pos_minmax$cluster_groups)), size = point_size)
    # negative_positions_sampleOrder <- geom_point(aes(x=interaction(sample_id, cluster_groups),
    #                                     y=position_neg_adjusted, color = as.character(cluster_groups)), size = point_size)

    # Order points by country
    p1 = p + negative_positions_sampleOrder +
        scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$cluster_groups), labels = chr_pos_minmax$FamCode)
        #scale_x_discrete(breaks = interaction(sample_id, cluster_groups), labels = FamCode)
    #print(p1)

    # Negative points added when sorting by genomic position
    # negative_positions <- geom_point(aes(x=reorder(sample_id, position_pos_adjusted), y=position_neg_adjusted,
    #                                      color = as.character(cluster_groups)), size = point_size)

    # Order points by positive genomics positions
    p2 = p +
        aes(x=reorder(sample_id, position_pos_adjusted)) +
        negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$position_pos_adjusted), labels = chr_pos_minmax$FamCode)
    #print(p2)

    # Order points by negative genomics positions
    p3 = p + aes(x=reorder(sample_id, -position_neg_adjusted)) + negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, -chr_pos_minmax$position_neg_adjusted), labels = chr_pos_minmax$FamCode)
    #print(p3)

    p4 = p + aes(x=reorder(sample_id, haplo_length)) + negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$haplo_length), labels = chr_pos_minmax$FamCode)
    #print(p4)

    return(list(p1,p2,p3,p4))
}

plot_haplotype <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group=NULL, otherCountry=F, col = "position_b37_adjusted"){
    pos <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted >= 0) %>%
        #slice(1:2) %>%
        slice(1) %>%
        rename(position_pos = position_b37, position_pos_adjusted = position_b37_adjusted,
               cM_pos_adjusted = cM_adjusted, snp_pos = SNP, cM_pos = cM)

    neg <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted <= 0) %>%
        arrange(sample_id, desc(position_b37_adjusted)) %>%
        #slice(1:2) %>%
        slice(1) %>%
        rename(position_neg = position_b37, position_neg_adjusted = position_b37_adjusted,
               cM_neg_adjusted = cM_adjusted, snp_neg = SNP, cM_neg = cM)

    chr_pos_minmax <- merge(pos, neg, by = intersect(names(pos), names(neg))) %>%
        mutate(haplo_length = position_pos-position_neg)

    # Define y axis step distance
    #tick_steps = 1000000; tick_steps_from = tick_steps; tick_steps_to = 15*tick_steps
    # Set the adjusted beginning and ending position of brca gene
    #brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
    # Number of samples
    num_samples = length(unique(chr_pos_minmax$sample_id))

    # Set up colors for plot
    if (!is.null(group)){
        label_cols <- colors_plot(palette = "sasha")[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]
        cols_groups <- colors_plot(palette = "Set1", color_names = names(sort(table(cluster_groups), decreasing = T)) )
        #label_cols <- colors_plot(palette = "Dark2")
    } else if (otherCountry){
        cols_groups = colors_plot(chr_pos_minmax, otherCountry = T)
        label_cols = colors_plot(chr_pos_minmax, otherCountry = T)
        chr_pos_minmax <- chr_pos_minmax %>% mutate(Country = fct_lump(Country, prop = 0.05))
    } else {
        cols_groups = colors_plot(brca_data, otherCountry = F)
        label_cols = colors_plot(brca_data, otherCountry = F)
    }

    # Set up genomic distance measure
    if (col %in% c("position_b37_adjusted", "physical_distance_from_brca", "physical_dist")) {
        # Define y axis step distance
        tick_steps = 1000000; tick_steps_from = tick_steps; tick_steps_to = 100*tick_steps
        tick_steps2 = -1000000; tick_steps_from2 = tick_steps2; tick_steps_to2 = 100*tick_steps2
        tick_steps3 = 1000000; tick_steps_from3 = tick_steps3; tick_steps_to3 = 100*tick_steps3
        tick_steps4 = -1000000; tick_steps_from4 = tick_steps4; tick_steps_to4 = 100*tick_steps4

        brca_pos = 0
        # Set the adjusted beginning and ending position of brca gene
        brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max

        p = chr_pos_minmax %>%
            # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
            # The physical position defines the y axis
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
            ylab(paste("Physical distance (bp) from", gene))

    }
    else if (col %in% c("position_b37", "physical_distance", "absolute_physical_dist")) {
        # Define y axis step distance
        brca_middle2 = round(brca_middle/10^6)*10^6 # Rounding brca_middle for nice output
        tick_steps = 1000000; tick_steps_from = tick_steps; tick_steps_to = 100*tick_steps
        tick_steps2 = -1000000; tick_steps_from2 = tick_steps2; tick_steps_to2 = 100*tick_steps2
        tick_steps3 = 1000000; tick_steps_from3 = brca_middle2+tick_steps3; tick_steps_to3 = brca_middle2+100*tick_steps3
        tick_steps4 = -1000000; tick_steps_from4 = brca_middle2+tick_steps4; tick_steps_to4 = brca_middle2+100*tick_steps4

        brca_pos = 0
        # Set the adjusted beginning and ending position of brca gene
        brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
        #brca_max = brca_stop; brca_min = brca_start

        # Dirty fix to incorporate different genomic measures
        # chr_pos_minmax <<- chr_pos_minmax %>% rename(pos_old=position_pos_adjusted, neg_old=position_neg_adjusted,
        #                           position_pos_adjusted=position_pos, position_neg_adjusted=position_neg)

        p = chr_pos_minmax %>%
            # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
            # The physical position defines the y axis
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
            ylab("Physical distance (bp)")

        # negative_positions <- geom_col(aes(x=interaction(sample_id, cluster_groups),
        #                                    y=position_neg, fill = as.character(cluster_groups)))
        #
        # negative_pos_reorder <- geom_col(aes(x=reorder(sample_id, position_pos), y=position_neg,
        #                                      fill = as.character(cluster_groups)))

    }
    else if (col %in% c("centimorgan", "genetic_distance", "absolute_genetic_dist")) {
        brca_pos = 0
        brca_max = 0.1; brca_min = -0.1

        brca_pos2 = round(brca_cM_middle)
        tick_steps = 1; tick_steps_from = tick_steps; tick_steps_to = 100*tick_steps
        tick_steps2 = -1; tick_steps_from2 = tick_steps2; tick_steps_to2 = 100*tick_steps2
        tick_steps3 = 1; tick_steps_from3 = brca_pos2+tick_steps3; tick_steps_to3 = brca_pos2+100*tick_steps3
        tick_steps4 = -1; tick_steps_from4 = brca_pos2+tick_steps4; tick_steps_to4 = brca_pos2+100*tick_steps4

        # Dirty fix to incorporate different genomic measures
        chr_pos_minmax <- chr_pos_minmax %>% rename(pos_old=position_pos_adjusted, neg_old=position_neg_adjusted,
                                                     position_pos_adjusted=cM_pos_adjusted, position_neg_adjusted=cM_neg_adjusted)

        p = chr_pos_minmax %>%
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
            ylab("Genetic distance (cM)")
    }
    else if (col %in% c("genetic_distance_from_gene", "genetic_dist")) {
        brca_pos = 0
        brca_max = 0.1; brca_min = -0.1

        # Dirty fix to incorporate different genomic measures
        chr_pos_minmax <- chr_pos_minmax %>% rename(pos_old=position_pos_adjusted, neg_old=position_neg_adjusted,
                                  position_pos_adjusted=cM_pos_adjusted, position_neg_adjusted=cM_neg_adjusted)

        tick_steps = 1; tick_steps_from = tick_steps; tick_steps_to = 100*tick_steps
        tick_steps2 = -1; tick_steps_from2 = tick_steps2; tick_steps_to2 = 100*tick_steps2
        tick_steps3 = 1; tick_steps_from3 = tick_steps3; tick_steps_to3 = 100*tick_steps3
        tick_steps4 = -1; tick_steps_from4 = tick_steps4; tick_steps_to4 = 100*tick_steps4

        p = chr_pos_minmax %>%
            ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
            ylab(paste("Genetic distance (cM) from", gene))
    }



    if (!is.null(group)){
        #p = chr_pos_minmax %>%
        #    ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
        p = p +
            geom_col(aes(fill = as.character(cluster_groups))) + #, width = 1) +
            ggtitle(paste("Mutation:", mut, "- Consensus: Group", group)) + # Adds plot title
            theme_minimal() + # theme
            theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = label_cols)) +
            scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$cluster_groups), labels = chr_pos_minmax$FamCode) +
            scale_fill_manual(values = cols_groups) +
            labs(fill="Groups")

        negative_positions <- geom_col(aes(x=interaction(sample_id, cluster_groups),
                                           y=position_neg_adjusted, fill = as.character(cluster_groups)))

        negative_pos_reorder <- geom_col(aes(x=reorder(sample_id, position_pos_adjusted), y=position_neg_adjusted,
                                             fill = as.character(cluster_groups)))

    }
    else { # Country
        p = chr_pos_minmax %>%
            ggplot(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_pos_adjusted)) +
            geom_col(aes(fill = Country)) +
            theme_minimal() + # theme
            ggtitle(paste("Mutation:", mut)) + # Adds plot title
            theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = label_cols[as.character(chr_pos_minmax$Country)])) +
            theme(legend.title = element_text(face="bold")) +
            scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$Country), labels = chr_pos_minmax$FamCode)+
            scale_fill_manual(values = cols_groups, breaks = names(cols_groups)) +
            labs(fill="Country")

        negative_positions <- geom_col(aes(x=interaction(sample_id, fct_infreq(Country)),
                                           y=position_neg_adjusted, fill = Country)) #, width=1)

        negative_pos_reorder <- geom_col(aes(x=reorder(sample_id, position_pos_adjusted), y=position_neg_adjusted,
                                             fill = Country))
    }

    # Common
    p = p +

        # Plots the area of the brca gene - plot before points to keep behind
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.1) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.1) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +

        # Defines the color of the bars
        scale_color_manual(values = cols_groups) +
        
        # Define y-ticks and their label names (sets BRCA name in y axis)
        scale_y_continuous(breaks = c(seq(from = tick_steps_from, to = tick_steps_to, by = tick_steps),
                                      seq(from = tick_steps_from2, to = tick_steps_to2, by = tick_steps2), brca_pos),
                           labels = c(seq(from = tick_steps_from3, to = tick_steps_to3, by = tick_steps3),
                                      seq(from = tick_steps_from4, to = tick_steps_to4, by = tick_steps4), brca_name)) +

        # Set x and y axis label names
        ylab(paste("Physical distance (bp) from", gene)) +
        xlab("Sample")

    #print(p)

    # Order points by group
    p1 = p + negative_positions
    # if (!is.null(group)){
    #     p1 = p1 + scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$cluster_groups), labels = chr_pos_minmax$FamCode)
    # } else {
    #     p1 = p1 + scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$Country), labels = chr_pos_minmax$FamCode)
    # }
    print(p1)

    # Order points by positive genomics positions
    p2 = p +
        aes(x=reorder(sample_id, position_pos_adjusted)) +
        negative_pos_reorder + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$position_pos_adjusted), labels = chr_pos_minmax$FamCode)
    #print(p2)

    # Order points by negative genomics positions
    p3 = p + aes(x=reorder(sample_id, -position_neg_adjusted)) + negative_pos_reorder + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, -chr_pos_minmax$position_neg_adjusted), labels = chr_pos_minmax$FamCode)
    #print(p3)

    p4 = p + aes(x=reorder(sample_id, haplo_length)) + negative_pos_reorder + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$haplo_length), labels = chr_pos_minmax$FamCode)
    #print(p4)

    return(list(p1,p2,p3,p4))
}

# Adds second label to breakpoint plots - describing country color coding
add_second_legend <- function(p, chr_pos, otherCountry = F){
    if (otherCountry){
        sample_country <- chr_pos %>% distinct(sample_id, Country)
        label_cols = colors_plot(sample_country, otherCountry = T)
        country_freq <- names(table(fct_lump(sample_country$Country, prop = 0.05)))
        countries <- sample_country %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > 0.05) %>% .$Country
        chr_pos <- chr_pos %>% mutate(Country = fct_other(Country, keep = countries))
    } else {
        label_cols <- colors_plot(palette = "sasha")[as.character(distinct(chr_pos, sample_id, Country)$Country)]
        country_freq <- names(sort(table(names(label_cols)), decreasing = T))
    }
    val_col <- label_cols[country_freq]
    names(val_col) <- val_col
    if(is.null(p)){
        return(val_col)
    }
    #print(val_col)
    p2 = chr_pos %>% distinct(Country, .keep_all = T) %>%
        ggplot(aes(x=FamCode, y = position_b37)) + geom_col(aes(fill=Country)) +
        scale_fill_manual(values = names(val_col), labels=paste(country_freq, "   ")) +
        theme(legend.position='bottom') + labs(fill="") #labs(fill="Sample \n labels ") +
        #guides(fill = guide_legend(default.unit = "cm", label.hjust = 2, keywidth = 2))
    #print(p2)
    legend <- gtable_filter(ggplotGrob(p2), 'guide-box')
    p3 = arrangeGrob(p, bottom=legend$grobs[[1]]$grobs[[1]])
    grid.arrange(p3)
    return(p3)
}

## ---- findNearestBreaks
nearestBreaksStatistics <- function(chr_pos, mut, country = NULL, group = NULL, input = NULL){
    # Extract (positive position) breakpoint closest to brca gene
    positive <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted > 0) %>%
        slice(which.min(position_b37_adjusted))
    # Extract (negative position) breakpoint closest to brca gene
    negative <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted < 0) %>%
        slice(which.max(position_b37_adjusted))

    # Merge breakpoints into one data frame
    haplotypeDistance <- merge(select(positive, sample_id, Country, positive_pos = position_b37_adjusted),
                               select(negative, sample_id, negative_pos = position_b37_adjusted))
    # Compute haplotype length for each sample
    haplotypeDistance <- haplotypeDistance %>% mutate(haplo_distance = positive_pos + abs(negative_pos))

    if (!is.null(country)){
        # Arrange samples according to plot i.e. most frequent country samples first
        ordered_haploDist <- haplotypeDistance %>% arrange(fct_infreq(Country), sample_id)

        # Compute mean haplotype length for all countries
        ordered_mean_haploDist <- haplotypeDistance %>%
            group_by(Country) %>%
            #summarize_each(funs(mean(., na.rm = TRUE)), haplo_distance) %>%
            #summarise_at(.vars=c("haplo_distance"), .funs=c(Mean="mean")) %>%
            summarise(mean_haplo_distance=mean(haplo_distance, na.rm=TRUE)) %>%
            mutate(mean_haplo_distance = round(mean_haplo_distance)) %>%
            arrange(desc(mean_haplo_distance))

        #print(haplotypeDistance$haplo_distance)
        #print(mean(haplotypeDistance$haplo_distance))
        total <- data.frame("Overall", round(mean(haplotypeDistance$haplo_distance)))
        names(total) <- c("Country","mean_haplo_distance")
        ordered_mean_haploDist <- rbind(ordered_mean_haploDist, total)
    }

    if (!is.null(group)){
        haplotypeDistance <- merge(haplotypeDistance, select(brca_data, sample_id=Onc_ID, cluster_groups), by = "sample_id")

        ordered_haploDist <- haplotypeDistance %>% arrange(cluster_groups, haplo_distance, sample_id)

        ordered_mean_haploDist <- haplotypeDistance %>%
            group_by(cluster_groups) %>%
            summarise(mean_haplo_distance=mean(haplo_distance, na.rm=TRUE)) %>%
            mutate(mean_haplo_distance = round(mean_haplo_distance, digits = 2))
            #arrange(desc(mean_haplo_distance))

        total <- data.frame("Overall", round(mean(haplotypeDistance$haplo_distance)))
        names(total) <- c("cluster_groups","mean_haplo_distance")
        ordered_mean_haploDist <- rbind(ordered_mean_haploDist, total)
        names(ordered_mean_haploDist) <- c("Cluster","Mean Distance")
    }

    # Write tables to file
    #consensus <- if (is.null(country)) "All_consensus" else paste0(country, "_consensus") # Decide consensus for filename
    if (!is.null(input) && !is.null(group)){
        fam = if (input$fam) "fam" else "individual"
        write.table(ordered_haploDist, file = paste0("cache/nearestBreaksStatistics/", gene,
                    "-", mut_name, "-nearestBreaksStatistics-distMethod-", input$distMethod, "-clustMethod-", input$clustMethod,
                    "-", nrow(positive), "_samples-", fam, "-cutoff_", input$cutoff, "-group_", group, ".txt"),
                    sep = "\t", col.names = T, quote = F, row.names = F)
        write.table(ordered_mean_haploDist, file = paste0("cache/nearestBreaksStatistics/", gene,
                    "-", mut, "-", "mean_nearestBreaksStatistics-distMethod-", input$distMethod, "-clustMethod-", input$clustMethod,
                    "-", nrow(positive), "_samples-", fam, "-cutoff_", input$cutoff, "-group_", group, ".txt"),
                    sep = "\t", col.names = T, quote = F, row.names = F)
    }

    #print("DONE")

    #return(haplotypeDistance)
    return(ordered_mean_haploDist)
}

###############################################################################
### Helper functions
###############################################################################
convertPhysicalToGeneticDistance <- function(){
    linearInterpolation = T
    
    if(!linearInterpolation){
        ### Nearest snp in centimorgan map approach ###
        getGeneticCoord <- function(row, gene){
            #print(row)
            #cM = morgan.brca1[which.min(abs(morgan.brca1$V4-as.integer(row[3]))),3]
            if (gene == "BRCA1"){
                cM = morgan.brca1[which.min(abs(morgan.brca1$position-as.integer(row[3]))),3]
            } else {
                cM = morgan.brca2[which.min(abs(morgan.brca2$position-as.integer(row[3]))),3]
            }
            return(cM)
        }
        # BRCA1
        #morgan.brca1 <- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase2/plink.chr17.GRCh37.map", header = F) # phase2 map
        morgan.brca1 <- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
        coords.brca1 = chr_coords_all %>% filter(Chr_numeric == 17) %>% arrange(position_b37)
        gene <<- "BRCA1"
        coords.brca1$cM <- apply(coords.brca1, 1, getGeneticCoord)
        #write.table(coords.brca1, file = "input/139_ouh_june_2017/brca1_snps_position_and_centiMorgan.txt", col.names = T, quote = F, sep = "\t", row.names = F)
    
        # BRCA2
        #morgan.brca2 <- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase2/plink.chr13.GRCh37.map", header = F) # phase2 map
        morgan.brca2 <- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
        coords.brca2 = chr_coords_all %>% filter(Chr_numeric == 13) %>% arrange(position_b37)
        gene <<- "BRCA2"
        coords.brca2$cM <- apply(coords.brca2, 1, getGeneticCoord)
        #write.table(coords.brca2, file = "input/139_ouh_june_2017/brca2_snps_position_and_centiMorgan.txt", col.names = T, quote = F, sep = "\t", row.names = F)
    
        coords.cm = rbind(coords.brca1, coords.brca2)
        write.table(coords.cm, file = "input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords.txt", col.names = T, quote = F, sep = "\t", row.names = F)
        
    } else {
        ### Linear interpolation approach ###
        
        # BRCA1
        morgan.brca1 <- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
        coords.brca1 = chr_coords_all %>% filter(Chr_numeric == 17) %>% arrange(position_b37)
        coords.brca1$cM <- approx(x = morgan.brca1[,1], y = morgan.brca1[,3], xout = coords.brca1$position_b37)$y
        
        # BRCA2
        morgan.brca2 <- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
        coords.brca2 = chr_coords_all %>% filter(Chr_numeric == 13) %>% arrange(position_b37)
        coords.brca2$cM <- approx(x = morgan.brca2[,1], y = morgan.brca2[,3], xout = coords.brca2$position_b37)$y
        
        coords.cm = rbind(coords.brca1, coords.brca2)
        write.table(coords.cm, file = "input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", col.names = T, quote = F, sep = "\t", row.names = F)
    }
    
    ### Hapmap phase 2 - Linear interpolation
    hapmap2 <- read.table("/Users/lars/Downloads/annotHapMap2U.txt", header = T)
    # BRCA1
    hapmap2.brca1 <- hapmap2 %>% filter(Chrom == "chr17")
    coords.brca1 = chr_coords_all %>% filter(Chr_numeric == 17) %>% arrange(position_b37)
    coords.brca1$cM <- approx(x = hapmap2.brca1$physical_position_build37, y = hapmap2.brca1$deCODE_genetic_map_position, xout = coords.brca1$position_b37)$y
    
    # BRCA2
    hapmap2.brca2 <- hapmap2 %>% filter(Chrom == "chr13")
    coords.brca2 = chr_coords_all %>% filter(Chr_numeric == 13) %>% arrange(position_b37)
    coords.brca2$cM <- approx(x = hapmap2.brca2$physical_position_build37, y = hapmap2.brca2$deCODE_genetic_map_position, xout = coords.brca2$position_b37)$y
    
    coords.cm = rbind(coords.brca1, coords.brca2)
    write.table(coords.cm, file = "input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation_hapmap-phase2.txt", col.names = T, quote = F, sep = "\t", row.names = F)

    ### Hapmap phase 3 - Linear interpolation
    hapmap3 <- read.table("/Users/lars/Downloads/annotHapMap3U.txt", header = T)
    # BRCA1
    hapmap3.brca1 <- hapmap2 %>% filter(Chrom == "chr17")
    coords.brca1 = chr_coords_all %>% filter(Chr_numeric == 17) %>% arrange(position_b37)
    coords.brca1$cM <- approx(x = hapmap3.brca1$physical_position_build37, y = hapmap3.brca1$deCODE_genetic_map_position, xout = coords.brca1$position_b37)$y
    
    # BRCA2
    hapmap3.brca2 <- hapmap3 %>% filter(Chrom == "chr13")
    coords.brca2 = chr_coords_all %>% filter(Chr_numeric == 13) %>% arrange(position_b37)
    coords.brca2$cM <- approx(x = hapmap3.brca2$physical_position_build37, y = hapmap3.brca2$deCODE_genetic_map_position, xout = coords.brca2$position_b37)$y
    
    coords.cm = rbind(coords.brca1, coords.brca2)
    chr_coords_all = coords.cm
    write.table(coords.cm, file = "input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation_hapmap-phase3.txt", col.names = T, quote = F, sep = "\t", row.names = F)
    
    
}

createBRCA_snp_coords_file <- function(){
    # Filter out the coordinate of SNPs needed i.e. on current chromosome (13 or 17)
    brca1_snp_coords <<- filter(chr_coords_all, Chr_numeric == 17) %>%
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(brca1_geno_plink[,-1])) %>%
        distinct(position_b37, .keep_all = T) %>%
        # sort snps according to location
        arrange(position_b37)
    dim(brca1_snp_coords)
    write.table(file = "Scripts/simulate-population/brca1_snp_coords.txt", x = brca1_snp_coords, col.names = F, row.names = F, sep = "\t", quote = F)
    
    # Filter out the coordinate of SNPs needed i.e. on current chromosome (13 or 17)
    brca2_snp_coords <<- filter(chr_coords_all, Chr_numeric == 13) %>%
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(brca2_geno_plink[,-1])) %>%
        distinct(position_b37, .keep_all = T) %>%
        # sort snps according to location
        arrange(position_b37)
    dim(brca2_snp_coords)
    write.table(file = "Scripts/simulate-population/brca2_snp_coords.txt", x = brca2_snp_coords, col.names = F, row.names = F, sep = "\t", quote = F)
}

# This function is used to test if countries affect a pattern or if what see could be several founder mutations
validateNoCountryPattern <- function(){
    # Create table with the same ratio as mutation c.7617+1G>A
    DenSample <- filter(brca2_pheno_merged, Country == "DENMARK", Mut1HGVS != mut) %>% sample_n(93)
    SpanSample <- filter(brca2_pheno_merged, Country == "SPAIN", Mut1HGVS != mut) %>% sample_n(20)
    otherSample <- filter(brca2_pheno_merged, Country != "SPAIN", Country != "DENAMRK", Mut1HGVS != mut) %>% sample_n(5)
    brca_data <- bind_rows(DenSample, SpanSample, otherSample)

    haplotype_mutation(mut, pop_ref, pop_factor = pop_factor)
    haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "SPAIN")
    haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = "DENMARK")
}


countHetAndHomo <- function(matched){
    het = 0
    homo = 0
    missing = 0
    #row=1
    #pos = 2
    het_list = vector()
    homo_list = vector()
    for (row in 1:nrow(matched)){
        het = 0
        homo = 0
        missing = 0
        for (pos in 2:ncol(matched)){
            #if (pos %% 500 == 0) print(paste("Process:", pos, "of", ncol(matched), "SNPs"))

            # Split string into one-character vector
            seq <- strsplit(paste0(matched[row, pos]), split = "")[[1]]

            if (seq[1] == "-" || seq[2] == "-"){
                missing = missing + 1
                next
            }

            # print(paste(seq[1], seq[2]))
            if (seq[1] == seq[2]) {
                homo = homo + 1
            } else {
                het = het + 1
            }


        }
        print(rownames(matched[row,]))
        print(paste("Het:", het))
        print(paste("Homo:", homo))
        print(paste("Missing:", missing))
        print("")
        het_list = c(het_list, het)
        homo_list = c(homo_list, homo)
    }
    print(mean(het_list))
    print(mean(homo_list))
}
#countHetAndHomo(brca1_geno_merged[100:200,])
#countHetAndHomo(brca2_geno_merged[100:200,])
#countHetAndHomo(brca2_geno_merged[100,])

# gene = BRCA1 or BRCA2, area = above or under
# Counts homozygous and heterozygous snps in area of 1.5 Mb from brca gene
countHetAndHomoArea <- function(gene = "BRCA1", area="above"){
    print(paste(gene, area))

    brca_geno_plink = rbind(brca1_geno_plink,brca2_geno_plink)

    if (gene == "BRCA1"){
        #brca_geno_plink = brca1_geno_plink
        print(dim(brca_geno_plink))
        chr_numeric = 17
    } else {
        #brca_geno_plink = brca2_geno_plink
        print(dim(brca_geno_plink))
        chr_numeric = 13
    }

    if (area == "above"){
        coords_area <- filter(chr_coords_all, Chr_numeric == chr_numeric) %>%
            filter(SNP %in% names(brca_geno_plink[,2:ncol(brca_geno_plink)])) %>%
            distinct(position_b37, .keep_all = T) %>% arrange(position_b37) %>%
            filter(position_b37 > brca_middle, position_b37 < (brca_middle+1000000))
        print(paste("snps in area:", nrow(coords_area)))
    } else { # Under
        coords_area <- filter(chr_coords_all, Chr_numeric == chr_numeric) %>%
            filter(SNP %in% names(brca_geno_plink[,2:ncol(brca_geno_plink)])) %>%
            distinct(position_b37, .keep_all = T) %>% arrange(position_b37) %>%
            filter(position_b37 < brca_middle, position_b37 > (brca_middle-1000000))
        print(paste("snps in area:", nrow(coords_area)))
    }

    index = match(coords_area$SNP,names(brca_geno_plink))
    mini_brca <- brca_geno_plink[,index]
    dim(mini_brca)
    # Heterozygous snps
    het_count=sum(apply(mini_brca, 1, function(row) sum(row==1)))
    # Homozygous snps
    hom_count=sum(apply(mini_brca, 1, function(row) sum(row==0)))
    hom_count_minor=sum(apply(mini_brca, 1, function(row) sum(row==2)))

    # Het proportion
    print(paste("Het count:", het_count))
    print(paste("Het proportion:", het_count/(hom_count+het_count+hom_count_minor)))

    # Hom proportion
    print(paste("Hom count:", (hom_count+hom_count_minor)))
    print(paste("Hom proportion:", (hom_count+hom_count_minor)/(hom_count+het_count+hom_count_minor)))
    print(paste("Hom count major:", hom_count))
    print(paste("Hom major proportion:", hom_count/(hom_count+het_count+hom_count_minor)))
    print(paste("Hom count minor:", hom_count_minor))
    print(paste("Hom minor proportion:", hom_count_minor/(hom_count+het_count+hom_count_minor)))

    ## Results for area 1.5 Mb:
    # "BRCA1 above"                                     # "BRCA1 under"
    # 26914 11444                                       # 26914 11444
    # "snps in area: 195"                               # "snps in area: 347"
    # "Het count: 1307046"                              # "Het count: 2495894"
    # "Het proportion: 0.249045106635952"               # "Het proportion: 0.267250430927499"
    # "Hom count: 3941184"                              # "Hom count: 6843264"
    # "Hom proportion: 0.750954893364048"               # "Hom proportion: 0.732749569072501"
    # "Hom count major: 3612914"                        # "Hom count major: 6232683"
    # "Hom major proportion: 0.688406186466675"         # "Hom major proportion: 0.667370977126632"
    # "Hom count minor: 328270"                         # "Hom count minor: 610581"
    # "Hom minor proportion: 0.062548706897373"         # "Hom minor proportion: 0.0653785919458692"

    # "BRCA2 above"                                     # "BRCA2 under"
    # 26914 11444                                       # 26914 11444
    # "snps in area: 519"                               # "snps in area: 661"
    # "Het count: 2298398"                              # "Het count: 2737465"
    # "Het proportion: 0.164543082562413"               # "Het proportion: 0.153875284047569"
    # "Hom count: 11669968"                             # "Hom count: 15052689"
    # "Hom proportion: 0.835456917437587"               # "Hom proportion: 0.846124715952431"
    # "Hom count major: 11084746"                       # "Hom count major: 14434973"
    # "Hom major proportion: 0.793560678464468"         # "Hom major proportion: 0.811402363352223"
    # "Hom count minor: 585222"                         # "Hom count minor: 617716"
    # "Hom minor proportion: 0.0418962389731197"        # "Hom minor proportion: 0.034722352600208"


    ## Results for 1 Mb area
    # "BRCA1 above"                                     # "BRCA1 under"
    # 26914 11444                                       # 26914 11444
    # "snps in area: 147"                               # "snps in area: 227"
    # "Het count: 1030814"                              # "Het count: 1658732"
    # "Het proportion: 0.260546189197236"               # "Het proportion: 0.271501427781555"
    # "Hom count: 2925544"                              # "Hom count: 4450746"
    # "Hom proportion: 0.739453810802764"               # "Hom proportion: 0.728498572218445"
    # "Hom count major: 2657888"                        # "Hom count major: 4055206"
    # "Hom major proportion: 0.671801692364543"         # "Hom major proportion: 0.663756543521394"
    # "Hom count minor: 267656"                         # "Hom count minor: 395540"
    # "Hom minor proportion: 0.0676521184382202"        # "Hom minor proportion: 0.0647420286970507"

    # "BRCA2 above"                                     # "BRCA2 under"
    # 26914 11444                                       # 26914 11444
    # "snps in area: 445"                               # "snps in area: 564"
    # "Het count: 1817680"                              # "Het count: 2010899"
    # "Het proportion: 0.151767636074287"               # "Het proportion: 0.132474688224168"
    # "Hom count: 10159050"                             # "Hom count: 13168597"
    # "Hom proportion: 0.848232363925713"               # "Hom proportion: 0.867525311775832"
    # "Hom count major: 9705684"                        # "Hom count major: 12718421"
    # "Hom major proportion: 0.810378458894874"         # "Hom major proportion: 0.837868464144001"
    # "Hom count minor: 453366"                         # "Hom count minor: 450176"
    # "Hom minor proportion: 0.037853905030839"         # "Hom minor proportion: 0.0296568476318318"


}
# countHetAndHomoArea("BRCA1", "above")
# countHetAndHomoArea("BRCA1", "under")
# countHetAndHomoArea("BRCA2", "above")
# countHetAndHomoArea("BRCA2", "under")

# Find combined family haplotypes for single mutation
findMutFamHaplotypes <- function(mut, gene){
    filename = paste0("cache/famHaplotypes-", mut, ".RData")
    if (file.exists(filename)){
        load(filename)
        return(famHaplotypes)
    }

    prepare_dataframe(mut, gene)

    #fam = "CIMF016397"
    #fam = "CIMF018512"
    #fam = "CIMF016994"
    #fam = "CIMF023250"
    famHaplotypes = data.frame()

    #famMembers_pheno <- filter(brca_data, FamCode == Fam)
    MultiplefamMembers_pheno <- brca_data %>% group_by(FamCode) %>% filter(n()>1)
    print(unique(MultiplefamMembers_pheno$FamCode))
    for (fam in unique(MultiplefamMembers_pheno$FamCode)){
        print(fam)
        famMembers_pheno = filter(MultiplefamMembers_pheno, FamCode == fam)
        index <- fmatch(famMembers_pheno$Onc_ID, matched$SNP)
        famMembers_geno <- matched[index, ]

        #famHaplotype = c(fam, famMembers_pheno$Onc_ID[1])
        #names(famHaplotype) = c("FamCode", "SNP")
        famHaplotype = c(fam)
        names(famHaplotype) = "FamCode"
        for (pos in 2:ncol(matched)){
        #for (pos in 2:7){
            counts <- table(unlist(strsplit(as.character(famMembers_geno[, pos]), split = "")))

            # If multiple candidates and one is '-' (missing), remove it
            if (names(counts[1]) == "-" && length(counts) > 1){
                counts = counts[2:length(counts)]
            }

            candidates = names(which(counts==max(counts)))

            new = paste(candidates, collapse = '')
            #if (nchar(new) == 1) new=paste0(new,new)
            names(new) <- names(matched[pos])
            famHaplotype <- c(famHaplotype, new)

        }
        # paste(famHaplotype[2:length(famHaplotype)], collapse = " ")
        famHaplotypes <- bind_rows(famHaplotypes, famHaplotype)
    }
    print(dim(famHaplotypes))
    print(famHaplotypes[,1:7])

    #write.table(famHaplotypes, file = paste0("cache/famHaplotypes-", mut_name, ".txt"), quote = F, sep = "\t", row.names = F)
    #save(famHaplotypes, file = filename)
    famHaplotypes <<- famHaplotypes
    return(famHaplotypes)
}

# Find combined family haplotypes for all mutations
findAllFamHaplotypes <- function(pheno_data, geno_data, gene){
    filename = paste0("cache/famHaplotypes-", gene, ".RData")
    # if (file.exists(filename)){
    #     load(filename)
    #     return(famHaplotypes)
    # }

    #fam = "CIMF016397"
    #fam = "CIMF018512"
    #fam = "CIMF016994"
    #fam = "CIMF023250"

    #fam = "CIMF000022"

    #famMembers_pheno <- filter(brca_data, FamCode == Fam)
    MultiplefamMembers_pheno <- pheno_data %>% group_by(FamCode, Mut1HGVS) %>% filter(n()>1)
    #MultiplefamMembers_pheno <- MultiplefamMembers_pheno[1:7,]
    print(unique(MultiplefamMembers_pheno$FamCode))
    print(length(unique(MultiplefamMembers_pheno$FamCode)))
    k=0
    start.time = Sys.time()

    #famHaplotypes = data.frame()
    famHaplotypes = matrix("", nrow=nrow(MultiplefamMembers_pheno %>% distinct(FamCode, Mut1HGVS)), ncol=ncol(geno_data)+1)

    for (fam in unique(MultiplefamMembers_pheno$FamCode)){
        k=k+1
        if (k %% 5 == 0){
            print(paste0(k," of ", length(unique(MultiplefamMembers_pheno$FamCode)), " families"))
            print(Sys.time()-start.time)
        }
        print(fam)
        famMembers_pheno = filter(MultiplefamMembers_pheno, FamCode == fam)
        index <- fmatch(famMembers_pheno$Onc_ID, geno_data$SNP)
        famMembers_geno <- geno_data[index, ]

        famHaplotype = c(fam, famMembers_pheno$Onc_ID[1])
        #names(famHaplotype) = c("FamCode", "SNP")
        #famHaplotype = c(fam)
        #names(famHaplotype) = "FamCode"
        for (pos in 2:ncol(geno_data)){
            #for (pos in 2:7){
            counts <- table(unlist(strsplit(as.character(famMembers_geno[, pos]), split = "")))

            # If multiple candidates and one is '-' (missing), remove it
            if (names(counts[1]) == "-" && length(counts) > 1){
                counts = counts[2:length(counts)]
            }

            candidates = names(counts[1])

            best = counts[1]
            # Check if more than one element has the highest frequency
            # If only one element, return that as consensus
            if (length(counts) == 1){
                new = paste0(candidates,candidates)
                names(new) <- names(geno_data[pos])
                famHaplotype <- c(famHaplotype, new)
                next
            }

            # If multiple candidates
            for (i in 2:length(counts)){
                if (counts[i] == best){
                    candidates = c(candidates, names(counts[i]))
                }
                if (counts[i] > best){
                    candidates = names(counts[i])
                    best = counts[i]
                }
                # print(candidates)
            }

            new = paste(candidates, collapse = '')
            if (nchar(new) == 1) new=paste0(new,new)
            names(new) <- names(geno_data[pos])
            famHaplotype <- c(famHaplotype, new)

        }
        # paste(famHaplotype[2:length(famHaplotype)], collapse = " ")
        #famHaplotypes <- bind_rows(famHaplotypes, famHaplotype)
        #famHaplotypes <- rbind(famHaplotypes, famHaplotype)
        famHaplotypes[k,] = famHaplotype
    }
    famHaplotypes = as.data.frame(famHaplotypes)
    names(famHaplotypes) = c("FamCode", "SNP", names(geno_data[2:ncol(geno_data)]))
    print(dim(famHaplotypes))
    print(famHaplotypes[1:10,1:10])

    write.table(famHaplotypes, file = paste0("cache/famHaplotypes-", gene, "-geno.txt"), quote = F, sep = "\t", row.names = F)
    save(famHaplotypes, file = filename)
    #famHaplotypes <<- famHaplotypes
    return(famHaplotypes)
}

# Find combined family haplotypes for all mutations (improved)
findAllFamHaplotypes2 <- function(pheno_data, geno_data, gene, filepath, multiCore=F){
    library(dplyr)
    library(fastmatch)
    computeFamConsensusHaplotype <- function(row, pheno_data, geno_data, MultiplefamMembers_pheno){
        print(row)
        fam = row[1]
        print(fam)
        mut = row[2]
        print(mut)
        famMembers_pheno = filter(MultiplefamMembers_pheno, FamCode == fam, Mut1HGVS == mut)
        famMembers_pheno <<- famMembers_pheno
        index <- fmatch(famMembers_pheno$Onc_ID, geno_data$SNP)
        famMembers_geno <- geno_data[index, ]#[,1:10]

        famHaplotype = c(fam, famMembers_pheno$Onc_ID[1])

        for (pos in 2:ncol(geno_data)){
            #for (pos in 2:7){
            counts <- table(unlist(strsplit(as.character(famMembers_geno[, pos]), split = "")))

            # If multiple candidates and one is '-' (missing), remove it
            if (names(counts[1]) == "-" && length(counts) > 1){
                counts = counts[2:length(counts)]
            }

            candidates = names(which(counts==max(counts)))

            new = paste(candidates, collapse = '')
            if (nchar(new) == 1) new=paste0(new,new)
            names(new) <- names(geno_data[pos])
            famHaplotype <- c(famHaplotype, new)

        }
        # paste(famHaplotype[2:length(famHaplotype)], collapse = " ")
        #famHaplotypes <- bind_rows(famHaplotypes, famHaplotype)
        #famHaplotypes <- rbind(famHaplotypes, famHaplotype)
        #famHaplotypes[k,] <<- famHaplotype
        return(famHaplotype)
    }

    
    print("Creating fam haplotypes...")
    #pheno_data <- brca1_pheno_merged %>% filter(FamCode %in% c("CIMF000851","CIMF000913"))
    MultiplefamMembers_pheno <<- pheno_data %>% group_by(FamCode, Mut1HGVS) %>% filter(n()>1)
    #MultiplefamMembers_pheno <- MultiplefamMembers_pheno[1:7,]
    start.time = Sys.time()

    #famHaplotypes = data.frame()
    #famHaplotypes = matrix("", nrow=nrow(MultiplefamMembers_pheno %>% distinct(FamCode, Mut1HGVS)), ncol=ncol(geno_data)+1)

    all_combinations <- MultiplefamMembers_pheno %>% group_by(FamCode, Mut1HGVS) %>% summarise()

    if (multiCore){
        # Multi-core version
        library(parallel)
        pheno_data <<- pheno_data
        geno_data <<- geno_data
        cl <- makeCluster(detectCores())
        #clusterExport(cl, list("pheno_data", "geno_data", "MultiplefamMembers_pheno"))
        fh <- parRapply(cl, all_combinations, computeFamConsensusHaplotype, pheno_data, geno_data, MultiplefamMembers_pheno)
        dim(fh) <- c(nrow=ncol(geno_data)+1, ncol=nrow(MultiplefamMembers_pheno %>% distinct(FamCode, Mut1HGVS)))
        famHaplotypes = as.data.frame(t(fh))
        stopCluster(cl)
    } else {
        # Single core version
        famHaplotypes = as.data.frame(t(apply(all_combinations, 1, computeFamConsensusHaplotype, pheno_data, geno_data, MultiplefamMembers_pheno)))
    }

    # Attach names to columns
    names(famHaplotypes) = c("FamCode", "SNP", names(geno_data[2:ncol(geno_data)]))
    print(dim(famHaplotypes))
    print(famHaplotypes[1:10,1:10])

    write.table(famHaplotypes, file = paste0(filepath, "famHaplotypes-", gene, "-geno.txt"), quote = F, sep = "\t", row.names = F)
    save(famHaplotypes, file = paste0(filepath, "famHaplotypes-", gene, "-geno.RData"))
    #famHaplotypes <<- famHaplotypes
    print(Sys.time()-start.time)
    return(famHaplotypes)
}

#famHaplotypes=findAllFamHaplotypes2(brca1_pheno_merged, brca1_geno_merged, "BRCA1")
#famHaplotypes=findAllFamHaplotypes2(brca2_pheno_merged, brca2_geno_merged, "BRCA2")
#famHaplotypes = findAllFamHaplotypes2(pheno_data = pujana_pheno, geno_data = pujana_geno, gene = "HMMR", filepath = "38_Pujana_OD7-BRCA1/")

computePopulationFreqs_genotypes <- function(){
    getFreqs <- function(col, n_samples){
        t = setNames(tabulate(col+1), 0:max(col)) / n_samples
        if (length(t) == 1) t=c(t, "1"=0, "2"=0)
        if (length(t) == 2) t=c(t, "2"=0)
        t
    }
    
    genotypes = rbind(brca1_geno_plink, brca2_geno_plink)
    #genotypes = genotypes[1:100,1:10]
    pop_freqs = as.data.frame(t(apply(genotypes[,-1], 2, getFreqs, n_samples=nrow(genotypes))))
    print(pop_freqs)
    
    # BRCA1
    chr_coords_b1 <- filter(chr_coords_all, Chr_numeric == 17) %>%
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(brca1_geno_plink[,-1])) %>%
        distinct(position_b37, .keep_all = T) %>%
        # sort snps according to location
        arrange(position_b37)
    
    pop_freqs_brca1 = pop_freqs[as.character(chr_coords_b1$SNP),]
    write.table(pop_freqs_brca1, file = "cache/pop_freqs_genotypes_brca1_ordered.txt", quote = F, row.names = F, sep = "\t", col.names = T)
    #save(pop_freqs_brca1, file = "cache/pop_freqs_minorAllele_brca1_ordered.RData")
    
    # BRCA2
    chr_coords_b2 <- filter(chr_coords_all, Chr_numeric == 13) %>%
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(brca2_geno_plink[,-1])) %>%
        distinct(position_b37, .keep_all = T) %>%
        # sort snps according to location
        arrange(position_b37)
    
    pop_freqs_brca2 = pop_freqs[as.character(chr_coords_b2$SNP),]
    write.table(pop_freqs_brca2, file = "cache/pop_freqs_genotypes_brca2_ordered.txt", quote = F, row.names = F, sep = "\t", col.names = T)
    #save(pop_freqs_brca2, file = "cache/pop_freqs_minorAllele_brca2_ordered.RData")
    
}

computePopulationFreqs_allele <- function(){
    getFreq <- function(col){
        #print(col)
        counts = table(col)
        if (is.na(counts["1"])){
            ones = 0
        } else {
            ones = counts["1"]
        }
        if (is.na(counts["2"])){
            twos = 0
        } else {
            twos = counts["2"]*2
        }
        minor_freq = (ones + twos) / (sum(counts)*2)
        return(minor_freq)
    }
    
    genotypes = rbind(brca1_geno_plink, brca2_geno_plink)
    #genotypes = genotypes[1:100,1:10]
    pop_freqs = apply(genotypes[,-1], 2, getFreq)
    pop_freqs2 = data.frame(as.list(pop_freqs))
    write.table(pop_freqs2, file = "cache/pop_freqs_minorAllele.txt", quote = F, row.names = F, sep = "\t", col.names = T)
    save(pop_freqs2, file = "cache/pop_freqs_minorAllele.RData")
    
    # BRCA1
    chr_coords_b1 <- filter(chr_coords_all, Chr_numeric == 17) %>%
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(brca1_geno_plink[,-1])) %>%
        distinct(position_b37, .keep_all = T) %>%
        # sort snps according to location
        arrange(position_b37)
    
    pop_freqs_brca1 = select(pop_freqs2, as.character(chr_coords_b1$SNP))
    write.table(pop_freqs_brca1, file = "cache/pop_freqs_minorAllele_brca1_ordered.txt", quote = F, row.names = F, sep = "\t", col.names = T)
    save(pop_freqs_brca1, file = "cache/pop_freqs_minorAllele_brca1_ordered.RData")
    
    # BRCA2
    chr_coords_b2 <- filter(chr_coords_all, Chr_numeric == 13) %>%
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(brca2_geno_plink[,-1])) %>%
        distinct(position_b37, .keep_all = T) %>%
        # sort snps according to location
        arrange(position_b37)
    
    pop_freqs_brca2 = select(pop_freqs2, as.character(chr_coords_b2$SNP))
    write.table(pop_freqs_brca2, file = "cache/pop_freqs_minorAllele_brca2_ordered.txt", quote = F, row.names = F, sep = "\t", col.names = T)
    save(pop_freqs_brca2, file = "cache/pop_freqs_minorAllele_brca2_ordered.RData")
}

# create_hapmap_files()
create_hapmap_files <- function(gene="BRCA1"){
    #load("cache/brca1_geno_merged.RData", .GlobalEnv)
    #load("cache/brca2_geno_merged.RData", .GlobalEnv)
    #brca_geno <- rbind(brca1_geno_merged, brca2_geno_merged)
    #chr_coords <- chr_coords_all 
    
    brca1_geno <- read.table("Prob_PCAScore/PCA/139_ouh_brca1_onco_geno.txt", header = T)
    brca2_geno <- read.table("Prob_PCAScore/PCA/139_ouh_brca2_onco_geno.txt", header = T)
    brca_geno <- rbind(brca1_geno, brca2_geno)
    chr_coords <- read.table("Prob_PCAScore/PCA/2318_snps_coords.txt", header = T)
    
    chr_coords_exist <- filter(chr_coords, SNP %in% names(brca_geno[,-1])) %>% distinct(position_b37, .keep_all = T) %>% arrange(Chr_numeric, position_b37)
    
    brca_geno_matched <- brca_geno[,c(1, fmatch(chr_coords_exist$SNP, names(brca_geno)))]
    #mut = "c.7617+1G>A"; gene = "BRCA2"
    #getBRCAinfo(gene, mut)
    #brca_geno_matched <- extractSamples(brca_data, brca_geno, chr_coords_exist)
    dim(brca_geno_matched)
    
    genotypes = {}
    j = 1
    remove_cols = {}
    for (i in 2:ncol(brca_geno_matched)){
    #for (i in 100:200){
        if      ('AC' %in% brca_geno_matched[,i]) genotypes[j] = "A/C"
        else if ('AG' %in% brca_geno_matched[,i]) genotypes[j] = "A/G"
        else if ('AT' %in% brca_geno_matched[,i]) genotypes[j] = "A/T"
        else if ('CG' %in% brca_geno_matched[,i]) genotypes[j] = "C/G"
        else if ('CT' %in% brca_geno_matched[,i]) genotypes[j] = "C/T"
        else if ('GT' %in% brca_geno_matched[,i]) genotypes[j] = "G/T"
        else if ('AA' %in% brca_geno_matched[,i]) genotypes[j] = "A/A"
        else if ('CC' %in% brca_geno_matched[,i]) genotypes[j] = "C/C"
        else if ('GG' %in% brca_geno_matched[,i]) genotypes[j] = "G/G"
        else if ('TT' %in% brca_geno_matched[,i]) genotypes[j] = "T/T"
        else if ('II' %in% brca_geno_matched[,i]) {
            remove_cols = c(remove_cols, i-1)
            j = j - 1
        }
        else if ('DD' %in% brca_geno_matched[,i]) {
            remove_cols = c(remove_cols, i-1)
            j = j - 1
        }
        else if ('ID' %in% brca_geno_matched[,i]) {
            remove_cols = c(remove_cols, i-1)
            j = j - 1
        }
        else if ('DI' %in% brca_geno_matched[,i]) {
            remove_cols = c(remove_cols, i-1)
            j = j - 1
        }
        else if ('--' %in% brca_geno_matched[,i]) {
            remove_cols = c(remove_cols, i-1)
            j = j - 1
        }
        else {
            print(table(brca_geno_matched[,i]))
            stop()
        }
        
        j = j + 1
        
        #print(table(brca_geno_matched[,i]))
    }
    length(genotypes)
    
    snp_data = t(brca_geno_matched[,-c(1, remove_cols+1)])
    dim(snp_data)
    
    if (length(remove_cols)>0){
        chr_coords_filtered = chr_coords_exist[-remove_cols, ]
    } else {
        chr_coords_filtered = chr_coords_exist
    }
    
    dots = rep('.', nrow(chr_coords_filtered))
    nas = rep(NA, nrow(chr_coords_filtered))
    
    hapmap = cbind(chr_coords_filtered[,1], genotypes, chr_coords_filtered[,2:3], dots, nas, nas, nas, nas, nas, nas, snp_data)
    names(hapmap) = c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID", "QCcode", brca_geno_matched$SNP)
    dim(hapmap)
    hapmap[1:10,1:20]
    
    #write.table(hapmap, file = "input/139_ouh_june_2017/brca1-2_samples_brca1-2_snps_hapmap_format.txt", quote = F,col.names = T, row.names = F, sep = "\t")
    write.table(hapmap, file = "input/139_ouh_june_2017/brca1-2_samples_2318snps.hapmap", quote = F,col.names = T, row.names = F, sep = "\t")
}

# create_ped_and_map_files()
create_ped_and_map_files <- function(matched_single, brca_data.single){
    # Create map and ped files for use with phasing tools like e.g. shapeit2

    # country = "SPAIN"
    # matched_single = filter(matched_single, SNP %in% subset(brca_data, Country == country)$Onc_ID)

    # MAP file
    chr_coords_exist <- filter(chr_coords, SNP %in% names(matched_single[,2:ncol(matched_single)])) %>% distinct(position_b37, .keep_all = T)
    # Chromosome
    chr = rep(brca_chr, nrow(chr_coords_exist))
    # SNP ID
    snp_id = as.character(chr_coords_exist$SNP)
    # SNP genetic position (cM)
    position_cM <- rep(0, nrow(chr_coords_exist))
    # SNP physical position (bp)
    position_bp <- chr_coords_exist$position_b37

    map <- cbind(chr, snp_id, position_cM, position_bp)

    #write.table(map, file = paste0("shapeit2/", gene, "-", mut_name, ".map"), quote = F, col.names = F, sep = "\t", row.names = F)
    write.table(map, file = paste0("shapeit2/", gene, "-", mut_name, "-", country, ".map"), quote = F, col.names = F, sep = "\t", row.names = F)


    # PED file
    split_genotypes <- function(row){
        geno <- unlist(strsplit(unlist(lapply(row, as.character)), split = ""))
        geno[geno=="-"] = 0
        return(geno)
    }
    alleles <- t(apply(matched_single[,fmatch(chr_coords_exist$SNP, names(matched_single))], 1, split_genotypes))
    rownames(alleles) <- matched_single$SNP

    pheno_info <- brca_data.single %>% filter(Onc_ID %in% rownames(alleles)) %>% select(FamCode, Onc_ID, sex)
    #pheno_info <- brca2_pheno_merged %>% filter(Onc_ID %in% rownames(alleles)) %>% select(FamCode, Onc_ID, sex)
    #pheno_info <- brca1_pheno_merged %>% filter(Onc_ID %in% rownames(alleles)) %>% select(FamCode, Onc_ID, sex)
    zeros <- rep(0, nrow(alleles))
    ped <- cbind(pheno_info[,1:2], farther_id=zeros, mom_id=zeros, sex=pheno_info[,3], phenotype=zeros, alleles)

    #write.table(ped, file = paste0("shapeit2/", gene, "-", mut_name, ".ped"), quote = F, col.names = F, sep = "\t", row.names = F)
    write.table(ped, file = paste0("shapeit2/", gene, "-", mut_name, "-", country, ".ped"), quote = F, col.names = F, sep = "\t", row.names = F)
}

# shapeitData()
shapeitData <- function(){
    # Read shapeit data
    #data <- t(read.table(file = "shapeit2/gwas.phased.haps.rehh"))
    data <- t(read.table(file = "shapeit2/BRCA2-c.7617+1G_A-noHapMap.phased.haps.rehh"))
    #data <- t(read.table(file = "shapeit2/BRCA2-c.7617+1G_A-SPAIN-noHapMap.phased.haps.rehh"))
    #data <- t(read.table(file = "shapeit2/BRCA2-c.7617+1G_A-SPAIN.phased.haps.rehh"))

    colnames(data) = names(matched[,2:ncol(matched)])

    snp_names = rep("", nrow(data))
    snp_names[1:nrow(data) %% 2 == 0] = matched$SNP
    snp_names[1:nrow(data) %% 2 == 1] = matched$SNP

    rownames(data) = 1:length(snp_names)

    data <- as.data.frame(cbind(SNP=snp_names, data))

    data[1:10,1:10]

    # Compute
    selectHaplo <- function(data){

        relative = 0
        print(data[1,1])
        index = 0

        while (relative == 0){
            #hap_min = c(0,0)
            #hap_max = c(0,0)
            hap = c(0,0)
            for (i in 1:2){
                # Compare ancestral haplotype to current haplotype
                condition <- (data[i, 2:ncol(data)] == haplotypes[2:length(haplotypes)])
                # Ignore snps where ancestral haplotype is ambigious
                condition[nchar(haplotypes[2:length(haplotypes)]) > 1] = T
                # Split vector according to positive or negative distance to brca gene
                cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
                cond_neg = condition[1:chr_pos.split]

                min_index=which(cond_pos==F)[1+index]
                max_index=which(cond_neg==F)[length(which(cond_neg==F))-index]
                print(min_index)
                print(max_index)

                pos.len = chr_pos.dist[chr_pos.split+min_index,4]
                neg.len = abs(chr_pos.dist[max_index,4])
                #print(pos.len)
                #print(neg.len)

                #hap_min[i] = min_index
                #hap_max[i] = max_index

                hap[i] = pos.len+neg.len
            }

            # if (hap_min[1] > hap_min[2] && hap_max[1] < hap_max[2]){
            #     relative = 1
            # }
            # if (hap_min[1] < hap_min[2] && hap_max[1] > hap_max[2]){
            #     relative = 2
            # }
            if (hap[1] > hap[2]){
                relative = 1
                print("Selected hap 1")
            }
            if (hap[1] < hap[2]){
                relative = 2
                print("Selected hap 2")
            }
            index = index + 1
        }

        return(data[relative,])
    }

    #haplotypes <<- findConsensus(filter(data, SNP %in% subset(brca_data, Country == "DENMARK")$Onc_ID))
    group = 1
    haplotypes <<- findConsensus(filter(data, SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID), cutoff = 0.7)
    #haplotypes <<- findConsensus(filter(data, SNP %in% subset(brca_data, Country == "SPAIN")$Onc_ID), cutoff = 0.7)
    #haplotypes <<- findConsensus(filter(data, SNP %in% subset(brca_data, Country == "DENMARK")$Onc_ID), cutoff = 0.5)

    chr_pos.dist <<- mapSNPsToCoords(chr_coords, data)
    chr_pos.split <<- nrow(filter(chr_pos.dist, brca_distance < 0))

    matched_plink <<- data %>% group_by(SNP) %>% do(selectHaplo(.))

    haplotypes <<- findConsensus(filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID), cutoff = 0.7)

    breaks <<- findHaplotypeBreaks_single_snp(matched_plink, haplotypes)
    chr_pos <<- mapSNPs_plink(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    p <- plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
    print(p)
    ggsave(plot = p, paste0("plots/phased_haplotypes-", brca_name, "-", mut_name, "-group_", group, ".png"), scale = 1.2, width = 16, height = 10)
    p2 <- plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
    print(p2)
    ggsave(plot = p2[[1]], paste0("plots/phased_haplotypes-", brca_name, "-", mut_name, "-group_", group, "-nearest_break-group_order.png"), scale = 1, width = 16, height = 10)
    ggsave(plot = p2[[2]], paste0("plots/phased_haplotypes-", brca_name, "-", mut_name, "-group_", group, "-nearest_break-positive_order.png"), scale = 1, width = 16, height = 10)
    ggsave(plot = p2[[3]], paste0("plots/phased_haplotypes-", brca_name, "-", mut_name, "-group_", group, "-nearest_break-negative_order.png"), scale = 1, width = 16, height = 10)
}

# source("Haplotype-projekt-ver2.R")
# readData()
# famHaplotypes = findAllFamHaplotypes(brca2_pheno_merged, brca2_geno_merged, "BRCA2")

prepare_famHaplotypes_for_simulations <- function(){
    # BRCA1
    mut = "c.68_69delAG"; gene = "BRCA1"
    readIndvAndFamData(gene, mut)
    brca_start <<- 41197695; brca_stop <<- 41276113
    brca_middle <<- (brca_stop - brca_start)/2 + brca_start
    brca_cM_middle <<- morgan.brca1[which.min(abs(morgan.brca1$position-brca_middle)),3]
    
    # Extract BRCA1 coords
    geno_merged = rbind(famHaplotypes_brca1_geno_plink, famHaplotypes_brca2_geno_plink)
    dim(geno_merged)
    index = match(chr_coords$SNP, names(geno_merged))
    geno_merged_brca1 = geno_merged[,c(2,index)]
    dim(geno_merged_brca1)
    
    # Extract one family from each mutation
    pheno_data_all = rbind(brca1_pheno_merged, brca2_pheno_merged) %>% filter(Onc_ID %in% geno_merged_brca1$SNP)
    pheno_data = pheno_data_all %>% distinct(Mut1HGVS, .keep_all = T)
    dim(pheno_data)
    dim(brca1_pheno_merged %>% filter(Onc_ID %in% geno_merged_brca1$SNP) %>% distinct(Mut1HGVS, .keep_all = T))
    dim(brca2_pheno_merged %>% filter(Onc_ID %in% geno_merged_brca1$SNP) %>% distinct(Mut1HGVS, .keep_all = T))
    
    # Keep only the single families in the geno data
    geno_merged_brca1_filtered = filter(geno_merged_brca1, SNP %in% pheno_data$Onc_ID)
    dim(geno_merged_brca1_filtered)
    
    ### FILES CONTAINING 1's
    # Save files to disk
    write.table(geno_merged_brca1_filtered, file = "Scripts/simulate-population/famHaplotypes_simulationBasis_brca1_geno.txt", quote = F, row.names = F, sep = "\t")
    write.table(pheno_data, file = "Scripts/simulate-population/famHaplotypes_simulationBasis_brca1_pheno.txt", quote = F, row.names = F, sep = "\t")
    
    # Convert 1's to 0's or 2's depended on most frequent base in population with that mutation
    #row = geno_merged_brca1_filtered[4,1:10]
    convert_rows <- function(row){
        library(dplyr)
        mut = pheno_data_all$Mut1HGVS[which(pheno_data_all$Onc_ID == as.character(row[1]))]
        pheno_pop = filter(pheno_data_all, Mut1HGVS == mut)
        geno_pop = subset(geno_merged_brca1, SNP %in% pheno_pop$Onc_ID)
        ancestral = findConsensus_plink(geno_pop, cutoff = 0.5)
        res = row
        for (i in 2:length(row)){
            if (row[i] == 1){
                if (ancestral[i-1] != 1){
                    res[i] = ancestral[i-1]
                } else {
                    res[i] = 2
                }
            }
        }
        #print(res)
        return(res)
    }
    cl = makeCluster(7)
    clusterExport(cl = cl, varlist = c("pheno_data_all", "geno_merged_brca1", "findConsensus_plink"))
    geno_converted = t(parApply(cl = cl, geno_merged_brca1_filtered, 1, convert_rows))
    
    # Save files to disk
    write.table(geno_converted, file = "Scripts/simulate-population/famHaplotypes_simulationBasis_brca1_geno_no-het.txt", quote = F, row.names = F, sep = "\t")
    #write.table(pheno_data, file = "Scripts/simulate-population/famHaplotypes_simulationBasis_brca1_pheno.txt", quote = F, row.names = F, sep = "\t")
    
    # BRCA1 snp names
    write.table(geno_merged_brca1_filtered[1,-1], file = "Scripts/simulate-population/brca1_snps_names.txt", quote = F, row.names = F, sep = "\t")
    
    # BRCA1 snp coords
    write.table(chr_coords, file = "Scripts/simulate-population/brca1_snp_coords.txt", quote = F, row.names = F, sep = "\t", col.names = F)
}

prepare_haplotypes_for_simulations <- function(){
    mut = "c.68_69delAG"; gene = "BRCA1"
    readIndvAndFamData(gene, mut)
    
    index <- fmatch(brca1_pheno_merged$Onc_ID, brca1_geno_plink$SNP)
    brca1_geno_sim <- brca1_geno_plink[index, c(1, na.omit(fmatch(chr_coords$SNP, colnames(matched))))]
    dim(brca1_geno_sim)
    
    write.table(brca1_geno_sim, file = "Scripts/simulate-population/haplotypes_brca1_geno_plink.txt", quote = F, row.names = F, sep = "\t")
    write.table(brca1_pheno_merged, file = "Scripts/simulate-population/haplotypes_brca1_pheno_plink.txt", quote = F, row.names = F, sep = "\t")
    
    geno = read.table2("Scripts/simulate-population/haplotypes_brca1_geno_plink.txt", header = T)
    dim(geno)
    pheno = read.table2("Scripts/simulate-population/haplotypes_brca1_pheno_plink.txt", header = T)
    dim(pheno)
}

