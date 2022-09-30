library(readxl)
library(dplyr)
library(forcats)
library(fastmatch)
library(ggplot2)
library(ggpubr)

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
haplotype_mutation <- function(mut, pop_ref = NULL, pop_factor = 1, country = NULL, distinct = F, allCountryConsensus = F, prepare = T, group = NULL, cutoff=0.5){
    "
	Haplotype samples with a given brca mutation
	"
    
    start.time = Sys.time()
    
    if (prepare){
        prepare_dataframe(mut, gene)
    }
    
    if (!is.null(pop_ref)){
        matched2 <- select_rare_snps(matched, pop_ref, pop_factor)
    }
    
    if (!is.null(country)){
        # If haplotype should be based on samples from specific country
        #haplotypes <<- findConsensus(filter(matched, SNP %in% subset(brca_data, Country == country)$Onc_ID))
        haplotypes <<- findConsensus_plink(filter(matched_plink, SNP %in% subset(brca_data, Country == country)$Onc_ID))
    } else if (!is.null(group)) {
        #brca_data <<- brca_data.filtered
        # haplotypes <<- findConsensus(filter(matched, SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID))
        haplotypes <<- findConsensus_plink(filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == group)$Onc_ID), cutoff)
    } else {
        # If the haplotype should be defined across all samples    
        #haplotypes <<- findConsensus(matched)
        haplotypes <<- findConsensus_plink(matched_plink)
    }
    #haplotypes <<- haplotypes
    
    # Find SNPs breaking the haplotype
    # breaks <<- findHaplotypeBreaks(matched, haplotypes)
    breaks <<- findHaplotypeBreaks_plink(matched_plink, haplotypes)
    
    # If no samples breaking haplotype, then stop analysis of current mutation
    if (length(breaks) == 0){
        print(paste("No haplotype breaks found in families for mut:", mut))
        return()
    }
    
    # Map found breaks (SNPs) to their genomic position
    # chr_pos <<- mapSNPs(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    chr_pos <<- mapSNPs_plink(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    
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
    haplotypeDistance <- findNearestBreaks(chr_pos, mut, country)
    
    # Visualize the breaks around the BRCA gene
    if (!is.null(country)){
        p <- plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, country)
        # p2 <- plot_nearest_brca_break_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, country)
    } else {
        p <- plot_haplotype_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
        p2 <- plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
    }
    
    # ggplot(haplotypeDistance, aes(x=sample_id, y=positive_pos)) + geom_bar()
    
    return(p)
    
    # Computing all haplotypes for different ancestral references (countries), higher than 10% of total data size
    if (allCountryConsensus == T){
        quantile = 0.1
        countries <- brca_data %>% group_by(Country) %>% summarise(n=n()/nrow(.)) %>% filter(n > quantile) %>% .$Country
        for (country in countries){
            haplotype_mutation(mut, pop_ref, pop_factor = pop_factor, country = country)
        }
    }
}

prepare_dataframe <- function(mut, gene){
    # Read pheno and geno data
    ## readData()
    
    # Get BRCA genes coordinate information
    getBRCAinfo(gene, mut)
    
    # fixed mut name for saving files on windows
    mut_name <- gsub(">", "_", mut)
    mut_name <<- gsub("\\?", "_", mut_name)
    
    # Computes general statistics of dataset
    computeStatistics()
    
    # Prints for checking output
    print(dim(brca_data)) # 8x67
    print(dim(haplotypes_brca)) # 15679x12468
    print(dim(chr_coords_all)) # 12467x3
    print(mut) # c.1016delA
    
    # Filter out the SNPs needed i.e. on current chromosome
    chr_coords <<- filter(chr_coords_all, Chr_numeric == brca_chr) %>% 
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(haplotypes_brca[,2:ncol(haplotypes_brca)])) %>% 
        distinct(position_b37, .keep_all = T) %>%
        # sort snps according to location
        arrange(position_b37)
    
    # Extract samples with given mutation and remove SNPs not related to actual BRCA gene
    #matched <<- extractSamples(brca_data, haplotypes_brca, chr_coords, distinct)
    matched <<- extractSamples(brca_data, haplotypes_brca, chr_coords, F)
    # matched <- matched[, c(1, na.omit(fmatch(chr_coords$SNP, colnames(matched))))]
}

## Read input data into global variables
readData <- function(){
    # Read phenotype data
    brca1_pheno <- read.csv("input/139_ouh_june_2017/B1_Onco_phenotype_distribution_311215.csv")
    brca2_pheno <- read.csv("input/139_ouh_june_2017/B2_Onco_phenotype_distribution_180816.csv")
    brca1_pheno_hisp <- read.csv("input/139_ouh_june_2017/B1_Hispanic_OncoArray_phenotypes_230817.csv")
    brca2_pheno_hisp <- read.csv("input/139_ouh_june_2017/B2_Hispanic_OncoArray_phenotypes_230817.csv")
    
    # Note there's 20 more columns in brca1_pheno compared to brca1_pheno_hisp.
    dim(brca2_pheno)
    dim(brca2_pheno_hisp)
    colnames(brca2_pheno[,1:47])
    colnames(brca2_pheno_hisp)
    # I think these are not important and therefore removed in the merge step.
    #brca1_pheno_merged <<- rbind(brca1_pheno[, 1:47], brca1_pheno_hisp) %>% arrange(FamCode)
    #brca2_pheno_merged <<- suppressWarnings(rbind(brca2_pheno[, 1:47], brca2_pheno_hisp) %>% arrange(FamCode))
    
    # Hisp is broken, thus they are excluded in the analyses
    brca1_pheno_merged <<- brca1_pheno
    brca2_pheno_merged <<- brca2_pheno
    
    # Read genotype data
    if(!file.exists("cache/brca1_geno_merged.RData") || !file.exists("cache/brca2_geno_merged.RData")){
        # Read and merge genotype data - then save R object
        brca1_geno <- read.table("input/139_ouh_june_2017/139_ouh_brca1_onco_geno.txt", header = T, stringsAsFactors = T)
        brca2_geno <- read.table("input/139_ouh_june_2017/139_ouh_brca2_onco_geno.txt", header = T, stringsAsFactors = T)
        brca1_geno_hisp <- read.table("input/139_ouh_june_2017/139_ouh_female_hisp_brca1_onco_geno.txt", header = T, stringsAsFactors = F)
        brca2_geno_hisp <- read.table("input/139_ouh_june_2017/139_ouh_female_hisp_brca2_onco_geno.txt", header = T, stringsAsFactors = F)

        brca1_geno_merged <- bind_rows(brca1_geno, brca1_geno_hisp) %>% arrange(SNP)
        brca2_geno_merged <- bind_rows(brca2_geno, brca2_geno_hisp) %>% arrange(SNP)
        
        # brca1_geno_merged <- brca1_geno
        # brca2_geno_merged <- brca2_geno
        
        save(brca1_geno_merged, file = "cache/brca1_geno_merged.RData")
        save(brca2_geno_merged, file = "cache/brca2_geno_merged.RData")
    }
    load("cache/brca1_geno_merged.RData", .GlobalEnv)
    load("cache/brca2_geno_merged.RData", .GlobalEnv)
    #load("cache/brca1_geno_merged_noFactors.RData", .GlobalEnv)
    #load("cache/brca2_geno_merged_noFactors.RData", .GlobalEnv)
    
    # Read genotype data in plink format
    if (!file.exists("cache/brca1_geno_plink.RData") || !file.exists("cache/brca2_geno_plink.RData")){
        brca1_geno_plink <- read.table("input/139_ouh_june_2017/139_ouh_brca1_onco_geno_plink_format.txt", header = T)
        save(brca1_geno_plink, file = "cache/brca1_geno_plink.RData")
        brca2_geno_plink <- read.table("input/139_ouh_june_2017/139_ouh_brca2_onco_geno_plink_format.txt", header = T)
        save(brca2_geno_plink, file = "cache/brca2_geno_plink.RData")
    }
    load("cache/brca1_geno_plink.RData", envir = .GlobalEnv)
    load("cache/brca2_geno_plink.RData", envir = .GlobalEnv)
    
    # Read SNP genomic coordinates
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed.txt", header = T, sep = ",")
}

getBRCAinfo <- function(brca_name, mut){
    # BRCA gene information
    if (brca_name == "BRCA1"){
        brca_name <<- "BRCA1"
        brca_start <<- 41197695; brca_stop <<- 41276113
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_chr <<- 17
        
        # Legacy names
        haplotypes_brca <<- brca1_geno_merged
        brca_data <<- filter(brca1_pheno_merged, Mut1HGVS %in% mut)
    } else { # BRCA2
        brca_name <<- "BRCA2"
        brca_start <<- 32890598; brca_stop <<- 32972907
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_chr <<- 13
        
        # Legacy names
        haplotypes_brca <<- brca2_geno_merged
        brca_data <<- filter(brca2_pheno_merged, Mut1HGVS %in% mut)
    }
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
    #MutPatientCount <- count(brca1_pheno_merged, Mut1HGVS)
    #MutPatientCount <- count(brca2_pheno_merged, Mut1HGVS)
    MutPatientCount <<- count(brca_data, Mut1HGVS)
    
    # Count the number of family members 
    #FamMemberCount <- brca1_pheno_merged %>% distinct(Onc_ID, FamCode) %>% count(FamCode) %>% count(FamMembers=n)
    #FamMemberCount <- brca2_pheno_merged %>% distinct(Onc_ID, FamCode) %>% count(FamCode) %>% count(FamMembers=n)
    FamMemberCount <<- brca_data %>% distinct(Onc_ID, FamCode) %>% count(FamCode) %>% count(FamMembers=n)
    
    # Count number of families with at least two members
    brca_data %>% count(FamCode) %>% filter(n>1)
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

## ---- extractSamples
extractSamples <- function(brca_data, haplotypes_brca, chr_coords, distinct = F){
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
    index <- fmatch(brca_data$Onc_ID, haplotypes_brca$SNP)
    matched <- haplotypes_brca[index, ]
    # as above, but slower
    # matched <- filter(haplotypes_brca, SNP %in% brca_data$Onc_ID)
    print(dim(matched))
    
	# Only keep SNPs (columns) that are related to current chromosome
    matched <- matched[, c(1, na.omit(fmatch(chr_coords$SNP, colnames(matched))))]
	print(dim(matched))
	
    return(matched)
}

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

## ---- findConsensus
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

findHaplotypeBreaks_plink <- function(matched_plink, consensus){
    #c <- c(1:9)
    #breaks <- t(apply(matched_plink[,2:ncol(matched_plink)], 1, function(row) abs(row-consensus)==2))
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

### Map found SNPs to their genomic position ###
mapSNPs_plink <- function(breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name){
    #mapSNPs <- function(){
    chr_pos = NULL
    # For each sample (patient) with given mutation
    for (i in 1:nrow(breaks)){
        # If no breaks for sample, then skip and continue
        if (length(names(which(breaks[i,]==T))) == 0) next
        # Extract rows with SNPs in brca1 chromosome (17) and where SNP is homozygote for that sample
        #chr_pos2 = filter(chr_coords, chr_coords$Chr_numeric == brca_chr, chr_coords$SNP %in% breaks[[i]])
        chr_pos2 = filter(chr_coords, chr_coords$SNP %in% names(which(breaks[i,]==T)))
        # print(dim(chr_pos2))
        # Removes SNPs in the brca1 gene
        chr_pos2 = filter(chr_pos2, as.numeric(as.character(position_b37)) < brca_start | as.numeric(as.character(position_b37)) > brca_stop)
        # Add sample_id column
        chr_pos2["sample_id"] = matched_plink[i,1]
        # Add position_b37_adjusted column, where brca1 = 0
        chr_pos2 = mutate(chr_pos2, position_b37_adjusted = as.numeric(as.character(position_b37)) - brca_middle)
        # Combine all samples into big data frame for plotting (likely not most efficient way)
        chr_pos = bind_rows(chr_pos, chr_pos2)
    }
    # Add columns Country and FamCode to data frame
    index <- fmatch(chr_pos$sample_id, brca_data$Onc_ID)
    # chr_pos[c("Country", "FamCode")] <- brca_data[index, ] %>% select(Country, FamCode)
    chr_pos[c("Country", "FamCode", "cluster_groups")] <- brca_data[index, ] %>% select(Country, FamCode, cluster_groups)
    
    # Add fake-break at max and min dist for samples with no breaks, so they show in plots
    for (sample in unique(matched_plink$SNP)[!(unique(matched_plink$SNP) %in% unique(chr_pos$sample_id))]){
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
        # Removes SNPs in the brca1 gene
        chr_pos2 = filter(chr_pos2, as.numeric(as.character(position_b37)) < brca_start | as.numeric(as.character(position_b37)) > brca_stop)
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

## ---- plot
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
    
    
    
    # Get number of colors according to number of groups
    n_palette_col = if (palette == "Set1") 9 else 8 # Dark2=8, Set1=9
    if (length(country_freq) <= n_palette_col){
        cols <- brewer.pal(n = max(length(country_freq), 3), name = palette)[1:(length(country_freq))]
    } else {
        cols <- colorRampPalette(brewer.pal(n = n_palette_col, name = palette))(length(country_freq))
    }
    
    #cols <- rainbow_hcl(length(country_freq))
    
    # Assign color to country or group
    if (is.null(color_names)){
        names(cols) <- names(country_freq)
    } else {
        names(cols) <- color_names
    }
    
    # Make it global (for convenience) and return
    cols <<- cols
    return(cols)
}

### Visual output of haplotypes and breaks ###
plot_haplotype <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, country = F){

    # Define y axis step distance
    tick_steps = 3000000
    # Set the adjusted beginning and ending position of brca gene 
    brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
    # Number of samples
    num_samples <<- length(unique(chr_pos$sample_id))
    
    cols = colors_plot(chr_pos, otherCountry = T)
    
    chr_pos2 = chr_pos %>% mutate(Country = fct_lump(Country, prop = 0.05))
    #lab_cols = cols[as.character(distinct(chr_pos, sample_id, Country)$Country)]
    #lab_cols = cols[lab_cols[lab_cols != "Other"]]
    # Plot countries together that has more than 5% samples of total datasets
    p = chr_pos2 %>%
    #p = chr_pos %>% mutate(Country = fct_lump(Country, prop = 0.05)) %>% # group_by(Country) %>%
        
        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        ggplot(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_b37_adjusted)) +
        #ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_b37_adjusted)) +
        
        # Plots the breakpoints around the brca gene
        geom_point(aes(color = Country), size = 0.8) +
        #geom_point(aes(color = cols[Country]), size = 0.8) +
        #geom_point(aes(color = as.character(cluster_groups)), size = 0.8) +
        #geom_point(size = 0.8, color = "darkblue") +
        
        # Defines the color of the dots
        #scale_color_brewer(palette="Set1") +
        scale_color_manual(values = cols) + 
        
        # Plots the area of the brca gene
        #geom_hline(yintercept = 0, color = "red", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.4) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.4) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +
        
        #geom_hline(aes(yintercept = 5*tick_steps), color = "white", size = 0.2) +
        #geom_hline(aes(yintercept = -5*tick_steps), color = "white", size = 0.2) +
        
        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", mut)) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title
        
        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols)) +  # vertical x-labels
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = distinct(chr_pos2, sample_id, Country)$Country)) +  # vertical x-labels
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = lab_cols)) +  # vertical x-labels
        
        # Define x-ticks and their label names
        scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country), labels = chr_pos$FamCode) +
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country)) +

        # Define y-ticks and their label names (sets BRCA name in y axis)
        scale_y_continuous(breaks = c(seq(from = tick_steps, to = 5*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -5*tick_steps, by = -tick_steps), 0),
                           labels = c(seq(from = tick_steps, to = 5*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -5*tick_steps, by = -tick_steps), brca_name)) +
        
        
        # Define visible y-axis
        #ylim(-15000000,15000000) + 
        
        # Set x and y axis label names
        ylab("Genomic position") + xlab("Sample")
    
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

### Produces three plots, ordered by country, by positve genomic points and by negative points - all breakpoints coloured by country.
plot_nearest_brca_break_country <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, country){
    
    # Extract (positive position) breakpoint closest to brca gene
    positive <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted > 0) %>%
        slice(which.min(position_b37_adjusted)) %>% 
        rename(position_pos = position_b37, position_pos_adjusted = position_b37_adjusted, snp_pos = SNP)
    # Extract (negative position) breakpoint closest to brca gene
    negative <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted < 0) %>%
        slice(which.max(position_b37_adjusted)) %>%
        rename(position_neg = position_b37, position_neg_adjusted = position_b37_adjusted, snp_neg = SNP)
    
    # chr_pos_minmax <- rbind(positive, negative)
    chr_pos_minmax <- merge(positive, negative, by = intersect(names(positive), names(negative)))
    
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
    
    #cols = colors_plot(brca_data, otherCountry = T)
    cols = colors_plot(chr_pos_minmax, otherCountry = T)
    
    # Plot countries together that has less than 5% samples of total datasets
    #p = chr_pos_minmax %>% #mutate(Country = fct_lump(Country, prop = 0.1)) %>% # group_by(Country) %>%
    p = chr_pos_minmax %>% mutate(Country = fct_lump(Country, prop = 0.05)) %>% #group_by(Country) %>% 
        
        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        ggplot(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_pos_adjusted)) +
        
        # Plots the breakpoints around the brca gene
        geom_point(aes(color = Country), size = point_size) +
        
        # OBS: Negative points need special care, depending on plot needed. See printing of plots below.
        #geom_point(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_neg_adjusted, color = Country), size = 0.8) +
        
        # Defines the color of the dots
        #scale_color_brewer(palette="Set1") +
        scale_color_manual(values = cols) + 
        
        # Plots the area of the brca gene
        #geom_hline(yintercept = 0, color = "darkblue", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.1) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.1) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +
        
        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", paste(mut, collapse = ","))) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title
        
        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)] )) +  # vertical x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols)) +  # vertical x-labels
        
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
        ylab("Genomic position") + xlab("Sample")
    
    # Decide consensus for filename 
    #consensus <- if (is.null(country)) "All_consensus" else paste0(country, "_consensus")
    
    # Negative points added when sorting by country
    negative_positions <- geom_point(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_neg_adjusted, 
                                         color = Country), size = point_size)
    
    # Order points by country
    p1 = p + negative_positions +
        scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$Country), labels = chr_pos_minmax$FamCode)
    p1
    
    # Save plots to SVG and PNG files
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-country-", num_samples, "_samples.png"))
    
    
    # Negative points added when sorting by genomic position
    negative_positions <- geom_point(aes(x=reorder(sample_id, position_pos_adjusted), y=position_neg_adjusted, 
                                         color = Country), size = point_size)
    # negative_positions <- geom_point(aes(x=sample_id, y=position_neg_adjusted, 
    #                                      color = Country), size = point_size)
    
    # Order points by positive genomics positions
    p2 = p + aes(x=reorder(sample_id, position_pos_adjusted)) + negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$position_pos_adjusted), labels = chr_pos_minmax$FamCode) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[as.character(distinct(chr_pos_minmax2, sample_id, Country)$Country)]))
    p2
    
    # Order points by positive genomics positions
    #print(p + aes(x=reorder(sample_id, position_pos_adjusted)) + negative_positions + xlab("Sample"))
    # Save plots to SVG and PNG files
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-positive_order-", num_samples, "_samples.png"))
    
    
    # Order points by negative genomics positions
    p3 = p + aes(x=reorder(sample_id, -position_neg_adjusted)) + negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$position_pos_adjusted), labels = chr_pos_minmax$FamCode) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = cols[as.character(distinct(chr_pos_minmax2, sample_id, Country)$Country)]))
    p3
    
    
    # Order points by negative genomics positions
    #print(p + aes(x=reorder(sample_id, -position_neg_adjusted)) + negative_positions + xlab("Sample"))
    # Save plots to SVG and PNG files
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-negative_order-", num_samples, "_samples.png"))
    
    return(list(p1,p2,p3))
}

### Visual output of haplotypes and breaks ###
plot_haplotype_group <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group){
    
    # Define y axis step distance
    tick_steps = 3000000
    # Set the adjusted beginning and ending position of brca gene 
    brca_max = (brca_stop-brca_start)/2; brca_min = -brca_max
    # Number of samples
    num_samples <<- length(unique(chr_pos$sample_id))
    
    # Set up colors for plot
    label_cols <- colors_plot(palette = "Dark2")[as.character(distinct(chr_pos, sample_id, Country)$Country)]
    cols_groups <- colors_plot(palette = "Set1", color_names = names(sort(table(cluster_groups), decreasing = T)) )
    #cols_groups <- colors_plot(palette = "Set1", color_names = 1:max(cluster_groups))
    
    # Plot countries together that has less than 5% samples of total datasets
    p = chr_pos %>% #mutate(Country = fct_lump(Country, prop = 0.05)) %>% # group_by(Country) %>%
        
        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        #ggplot(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_b37_adjusted)) +
        ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_b37_adjusted)) +
        #ggplot(aes(x=interaction(FamCode, cluster_groups), y=position_b37_adjusted)) +
        
        # Plots the breakpoints around the brca gene
        #geom_point(aes(color = Country), size = 0.8) +
        geom_point(aes(color = as.character(cluster_groups)), size = 0.8) +
        #geom_point(size = 0.8, color = "darkblue") +
        
        # Defines the color of the dots
        #scale_color_brewer(palette="Set1") +
        scale_color_manual(values = cols_groups) + 
        
        # Plots the area of the brca gene
        #geom_hline(yintercept = 0, color = "red", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.4) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.4) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +
        
        #geom_hline(aes(yintercept = 5*tick_steps), color = "white", size = 0.2) +
        #geom_hline(aes(yintercept = -5*tick_steps), color = "white", size = 0.2) +
        
        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", mut, "- Group consensus:", group)) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title
        
        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = label_cols)) +  # vertical x-labels
        
        # Define x-ticks and their label names
        scale_x_discrete(breaks = interaction(chr_pos$sample_id, chr_pos$cluster_groups), 
                         labels = chr_pos$FamCode) +
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country)) +
        
        # Define y-ticks and their label names (sets BRCA name in y axis)
        scale_y_continuous(breaks = c(seq(from = tick_steps, to = 5*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -5*tick_steps, by = -tick_steps), 0),
                           labels = c(seq(from = tick_steps, to = 5*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -5*tick_steps, by = -tick_steps), brca_name)) +
        
        
        # Define visible y-axis
        #ylim(-15000000,15000000) + 
        
        # Legend title
        labs(colour="Groups") +
        
        # Set x and y axis label names
        ylab("Genomic position") + xlab("Sample")
    
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

plot_nearest_brca_break_group <- function(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group){
    
    # Extract (positive position) breakpoint closest to brca gene
    positive <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted > 0) %>%
        slice(which.min(position_b37_adjusted)) %>% 
        rename(position_pos = position_b37, position_pos_adjusted = position_b37_adjusted, snp_pos = SNP)
    # Extract (negative position) breakpoint closest to brca gene
    negative <- chr_pos %>%
        group_by(sample_id) %>%
        filter(position_b37_adjusted < 0) %>%
        slice(which.max(position_b37_adjusted)) %>%
        rename(position_neg = position_b37, position_neg_adjusted = position_b37_adjusted, snp_neg = SNP)
    
    # chr_pos_minmax <- rbind(positive, negative)
    chr_pos_minmax <- merge(positive, negative, by = intersect(names(positive), names(negative))) %>% mutate(haplo_length = position_pos-position_neg)
    
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
    
    # Set up colors for plot
    label_cols <- colors_plot(palette = "Dark2")[as.character(distinct(chr_pos_minmax, sample_id, Country)$Country)]
    cols_groups <- colors_plot(palette = "Set1", color_names = names(sort(table(cluster_groups), decreasing = T)) )
    
    # Plot countries together that has less than 5% samples of total datasets
    #p = chr_pos_minmax %>% #mutate(Country = fct_lump(Country, prop = 0.1)) %>% # group_by(Country) %>%
    p = chr_pos_minmax %>% # mutate(Country = fct_lump(Country, prop = 0.05)) %>% #group_by(Country) %>% 
        
        # Set up the ggplot, with x sorted by most frequent country, and sample_id as tie breaker
        # The genomic position defines the y axis
        ggplot(aes(x=interaction(sample_id, cluster_groups), y=position_pos_adjusted)) +
        
        # Plots the breakpoints around the brca gene
        geom_point(aes(color = as.character(cluster_groups)), size = point_size) +
        #geom_point(aes(color = Country), size = point_size) +
        
        # OBS: Negative points need special care, depending on plot needed. See printing of plots below.
        #geom_point(aes(x=interaction(sample_id, fct_infreq(Country)), y=position_neg_adjusted, color = Country), size = 0.8) +
        
        
        # Defines the color of the dots
        # scale_color_brewer(palette="Set1") +
        scale_color_manual(values = cols_groups) + 
        
        # Plots the area of the brca gene
        #geom_hline(yintercept = 0, color = "darkblue", size = 1) +
        geom_hline(aes(yintercept = brca_max), color = "darkblue", size = 0.1) +
        geom_hline(aes(yintercept = brca_min), color = "darkblue", size = 0.1) +
        geom_rect(aes(xmin = 0, xmax = num_samples+1, ymin=brca_min, ymax=brca_max), fill="darkblue") +
        
        # Inserts the plot title and centers it
        ggtitle(paste("Mutation:", mut, "- Group consensus:", group)) + # Adds plot title
        theme(plot.title = element_text(hjust = 0.5)) + # Centering of plot title
        
        # Change the angle of the x-labels
        # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # 45 degree x-labels
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = label_cols)) +  # vertical x-labels
        
        # Define x-ticks and their label names
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country), labels = chr_pos$sample_id) +
        #scale_x_discrete(breaks=interaction(chr_pos$sample_id, chr_pos$Country)) +
        #scale_x_discrete(breaks = interaction(chr_pos$sample_id, chr_pos$cluster_groups), 
        #                 labels = chr_pos$FamCode) +
        
        # Define y-ticks and their label names (sets BRCA name in y axis)
        scale_y_continuous(breaks = c(seq(from = tick_steps, to = 15*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -15*tick_steps, by = -tick_steps), 0),
                           labels = c(seq(from = tick_steps, to = 15*tick_steps, by = tick_steps),
                                      seq(from = -tick_steps, to = -15*tick_steps, by = -tick_steps), brca_name)) +
        
        
        # Define visible y-axis
        #ylim(-1000000,1000000) +
        
        # Legend title
        labs(colour="Groups") +
        
        # Set x and y axis label names
        ylab("Genomic position") + xlab("Sample")
    
    # Decide consensus for filename 
    #consensus <- if (is.null(country)) "All_consensus" else paste0(country, "_consensus")
    
    # nearest_break = "NearestBreak_country"
    # 
    # filename3 = paste0("cache/nearestBreakpointPlots/", gene, "-", mut, "-", nearest_break, "-", 
    #                    "distMethod-", input$distMethod, "-clustMethod-", input$clustMethod, "-",
    #                    fam, "-k_", k, "-group_", i, ".png")
    
    # Negative points added when sorting by country
    negative_positions <- geom_point(aes(x=interaction(sample_id, cluster_groups),
                                        y=position_neg_adjusted, color = as.character(cluster_groups)), size = point_size)
    # negative_positions <- geom_point(aes(x=interaction(sample_id, cluster_groups), 
    #                                      y=position_neg_adjusted, color = Country), size = point_size)
    # Order points by country
    p1 = p + negative_positions + 
         scale_x_discrete(breaks = interaction(chr_pos_minmax$sample_id, chr_pos_minmax$cluster_groups), labels = chr_pos_minmax$FamCode)
    #print(p1)
    # Save plots to PNG files
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-country-", num_samples, "_samples.png"))
    
    
    # Negative points added when sorting by genomic position
    negative_positions <- geom_point(aes(x=reorder(sample_id, position_pos_adjusted), y=position_neg_adjusted, 
                                         color = as.character(cluster_groups)), size = point_size)
    # negative_positions <- geom_point(aes(x=sample_id, y=position_neg_adjusted, 
    #                                      color = Country), size = point_size)
    
    # Order points by positive genomics positions
    p2 = p + 
        aes(x=reorder(sample_id, position_pos_adjusted)) + 
        negative_positions + xlab("Sample") + 
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, chr_pos_minmax$position_pos_adjusted), labels = chr_pos_minmax$FamCode)
    #print(p2)
    
    # Save plots to PNG files
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-positive_order-", num_samples, "_samples.png"))
    
    # Order points by negative genomics positions
    p3 = p + aes(x=reorder(sample_id, -position_neg_adjusted)) + negative_positions + xlab("Sample") +
        scale_x_discrete(breaks = reorder(chr_pos_minmax$sample_id, -chr_pos_minmax$position_neg_adjusted), labels = chr_pos_minmax$FamCode)
    #print(p3)
    # Save plots PNG files
    #ggsave(paste0("plots/", brca_name, "-", mut_name, "-", consensus, "-", "nearest_break-negative_order-", num_samples, "_samples.png"))
    
    p4 = p + aes(x=reorder(sample_id, haplo_length)) + negative_positions + xlab("Sample")
    
    return(list(p1,p2,p3,p4))
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
            #mutate(mean_haplo_distance = round(mean_haplo_distance, digits = 1)) %>%
            arrange(desc(mean_haplo_distance))
        
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

# Find combined family haplotypes for single mutation
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

findAllFamHaplotypes2 <- function(pheno_data, geno_data, gene){
    computeFamConsensusHaplotype <- function(row){
        library(dplyr)
        library(fastmatch)
        print(row)
        fam = row[1]
        print(fam)
        mut = row[2]
        print(mut)
        famMembers_pheno = filter(MultiplefamMembers_pheno, FamCode == fam, Mut1HGVS == mut)
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
    
    filename = paste0("cache/famHaplotypes-", gene, ".RData")
    # if (file.exists(filename)){
    #     load(filename)
    #     return(famHaplotypes)
    # }
    
    print("Creating fam haplotypes...")
    #pheno_data <- brca1_pheno_merged %>% filter(FamCode %in% c("CIMF000851","CIMF000913"))
    MultiplefamMembers_pheno <<- pheno_data %>% group_by(FamCode, Mut1HGVS) %>% filter(n()>1)
    #MultiplefamMembers_pheno <- MultiplefamMembers_pheno[1:7,]
    start.time = Sys.time()
    
    #famHaplotypes = data.frame()
    #famHaplotypes = matrix("", nrow=nrow(MultiplefamMembers_pheno %>% distinct(FamCode, Mut1HGVS)), ncol=ncol(geno_data)+1)
    
    all_combinations <- MultiplefamMembers_pheno %>% group_by(FamCode, Mut1HGVS) %>% summarise()
    
    # Multi-core version
    library(parallel)
    pheno_data <<- pheno_data
    geno_data <<- geno_data
    cl <- makeCluster(detectCores())
    clusterExport(cl, list("pheno_data", "geno_data", "MultiplefamMembers_pheno"))
    fh <- parRapply(cl, all_combinations, computeFamConsensusHaplotype)
    dim(fh) <- c(nrow=ncol(geno_data)+1, ncol=nrow(MultiplefamMembers_pheno %>% distinct(FamCode, Mut1HGVS)))
    famHaplotypes = as.data.frame(t(fh))
    stopCluster(cl)
    
    # Single cores
    #famHaplotypes = as.data.frame(t(apply(all_combinations, 1, computeFamConsensusHaplotype)))
    
    # Attach names to columns
    names(famHaplotypes) = c("FamCode", "SNP", names(geno_data[2:ncol(geno_data)]))
    print(dim(famHaplotypes))
    print(famHaplotypes[1:10,1:10])
    
    write.table(famHaplotypes, file = paste0("cache/famHaplotypes-", gene, "-geno.txt"), quote = F, sep = "\t", row.names = F)
    save(famHaplotypes, file = filename)
    #famHaplotypes <<- famHaplotypes
    print(Sys.time()-start.time)
    return(famHaplotypes)
}

#famHaplotypes=findAllFamHaplotypes(brca1_pheno_merged, brca1_geno_merged, "BRCA1")
#famHaplotypes=findAllFamHaplotypes(brca2_pheno_merged, brca2_geno_merged, "BRCA2")

# create_ped_and_map_files()
create_ped_and_map_files <- function(){
    # Create map and ped files for use with phasing tools like e.g. shapeit2
    
    country = "SPAIN"
    matched_single = filter(matched_single, SNP %in% subset(brca_data, Country == country)$Onc_ID)
    
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
    p <- plot_haplotype_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
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