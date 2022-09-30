#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)<12){
    stop("Should specify 12 arguments in the form: Rscript --vanilla prepare_simulated_data.R n_samples n_simulations gene gen_start gen_stop gen_step seed genealogy workdir option ancestral_method cs_correction")
} else {
    #string_id = args[1]
    n_samples = as.integer(args[1])
    n_simulations = as.integer(args[2])
    gene = args[3]
    gen_start = as.integer(args[4])
    gen_stop = as.integer(args[5])
    gen_step = as.integer(args[6])
    seed = as.integer(args[7])
    genealogy = args[8]
    setwd(args[9])
    option = args[10]
    ancestral_method = args[11]
    cs_correction = args[12]
}

source("../../Scripts/Haplotype-projekt-ver2.R")
source("../../Scripts/Clustering.R")
## Read Data
chr_coords_all <<- read.table("../../input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
morgan.brca1 <<- read.table("../../input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
morgan.brca2 <<- read.table("../../input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
alleleFreqs_brca1 <<- as.numeric(read.table2("../../cache/pop_freqs_minorAllele_brca1_ordered.txt", header = T))

print(c(n_samples, n_simulations, gene, gen_start, gen_stop, gen_step, seed, genealogy, args[9], option, ancestral_method, cs_correction))

dir.create("classification_data/", showWarnings = F)

writeSimulatedDistMatrices <- function(samples, simulations, in_gene, gen_start, gen_stop, gen_step, seed, genealogy){
    dists = list()
    l_index=1
    for (generations in seq(gen_start,gen_stop,gen_step)){
    #for (generations in seq(20,500,20)){
        #samples = 100; simulations = 10#; generations = 10
        filename_prefix = paste0("simulated_population/",in_gene,"_simulated-",genealogy,
                                 "-",samples,"_samples-",simulations,"_simulations-","generations_",
                                 generations,"-seed_",seed)
        brca1_simulated_pheno_merged <<- as.data.frame(data.table::fread(file = paste0(filename_prefix, "-pheno.txt"), header = T))
        brca1_simulated_geno_plink <<- as.data.frame(data.table::fread(file = paste0(filename_prefix, "-geno.txt"), header = T))

        for (i in 1:simulations){
            #for (i in 1:1){
            mut = paste0("mutSim_BRCA1_",generations,"_",i); gene = "BRCA1_simulated"

            mut <<- mut; gene <<- gene
            ## Prepare data frames
            prepare_dataframe(mut, gene)

            dst <- firstBreakDist(fam=F)
            dists[[l_index]] = as.vector(dst)
            l_index = l_index + 1
        }
        
    }
    dists_df = do.call(rbind, dists)
    
    # write training data
    write.table(x = dists_df, file = paste0("classification_data/Distances-",gene,"-",
                                            genealogy,"-",samples,"_samples-",simulations,"_simulations-generations_",
                                            gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,".txt"),
                                            quote = F, sep = "\t", row.names = F, col.names = F)
    # write results
    age_results = rep(seq(gen_start,gen_stop,gen_step), each = simulations)
    write.table(x = age_results, file = paste0("classification_data/results-generations_",gen_start,"_",gen_stop,"_step",gen_step,".txt"), row.names = F, col.names = F)
}

writeSimulatedHaploLength <- function(samples, simulations, in_gene, gen_start, gen_stop, gen_step, seed, genealogy){
    # Writes length of haplotypes in cM to simulated founder
    haplo_length_cM = list()
    l_index=1
    for (generations in seq(gen_start,gen_stop,gen_step)){
    #for (generations in seq(20,500,20)){
        #samples = 100; simulations = 10#; generations = 10
        filename_prefix = paste0("simulated_population/",in_gene,"_simulated-",genealogy,
                                 "-",samples,"_samples-",simulations,"_simulations-","generations_",
                                 generations,"-seed_",seed)
        brca1_simulated_pheno_merged <<- as.data.frame(data.table::fread(file = paste0(filename_prefix, "-pheno.txt"), header = T))
        brca1_simulated_geno_plink <<- as.data.frame(data.table::fread(file = paste0(filename_prefix, "-geno.txt"), header = T))
        brca1_simulated_founder <<- as.data.frame(data.table::fread(file = paste0(filename_prefix, "-founder.txt"), header = T))[-1]

        for (i in 1:simulations){
            mut = paste0("mutSim_BRCA1_",generations,"_",i); gene = "BRCA1_simulated"

            mut <<- mut; gene <<- gene
            ## Prepare data frames
            prepare_dataframe(mut, gene)

            dst <- firstBreakDist(fam=F)

            hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = 1, doHapPlot = F)
            #matched_plink2 <- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == age_group)$Onc_ID)
            haplotypes <- findConsensus_plink(filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID))
            breaks <- findHaplotypeBreaks_plink(matched_plink, haplotypes)
            #names(breaks) = names(brca1_simulated_founder)
            chr_pos <- mapNearestSNPs_plink(matched_plink, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)

            haplo_length_cM[[l_index]] = chr_pos$haplo_length_cM
            l_index = l_index + 1
        }
    }
    hap_lengths_df = do.call(rbind, haplo_length_cM)
    # write training data
    write.table(x = hap_lengths_df, file = paste0("classification_data/Haplo_lengths_cM-",gene,"-",
                                            genealogy,"-",samples,"_samples-",simulations,"_simulations-generations_",
                                            gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,".txt"),
                                            quote = F, sep = "\t", row.names = F, col.names = F)
    
    # write results
    age_results = rep(seq(gen_start,gen_stop,gen_step), each = simulations)
    write.table(x = age_results, file = paste0("classification_data/results-generations_",gen_start,"_",gen_stop,"_step",gen_step,".txt"), row.names = F, col.names = F)
}

average_dendrogram_heights <- function(){

    heights_mean = {}
    g_index = 1
    generations <- seq(10,100,10)
    #generations <- seq(10,20,10)
    for (generation in generations){
        print(paste("Generation:", generation))
        l_index=1
        heights = {}
        simulations = 10#; generations = 10
        brca1_simulated_pheno_merged <<- read.table(file = paste0("Scripts/simulate-population/BRCA1_simulated_starGenealogy_200-samples_",simulations,"-muts_",generation,"-generations_pheno.txt"), header = T)
        brca1_simulated_geno_plink <<- read.table(file = paste0("Scripts/simulate-population/BRCA1_simulated_starGenealogy_200-samples_",simulations,"-muts_",generation,"-generations_geno.txt"), header = T)

        for (i in 1:simulations){
            #for (i in 1:2){
            mut = paste0("mutSim_BRCA1_",i); gene = "BRCA1_simulated"

            mut <<- mut; gene <<- gene
            ## Prepare data frames
            prepare_dataframe(mut, gene)
            ## Get matched object in plink format
            get_matched_plink()

            dst <- firstBreakDist(fam=F)
            hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = 1, doHapPlot = F)
            heights[l_index] = hc$height[length(hc$height)]
            l_index = l_index + 1
        }
        heights_mean[g_index] = mean(heights)
        g_index = g_index + 1
    }
    dists_df = data.frame(generation = generations, height = heights_mean)
    write.table(x = dists_df, file = paste0("classify-simulated-population-age/average_heights-generations_10-100.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
}

writeCompareAgeMethodsFiles <- function(samples, simulations, in_gene, gen_start, gen_stop, gen_step, seed, genealogy){
    source("../../Scripts/Mutation-Age.R")
    source("../../Scripts/Mutation-age-tripleA.R")
    source("../../Scripts/Gandolfo_Speed_Mutation_Age_Estimation.R")
    source("../../Scripts/maxLikelihood.R")
    dir.create("compare-age-methods/", showWarnings = F)
    # compareAgeMethod_generateFiles(gen_start, gen_stop, gen_step, simulations, 
    #                                samples, seed, in_gene, genealogy, ancestral = ancestral_method, 
    #                                cs_correction = cs_correction)
    compareAgeMethod_generateFiles_parallel(gen_start = gen_start, gen_stop = gen_stop, gen_step = gen_step, 
                                            sims = simulations, samples = samples, seed = seed, in_gene = in_gene, 
                                            genealogy = genealogy, ancestral = ancestral_method, minSampleBB = 5, 
                                            cs_correction = cs_correction, pre_path = "")
}

# Call function(s)
if (option == "all"){
    writeSimulatedDistMatrices(n_samples, n_simulations, gene, gen_start, gen_stop, gen_step, seed, genealogy)
    writeSimulatedDistToFounder(n_samples, n_simulations, gene, gen_start, gen_stop, gen_step, seed, genealogy)
    writeCompareAgeMethodsFiles(n_samples, n_simulations, gene, gen_start, gen_stop, gen_step, seed, genealogy)
} else if (option == "dist_mat"){
    writeSimulatedDistMatrices(n_samples, n_simulations, gene, gen_start, gen_stop, gen_step, seed, genealogy)
} else if (option == "haplo_length"){
    writeSimulatedDistToFounder(n_samples, n_simulations, gene, gen_start, gen_stop, gen_step, seed, genealogy)
} else { # option == "compare_age_methods"
    writeCompareAgeMethodsFiles(n_samples, n_simulations, gene, gen_start, gen_stop, gen_step, seed, genealogy)
}