#library(tidyverse)
library(dplyr)
library(parallel)
#library(snow)

plot_ancestral_reconstruction_errors <- function(){
    
    if (manual){
        gen_start = 10
        gen_stop = 510#2510#10010
        gen_step = 10#20#10#20
        sims = 20#20#10
        samples = 100
        in_gene = "BRCA1"
        genealogy="starGenealogy"
        seed = 42
        ancestral = "branchBoundIndep"
        #ancestral = "branchBound"
        #ancestral = "mostFreqBase"
        #sim_type = "knownBreaks"
        sim_type = "famHaplotypes"
        minSampleBB=5
    }
    
    pre_path = paste0("classify-simulated-population-age/", in_gene,"_simulated-",sim_type,"-",genealogy,"-",samples,"_samples-",sims,
                      "_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/")
    
    gen = 360
    
    filename = paste0(pre_path, "simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen,"-seed_",seed)
    brca1_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
    brca1_simulated_geno_plink <<- read.table2(file = paste0(filename, "-geno.txt"), header = T)
    
    for (i in 1:sims){
        mut = paste0("mutSim_BRCA1_",gen,"_",i); gene = "BRCA1_simulated"
        mut <<- mut; gene <<- gene
        readIndvAndFamData(gene, mut)
        # dst <- firstBreakDist(fam=F)
        # hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = 1, doHapPlot = F)
        # brca.mut.cM <<- brca_cM_middle
        
        #mut_location = nrow(filter(chr_coords, cM < brca.mut.cM))
        
        # Simulated founder
        haplotypes <- read.table(paste0(filename, "-founder.txt"), header = T)
        sim_founder <<- as.integer(haplotypes[i,-1])
        if (ancestral == "mostFreqBase"){
            #matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
            matched_plink2 = matched_plink
            haplotypes2 <<- findConsensus_plink(matched_plink2, cutoff = 0.5)
        } else { # star_pop_freqs
            #matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
            matched_plink2 = matched_plink
            if (ancestral == "branchBoundIndep"){
                res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = T, min_samples = minSampleBB)
            } else { #branchBound
                res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = F, min_samples = minSampleBB)
            }
            haplotypes2 <<- res[[1]]
            breaks <<- res[[2]]
            chr_pos <<- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
            print(mutationAge_method8_tweedie(chr_pos, brca_cM_middle))
        }
        
        get_score <- function(row){
            if (row[1] == row[2]){
                val = 0
            } else if (row[1] == 0 && row[2] == 2){
                val = -3 # -2
            } else if (row[1] == 2 && row[2] == 0){
                val = 3 # 3
            } else if ((row[1] == 0 && row[2] == 1)){# || (row[1] == 2 && row[2] == 1)){
                val = -1
            } else if (row[1] == 2 && row[2] == 1){
                val = -2
            } else if ((row[1] == 1 && row[2] == 0)){# || (row[1] == 1 && row[2] == 2)){
                val = 1
            } else if (row[1] == 1 && row[2] == 2){
                val = 2
            }
            names(val) = "combined"
            return(c(row, val))
        }
        scores = t(apply(data.frame(sim_founder,haplotypes2), 1, get_score))
        scores[1:10,]
        
        df = data.frame(x=rep(chr_coords$cM,3), y = c(scores[,1], scores[,2], scores[,3]), group = rep(c("sim_founder", "reconstructed", "combined"), each = length(sim_founder)))
        df2 = filter(df, y==0, group == "reconstructed")
        
        #df = filter(df,group == "sim_founder")
        #df = filter(df,group == "reconstructed")
        df = filter(df,group == "combined")
        p = ggplot(df, aes(x=x, y=y, color=group)) +
            geom_point() +
            #geom_vline(xintercept = mut_location, color = "blue", size=0.8, linetype = 2) +
            geom_vline(xintercept = brca_cM_middle, color = "blue", size=0.8, linetype = 2) +
            xlim(c(min(df2$x)-2,max(df2$x)+2)) +
            theme_bw()
        print(p)
    }
}

optimal_parameters_ancestral <- function(){
    compute_age <- function(par, gen, ancestral_method){
        #par = round(par)
        print(paste(par, gen, ancestral_method))
        chr_pos = getAncestralHaplotype(ancestral = ancestral_method, matched_plink, brca_data, minSampleBB = par)
        print(dim(chr_pos))
        #age8 = mutationAge_method8_tweedie(chr_pos, brca_cM_middle)
        chr_pos <- chr_pos %>% mutate(positive_cM_adj = positive_cM - brca.mut.cM, negative_cM_adj = abs(negative_cM - brca.mut.cM))
        arm.lengths = c(chr_pos$positive_cM_adj, chr_pos$negative_cM_adj)
        age = 100/mean(arm.lengths)
        out = abs(gen-age[1])
        return(out)
    }
    optim(par = 5, fn = compute_age, gen = gen, ancestral_method = ancestral, lower = 2, upper = samples*0.25, method = "L-BFGS-B")
    best = c()
    bestVal = c()
    for (i in 1:sims){
        mut = paste0("mutSim_BRCA1_",gen,"_",i); gene = "BRCA1_simulated"
        mut <<- mut; gene <<- gene
        readIndvAndFamData(gene, mut)
        minVal = gen
        s = 2
        for (minSamples in seq(2,round(samples*0.3),2)){
            val = compute_age(minSamples, gen, ancestral)
            #gridSearch(fun = compute_age, levels = list(2:round(samples*0.3)), sim=i, gen = gen, ancestral_method = ancestral)
            #gridSearch(fun = compute_age, lower = 2, upper = round(samples*0.3), sim=i, gen = gen, ancestral_method = ancestral, n=20)
            if (val<minVal){
                minVal=val
                s = minSamples
            }
        }
        print(s)
        best[i] = s
        bestVal[i] = minVal
    }
    
    if (manual){
        gen_start = 10
        gen_stop = 510#510#10010
        gen_step = 50#10#20
        sims = 20#20#10
        samples = 100
        in_gene = "BRCA1"
        genealogy="starGenealogy"
        #genealogy="correlatedGenealogy"
        seed = 42
        ancestral = "branchBoundIndep"
        ancestral = "simulatedFounder"
        #cs_correction="gandolfo"
        cs_correction="None"
        cs_correction="adjustHaploFreqs"
        #sim_type = "popFreqs"
        sim_type = "famHaplotypes"
        #sim_type = "haplotypes"
        #sim_type = "alleleFreqs"
        # sim_type = "knownBreaks"
    }
    
    pre_path = paste0("classify-simulated-population-age/", in_gene,"_simulated-",sim_type,"-",genealogy,"-",samples,"_samples-",sims,
                      "_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/")
    for (gen in seq(110, 111,2)){
        filename = paste0(pre_path, "simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen,"-seed_",seed)
        brca1_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
        brca1_simulated_geno_plink <<- read.table2(file = paste0(filename, "-geno.txt"), header = T)
        
        for (i in 1:sims){
            optim(par = 5, fn = compute_age, gen = gen, sim = i, ancestral_method = ancestral, lower = 2, upper = samples*0.25, method = "Brent")
        }
    }
}


compute_generations <- function(gen, gen_start, gen_stop, gen_step, sims, samples, genealogy, 
                                ancestral, sim_type, in_gene = "BRCA1", seed = 42, manual=F){
    if (manual){
        print("WHAT!")
        gen_start = 10
        gen_stop = 510#2510#10010
        gen_step = 10#20#10#20
        sims = 20#20#10
        samples = 100
        in_gene = "BRCA1"
        genealogy="starGenealogy"
        seed = 42
        ancestral = "branchBoundIndep"
        #ancestral = "branchBound"
        #ancestral = "mostFreqBase"
        #sim_type = "knownBreaks"
        sim_type = "famHaplotypes"
        #minSampleBB=5
    }
    
    if (Sys.info()["sysname"] == "Linux"){
        setwd("/work/sduvarcall/haplotype-project")
    } else {
        setwd("~/Documents/PhD/haplotype-project")
    }
    
    source("Scripts/Haplotype-projekt-ver2.R")
    source("Scripts/Clustering.R")
    source("Scripts/ValidateAncestralReconstruction.R")
    
    ## Read Data
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    morgan.brca1 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    alleleFreqs_brca1 <<- as.numeric(read.table2("cache/pop_freqs_minorAllele_brca1_ordered.txt", header = T))
    pre_path = paste0("classify-simulated-population-age/", in_gene,"_simulated-",sim_type,"-",genealogy,"-",samples,"_samples-",sims,
                      "_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/")
    filename <<- paste0(pre_path, "simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen,"-seed_",seed)
    brca1_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
    brca1_simulated_geno_plink <<- read.table2(file = paste0(filename, "-geno.txt"), header = T)
    print(filename)
    
    
    #par = expand.grid(sims = 1:sims, minSamples = 2:(samples*0.5))
    # cl2 = parallel::makeCluster(spec=6)
    # parallel::clusterExport(cl2, c("brca1_simulated_pheno_merged", "brca1_simulated_geno_plink", "morgan.brca1", "chr_coords_all"))
    # data = do.call(rbind, parLapply(cl = cl2, X = 1:sims, fun = compute_simulations, gen = gen, minSamples=2:(samples*0.5), ancestral = ancestral, filename = filename))
    # data = do.call(rbind, parLapply(cl = cl2, X = 1:2, fun = compute_simulations, gen = gen, minSamples=2:50, ancestral = ancestral, filename = filename))
    # data = do.call(rbind, parLapply(cl = cl2, X = 50, fun = compute_simulations, gen = gen, minSamples=2:3, ancestral = ancestral, filename = filename))
    # data = do.call(rbind, lapply(X = 1:2, FUN = compute_simulations, gen = gen, minSamples=2:3, ancestral = ancestral, filename = filename))
    data = do.call(rbind, lapply(X = 1:sims, FUN = compute_simulations, gen = gen, minSamples=2:(samples*0.5), ancestral = ancestral, filename = filename))
}

compute_simulations <- function(sim, gen, minSamples, ancestral, filename){
    #for (i in 1:sims){
    source("Scripts/Haplotype-projekt-ver2.R")
    source("Scripts/Clustering.R")
    source("Scripts/ValidateAncestralReconstruction.R")
    mut = paste0("mutSim_BRCA1_",gen,"_",sim); gene = "BRCA1_simulated"
    mut <<- mut; gene <<- gene
    readIndvAndFamData(gene, mut)
    # dst <- firstBreakDist(fam=F)
    # hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = 1, doHapPlot = F)
    # brca.mut.cM <<- brca_cM_middle
    
    #mut_location = nrow(filter(chr_coords, cM < brca.mut.cM))
    
    # Simulated founder
    haplotypes <- read.table(paste0(filename, "-founder.txt"), header = T)
    sim_founder <<- as.integer(haplotypes[sim,-1])
    
    data = do.call(rbind, lapply(minSamples, FUN = compute_minSample, ancestral, gen, sim))
    return(data)
}

compute_minSample <- function(n_samples, ancestral, gen, sim){
    if (n_samples>=nrow(matched_plink)) return()
    matched_plink_subset = matched_plink
    if (ancestral == "mostFreqBase"){
        haplotypes2 <<- findConsensus_plink(matched_plink_subset, cutoff = 0.5)
        # breaks <<- findHaplotypeBreaks_plink(matched_plink_subset, haplotypes2)
        # chr_pos <- mapNearestSNPs_plink(matched_plink_subset, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    } else {
        if (ancestral == "branchBoundIndep"){
            res = branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = 0.5, indep_sides = T, min_samples = n_samples, verbose = "left")
        } else if (ancestral == "branchBound"){ #branchBound
            res = branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = 0.5, indep_sides = F, min_samples = n_samples)
        } else if (ancestral == "branchBoundIndep_alleleFreq"){
            res = branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = 0.5, indep_sides = T, 
                                             min_samples = n_samples, alleleFreqs = alleleFreqs, verbose = "right")
        } else if (ancestral == "branchBound_alleleFreq"){
            res = branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = 0.5, indep_sides = F, 
                                             min_samples = n_samples, alleleFreqs = alleleFreqs)
        } else {
            stop("Ancestral method does not exist")
        }
        haplotypes2 <- res[[1]]
        breaks <- res[[2]]
        break_left <- ifelse(is.na(res[[3]]), 1, res[[3]])
        break_right <- ifelse(is.na(res[[4]]), nrow(chr_coords), res[[4]])
        chr_pos <- mapNearestSNPs_plink(matched_plink_subset, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    }
    # sim_founder[2072:2076]
    # haplotypes2[2072:2076]
    #print(break_left)
    sim_founder.limited = sim_founder[break_left:break_right]
    haplotypes2.limited = haplotypes2[break_left:break_right]
    diff = abs(sim_founder.limited-haplotypes2.limited)==2
    if (sum(diff)==0){
        diff_coords = data.frame(position_b37=-1, cM=-1)
        error_type = NA
    } else {
        diff_coords = filter(chr_coords, SNP %in% names(diff[diff==1]))
        index=which(diff==1)
        error_type = paste0(sim_founder.limited[index], "->", haplotypes2.limited[index])
    }
    #error_df = data.frame(actual_age = gen, sim = i, error_pos = diff_coords$position_b37, error_cM = diff_coords$cM)
    
    chr_pos <- chr_pos %>% mutate(positive_cM_adj = positive_cM - brca.mut.cM, negative_cM_adj = abs(negative_cM - brca.mut.cM))
    arm.lengths = c(chr_pos$positive_cM_adj, chr_pos$negative_cM_adj)
    age = 100/mean(arm.lengths)
    length = chr_coords[break_right, 3]-chr_coords[break_left, 3]
    length_cM = chr_coords[break_right, "cM"]-chr_coords[break_left, "cM"]
    ancestral_info = data.frame(actual_age=gen, sim=sim, minSamples = n_samples, predicted_age = age,
                                ancestral_length = length, ancestral_length_cM = length_cM, 
                                error_type = error_type, error_pos = diff_coords$position_b37, error_cM = diff_coords$cM)
    #print(ancestral_info)
    return(ancestral_info)
}

find_parameters <- function(){
    if (Sys.info()["sysname"] == "Linux"){
        setwd("/work/sduvarcall/haplotype-project")
        cores=24
    } else {
        setwd("~/Documents/PhD/haplotype-project")
        cores=1
    }
    
    # gen_start = 10; gen_stop = 510; gen_step = 10; sims = 50; samples = 100; genealogy = "starGenealogy"; ancestral = "branchBoundIndep"; sim_type = "famHaplotypes"
    # gen_start = 10; gen_stop = 510; gen_step = 10; sims = 50; samples = 20; genealogy = "starGenealogy"; ancestral = "branchBoundIndep"; sim_type = "famHaplotypes"
    gen_start = 10; gen_stop = 510; gen_step = 10; sims = 50; samples = 100; genealogy = "starGenealogy"; ancestral = "branchBoundIndep"; sim_type = "haplotypes"
    #gen_start = 10; gen_stop = 2510; gen_step = 50; sims = 20; samples = 100; genealogy = "starGenealogy"; ancestral = "branchBoundIndep"; sim_type = "famHaplotypes"
    # gen_start = 10; gen_stop = 2010; gen_step = 50; sims = 50; samples = 100; genealogy = "starGenealogy"; ancestral = "branchBoundIndep"; sim_type = "famHaplotypes"
    
    # gen_start = 10; gen_stop = 510; gen_step = 10; sims = 50; samples = 100; genealogy = "starGenealogy"; ancestral = "branchBoundIndep"; sim_type = "famHaplotypesNoHet"
    # gen_start = 10; gen_stop = 510; gen_step = 10; sims = 20; samples = 100; genealogy = "starGenealogy"; ancestral = "branchBoundIndep_alleleFreq"; sim_type = "famHaplotypesNoHet"
    # gen_start = 10; gen_stop = 510; gen_step = 10; sims = 20; samples = 100; genealogy = "starGenealogy"; ancestral = "branchBoundIndep"; sim_type = "famHaplotypesNoHet"
    # gen_start = 10; gen_stop = 510; gen_step = 10; sims = 20; samples = 100; genealogy = "correlatedGenealogy"; ancestral = "branchBoundIndep_alleleFreq"; sim_type = "famHaplotypesNoHet"
    # gen_start = 10; gen_stop = 510; gen_step = 10; sims = 20; samples = 100; genealogy = "correlatedGenealogy"; ancestral = "branchBoundIndep"; sim_type = "famHaplotypesNoHet"
    
    cl = makeCluster(spec=cores)
    data = do.call(rbind, parLapply(cl = cl, X = seq(gen_start, gen_stop, gen_step), fun = compute_generations,
                                    gen_start, gen_stop, gen_step,
                                    sims, samples, genealogy, ancestral, sim_type))
    #data = do.call(rbind, parLapply(cl = cl, X = seq(10,20,10), fun = compute_generations))
    # data = do.call(rbind, lapply(X = seq(110,120,gen_step), FUN = compute_generations,
    #                              gen_start, gen_stop, gen_step,
    #                              sims, samples, genealogy, ancestral, sim_type))
    dir.create("validateAncestral", showWarnings = F)
    write.table(data, file = paste0("validateAncestral/ancestralReconstruction_",gen_start,"-",gen_stop,"-step",gen_step,"_", sims,"sims_",samples,"samples_",genealogy,"_",ancestral,"_ancestralMethod_",sim_type,".txt"), 
                quote = F, row.names = F, sep = "\t")
}

plot_parameters <- function(){
    setwd("~/Documents/PhD/haplotype-project")
    
    get_info <- function(){
        # chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
        morgan.brca1 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
        brca_name <<- "BRCA1"
        brca_start <<- 41197695; brca_stop <<- 41276113
        brca_middle <<- (brca_stop - brca_start)/2 + brca_start
        brca_cM_middle <<- morgan.brca1[which.min(abs(morgan.brca1$position-brca_middle)),3]
        brca_chr <<- 17
    }
    get_info()
    
    # data = read.table("validateAncestral/ancestralReconstruction_10-510-step10_20sims_100samples.txt", header = T)
    #data = read.table2("validateAncestral/ancestralReconstruction_10-2510-step50_20sims_100samplesbranchBoundIndep_ancestralMethod.txt", header = T)
    # data = read.table2("validateAncestral/ancestralReconstruction_10-510-step10_20sims_100samples_correlatedGenealogy_branchBoundIndep_ancestralMethod.txt", header = T)
    # data = read.table2("validateAncestral/ancestralReconstruction_10-510-step10_20sims_100samples_correlatedGenealogy_branchBoundIndep_alleleFreq_ancestralMethod.txt", header = T)
    # data = read.table2("validateAncestral/ancestralReconstruction_10-510-step10_20sims_100samples_starGenealogy_branchBoundIndep_alleleFreq_ancestralMethod.txt", header = T)
    data = read.table2("validateAncestral/ancestralReconstruction_10-510-step10_50sims_100samples_starGenealogy_branchBoundIndep_ancestralMethod_famHaplotypes.txt", header = T)
    data = read.table2("validateAncestral/ancestralReconstruction_10-510-step10_50sims_100samples_starGenealogy_branchBoundIndep_ancestralMethod_haplotypes.txt", header = T)
    data = read.table2("validateAncestral/ancestralReconstruction_10-510-step10_50sims_100samples_starGenealogy_branchBoundIndep_ancestralMethod_famHaplotypesNoHet.txt", header = T)
    data = read.table2("validateAncestral/ancestralReconstruction_10-510-step10_50sims_20samples_starGenealogy_branchBoundIndep_ancestralMethod_famHaplotypes.txt", header = T)
    data2 = filter(data, error_pos>-1)
    
    dim(data2)
    count(data2,error_type)
    data2 %>% filter(minSamples>3) %>% count(error_type)
    count(data2, error_type, minSamples) %>% filter(minSamples>3) %>% ggplot(aes(x=minSamples, y=n)) + geom_col(fill="blue")
    filter(data2, minSamples>3) %>% count(actual_age) %>% ggplot(aes(x=actual_age, y=n)) + geom_col(fill="blue")
    #data2 %>% filter(actual_age>100 & actual_age<170) %>% filter(minSamples>2) %>% count(actual_age, sim, minSamples) %>% View()
    data2 %>% filter(minSamples>3) %>% count(actual_age, sim, minSamples) %>% filter(n>3) %>% View()
    
    #data3 = data %>% filter(actual_age %in% seq(10,510,100)) %>% mutate(age_label = paste(actual_age, "generations"))
    # data3 = data %>% mutate(age_label = paste(actual_age, "generations")) 
    # #d=cut(data3$actual_age, breaks = c(0, seq(60,510,50)), labels = paste0("=<",seq(60,510,50), " generations"))
    # d=cut(data3$actual_age, breaks = c(0, seq(60,510,50)), labels = paste0(c(10, seq(70,470,50)), "-", seq(60,510,50), " generations"))
    # data3$group_label = d
    
    data3 = data %>% mutate(age_label = paste(actual_age, "generations")) %>% filter(actual_age<510)
    d=cut(data3$actual_age, breaks = c(0, seq(50,500,50)), labels = paste0(seq(10,470,50), "-", seq(50,500,50), " generations"))
    data3$group_label = d
    
    # data3 = data %>% mutate(age_label = paste(actual_age, "generations")) %>% filter(actual_age<2510)
    # d=cut(data3$actual_age, breaks = c(0, seq(250,2500,250)), labels = paste0(seq(10,2500,250), "-", seq(250,2500,250), " generations"))
    # data3$group_label = d
    
    # data3$group_label = paste(data3$actual_age, "generations")
    
    neg = data3 %>% filter(error_pos>0, error_cM < brca_cM_middle) %>% group_by(age_label, minSamples, sim) %>% slice(which.min(error_cM))
    pos = data3 %>% filter(error_pos>0, error_cM > brca_cM_middle) %>% group_by(age_label, minSamples, sim) %>% slice(which.max(error_cM))
    data4 = bind_rows(neg, pos)
    data4 %>% mutate(dist=abs(brca_cM_middle-error_cM)) %>% filter(minSamples>2) %>%
        group_by(minSamples, group_label) %>% summarise(n=n(), mean_dist=mean(dist)) %>%
        ggplot(aes(x=minSamples, y = mean_dist, group=group_label, color=group_label)) +
        geom_line(size=1.1) +
        scale_y_continuous(breaks = seq(0,100,1)) +
        #scale_x_continuous(breaks = seq(0,50,1)) +
        #scale_color_brewer(palette = "Dark2") +
        scale_color_manual(values = colors_plot(color_names = unique(data4$group_label), palette = "sasha")) +
        theme_classic() +
        theme(panel.grid.major.y = element_line(size = 0.6, colour = "lightgrey", linetype = 2))
    
    
    # Plot errors per group
    data3 %>% filter(error_pos>0, minSamples>2) %>% group_by(group_label, minSamples) %>% summarise(n=n()) %>%
        ggplot(aes(x=minSamples, y = n, group = group_label, color = group_label)) +
        geom_point() +
        geom_line(size=1) +
        theme_bw() +
        scale_x_continuous(breaks = seq(0,50,5)) +
        scale_y_continuous(breaks = seq(0,5000,20)) +
        scale_color_manual(values = colors_plot(color_names = unique(data4$group_label), palette = "sasha")) +
        ylab("Number of errors")
    
    # Plot errors over all generations
    data %>% filter(error_pos>0, minSamples>2) %>% mutate(dist=abs(brca_cM_middle-error_cM)) %>% 
        group_by(minSamples) %>% summarise(n=n(), mean_dist=mean(dist)) %>% 
        ggplot(aes(x=minSamples, y = n)) +
        #geom_line(color="blue", size=1.1) +
        geom_col(fill="blue") +
        scale_x_continuous(breaks = seq(0,50,5)) +
        scale_y_continuous(breaks = seq(0,5000,250)) +
        theme_bw() +
        ylab("Number of errors")
    
    # Plot location of errors
    data %>% filter(error_pos>0, minSamples>2) %>% mutate(dist=abs(brca_cM_middle-error_cM)) %>% #filter(actual_age==160 & sim == 17) %>%
        ggplot(aes(x=error_cM)) +
        #geom_line(color="blue", size=1.1) +
        geom_histogram(binwidth = 0.05) +
        geom_vline(xintercept = brca_cM_middle, color = "black", size=0.8, linetype = 2) +
        theme_bw() +
        xlim(c(64.5,69.5)) +
        #xlim(c(66.5,67.5)) +
        ylab("Number of errors") +
        xlab("Location of errors in cM")
    
    # Predicted age versus minSamples - extracted subset
    data5 = data %>% filter(actual_age %in% seq(10,510,10)) %>% mutate(age_label = paste(actual_age, "generations")) %>%
    # data5 = data %>% filter(actual_age %in% seq(10,510,50)) %>% mutate(age_label = paste(actual_age, "generations")) %>%
    # data5 = data %>% filter(actual_age %in% seq(60,140,20)) %>% mutate(age_label = paste(actual_age, "generations")) %>% 
        distinct(actual_age, sim, minSamples, .keep_all = T)
    data5_2 = data5 %>% filter(minSamples>2) %>% group_by(age_label, minSamples) %>% #mutate(age_diff=abs(predicted_age-actual_age))
        summarise(predicted_age=mean(predicted_age), actual_age=mean(actual_age)) %>% mutate(age_diff=abs(predicted_age-actual_age))
    
    data5_mini <- data5_2 %>% group_by(actual_age) %>% slice(which.min(age_diff))
    model = lm(predicted_age~minSamples, data5_mini)
    
    p = data5_2 %>% 
        ggplot(aes(x=minSamples, y = predicted_age, group=age_label, color=reorder(age_label, actual_age))) +
        geom_line(size=1.1) +
        geom_abline(slope=model$coefficients[2], intercept = model$coefficients[1], linetype=2, color="darkblue") +
        scale_y_continuous(breaks = seq(10,3010,10)) +
        # scale_y_continuous(breaks = seq(10,3010,50)) +
        scale_x_continuous(breaks = seq(0,50,5)) +
        # scale_color_brewer(palette = "Dark2") +
        scale_color_manual(values = colors_plot(color_names = unique(data5$age_label), palette = "sasha")) +
        theme_classic() +
        theme(text = element_text(size=14)) +
        theme(legend.title = element_text(face = "bold")) +
        theme(panel.grid.major.y = element_line(size = 0.6, colour = "lightgrey", linetype = 2)) +
        labs(y="Estimated age", x="Min samples", color="True age")
    print(p)
    ggsave(filename = "validateAncestral/optimal_minSamples_10-510-step10_50sims_100samples_starGenealogy_branchBoundIndep_ancestralMethod_famHaplotypes.pdf",
           width = 9, height = 9, scale = 0.9)
    
    # Predicted age versus minSamples - group intervals
    data6 = data3 %>% distinct(actual_age, sim, minSamples, .keep_all = T) %>% filter(minSamples>2) %>%
        group_by(group_label, minSamples) %>% summarise(predicted_age=mean(predicted_age), actual_age=mean(actual_age)) %>%
        mutate(age_diff=abs(predicted_age-actual_age))
    
    data6_mini <- data6 %>% group_by(actual_age) %>% slice(which.min(age_diff))
    model = lm(predicted_age~minSamples, data6_mini)
        
    data6 %>% ggplot(aes(x=minSamples, y = predicted_age, group=group_label, color=group_label)) +
        geom_line(size=1.1) +
        geom_abline(slope=model$coefficients[2], intercept = model$coefficients[1], linetype=2, color="darkblue") +
        scale_y_continuous(breaks = seq(30,3010,50)) +
        scale_x_continuous(breaks = seq(0,50,5)) +
        #scale_color_brewer(palette = "Dark2") +
        scale_color_manual(values = colors_plot(color_names = unique(data3$group_label), palette = "sasha")) +
        theme_classic() +
        theme(panel.grid.major.y = element_line(size = 0.6, colour = "lightgrey", linetype = 2)) +
        ylab("predicted age")
    
    # Ranking minSamples parameter
    data7 = data3 %>% distinct(actual_age, sim, minSamples, .keep_all = T) %>%
        group_by(age_label, minSamples) %>% 
        summarise(predicted_age=mean(predicted_age), actual_age=mean(actual_age)) %>% 
        mutate(age_diff=abs(predicted_age-actual_age)) %>%
        arrange(actual_age, age_diff) %>% 
        mutate(rank = rank(age_diff)) %>% 
        group_by(minSamples) %>%
        summarise(predicted_age=mean(predicted_age), rank = mean(rank))
    
    data7 %>% 
        ggplot(aes(x=minSamples, y=rank)) +
        geom_line(color="blue") +
        geom_point(color="blue") +
        scale_x_continuous(breaks = seq(0,50,5)) +
        theme_bw()
    
    #View(data7)
    best_minSamples = data7 %>% slice(which.min(rank)) %>% .$minSamples
    
    # Plot actual age vs predicted age for each minSample
    d3 = data3 %>% filter(minSamples %in% seq(5,50,5)) %>% #filter(minSamples==best_minSamples) %>%
        group_by(group_label, minSamples) %>% summarise(actual_age=mean(actual_age), predicted_age=mean(predicted_age))
    #d3$minSamples = factor(d3$minSamples)
    d3 %>%
        ggplot(aes(x=actual_age, y=predicted_age, group=factor(minSamples), color=factor(minSamples))) +
        geom_line() +
        geom_point() + 
        scale_x_continuous(breaks = seq(10,1000,50)) +
        scale_y_continuous(breaks = seq(10,1000,50)) +
        theme_bw() +
        scale_color_manual(values = colors_plot(color_names = unique(d3$minSamples), palette = "sasha")) +
        geom_abline(slope=1, intercept = 0, color="darkblue")
        
    
    
    # Ancestral length versus minSamples
    # data5, gen = 10,100,...,500
    data5 %>% group_by(age_label, minSamples) %>% #filter(actual_age>10) %>%
        summarise(ancestral_length=mean(ancestral_length), ancestral_length_cM=mean(ancestral_length_cM)) %>%
        ggplot(aes(x=minSamples, y = ancestral_length_cM, group=age_label, color=age_label)) +
        geom_line(size=1.1) +
        scale_y_continuous(breaks = seq(0,1010,1)) +
        scale_x_continuous(breaks = seq(0,50,5)) +
        #scale_color_brewer(palette = "Dark2") +
        scale_color_manual(values = colors_plot(color_names = unique(data5$age_label), palette = "sasha")) +
        theme_classic() +
        theme(panel.grid.major.y = element_line(size = 0.6, colour = "lightgrey", linetype = 2)) +
        ylab("Ancestral reconstruction length in cM")
    
    data3 %>% distinct(actual_age, sim, minSamples, .keep_all = T) %>% 
        group_by(group_label, minSamples) %>% summarise(ancestral_length=mean(ancestral_length), ancestral_length_cM=mean(ancestral_length_cM)) %>%
        ggplot(aes(x=minSamples, y = ancestral_length_cM, group=group_label, color=group_label)) +
        geom_line(size=1.1) +
        scale_y_continuous(breaks = seq(0,1010,1)) +
        scale_x_continuous(breaks = seq(0,50,5)) +
        #scale_color_brewer(palette = "Dark2") +
        scale_color_manual(values = colors_plot(color_names = unique(data3$group_label), palette = "sasha")) +
        theme_classic() +
        theme(panel.grid.major.y = element_line(size = 0.6, colour = "lightgrey", linetype = 2)) +
        ylab("Ancestral reconstruction length in cM")
}    
