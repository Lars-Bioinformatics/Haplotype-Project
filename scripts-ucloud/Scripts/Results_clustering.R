# setwd("~/Documents/PhD/haplotype-project")

# library(fossil)
library(clues)
library(foreach)
library(doMC)
source("Scripts/Haplotype-projekt-ver2.R")
source("Scripts/Clustering.R")
source("Scripts/Mutation-Age.R")
source("Scripts/loadShinyData.R")
source("Scripts/Gandolfo_Speed_Mutation_Age_Estimation.R")
source("Scripts/Mutation-age-tripleA.R")
source("Scripts/maxLikelihood.R")

# Result dir
dir.create("compare_cluster_methods", showWarnings = F)

comb2 <- function(n){
    if (n<2) return(0)
    n * (n-1) / 2
}

comb <- function(n, x) {
    print(n)
    if (n<x) return(0)
    factorial(n) / (factorial(n-x) * factorial(x))
}

cluster <- function(pheno, geno, pops, clusterMethods, groups, labels){
    # brca_data.filtered = filter(pheno, Country %in% paste0("simuland", n:(k+n-1)))
    brca_data.filtered = filter(pheno, Country %in% paste0("simuland", pops))
    brca_data.filtered <<- brca_data.filtered
    dim(brca_data.filtered)
    matched_plink_single = filter(geno, SNP %in% brca_data.filtered$Onc_ID)
    dim(matched_plink_single)
    
    dst = firstBreakDist(matched_custom = matched_plink_single)
    length(dst)
    
    res = data.frame(clusterMethod = clusterMethods, Founders = groups, Fmeasure=0, randIndex = 0, adjRandIndex = 0, parameters = 0, failed_count = 0)
    rownames(res) = clusterMethods
    for (m in clusterMethods){
        k = groups
        # randIndex_best = 0
        # F_measure_best = 0
        while (k < 20){
            if (m == "DBSCAN"){
                # kNNdistplot(dst)
                minPts = 4
                eps = 0.6
                db <- fpc::dbscan(data = dst, MinPts = minPts, eps = eps, method = "dist")
                cluster_groups <- db$cluster
                brca_data.filtered$cluster_groups <- db$cluster
                brca_data <- brca_data.filtered
            } else if (m == "HDBSCAN") {
                minPts = 4
                hdb <- dbscan::hdbscan(x = dst, minPts = minPts)
                cluster_groups <- hdb$cluster
                brca_data.filtered$cluster_groups <- hdb$cluster
                brca_data <- brca_data.filtered
            } else if (m == "k-means"){
                km <- kmeans(dst, k)
                cluster_groups <- km$cluster
                brca_data.filtered$cluster_groups <- km$cluster
                brca_data <- brca_data.filtered
            } else {
                hc <- hiearchical_clustering(dst, clustMethod = m, k = k, doHapPlot = F)
            }
            
            cls = brca_data[,c("Country", "cluster_groups")] %>% count(Country,cluster_groups)
            big_clusters = cls %>% group_by(Country) %>% slice(which.max(n)) %>% group_by(cluster_groups) %>% slice(which.max(n)) %>% arrange(cluster_groups)
            totals = cls %>% group_by(cluster_groups) %>% summarise(n=sum(n)) %>% arrange(cluster_groups)
            
            if (k < 20 && nrow(big_clusters)<groups && !(m %in% c("DBSCAN", "HDBSCAN"))){
                # cluster(k+1, m, groups, labels,eps,minPts)
                k = k + 1
                # print(cls)
                # if (randIndex > randIndex_best && F_measure > F_measure_best){
                #     randIndex_best = randIndex
                #     F_measure_best = F_measure
                # }

            } else {
                # print(paste("Method:", m, "k:", k))
                # print(big_clusters)
                # print(cls)
                
                RI = adjustedRand(labels, cluster_groups)
                
                # Rand Index
                # randIndex = rand.index(labels, cluster_groups)
                randIndex = RI[1]
                
                # Adjusted Rand Index
                ARI = RI[2]
                
                # Manual
                # contingency.table = xtabs(n~Country+cluster_groups, cls)
                # row_col = (sum(apply(contingency.table, 1, function(n) comb2(sum(n)))) * sum(apply(contingency.table, 2, function(n) comb2(sum(n)))))
                # nominator = sum(sapply(contingency.table, function(n) comb2(n))) - row_col / comb2(sum(contingency.table))
                # denominator = 1/2 * (sum(apply(contingency.table, 1, function(n) comb2(sum(n)))) + sum(apply(contingency.table, 2, function(n) comb2(sum(n))))) - row_col / comb2(sum(contingency.table))
                # ARI = nominator / denominator
                
                # print(paste("Rand Index:", randIndex))
                # F-score
                totals = filter(totals, cluster_groups %in% big_clusters$cluster_groups)
                precision = big_clusters$n/totals$n
                recall = big_clusters$n/table(labels)[1]
                
                # Old - meaningfull?
                # precision = mean(big_clusters$n/totals$n)
                # recall = mean(big_clusters$n/table(labels)[1])
                
                # micro-F_measure. Below two formulas shouls be equal
                # precision = sum(big_clusters$n)/sum(totals$n)
                # recall = sum(big_clusters$n)/length(labels)
                
                F_measure = 2 * precision * recall / (precision + recall)
                # Averaged F_measure
                F_measure = sum(F_measure) / groups
                # print(paste("F-measure:", F_measure))
                
                if (m == "DBSCAN"){
                    res[m, 3:6] = c(F_measure, randIndex, ARI, paste0("minPts=",minPts,", eps=", eps))
                } else if (m == "HDBSCAN") {
                    res[m, 3:6] = c(F_measure, randIndex, ARI, paste0("minPts=",minPts))
                } else {
                    res[m, 3:6] = c(F_measure, randIndex, ARI, paste0("k=",k))
                }
                
                break()
            }
        }
        if (k == 20){
            # print(paste("Method:", m, "k:", k))
            # print("CLUSTERING FAILED")
            # res[m,3:7] = c(randIndex_best, F_measure_best, "FAILED", 1)
            res[m,7] = 1
        }
    }
    return(res)
}

# mut = "c.7617+1G>A"; gene = "BRCA2"
# readIndvAndFamData(gene, mut)
# dst <- firstBreakDist(fam=F)
# clusterMethods = c("ward.D2")#, "average", "single", "centroid", "k-means")


get_data <- function(gen = 200){
    # geno = read.table2(paste0("classify-simulated-population-age/BRCA1_simulated-famHaplotypesNoHet-starGenealogy-100_samples-50_simulations-generations_10_2010_step50-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-50_simulations-generations_",gen,"-seed_42-geno.txt"), header = T)
    # pheno = read.table2(paste0("classify-simulated-population-age/BRCA1_simulated-famHaplotypesNoHet-starGenealogy-100_samples-50_simulations-generations_10_2010_step50-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-50_simulations-generations_",gen,"-seed_42-pheno.txt"), header = T)
    # geno <<- read.table2(paste0("classify-simulated-population-age/BRCA1_simulated-famHaplotypes-starGenealogy-100_samples-50_simulations-generations_10_510_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-50_simulations-generations_",gen,"-seed_42-geno.txt"), header = T)
    # pheno <<- read.table2(paste0("classify-simulated-population-age/BRCA1_simulated-famHaplotypes-starGenealogy-100_samples-50_simulations-generations_10_510_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-50_simulations-generations_",gen,"-seed_42-pheno.txt"), header = T)
    # geno <<- read.table2(paste0("classify-simulated-population-age/BRCA1_simulated-haplotypes-starGenealogy-100_samples-50_simulations-generations_10_510_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-50_simulations-generations_",gen,"-seed_42-geno.txt"), header = T)
    # pheno <<- read.table2(paste0("classify-simulated-population-age/BRCA1_simulated-haplotypes-starGenealogy-100_samples-50_simulations-generations_10_510_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-50_simulations-generations_",gen,"-seed_42-pheno.txt"), header = T)
    geno <<- read.table2(paste0("classify-simulated-population-age/BRCA1_simulated-famHaplotypes-starGenealogy-20_samples-50_simulations-generations_10_510_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-20_samples-50_simulations-generations_",gen,"-seed_42-geno.txt"), header = T)
    pheno <<- read.table2(paste0("classify-simulated-population-age/BRCA1_simulated-famHaplotypes-starGenealogy-20_samples-50_simulations-generations_10_510_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-20_samples-50_simulations-generations_",gen,"-seed_42-pheno.txt"), header = T)
    
    # Coords
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    # Read Genetic coordinates (centimorgan)
    morgan.brca1 <<- read.table2("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    
    # BRCA info
    brca_name <<- "BRCA1"
    brca_start <<- 41197695; brca_stop <<- 41276113
    brca_middle <<- (brca_stop - brca_start)/2 + brca_start
    brca_cM_middle <<- morgan.brca1[which.min(abs(morgan.brca1$position-brca_middle)),3]
    brca_chr <<- 17
    
    # Chr_coords
    chr_coords <<- filter(chr_coords_all, Chr_numeric == brca_chr) %>% 
        # Filter away duplicate snps i.e. snps with multiple names listed multiple times
        filter(SNP %in% names(geno[,-1])) %>% 
        distinct(position_b37, .keep_all = T) %>% 
        # sort snps according to location
        arrange(position_b37)
}

run_cluster_results <- function(){
    if (Sys.info()["sysname"] == "Linux"){
        # setwd("/work/sduvarcall/haplotype-project")
        cores=24
    } else {
        # setwd("~/Documents/PhD/haplotype-project")
        cores=6
    }
    doMC::registerDoMC(cores = cores)
    for (age in c(50,100,200)){
    # for (age in c(50,100)){
        print(paste("Age:", age))
        get_data(age)
        for (k in 2:4){
        # for (k in 2:3){
            print(paste("k:", k))
            set.seed(42)
            clusterMethods = c("ward.D2", "average", "single", "centroid", "complete", "k-means", "HDBSCAN", "DBSCAN")
            sims = 1000
            samples = 20
            
            pops_combinations = list()
            for (i in 1:sims){
                pops = sort(sample(50,k))
                while(list(pops) %in% pops_combinations){
                    pops = sort(sample(50,k))
                }
                pops_combinations[[i]] = pops
            }
            
            res = foreach(i=1:sims, .combine = rbind) %dopar% {
                print(paste("Sim:", i))
                pops = pops_combinations[[i]]
                # cluster(pheno, geno, pops, clusterMethods, groups = k, labels = rep(1:k, each = 100), eps=1.5, minPts=30)
                cluster(pheno, geno, pops, clusterMethods, groups = k, labels = rep(1:k, each = samples))
            }
            # Save complete table
            write.table(res, file = paste0("compare_cluster_methods/Results_Clustering_Full_Table_age_",age,"_k_",k,"_sims_",sims,".txt"), quote = F, row.names = F, sep = "\t")
            
            # age=200; k=4;
            # res = read.table(paste0("compare_cluster_methods/Results_Clustering_Full_Table_age_",age,"_k_",k,"_sims_1000.txt"), header = T, sep = "\t")
            
            # Summary table
            failed = res %>% group_by(clusterMethod, Founders) %>% dplyr::summarise(failed_count = sum(failed_count))
            s = res %>% filter(randIndex>0) %>% group_by(clusterMethod, Founders) %>% dplyr::summarise(Fmeasure = mean(as.numeric(Fmeasure)), randIndex = mean(as.numeric(randIndex)), adjRandIndex = mean(as.numeric(adjRandIndex)))
            res2 = merge(s, failed, by = c("clusterMethod", "Founders"))
            write.table(res2, file = paste0("compare_cluster_methods/Results_Clustering_Summary_age_",age,"_k_",k,"_sims_",sims,".txt"), quote = F, row.names = F, sep = "\t")
        
            # population combinations
            write.table(do.call(rbind, pops_combinations), file = paste0("compare_cluster_methods/population_combinations_k_",k,"_sims_",sims,".txt"), quote = F, row.names = F, col.names = F, sep = "\t")
        }
    }
}

create_pretty_table <- function(){
    data1 <- do.call(rbind, lapply(list.files(path = "compare_cluster_methods/", full.names = T, pattern = "Summary_age_50"), read.table, header = T))
    data2 <- do.call(rbind, lapply(list.files(path = "compare_cluster_methods/", full.names = T, pattern = "Summary_age_100"), read.table, header = T))
    data3 <- do.call(rbind, lapply(list.files(path = "compare_cluster_methods/", full.names = T, pattern = "Summary_age_200"), read.table, header = T))
    # data1 <- do.call(rbind, lapply(list.files(path = "compare_cluster_methods/famHaplotypes_1000_simulations_adjRandIndex/", full.names = T, pattern = "Summary_age_50"), read.table, header = T))
    # data2 <- do.call(rbind, lapply(list.files(path = "compare_cluster_methods/famHaplotypes_1000_simulations_adjRandIndex/", full.names = T, pattern = "Summary_age_100"), read.table, header = T))
    # data3 <- do.call(rbind, lapply(list.files(path = "compare_cluster_methods/famHaplotypes_1000_simulations_adjRandIndex/", full.names = T, pattern = "Summary_age_200"), read.table, header = T))
    data <- rbind(data1,data2,data3)
    # data <- rbind(data1,data2)
    names(data) = c("Method", "Founders", "F-measure", "Rand Index", "Adjusted Rand Index", "Failures")
    data$Age = rep(c(50,100,200), each = 24)
    # data$Age = rep(c(50,100), each = 24)
    data_ordered = data[,c(1,2,7,3:6)] %>% arrange(Method, Founders, Age) %>% filter(Founders<4)
    data_ordered$`Rand Index` = format(round(data_ordered$`Rand Index`, 4), digits=4)
    data_ordered$`F-measure` = format(round(data_ordered$`F-measure`, 4), digits=4)
    data_ordered$`Adjusted Rand Index` = format(round(data_ordered$`Adjusted Rand Index`, 4), digits=4)
    data_ordered$Failures = format(data_ordered$Failures/10, digits = 1)
    # Remove some method names, for prettier output
    data_ordered$Method = as.character(data_ordered$Method)
    data_ordered$Method[2:6 + rep(seq(0,42,6), each=5)] = ""
    # skip Rand index
    data_ordered = data_ordered[, -5]
    # Skip age 200
    data_ordered = filter(data_ordered, Age < 200)
    write.table(data_ordered, file = "compare_cluster_methods/Pretty_table_summary_latex.txt", quote = F, row.names = F, sep = " & ", eol = " \\\\\n")
}






