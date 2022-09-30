#setwd("~/Documents/PhD/haplotype-project")

#library(foreach)
#library(doParallel)

source("Scripts/Haplotype-projekt-ver2.R")
source("Scripts/Clustering.R")

## Read input data into global variables
#readData()
#readFamData()

readBRCA1Data <- function(){
    brca1_pheno_merged <<- read.csv("input/139_ouh_june_2017/B1_Onco_phenotype_distribution_311215.csv")
    load("cache/brca1_geno_merged.RData", .GlobalEnv)
    load("cache/brca1_geno_plink.RData", envir = .GlobalEnv)
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed.txt", header = T, sep = ",")
    load("cache/famHaplotypes-BRCA1-geno.RData", envir = .GlobalEnv)
    load("cache/famHaplotypes-BRCA1-geno_plink_format.RData", envir = .GlobalEnv)
}

precompute <- function(gene, mut, dst, distMethod, clustMethod, fam=F, k, cutoff){

    # Re-creates shiny input object
    shinyInput = list()
    shinyInput$distMethod = distMethod
    shinyInput$clustMethod = clustMethod
    shinyInput$fam = fam
    shinyInput$cutoff = cutoff

    # Compute data
    # dst <- firstBreakDist(fam)
    hc <- hiearchical_clustering(dst, clustMethod = clustMethod, k = k, doHapPlot = F)
    brca_data <<- brca_data.filtered
    fam = if (fam) "fam" else "individual"

    for (i in 1:k){
        dir.create(paste0("cache/breakpointPlots/", mut_name, "/"), showWarnings = F)
        filename2 = paste0("cache/breakpointPlots/", mut_name, "/", gene, "-", mut_name, "-BreakpointPlot-",
                           "distMethod-", distMethod, "-clustMethod-", clustMethod, "-",
                           fam, "-k_", k, "-cutoff_", cutoff, "-group_", i, ".png")
        dir.create(paste0("cache/nearestBreakpointPlots/", mut_name, "/"), showWarnings = F)
        filename3 = paste0("cache/nearestBreakpointPlots/", mut_name, "/", gene, "-", mut_name, "-NearestBreakpointPlot-",
                           "distMethod-", distMethod, "-clustMethod-", clustMethod, "-",
                           fam, "-k_", k, "-cutoff_", cutoff, "-group_", i, ".png")

        if (!file.exists(filename2)){
            chr_pos = haplotype_mutation(mut, prepare = F, group = i, cutoff = cutoff)
            nearestBreaksStatistics(chr_pos, mut, country = NULL, group = i, input = shinyInput)

            p = plot_haplotype_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = i)
            ggsave(filename2, scale = 1, width = 17.5, height = 10)

            nearestBreakPlots = plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name,
                                                              brca_start, brca_stop, mut, group = i)
            p2 = ggarrange(nearestBreakPlots[[1]], nearestBreakPlots[[4]], nearestBreakPlots[[2]], nearestBreakPlots[[3]], ncol=2, nrow=2)
            ggsave(filename3, scale = 1, width = 17.5, height = 11.429)

        }
    }
    return("Success")
}

indexScores <- c("mcclain", "cindex", "silhouette", "dunn", "kl", "ch",
                 "hartigan", "ball", "ptbiserial", "gap", "gamma", "gplus",
                 "tau", "sdindex", "sdbw", "db", "hubert", "dindex")
#indexScores <- c("duda", "pseudot2", "hubert", "dindex")

# Num breakpoint plots per mutation:
# 3 * 2 * ((10+1)*10/2-1) * 2 = 648
# Num nearest breakpoint plots per mutation:
# 3 * 2 * ((10+1)*10/2-1) * 2 = 648
# Num cluster estimation plots
# 3 * 2 * 18 = 108
# Total plots per mutation
# 648 * 2 + 108 = 1404
runBRCA1 <- function(){
    gene <<- "BRCA1"
    # muts = c("c.-200-?_80+?del", "c.1687C>T", "c.181T>G", "c.211A>G")
    # muts = c("c.-200-?_80+?del", "c.1687C>T", "c.181T>G", "c.211A>G", "c.2475delC",
    #          "c.2681_2682delAA", "c.3319G>T", "c.3331_3334delCAAG", "c.3481_3491del11",
    #          "c.3700_3704del5", "c.3756_3759delGTCT")

    # Computed: touch c.-200-__80+_del c.181T_G c.1687C_T c.211A_G
    muts = c("c.-200-?_80+?del", "c.2681_2682delAA", "c.5123C>A", "c.3319G>T",
             "c.427G>T", "c.3700_3704del5", "c.5333-36_5406+400del510",
             "c.3331_3334delCAAG", "c.5503C>T", "c.3481_3491del11", "c.4035delA",
             "c.4186-?_4357+?dup", "c.1687C>T", "c.2475delC", "c.4327C>T",
             "c.3756_3759delGTCT", "c.4065_4068delTCAA", "c.181T>G", "c.68_69delAG",
             "c.5266dupC")
    distMethod = "firstBreakDist"
    # foreach(i=1:length(muts), .combine = "c") %dopar% {
    #     mut <<- muts[i]
    for (mut in muts){
        mut <<- mut
        mut_name <- gsub(">", "_", mut)
        mut_name <- gsub("\\?", "_", mut_name)
        if (file.exists(paste0("cache/computed/", mut_name))) next
        system(paste0("touch cache/computed/", mut_name), intern = TRUE)

        # init data
        readIndvAndFamData(gene, mut)

        for (clustMethod in c("ward.D", "complete", "ward.D2")){
            for (fam in c(F, T)){
                fam_name <- if (fam) "fam" else "individual"
                filename = paste0("cache/distMatrices/", gene, "-", mut_name, "-distMethod-",
                                  distMethod, "-", fam_name, ".rds")
                dst <- firstBreakDist(fam, filename)
                # numCluster_estimation(max.clust = 12, distMethod, clustMethod, fam_name,
                #                       dist = dst, req_members = 2, index = indexScores)
                for (k in 1:10){
                    for (cutoff in c(0.5)){ # 0.6, 0.7, 0.8
                        print(paste("Running:", gene, mut, distMethod, clustMethod, fam_name, k, cutoff))
                        precompute(gene, mut, dst, distMethod, clustMethod, fam, k, cutoff)
                    }
                }
            }
        }
        system(paste0("touch cache/computed/", mut_name, "-completed"), intern = TRUE)
    }
}

# Parallelization now done using helper bash script
readBRCA1Data()
runBRCA1()



#setup parallel backend to use many processors
#cores=detectCores()
#cl <- makeCluster(cores[1]-1) #not to overload your computer
# clusterExport(cl, list("brca_data.single", "matched_single", "matched_single_plink", "chr_coords", "brca_middle",
#                       "indexScores", "gene", "mut_name", "dst"))
#clusterExport(cl, list("brca1_geno_merged", "brca1_geno_plink", "brca1_pheno_merged", "chr_coords_all", "famHaplotypes_brca1_geno","famHaplotypes_brca1_geno_plink"))
#registerDoParallel(cl)
#
#runBRCA1()
#
#stop cluster
#stopCluster(cl)
