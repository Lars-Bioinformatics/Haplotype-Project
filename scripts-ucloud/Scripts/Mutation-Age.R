#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(parallel)
library(boot)
library(tweedie)
library(statmod)
library(foreach)
library(doMC)

getData <- function(){
    mut = "c.7617+1G>A"; gene = "BRCA2"
    readIndvAndFamData(gene, mut)
    matched_plink_fam = matched_plink_fam[,-1]

    #dst <- firstBreakDist(fam=T)
    dst <- firstBreakDist(fam=F)
    k = 2
    hc <- hiearchical_clustering(dst, clustMethod = "complete", k = k, doHapPlot = F)

    #
    # matched_plink_fam = filter(matched_plink_fam, SNP %in% subset(brca_data, Country == "DENMARK")$Onc_ID)
    #
    # haplotypes <<- findConsensus_plink(matched_plink_fam)
    # breaks <<- findHaplotypeBreaks_plink(matched_plink_fam, haplotypes)
    # chr_pos <<- mapNearestSNPs_plink(matched_plink_fam, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    # TESTING
    #chr_pos <<- subset(chr_pos, Country=="DENMARK")
    #chr_pos <<- subset(chr_pos, Country=="SPAIN")
}

getAncestralHaplotype <- function(ancestral, matched_plink, brca_data, filename = NULL, minSampleBB = 5){
    if (ancestral == "simulatedFounder"){
        haplotypes <- read.table(paste0(filename, "-founder.txt"), header = T)
        #haplotypes2 <<- as.integer(haplotypes[haplotypes$SNP == mut,-1])
        haplotypes2 <<- as.integer(haplotypes[i,-1])
        #matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
        matched_plink2 = matched_plink
        breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
        chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    } else if (ancestral == "knownBreaks"){
        haplotypes2 <<- rep(0, ncol(matched_plink)-1)
        #matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
        matched_plink2 = matched_plink
        breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
        chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    } else if (ancestral == "mostFreqBase"){
        #matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
        matched_plink2 = matched_plink
        haplotypes2 <<- findConsensus_plink(matched_plink2, cutoff = 0.5)
        breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
        chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    } else {
        #matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
        matched_plink2 = matched_plink
        if (ancestral == "branchBoundIndep"){
            res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = T, min_samples = minSampleBB)
        } else { #branchBound
            res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = F, min_samples = minSampleBB)
        }
        haplotypes2 <<- res[[1]]
        breaks <<- res[[2]]
        chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    }
    return(chr_pos)
}

adjustHaploFreqs <- function(chr_pos, matched_plink, chr_coords, epsilon=0.01){
    get_snp_correction <- function(index, side){
        t=1; k=0
        while (t > epsilon){
            k = k + 1
            if (side == "positive_cM"){
                #t = suppressMessages(nrow(plyr::match_df(plyr::count(hapFreqs[,(index-k):index])[,-length((index-k):index)-1], matched_plink[,-1][row_index,(index-k):index]))/nrow(hapFreqs))
                t = suppressMessages(nrow(plyr::match_df(hapFreqs[,(index-k):index], matched_plink[,-1][row_index,(index-k):index]))/nrow(hapFreqs))
            } else {
                t = suppressMessages(nrow(plyr::match_df(hapFreqs[,index:(index+k)], matched_plink[,-1][row_index,index:(index+k)]))/nrow(hapFreqs))
            }
            #print(t)
        }
        #print(k)
        return(k)
    }
    # geno_merged = rbind(famHaplotypes_brca1_geno_plink, famHaplotypes_brca2_geno_plink)
    # dim(geno_merged)
    # index = match(chr_coords$SNP, names(geno_merged))
    # geno_merged_brca1 = geno_merged[,index]
    # dim(geno_merged_brca1)
    # geno_merged_brca1[1:10,1:10]
    
    for (row_index in 1:nrow(chr_pos)){
        print(row_index)
        index = match(chr_pos[row_index, "positive_cM"], chr_coords$cM)
        # c=plyr::count(geno_merged_brca1[,(index-3):index])
        # plyr::match_df(c, matched_plink[,-1][row_index,(index-3):index])$freq/nrow(geno_merged_brca1)
        k = get_snp_correction(index, "positive_cM")
        chr_pos[row_index, "positive_cM"] = max(chr_coords$cM[index-k], brca_cM_middle)
        index = match(chr_pos[row_index, "negative_cM"], chr_coords$cM)
        k = get_snp_correction(index, "negative_cM")
        chr_pos[row_index, "negative_cM"] = min(chr_coords$cM[index+k], brca_cM_middle)
    }
    return(chr_pos)
}

adjustHapLengthsGandolfoAverage <- function(chr_pos, chr_coords, popFreqs, haplotypes2, epsilon=0.01){
    #chr_pos = mutate(chr_pos, org_positive_cM = positive_cM, org_negative_cM = negative_cM)
    ancestral_freqs = sapply(1:nrow(chr_coords), function(i) popFreqs[i, as.character(haplotypes2[i])])
    p = median(ancestral_freqs)^2 + (1-median(ancestral_freqs))^2
    length.of.chromosome = max(chr_coords[,4])-min(chr_coords[,4])
    markers.on.chromosome = length(chr_coords[,4])
    phi = (length.of.chromosome/100)/markers.on.chromosome
    loci = log(epsilon)/log(p)
    cs.correction = loci*phi
    chr_pos$positive_cM <- max(chr_pos$positive_cM-cs.correction, min(chr_coords$cM[chr_coords$cM>brca_cM_middle]))
    chr_pos$negative_cM <- min(chr_pos$negative_cM-cs.correction, max(chr_coords$cM[chr_coords$cM<brca_cM_middle]))
    return(chr_pos)
}

adjustHapLengthsGandolfo <- function(chr_pos, chr_coords, popFreqs, haplotypes2, epsilon=0.01){
    chr_pos = mutate(chr_pos, org_positive_cM = positive_cM, org_negative_cM = negative_cM)
    ancestral_freqs = sapply(1:nrow(chr_coords), function(i) popFreqs[i, as.character(haplotypes2[i])])
    p = median(ancestral_freqs)^2 + (1-median(ancestral_freqs))^2
    k = round(log(epsilon)/log(p), digits = 0)
    for (row_index in 1:nrow(chr_pos)){
        index = match(chr_pos[row_index, "positive_cM"], chr_coords$cM)
        chr_pos[row_index, "positive_cM"] = max(chr_coords$cM[index-k], min(chr_coords$cM[chr_coords$cM>brca_cM_middle]))
        index = match(chr_pos[row_index, "negative_cM"], chr_coords$cM)
        chr_pos[row_index, "negative_cM"] = min(chr_coords$cM[index+k], max(chr_coords$cM[chr_coords$cM<brca_cM_middle]))
    }
    
    # for (side in c("positive_cM", "negative_cM")){
    #     # mi=which(chr_coords$cM == min(chr_pos[,side]))
    #     # ma=which(chr_coords$cM == max(chr_pos[,side]))
    #     # 
    #     # ancestral_freqs <- c()
    #     # j = 1
    #     # for (i in mi:ma){
    #     #     ancestral_freqs[j] = popFreqs[i, as.character(haplotypes2[i])]
    #     #     j = j+1
    #     # }
    #     ancestral_freqs = sapply(1:nrow(chr_coords), function(i) popFreqs[i, as.character(haplotypes2[i])])
    #     p = median(ancestral_freqs)^2 + (1-median(ancestral_freqs))^2
    #     k = round(log(epsilon)/log(p), digits = 0)
    #     for (row_index in 1:nrow(chr_pos)){
    #         index = match(chr_pos[row_index, side], chr_coords$cM)
    #         if (side == "positive_cM"){
    #             chr_pos[row_index, side] = chr_coords$cM[index-k]
    #         } else {
    #             chr_pos[row_index, side] = chr_coords$cM[index+k]
    #         }
    #     }
    # }
    return(chr_pos)
}

adjustHapLengths <- function(chr_pos, chr_coords, popFreqs, matched_plink, epsilon=0.01){
    chr_pos = mutate(chr_pos, org_positive_cM = positive_cM, org_negative_cM = negative_cM)
    
    for (side in c("positive_cM", "negative_cM")){
        for (row_index in 1:nrow(chr_pos)){
            index = match(chr_pos[row_index, side], chr_coords$cM)
            if (matched_plink[row_index, index+1] == 5){
                prob = 0.2
            } else {
                prob = popFreqs[index, as.character(matched_plink[row_index, index+1])]
            }
            while(prob > epsilon){
                index = if (side=="positive_cM") index-1 else index+1
                genotype = as.character(matched_plink[row_index, index+1])
                if (matched_plink[row_index, index+1] == 5){
                    prob = prob * 0.2
                } else {
                    prob = prob * popFreqs[index, genotype]
                }
            }
            if (side=="positive_cM"){
                chr_pos[row_index, side] = max(chr_coords$cM[index], min(chr_coords$cM[chr_coords$cM>brca_cM_middle]))
            } else {
                chr_pos[row_index, side] = min(chr_coords$cM[index], max(chr_coords$cM[chr_coords$cM<brca_cM_middle]))
            }
        }
    }
    
    return(chr_pos)
}

# Mutation age - Method 1
mutationAge_geneticDist <- function(chr_pos, i=NULL){

    if (!is.null(i)){
        chr_pos = chr_pos[i,]
    }
    
    # Version 1
    haplotype_unbroken = median(abs(c(chr_pos$positive_cM, chr_pos$negative_cM) - brca.mut.cM)) / 100
    age = log(0.5) / log(1-haplotype_unbroken)
    #print(paste("Age1:", age))

    return(age)

    # haplotype_unbroken_pos = (median(chr_pos$positive_cM) - brca.mut.cM) / 100
    # haplotype_unbroken_neg = (brca.mut.cM - median(chr_pos$negative_cM)) / 100
    # age.neg = log(0.5) / log(1-haplotype_unbroken_neg)
    # age.neg
    # age.pos = log(0.5) / log(1-haplotype_unbroken_pos)
    # age.pos
    #
    # print(paste("V1:", mean(c(age.pos, age.neg))))

    # Version 2
    # haplo_unbroken = (median(chr_pos$positive_cM) - median(chr_pos$negative_cM)) / 100
    # age = log(0.5) / log(1-haplo_unbroken)
    # print(paste("V2:", age))
}

mutationAge_method1_mean <- function(chr_pos, brca.mut.cM){
    
    # Version 1
    haplotype_unbroken = mean(abs(c(chr_pos$positive_cM, chr_pos$negative_cM) - brca.mut.cM)) / 100
    haplotype_unbroken
    
    age = log(0.5) / log(1-haplotype_unbroken)
    #print(paste("Age1:", age))
    
    return(age)
}

mutationAge_method1_median_with_bootstrap <- function(chr_pos){
    set.seed(626)
    #cl = makeCluster(detectCores())
    #clusterExport(cl, c("marker_before_break", "chr_coords", "brca.mut.cM"))
    #bootcorr <- boot(data = chr_pos, statistic = mutationAge_geneticDist, R=100, cl = cl, parallel = "multicore", ncpus = detectCores())
    bootcorr <- boot(data = chr_pos, statistic = mutationAge_geneticDist, R=100)
    bootcorr
    #stopCluster(cl)
    #plot(bootcorr)
    #ci=boot.ci(boot.out = bootcorr, type = c("norm", "basic", "perc"))
    if (sd(bootcorr$t)==0){
        age = c(bootcorr$t0, NA, NA)
    } else {
        ci=boot.ci(boot.out = bootcorr, type = c("norm"))
        age = c(ci$t0, ci$normal[2], ci$normal[3])
    }
    print(age)
    return(age)
}

# Mutation age - Method 2
mutationAge_physicalDist <- function(chr_pos){
    #age <- 100/(mean(chr_pos$haplo_length)/1000000)
    #age <- 100/(median(chr_pos$haplo_length)/1000000)
    age <- 100/(mean(abs(c(chr_pos$positive_pos, chr_pos$negative_pos) - brca_middle))/1000000)
    #print(paste("Age2:", age))
    return(age)
}

# right_prob brca.mut.cM - 1.165896 = 29.03979
# left_prob 46.17665 - brca.mut.cM = 15.97097

# Mutation age - Method 3
mutationAge_method3 <- function(chr_pos, i=NULL){
    #age <- 100/median(chr_pos$haplo_length_cM)
    if (is.null(i)){
        arm.lengths = abs(c(chr_pos$positive_cM, chr_pos$negative_cM) - brca.mut.cM)
    } else {
        arm.lengths = chr_pos[i]
    }
    age <- 100/mean(arm.lengths)
    
    # Normal dist CI
    # https://www.cyclismo.org/tutorial/R/confidence.html
    # CI_min = 100/(mean(arm.lengths)+qnorm(0.975)*sd(arm.lengths)/sqrt(length(arm.lengths)))
    # CI_max = 100/(mean(arm.lengths)-qnorm(0.975)*sd(arm.lengths)/sqrt(length(arm.lengths)))
    # 
    # # Exponential dist CI
    # # https://stats.stackexchange.com/questions/155805/confidence-interval-for-exponential-distribution
    # CI_min = 100/(mean(arm.lengths)+qnorm(0.975)*sqrt(mean(arm.lengths)^2/length(arm.lengths)))
    # CI_max = 100/(mean(arm.lengths)-qnorm(0.975)*sqrt(mean(arm.lengths)^2/length(arm.lengths)))
    
    # return(c(age, CI_min, CI_max))
    return(age)
}

mutationAge_method3_with_bootstrap <- function(chr_pos){
    set.seed(626)
    #cl = makeCluster(detectCores())
    #clusterExport(cl, c("marker_before_break", "chr_coords", "brca.mut.cM"))
    #bootcorr <- boot(data = chr_pos, statistic = mutationAge_method3, R=100, cl = cl, parallel = "multicore", ncpus = detectCores())
    arm.lengths = abs(c(chr_pos$positive_cM, chr_pos$negative_cM) - brca.mut.cM)
    bootcorr <- boot(data = arm.lengths, statistic = mutationAge_method3, R=100)
    bootcorr
    #stopCluster(cl)
    #plot(bootcorr)
    #ci=boot.ci(boot.out = bootcorr, type = c("norm", "basic", "perc"))
    if (sd(bootcorr$t)==0){
        age = c(bootcorr$t0, NA, NA)
    } else {
        ci=boot.ci(boot.out = bootcorr, type = c("norm"))
        age = c(ci$t0, ci$normal[2], ci$normal[3])
    }
    print(age)
    return(age)
}

mutationAge_method3_median <- function(chr_pos){
    #age <- 100/median(chr_pos$haplo_length_cM)
    age <- 100/median(abs(c(chr_pos$positive_cM, chr_pos$negative_cM) - brca.mut.cM))
    return(age)
}


# Mutation age - Method 4 upper
mutationAge_method4upper <- function(chr_pos, brca.mut.cM){

    haplotype_unbroken = median(abs(chr_pos$positive_cM - brca.mut.cM)) / 100
    haplotype_unbroken

    age = log(0.5) / log(1-haplotype_unbroken)

    return(age)
}

# Mutation age - Method 4 under
mutationAge_method4under <- function(chr_pos, brca.mut.cM){

    haplotype_unbroken = median(abs(chr_pos$negative_cM - brca.mut.cM)) / 100
    haplotype_unbroken

    age = log(0.5) / log(1-haplotype_unbroken)

    return(age)
}

getNumAbove <- function(num, above=0.1){
    if (num==0){
        #stop("getNumAbove: Value is zero")
        return(c(308, 0.1))
    }
    pre_zeros = 0
    while(num < above){
        num = num * 10
        pre_zeros = pre_zeros + 1
    }
    return(c(pre_zeros, num))
}

powerNorm <- function(num, power, above=0.1){
    if (num==0){stop("powerNorm: Value is zero")}
    if (power==0){return(c(0,1))}
    if (power==1){return(c(0,num))}
    pre_zeros = 0
    res = num
    for (i in 2:power){
        res = res * num
        while (res < above){
            res = res * 10
            pre_zeros = pre_zeros + 1
        }
    }
    return(c(pre_zeros, res))
    
}

productNorm <- function(num, above=0.1){
    res = 1
    pre_zeros = 0
    for (i in num){
        res = res * i
        while(res < above){
            res = res * 10
            pre_zeros = pre_zeros + 1
        }
    }
    return(c(pre_zeros, res))
}

getUnbrokenMarker <- function(chr_pos, matched_plink2, chr_coords, haplotypes){
    matched = matched_plink2[,-1]
    unbroken_pos = c()
    for (i in 1:nrow(chr_pos)){
        #index = which(names(matched) == chr_coords[match(chr_pos$positive_pos[i], chr_coords$position_b37),1])
        index = match(chr_pos$positive_pos[i], chr_coords$position_b37)
        j = 1
        repeat{
            index2 = index-j
            cmp = matched[i, index2] - haplotypes[index2]
            if (haplotypes[index2] != 1 && cmp == 0){
                unbroken_cM = chr_coords$cM[index2]
                if (unbroken_cM - chr_coords$cM[index] != 0) break
            }
            j = j + 1
        }
        unbroken_pos[i] = unbroken_cM
    }
    print(unbroken_pos)
    chr_pos["unbroken_pos_cM"] = unbroken_pos
    
    unbroken_neg = c()
    for (i in 1:nrow(chr_pos)){
        index = match(chr_pos$negative_pos[i], chr_coords$position_b37)
        j = 1
        repeat{
            index2 = index+j
            cmp = matched[i, index2] - haplotypes[index2]
            print(index2)
            print(haplotypes[index2])
            if (haplotypes[index2] != 1 && cmp == 0){
                unbroken_cM = chr_coords$cM[index2]
                if (unbroken_cM - chr_coords$cM[index] != 0) break
            }
            j = j + 1
        }
        unbroken_neg[i] = unbroken_cM
    }
    print(unbroken_neg)
    chr_pos["unbroken_neg_cM"] = unbroken_neg
    
    return(chr_pos)
}

computeGenerationLikelihood_chanceSharing <- function(chr_pos, gen){
    print(paste("generation:", gen))
    prob_left = sapply(chr_pos$sample_id, g_chanceSharing, n=gen, right_side=F, chr_pos = chr_pos)
    prob_right = sapply(chr_pos$sample_id, g_chanceSharing, n=gen, right_side=T, chr_pos = chr_pos)
    
    likelihood = prod(prob_left, prob_right)
    likelihood_norm = getNumAbove(likelihood)
    
    res = c(gen, likelihood_norm[1], likelihood_norm[2])
    names(res) = c("Generations", "zeros", "likelihood")
    return(res)
}

computeGenerationLikelihood <- function(chr_pos, gen){
    # Left side
    print(paste("generation:", gen))
    #prob_left = sapply(chr_pos$negative_cM, f, n=gen, right_side=F)
    prob_left = sapply(chr_pos$sample_id, g, n=gen, right_side=F, chr_pos = chr_pos)
    #prob_left = sapply(chr_pos$sample_id, g_chanceSharing, n=gen, right_side=F, chr_pos = chr_pos)
    prob_left_norm = sapply(prob_left, getNumAbove)
    
    # y=nrow(chr_pos)
    # fv = sapply(sapply(chr_pos$negative_cM, f, n=gen, right_side=F), powerNorm, power=y-1, above=0.0001)
    # Sv = sapply(sapply(chr_pos$negative_cM, S, n=gen), getNumAbove)
    # prob_left_2 = sapply(Sv[2,] * fv[2,], getNumAbove)
    #chr_pos[,paste0("left_",gen,"_zero")] = prob_left_norm[1,]
    #chr_pos[,paste0("left_",gen,"_prob")] = prob_left_norm[2,]
    # OLD
    ##zeros_left = sum(prob_left_norm[1,])
    ##likelihood_left = prod(prob_left_norm[2,])
    
    likelihood_left = productNorm(prob_left_norm[2,])
    zeros_left = sum(prob_left_norm[1,]) + likelihood_left[1]
    likelihood_left = likelihood_left[2]
    
    # Right side
    #prob_right = sapply(chr_pos$positive_cM, f, n=gen, right_side=T)
    prob_right = sapply(chr_pos$sample_id, g, n=gen, right_side=T, chr_pos = chr_pos)
    #prob_right = sapply(chr_pos$sample_id, g_chanceSharing, n=gen, right_side=T, chr_pos = chr_pos)
    prob_right_norm = sapply(prob_right, getNumAbove)
    # chr_pos[,paste0("right_",gen,"_zero")] = prob_right_norm[1,]
    # chr_pos[,paste0("right_",gen,"_prob")] = prob_right_norm[2,]
    #chr_pos <<- chr_pos
    #zeros_right = sum(prob_right_norm[1,])
    #likelihood_right = prod(prob_right_norm[2,])
    
    likelihood_right = productNorm(prob_right_norm[2,])
    zeros_right = sum(prob_right_norm[1,]) + likelihood_right[1]
    likelihood_right = likelihood_right[2]
    
    # Combine 
    likelihood = prod(likelihood_left, likelihood_right)
    likelihood_norm = getNumAbove(likelihood)
    zeros_total = zeros_left + zeros_right + likelihood_norm[1]
    likelihood_final = likelihood_norm[2]
    #print(paste("gen:", gen, "zeros:", zeros_total, "likelihood:", likelihood_final))
    
    res = c(gen, zeros_total, likelihood_final)
    names(res) = c("Generations", "zeros", "likelihood")
    return(res)
}

computeGenerationLikelihood_newMethod5 <- function(chr_pos, gen){
    # Left side
    print(paste("generation:", gen))
    #prob_left = sapply(chr_pos$negative_cM, f, n=gen, right_side=F)
    prob_left = sapply(chr_pos$sample_id, g_new, n=gen, right_side=F, chr_pos = chr_pos)
    prob_left_norm = sapply(prob_left[2,], getNumAbove)
    
    # y=nrow(chr_pos)
    # fv = sapply(sapply(chr_pos$negative_cM, f, n=gen, right_side=F), powerNorm, power=y-1, above=0.0001)
    # Sv = sapply(sapply(chr_pos$negative_cM, S, n=gen), getNumAbove)
    # prob_left_2 = sapply(Sv[2,] * fv[2,], getNumAbove)
    #chr_pos[,paste0("left_",gen,"_zero")] = prob_left_norm[1,]
    #chr_pos[,paste0("left_",gen,"_prob")] = prob_left_norm[2,]
    # OLD
    ##zeros_left = sum(prob_left_norm[1,])
    ##likelihood_left = prod(prob_left_norm[2,])
    
    likelihood_left = productNorm(prob_left_norm[2,])
    zeros_left = sum(prob_left[1,]) + sum(prob_left_norm[1,]) + likelihood_left[1]
    likelihood_left = likelihood_left[2]
    
    # Right side
    #prob_right = sapply(chr_pos$positive_cM, f, n=gen, right_side=T)
    prob_right = sapply(chr_pos$sample_id, g_new, n=gen, right_side=T, chr_pos = chr_pos)
    prob_right_norm = sapply(prob_right[2,], getNumAbove)
    # chr_pos[,paste0("right_",gen,"_zero")] = prob_right_norm[1,]
    # chr_pos[,paste0("right_",gen,"_prob")] = prob_right_norm[2,]
    #chr_pos <<- chr_pos
    #zeros_right = sum(prob_right_norm[1,])
    #likelihood_right = prod(prob_right_norm[2,])
    
    likelihood_right = productNorm(prob_right_norm[2,])
    zeros_right = sum(prob_right[1,]) + sum(prob_right_norm[1,]) + likelihood_right[1]
    likelihood_right = likelihood_right[2]
    
    # Combine 
    likelihood = prod(likelihood_left, likelihood_right)
    likelihood_norm = getNumAbove(likelihood)
    zeros_total = zeros_left + zeros_right + likelihood_norm[1]
    likelihood_final = likelihood_norm[2]
    #print(paste("gen:", gen, "zeros:", zeros_total, "likelihood:", likelihood_final))
    
    res = c(gen, zeros_total, likelihood_final)
    names(res) = c("Generations", "zeros", "likelihood")
    return(res)
}

computeGenerationLikelihoodTest <- function(chr_pos, gen){
    print(paste("generation:", gen))
    prob_left = sapply(chr_pos$negative_cM, f, n=gen)#^(nrow(chr_pos))
    #prob_left = sapply(chr_pos$negative_cM, f, n=gen)
    prob_left_norm = sapply(prob_left, getNumAbove)
    zeros_left = sum(prob_left_norm[1,])
    prob_left_norm = sapply(prob_left_norm[2,]^(nrow(chr_pos)), getNumAbove)
    zeros_left = zeros_left+sum(prob_left_norm[1,])
    prob_left_norm = prob_left_norm[2,]
    likelihood_left = nrow(chr_pos)*sapply(chr_pos$negative_cM, S, n=gen)*prob_left_norm
    likelihood_left = sapply(likelihood_left, getNumAbove)
    zeros_left = zeros_left+sum(likelihood_left[1,])
    likelihood_left = likelihood_left[2,]
    #zeros_left = sum(prob_left_norm[1,])
    #likelihood_left = prod(prob_left_norm[2,])
    
    #prob_right = nrow(chr_pos)*sapply(chr_pos$positive_cM, S, n=gen)*sapply(chr_pos$positive_cM, f, n=gen)^(nrow(chr_pos))
    #prob_right = sapply(chr_pos$positive_cM, f, n=gen)
    prob_right = sapply(chr_pos$positive_cM, f, n=gen)
    prob_right_norm = sapply(prob_right, getNumAbove)
    zeros_right = sum(prob_right_norm[1,])
    prob_right_norm = sapply(prob_right_norm[2,]^(nrow(chr_pos)), getNumAbove)
    zeros_right = zeros_right+sum(prob_right_norm[1,])
    prob_right_norm = prob_right_norm[2,]
    likelihood_right = nrow(chr_pos)*sapply(chr_pos$positive_cM, S, n=gen)*prob_right_norm
    likelihood_right = sapply(likelihood_right, getNumAbove)
    zeros_right = zeros_right+sum(likelihood_right[1,])
    likelihood_right = likelihood_right[2,]
    #likelihood_right = prod(prob_right_norm[2,])
    
    likelihood = prod(likelihood_left, likelihood_right)
    likelihood_norm = getNumAbove(likelihood)
    zeros_total = zeros_left + zeros_right + likelihood_norm[1]
    likelihood_final = likelihood_norm[2]
    #print(paste("gen:", gen, "zeros:", zeros_total, "likelihood:", likelihood_final))
    
    # prob_left = (1-(abs(chr_pos$negative_cM - brca.mut.cM)/100))^gen
    # likelihood_left = prod(prob_left)
    # prob_right = (1-(abs(chr_pos$negative_cM - brca.mut.cM)/100))^gen
    # likelihood_right = prod(prob_right)
    # likelihood = prod(likelihood_left, likelihood_right)
    # print(paste("gen:", gen, "likelihood:", likelihood))
    
    res = c(gen, zeros_total, likelihood_final)
    names(res) = c("Generations", "zeros", "likelihood")
    return(res)
}



normalizeAgeEstimations <- function(ageMatRow, numZeros){
    #numZeros = min(ageMat$zeros)
    # for (row in 1:nrow(ageMat)){
    #     
    # }
    zeros = ageMatRow[2]
    likelihood = ageMatRow[3]
    while(zeros > numZeros){
        likelihood = likelihood / 10
        zeros = zeros - 1
    }
    return(c(ageMatRow[1], zeros, likelihood))
}

g <- function(sid,n,right_side,chr_pos){
    #print(sid)
    j = which(sid == chr_pos$sample_id)
    i = 1
    x = if (right_side) chr_pos[j, "positive_cM"] else chr_pos[j, "negative_cM"]
    
    if (x == max(chr_coords$cM) || x == min(chr_coords$cM)){
        prob = (1-abs(x-brca.mut.cM)/100)^n
        return(prob)
    }
    
    index = match(x, chr_coords$cM)
    count = 0
    repeat{
        index2 = if (right_side) index-i else index+i
        #index2 = if (right_side) index+i else index-i
        cmp = matched_plink2[j, index2] - haplotypes2[index2]
        if (haplotypes2[index2] != 1 && cmp == 0){
            count = count + 1
            if (count >= 1){
                x2 = chr_coords$cM[index2]
                # Check for genetic distance above 0
                if (x2 - x != 0) break
            }
        }
        i = i + 1
    }
    
    #print((1-abs(x2-brca.mut.cM)/100)^n * (1-(1-abs(x-x2)/100)^n))
    #print(1-((1-abs(x2-brca.mut.cM)/100)^n))
    prob = (1-abs(x2-brca.mut.cM)/100)^n * (1-(1-abs(x-x2)/100)^n) #+ 10*abs(x2-brca.mut.cM)/100
    #prob = (1-abs(x2-brca.mut.cM)/100)^(n)
    #prob = 1-(1-abs(x-x2)/100)^n
    #prob = (1-abs(x2-brca.mut.cM)/100)^n * (n*(abs(x-x2)/100)*(1-abs(x-x2)/100)^(n-1))
    
    # print(x); print(x2); print(n)
    #prob = S(x,n)*abs(S(x2,n) - S(x,n))
    return(prob)
}

g_chanceSharing <- function(sid,n,right_side,chr_pos){
    j = which(sid == chr_pos$sample_id)
    i = 1
    x = if (right_side) chr_pos[j, "positive_cM"] else chr_pos[j, "negative_cM"]
    
    if (x == max(chr_coords$cM) || x == min(chr_coords$cM)){
        prob = (1-abs(x-brca.mut.cM)/100)^n
        return(prob)
    }
    
    final_prob = 0
    
    brca_index = filter(chr_coords, cM < brca.mut.cM) %>% nrow()
    if (right_side) brca_index = brca_index + 5 else brca_index = brca_index - 4
    if (right_side) brca_index2 = brca_index + 1 else brca_index2 = brca_index - 1
    index = match(x, chr_coords$cM)
    for (k in brca_index2:index){
    #for (k in (index-20):index){
    #for (k in brca_index2:(brca_index2-1)){
        x3 = chr_coords$cM[k]
        count = 0
        repeat{
            index2 = if (right_side) k-i else k+i
            if (right_side && index2 < brca_index){
                x2 = chr_coords$cM[index2 + 1]
                break
            }
            if (!right_side && index2 > brca_index){
                x2 = chr_coords$cM[index2 - 1]
                break
            }
            #index2 = if (right_side) index+i else index-i
            cmp = matched_plink2[j, index2] - haplotypes2[index2]
            if (haplotypes2[index2] != 1 && cmp == 0){
                count = count + 1
                if (count >= 1){
                    x2 = chr_coords$cM[index2]
                    # Check for genetic distance above 0
                    if (x2 - x3 != 0) break
                }
            }
            i = i + 1
        }
        
        #print((1-abs(x2-brca.mut.cM)/100)^n * (1-(1-abs(x-x2)/100)^n))
        #print(1-((1-abs(x2-brca.mut.cM)/100)^n))
        prob = (1-abs(x2-brca.mut.cM)/100)^n * (1-(1-abs(x3-x2)/100)^n) #+ 10*abs(x2-brca.mut.cM)/100
        #prob = (1-abs(x2-brca.mut.cM)/100)^(n)
        #prob = 1-(1-abs(x-x2)/100)^n
        #prob = (1-abs(x2-brca.mut.cM)/100)^n * (n*(abs(x-x2)/100)*(1-abs(x-x2)/100)^(n-1))
        
        #freqs = prod(sapply((k+1):(index+1), function(i) popFreqs[i, as.character(matched_plink2[matched_plink2$SNP==sid, i])]))
        freqs = ancestral_probs[(k+1):(index+1)]
        res = prod(c(freqs,prob))
        
        
        final_prob = final_prob + res
    }
    
    # print(x); print(x2); print(n)
    #prob = S(x,n)*abs(S(x2,n) - S(x,n))
    return(final_prob)
}

g_new <- function(sid,n,right_side,chr_pos){
    #print(sid)
    j = which(sid == chr_pos$sample_id)
    i = 1
    x = if (right_side) chr_pos[j, "positive_cM"] else chr_pos[j, "negative_cM"]
    
    if (x == max(chr_coords$cM) || x == min(chr_coords$cM)){
        if (right_side){
            markers = chr_coords$cM[brca.mut.cM <= chr_coords$cM & chr_coords$cM <= x]
        } else {
            markers = chr_coords$cM[brca.mut.cM >= chr_coords$cM & chr_coords$cM >= x]
        }
        intervals = (markers[2:(length(markers))]-markers[1:(length(markers)-1)])/100
        # pNorm = sapply((1-intervals), powerNorm, power = n)
        # l = productNorm(pNorm[2,])
        # zeros = sum(pNorm[1,]) + l[1]
        # prob = l[2]
        l = productNorm((1-intervals)^n)
        zeros = l[1]
        prob = l[2]
        return(c(zeros,prob))
    }
    
    index = match(x, chr_coords$cM)
    count = 0
    repeat{
        index2 = if (right_side) index-i else index+i
        #index2 = if (right_side) index+i else index-i
        cmp = matched_plink2[j, index2] - haplotypes2[index2]
        if (haplotypes2[index2] != 1 && cmp == 0){
            count = count + 1
            if (count >= 1){
                x2 = chr_coords$cM[index2]
                # Check for genetic distance above 0
                if (x2 - x != 0) break
            }
        }
        i = i + 1
    }
    
    if (right_side){
        markers = chr_coords$cM[brca.mut.cM <= chr_coords$cM & chr_coords$cM <= x2]
    } else {
        markers = chr_coords$cM[brca.mut.cM >= chr_coords$cM & chr_coords$cM >= x2]
    }
    if (length(markers) == 0){
        #return(c(0, (1-(1-abs(x-brca.mut.cM)/100)^n)))
        return(c(0, (n*(abs(x-brca.mut.cM)/100)*(1-abs(x-brca.mut.cM)/100)^(n-1))))
    }
    intervals = (markers[2:(length(markers))]-markers[1:(length(markers)-1)])/100
    #pNorm = sapply((1-intervals), powerNorm, power = n)
    #l = productNorm(pNorm[2,])
    #zeros = sum(pNorm[1,]) + l[1]
    #l = productNorm((1-intervals)^n)
    l = productNorm((1-intervals)^n)
    zeros = l[1]
    #prob = l[2] * (n*(abs(x-x2)/100)*(1-abs(x-x2)/100)^(n-1))
    prob = l[2] * (1-(1-abs(x-x2)/100)^n)
    
    #print((1-abs(x2-brca.mut.cM)/100)^n * (1-(1-abs(x-x2)/100)^n))
    #print(1-((1-abs(x2-brca.mut.cM)/100)^n))
    #prob = (1-abs(x2-brca.mut.cM)/100)^n * (1-(1-abs(x-x2)/100)^n) #+ 10*abs(x2-brca.mut.cM)/100
    #prob = (1-abs(x2-brca.mut.cM)/100)^(n)
    #prob = 1-(1-abs(x-x2)/100)^n
    
    # print(x); print(x2); print(n)
    #prob = S(x,n)*abs(S(x2,n) - S(x,n))
    return(c(zeros, prob))
}

#testSimAge()
#mutationAge_method5(chr_pos, brca.mut.cM, start=10, stop=200, step=20)

mutationAge_method5 <- function(chr_pos, brca.mut.cM, start=10, stop=200, step=10){
    
    generations = seq(start,stop,step)
    #generations = seq(50,150,10)
    #generations = seq(10,100,1)
    ageMat = as.data.frame(t(sapply(generations, computeGenerationLikelihood, chr_pos=chr_pos)))
    #ageMat = as.data.frame(t(sapply(generations, computeGenerationLikelihood_chanceSharing, chr_pos=chr_pos)))
    #print(ageMat)
    
    age <- ageMat %>% filter(zeros == min(zeros)) %>% slice(which.max(likelihood)) %>% .[,1]
    
    ageMat_norm = as.data.frame(t(apply(ageMat, 1, normalizeAgeEstimations, numZeros = min(ageMat$zeros))))
    print(ageMat_norm)
    print(ageMat_norm[which.max(ageMat_norm$likelihood),])
    plot(ageMat_norm[,1], ageMat_norm[,3], col = 'red', type = "o")
    
    #age = ageMat_norm[which.max(ageMat_norm$likelihood),1]
    return(age)
}

mutationAge_newMethod5 <- function(chr_pos, brca.mut.cM, start=10, stop=200, step=10){
    
    generations = seq(start,stop,step)
    #generations = seq(50,150,10)
    #generations = seq(10,100,1)
    ageMat = as.data.frame(t(sapply(generations, computeGenerationLikelihood_newMethod5, chr_pos=chr_pos)))
    #print(ageMat)
    
    age <- ageMat %>% filter(zeros == min(zeros)) %>% slice(which.max(likelihood)) %>% .[,1]
    
    ageMat_norm = as.data.frame(t(apply(ageMat, 1, normalizeAgeEstimations, numZeros = min(ageMat$zeros))))
    print(ageMat_norm)
    print(ageMat_norm[which.max(ageMat_norm$likelihood),])
    plot(ageMat_norm[,1], ageMat_norm[,3], col = 'red', type = "o")
    
    #age = ageMat_norm[which.max(ageMat_norm$likelihood),1]
    return(age)
}

#system.time(mutationAge_method5(chr_pos, brca.mut.cM, start = 10, stop = 5000, step = 1))

mutationAge_method6 <- function(chr_pos, brca_cM_middle, confidence.coefficient = 0.95, chance.sharing.correction = F, epsilon = 0.01,
                                chr_coords = chr_coords, popFreqs = popFreqs, haplotypes2 = haplotypes2){
    age = gandolfo_age(l.lengths = chr_pos$negative_cM,
                       r.lengths = chr_pos$positive_cM,
                       confidence.coefficient = confidence.coefficient,
                       chance.sharing.correction = chance.sharing.correction,
                       e = epsilon,
                       chr_coords = chr_coords, 
                       popFreqs = popFreqs, 
                       haplotypes2 = haplotypes2)
    print(paste("Gandolfo:", age))
    return(age)
}

## Method 7 - ...
f <- function(x,n,right_side) {
    i = 1
    index = match(x, chr_coords$cM)
    
    if (x == max(chr_coords$cM) || x == min(chr_coords$cM)){
        prob = (1-abs(x-brca.mut.cM)/100)^n
        return(prob)
    }
    
    repeat{
        x2 = if (right_side) chr_coords$cM[index-i] else chr_coords$cM[index+i]
        if ((x-x2) != 0) break
        i = i + 1
    }
    # print(S(x2, n) - S(x,n))
    # if ((S(x2, n) - S(x,n))==0){
    #     print(paste(x,n,right_side))
    # }
    return(abs(S(x2, n) - S(x,n)))
}

S <- function(x,n) {
    # Kosambi map: http://www.crypticlineage.net/lib/Kosambi.pdf
    recombFrac <- function(x){1/2 * ((exp(4*x/100)-1)/(exp(4*x/100)+1))}
    return((1-recombFrac(abs(x - brca.mut.cM)))^n)
    
    # return((1-(abs(x - brca.mut.cM)/100))^n)
}

# Compute f with no grouping
computeMethod7Likelihood <- function(chr_pos, gen){
    # Left side
    print(paste("generation:", gen))
    prob_left = sapply(chr_pos$negative_cM, f, n=gen, right_side=F)
    prob_left_norm = sapply(prob_left, getNumAbove)
    
    likelihood_left = productNorm(prob_left_norm[2,])
    zeros_left = sum(prob_left_norm[1,]) + likelihood_left[1]
    likelihood_left = likelihood_left[2]
    
    # Right side
    prob_right = sapply(chr_pos$positive_cM, f, n=gen, right_side=T)
    prob_right_norm = sapply(prob_right, getNumAbove)

    likelihood_right = productNorm(prob_right_norm[2,])
    zeros_right = sum(prob_right_norm[1,]) + likelihood_right[1]
    likelihood_right = likelihood_right[2]
    
    # Combine 
    likelihood = prod(likelihood_left, likelihood_right)
    likelihood_norm = getNumAbove(likelihood)
    zeros_total = zeros_left + zeros_right + likelihood_norm[1]
    likelihood_final = likelihood_norm[2]

    res = c(gen, zeros_total, likelihood_final)
    names(res) = c("Generations", "zeros", "likelihood")
    return(res)
}

mutationAge_method7 <- function(chr_pos, brca.mut.cM, start=10, stop=200, step=10){
    
    generations = seq(start,stop,step)
    #generations = seq(50,150,10)
    #generations = seq(10,100,1)
    ageMat = as.data.frame(t(sapply(generations, computeMethod7Likelihood, chr_pos=chr_pos)))
    #print(ageMat)
    
    age <- ageMat %>% filter(zeros == min(zeros)) %>% slice(which.max(likelihood)) %>% .[,1]
    
    ageMat_norm = as.data.frame(t(apply(ageMat, 1, normalizeAgeEstimations, numZeros = min(ageMat$zeros))))
    print(ageMat_norm)
    print(ageMat_norm[which.max(ageMat_norm$likelihood),])
    plot(ageMat_norm[,1], ageMat_norm[,3], col = 'red', type = "o")
    
    #age = ageMat_norm[which.max(ageMat_norm$likelihood),1]
    return(age)
}

mutationAge_method7_group <- function(chr_pos, brca.mut.cM, start=10, stop=200, step=10){
    
    longest_hap = chr_pos %>% count(positive_pos) %>% filter(n>1) %>% summarise(max=max(positive_pos))
    group1_right = chr_pos %>% filter(positive_pos >= longest_hap$max)
    group1_right$positive_pos = longest_hap$max
    group2_right = chr_pos %>% filter(positive_pos < longest_hap$max)
    
    longest_hap = chr_pos %>% count(negative_pos) %>% filter(n>1) %>% summarise(min=min(negative_pos))
    group1_left = chr_pos %>% filter(negative_pos <= longest_hap$min)
    group1_left$negative_pos = longest_hap$min
    group2_left = chr_pos %>% filter(negative_pos > longest_hap$min)
    
    generations = seq(start,stop,step)
    #generations = seq(50,150,10)
    #generations = seq(10,100,1)
    ageMat = as.data.frame(t(sapply(generations, computeMethod7GroupLikelihood, group1_right = group1_right, 
                                    group2_right = group2_right, group1_left = group1_left, group2_left = group2_left)))
    #print(ageMat)
    
    age <- ageMat %>% filter(zeros == min(zeros)) %>% slice(which.max(likelihood)) %>% .[,1]
    
    ageMat_norm = as.data.frame(t(apply(ageMat, 1, normalizeAgeEstimations, numZeros = min(ageMat$zeros))))
    print(ageMat_norm)
    print(ageMat_norm[which.max(ageMat_norm$likelihood),])
    plot(ageMat_norm[,1], ageMat_norm[,3], col = 'red', type = "o")
    
    #age = ageMat_norm[which.max(ageMat_norm$likelihood),1]
    return(age)
}

mutationAge_method8_tweedie <- function(chr_pos, brca.mut.cM, phi = 1/2){
    # library(tweedie)
    # library(statmod)
    
    chr_pos <- chr_pos %>% mutate(positive_cM_adj = positive_cM - brca.mut.cM, negative_cM_adj = abs(negative_cM - brca.mut.cM))
    arm.lengths = c(chr_pos$positive_cM_adj, chr_pos$negative_cM_adj)
    
    ll.tweedie <- function(par, y, phi=1/2) {
        # print(par)
        #out <- sum(log(dtweedie(y = y, mu = par[1], phi = exp(par[2]), power = par[3])))
        # out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = 1/2, power = 2)))
        out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = par[2], power = par[3])))
        # out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = phi, power = par[2])))
        # out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = 1/2, power = par[2])))
        #out <- sum(log(tweedie.dev(y = y, mu = 2/par[1], power = par[2])))
        #out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = exp(par[2]), power = par[3])))
    
        # phi given from definition, see Dunn and Smyth (2005)
        # out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = 1/((2-par[2])/(par[2]-1)), power = par[2])))
        # print(out)
        return(out)
    }

    # tt.tweedie = optim(par = c(2/mean(arm.lengths), 1.5), fn = ll.tweedie, y = arm.lengths, phi = phi, control = list(fnscale = -1), hessian = TRUE,
    # lower = c(0.01,1.01), upper = c(100,2.00), method = "L-BFGS-B" )
    tt.tweedie = optim(par = c(2/mean(arm.lengths), 1/2, 2), fn = ll.tweedie, y = arm.lengths, control = list(fnscale = -1), hessian = TRUE,
                       lower = c(0.01,0.02,1.02), upper = c(100000,2.9,2.00), method = "L-BFGS-B")
    # tt.tweedie2 = nlminb(start = c(2/mean(arm.lengths), 1/2, 1.5), objective = ll.tweedie, y = arm.lengths, hessian = TRUE,
    #                     lower = c(0.01,0.01,1.01), upper = c(100,100,1.99))
    # tt.tweedie = optim(par = c(2/mean(arm.lengths)), fn = ll.tweedie, y = arm.lengths, control = list(fnscale = -1), hessian = TRUE,
    # lower = c(0.01), upper = c(100), method = "L-BFGS-B" )
    
    std.error.tw = sqrt(abs(diag(solve(-tt.tweedie$hessian))))
    
    # df = data.frame("Ic_Min" = 100/(2/(tt.tweedie$par[1] - qnorm(0.975) *std.error.tw[1])),
    #                 "Age_estimate" = 100/(2/tt.tweedie$par[1]),
    #                 "Ic_Max" = 100/(2/(tt.tweedie$par[1] + qnorm(0.975) *std.error.tw[1])))
    # print(df)
    
    age = 100/(2/tt.tweedie$par[1])
    if (tt.tweedie$par[3]==2){
        n = nrow(chr_pos)
        Ic_Min = age * qgamma(shape=2*n,scale=1/(2*n),((1-0.95)/2))
        Ic_Max = age * qgamma(shape=2*n,scale=1/(2*n),(0.95+(1-0.95)/2))
    } else {
        Ic_Min = 100/(2/(tt.tweedie$par[1] - qnorm(0.975) * std.error.tw[1]))
        Ic_Max = 100/(2/(tt.tweedie$par[1] + qnorm(0.975) * std.error.tw[1]))
        # Ic_Min = 100/(2/(tt.tweedie$par[1] - 1.96 *std.error.tw[1]))
        # Ic_Max = 100/(2/(tt.tweedie$par[1] + 1.96 *std.error.tw[1]))
        100/(2/qtweedie(p = ((1-0.95)/2), power = tt.tweedie$par[3], mu = 2/tt.tweedie$par[1], phi = tt.tweedie$par[2]) * std.error.tw[1])
        qtweedie(p = (0.95+(1-0.95)/2), power = tt.tweedie$par[3], mu = 2/tt.tweedie$par[1], phi = tt.tweedie$par[2])
    }
    
    #age_alt = 100/(2/tt.tweedie$par[1])^tt.tweedie$par[2]
    
    return(c(round(age, 1), round(Ic_Min, 1), round(Ic_Max, 1), tt.tweedie$par))
}

mutationAge_tweedie <- function(){
    library(tweedie)
    library(statmod)
    
    get_sim_data <- function(gen = 110, i = 1){
        #filename = paste0("classify-simulated-population-age/BRCA1_simulated-knownBreaks-starGenealogy-100_samples-10_simulations-generations_10_2510_step50-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-10_simulations-generations_",gen,"-seed_42")
        filename = paste0("classify-simulated-population-age/BRCA1_simulated-knownBreaks-starGenealogy-100_samples-20_simulations-generations_10_3010_step50-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_",gen,"-seed_42")
        brca1_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
        brca1_simulated_geno_plink <<- read.table(file = paste0(filename, "-geno.txt"), header = T)
        
        # Extract lengths
        mut = paste0("mutSim_BRCA1_",gen,"_",i); gene = "BRCA1_simulated"
        mut <<- mut; gene <<- gene
        readIndvAndFamData(gene, mut)
        dst <- firstBreakDist(fam=F)
        hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = 1, doHapPlot = F)
        brca.mut.cM <<- brca_cM_middle
        haplotypes <- read.table(paste0(filename, "-founder.txt"), header = T)
        #haplotypes2 <<- as.integer(haplotypes[haplotypes$SNP == mut,-1])
        haplotypes2 <<- as.integer(haplotypes[i,-1])
        matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
        breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
        chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
        
        # Adjust lengths relative to mutation
        chr_pos <- chr_pos %>% mutate(positive_cM_adj = positive_cM - brca.mut.cM, negative_cM_adj = abs(negative_cM - brca.mut.cM))
        return(chr_pos)
    }
    
    write_input_data <- function(){
        # write input files
        data = list()
        i = 1
        for (gen in c(110,510,1010,1510,2010)){
            for (sim in 1:10){
                chr_pos = get_sim_data(gen = gen, i = sim)
                arm.lengths = c(chr_pos$negative_cM_adj, chr_pos$positive_cM_adj)
                data[[i]] = data.frame(Age = gen, sim_num = sim, side = c(rep("left", length(arm.lengths)/2),rep("right", length(arm.lengths)/2)), arm.lengths = arm.lengths)
                i = i + 1
            }
        }
        df = do.call(rbind,data)
        dim(df)
        write.table(df, file = "Jacob Hjelmborg/input_110-2010-step500_10sims.txt", row.names = F, quote = F, sep = "\t")
    }
    
    for (i in 1:10){
        chr_pos = get_sim_data(gen = 2010, i = i)
        arm.lengths = c(chr_pos$positive_cM_adj, chr_pos$negative_cM_adj)
        
        # method 3
        print(100/mean(arm.lengths))
        #method 5
        age5=method5_fast(chr_pos)
        # method 6
        print(mutationAge_method6(chr_pos = chr_pos, brca_cM_middle = brca.mut.cM))
        
        # out <- tweedie.profile(arm.lengths ~ 1, p.vec=seq(1.5, 2.5, by=0.2))
        # out$p.max
        # out$ci
        # 
        # summary(glm(arm.lengths ~ 1, family = tweedie(var.power = 2, link.power = 0)))
        
        ## Fitting Tweedie distributions
        ll <- function(par, y) {
            print(par)
            #out <- sum(log(dtweedie(y = y, mu = par[1], phi = exp(par[2]), power = par[3])))
            out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = 1/2, power = par[2])))
            #out <- sum(log(dtweedie.saddle(y = y, mu = 2/par[1], phi = 1/2, power = par[2])))
            #out <- sum(log(tweedie.dev(y = y, mu = 2/par[1], power = par[2])))
            #out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = exp(par[2]), power = par[3])))
            return(out)
        }
        ll.gamma <- function(par, y) {
            print(par)
            out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = exp(par[2]), power = 2)))
            #out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = 1/2, power = 2)))
            return(out)
        }
        
        ll(par = c(2/mean(arm.lengths), 1.5), y = arm.lengths)
        # tt.tweedie = optim(par = c(2/mean(arm.lengths), 1/2, 1.5), fn = ll, y = arm.lengths, control = list(fnscale = -1), hessian = TRUE,
        #                    lower = c(0.01,0.01,1.02), upper = c(1000,100,2), method = "L-BFGS-B" )
        tt.tweedie = optim(par = c(2/mean(arm.lengths), 1.5), fn = ll, y = arm.lengths, control = list(fnscale = -1), hessian = TRUE,
                           lower = c(0.01,1.01), upper = c(100,1.99), method = "L-BFGS-B" )
        tt.gamma = optim(par = c(2/mean(arm.lengths),1/2), fn = ll.gamma, y = arm.lengths, control = list(fnscale = -1), hessian = TRUE,
                         lower = c(0.01,0.01), upper = c(1000, 10), method = "L-BFGS-B" ) ##method = "Brent", 
        
        std.error.tw = sqrt(diag(solve(-tt.tweedie$hessian)))
        std.error.gm = sqrt(diag(solve(-tt.gamma$hessian)))
        
        ll(par = tt.tweedie$par, y = arm.lengths) ### Tweedie model
        ll.gamma(par = tt.gamma$par, y = arm.lengths) ### Gamma model
        -2*(ll(par = tt.tweedie$par, y = arm.lengths) - ll.gamma(par = tt.gamma$par, y = arm.lengths))
        
        data.frame("Ic_Min"  = tt.tweedie$par - qnorm(0.975) *std.error.tw, "Estimates" = tt.tweedie$par, "Ic_Max" = tt.tweedie$par + qnorm(0.975) *std.error.tw)
        data.frame("Ic_Min"  = tt.gamma$par - qnorm(0.975) *std.error.gm, "Estimates" = tt.gamma$par, "Ic_Max" = tt.gamma$par + qnorm(0.975) *std.error.gm)
        
        #print(paste("Gamma:", 100*length(arm.lengths)/sum(arm.lengths))) ## Gamma age
        #print(paste("Gamma:", 100/(tt.gamma$par[1]))) ## Gamma age
        print(paste("Gamma:", 100/(2/tt.gamma$par[1]))) ## Gamma age
        print(paste("Tweedie:", 100/(2/tt.tweedie$par[1]))) ## Tweedie age
        #print(paste("Tweedie:", 100/(tt.tweedie$par[1]))) ## Tweedie age
        
        #print(paste("Tweedie-method5-mean:", mean(c(100/(2/tt.tweedie$par[1]), age5))))
        print(paste("Tweedie-method5-mean:", mean(c(100/(tt.tweedie$par[1]), age5))))
        
        df = data.frame("Ic_Min" = 100/(2/(tt.tweedie$par[1] - qnorm(0.975) *std.error.tw[1])),
                   "Age_estimate" = 100/(2/tt.tweedie$par[1]),
                   "Ic_Max" = 100/(2/(tt.tweedie$par[1] + qnorm(0.975) *std.error.tw[1])))
        print(df)
    }
    
    ## Example
    # Generate some fictitious data
    test.data <- rgamma(n=200, scale=1, shape=1)
    # The gamma is a Tweedie distribution with power=2;
    # let's see if p=2 is suggested by tweedie.profile:
    out <- tweedie.profile( test.data ~ 1, p.vec=seq(1.5, 2.5, by=0.2) )
    out$p.max
    out$ci
}

mutationAge_method9_rawGamma <- function(chr_pos, brca.mut.cM){
    chr_pos <- chr_pos %>% mutate(positive_cM_adj = positive_cM - brca.mut.cM, negative_cM_adj = abs(negative_cM - brca.mut.cM))
    n = length(chr_pos$positive_cM_adj)
    sum.lengths = sum(c(chr_pos$positive_cM_adj/100, chr_pos$negative_cM_adj/100))
    i.tau.hat <- 2*n/sum.lengths
    cc = 0.95
    g_l <- qgamma(shape=2*n,scale=1/(2*n),((1-cc)/2))
    g_u <- qgamma(shape=2*n,scale=1/(2*n),(cc+(1-cc)/2))		
    i.l = g_l*i.tau.hat
    i.u = g_u*i.tau.hat
    return(c(round(i.tau.hat, 1), round(i.l, 1), round(i.u, 1)))
}

computeMutationAge <- function(gene, mut, isFam=F, filename1 = "", filename2 = "", filename3 = "", 
                               filename4_1 = "", filename4_2 = "", filename5 = "", filename6 = "",
                               filename6_cor = "", filename8 = "", cluster = T, method = 1, cutoff = 0.5,
                               ancestral_method = "branchBoundIndep", minSamples = 5, progress = 0){
    #filename1 <<- filename1

    # if (gene == "BRCA1"){
    #     #morgan.brca1 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_Phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    #     brca.mut.cM <<- morgan.brca1[which.min(abs(morgan.brca1$position-brca_middle)),3]
    # } else {
    #     #morgan.brca2 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_Phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
    #     brca.mut.cM <<- morgan.brca2[which.min(abs(morgan.brca2$position-brca_middle)),3]
    # }

    brca.mut.cM = brca_cM_middle
    print(brca.mut.cM)

    num_groups = max(brca_data$cluster_groups)

    age_table1 = NULL
    age_table2 = NULL
    age_table3 = NULL
    age_table4_1 = NULL
    age_table4_2 = NULL
    age_table5 = NULL
    age_table6 = NULL
    age_table6_cor = NULL
    age_table8 = NULL
    if (cluster){
        for (consensus_group in 1:num_groups){
            #print(paste("consensus group:", consensus_group))
            #consensus_group = 1

            if (isFam){
                matched_plink = matched_plink_fam[,-1]
            } else {
                matched_plink = matched_plink_single
            }
            
            if (ancestral_method == "simulatedFounder"){
                print("FOUNDER")
                #haplotypes2 <<- rep(0, ncol(matched_plink)-1)
                haplotypes = read.table(paste0(filename_sim, "-founder.txt"), header = T)
                haplotypes2 <- as.integer(haplotypes[haplotypes$SNP == mut,-1])
            } else if (ancestral_method == "mostFreqBase"){
                haplotypes2 <<- findConsensus_plink(filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == consensus_group)$Onc_ID))
            } else if (ancestral_method == "branchBoundIndep"){
                matched_plink_subset = matched_plink[matched_plink$SNP %in% subset(brca_data, cluster_groups == consensus_group)$Onc_ID,]
                haplotypes2 <<- branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = cutoff, 
                                                           indep_sides = TRUE, min_samples = minSamples)[[1]]
            } else if (ancestral_method == "branchBound"){
                matched_plink_subset = matched_plink[matched_plink$SNP %in% subset(brca_data, cluster_groups == consensus_group)$Onc_ID,]
                haplotypes2 <<- branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = cutoff, 
                                                           indep_sides = FALSE, min_samples = minSamples)[[1]]
            } else if (ancestral_method == "branchBoundIndep_alleleFreq"){
                matched_plink_subset = matched_plink[matched_plink$SNP %in% subset(brca_data, cluster_groups == consensus_group)$Onc_ID,]
                haplotypes2 <<- branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = cutoff, 
                                                         indep_sides = T, min_samples = minSamples, alleleFreqs = alleleFreqs)[[1]]
            } else if (ancestral_method == "branchBound_alleleFreq"){
                matched_plink_subset = matched_plink[matched_plink$SNP %in% subset(brca_data, cluster_groups == consensus_group)$Onc_ID,]
                haplotypes2 <<- branch_and_bound_ancestral(matched_plink, matched_plink_subset, cutoff = cutoff, 
                                                         indep_sides = F, min_samples = minSamples, alleleFreqs = alleleFreqs)[[1]]
            } else {
                stop("Ancestral method does not exist")
            }

            age_temp1 = NULL
            age_temp2 = NULL
            age_temp3 = NULL
            age_temp4_1 = NULL
            age_temp4_2 = NULL
            age_temp5 = NULL
            age_temp6 = NULL
            age_temp6_cor = NULL
            age_temp8 = NULL
            for (age_group in 1:num_groups){
                #print(paste("age group:", age_group))
                # age_group = 1
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == age_group)$Onc_ID)
                breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
                chr_pos <<- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
                #chr_pos <- filter(chr_pos, sample_id %in% subset(brca_data, cluster_groups == age_group)$Onc_ID)
                #dim(chr_pos)
                
                if (is.null(chr_pos)){
                    age_temp1 = c(age_temp1, 0)
                    age_temp2 = c(age_temp2, 0)
                    age_temp3 = c(age_temp3, 0)
                    age_temp4_1 = c(age_temp4_1, 0)
                    age_temp4_2 = c(age_temp4_2, 0)
                    age_temp5 = c(age_temp5, 0)
                    age_temp6 = c(age_temp6, 0)
                    age_temp6_cor = c(age_temp6_cor, 0)
                    age_temp8 = c(age_temp8, 0)
                } else {

                    # Chance sharing correction
                    if (grepl("BRCA1", gene) && exists("popFreqs_brca1")){
                        chr_pos <<- adjustHapLengths(chr_pos, chr_coords, popFreqs = popFreqs_brca1, matched_plink2, epsilon = 0.0001)
                    } else if (grepl("BRCA2", gene) && exists("popFreqs_brca2")){
                        chr_pos <<- adjustHapLengths(chr_pos, chr_coords, popFreqs = popFreqs_brca2, matched_plink2, epsilon = 0.0001)
                    }
                    
                    # Method 1
                    if (nchar(filename1) > 0 || method == 1){
                        #age1 <- max(mutationAge_geneticDist(chr_pos, brca.mut.cM), -1)
                        #age_temp1 = c(age_temp1, round(age1, digits = 2))
                        age1 <- round(mutationAge_method1_median_with_bootstrap(chr_pos), digits = 1)
                        age_temp1 = c(age_temp1, paste0(age1[1], " (CI ", max(age1[2],0), ",", age1[3], ")"))
                        #print(age1)
                    }
    
                    # Method 2
                    if (nchar(filename2) > 0 || method == 2){
                        age2 <- max(mutationAge_physicalDist(chr_pos), -1)
                        age_temp2 = c(age_temp2, round(age2, digits = 1))
                        #print(age2)
                    }
    
                    # Method 3
                    if (nchar(filename3) > 0 || method == 3){
                        #age3 <- max(mutationAge_method3(chr_pos), -1)
                        #age_temp3 = c(age_temp3, round(age3, digits = 2))
                        age3 <- round(mutationAge_method3_with_bootstrap(chr_pos), digits = 1)
                        age_temp3 = c(age_temp3, paste0(age3[1], " (CI ", max(age3[2],0), ",", age3[3], ")"))
                        #print(age3)
                    }
    
                    # Method 4_1
                    if (nchar(filename4_1) > 0 || method == "4_1"){
                        age4_1 <- max(mutationAge_method4upper(chr_pos, brca.mut.cM), -1)
                        age_temp4_1 = c(age_temp4_1, round(age4_1, digits = 1))
                        #print(age4_1)
                    }
    
                    # Method 4_2
                    if (nchar(filename4_2) > 0 || method == "4_2"){
                        age4_2 <- max(mutationAge_method4under(chr_pos, brca.mut.cM), -1)
                        age_temp4_2 = c(age_temp4_2, round(age4_2, digits = 1))
                        #print(age4_2)
                    }
                    
                    # Method 5
                    if (nchar(filename5) > 0 || method == 5){
                        # age5 <- mutationAge_method5(chr_pos, brca.mut.cM, start = 60, stop = 5000, step = 100)
                        # age5 <- mutationAge_method5(chr_pos, brca.mut.cM, start = max(age5-50,1), stop = age5+50, step = 5)
                        # age5 <- mutationAge_method5(chr_pos, brca.mut.cM, start = max(age5-5,1), stop = age5+5, step = 1)
                        age5 <- round(method5_fast_with_bootstrap(chr_pos), digits = 1)
                        age_temp5 = c(age_temp5, paste0(age5[1], " (CI ", max(age5[2],0), ",", age5[3], ")"))
                    }
                    
                    # Method 6
                    if (nchar(filename6) > 0 || method == 6){
                        age6 <- mutationAge_method6(chr_pos, brca.mut.cM)
                        # confidence intervals
                        age_temp6 = c(age_temp6, paste0(age6[1], " (CI ", age6[2], ",", age6[3], ")"))
                        age_temp6_cor = c(age_temp6_cor, paste0(age6[4], " (CI ", age6[5], ",", age6[6], ")"))
                    }
                    
                    # Method 8
                    if (nchar(filename8) > 0 || method == 8){
                        age8 = mutationAge_method8_tweedie(chr_pos, brca.mut.cM)
                        age_temp8 = c(age_temp8, paste0(age8[1], " (CI ", age8[2], ",", age8[3], ")"))
                    }
                }
                if (progress>0){
                    incProgress(amount = 1/progress)
                }
            }

            if (nchar(filename1) > 0 || method == 1) age_table1 = rbind(age_table1, c(paste("Group", consensus_group), age_temp1))
            if (nchar(filename2) > 0 || method == 2) age_table2 = rbind(age_table2, c(paste("Group", consensus_group), age_temp2))
            if (nchar(filename3) > 0 || method == 3) age_table3 = rbind(age_table3, c(paste("Group", consensus_group), age_temp3))
            if (nchar(filename4_1) > 0 || method == "4_1") age_table4_1 = rbind(age_table4_1, c(paste("Group", consensus_group), age_temp4_1))
            if (nchar(filename4_2) > 0 || method == "4_2") age_table4_2 = rbind(age_table4_2, c(paste("Group", consensus_group), age_temp4_2))
            if (nchar(filename5) > 0 || method == 5) age_table5 = rbind(age_table5, c(paste("Group", consensus_group), age_temp5))
            if (nchar(filename6) > 0 || method == 6){
                age_table6 = rbind(age_table6, c(paste("Group", consensus_group), age_temp6))
                age_table6_cor = rbind(age_table6_cor, c(paste("Group", consensus_group), age_temp6_cor))
            }
            if (nchar(filename8) > 0 || method == 8) age_table8 = rbind(age_table8, c(paste("Group", consensus_group), age_temp8))
        }

        #print(age_table1)
        #print(age_table2)

    } else { # Country
        for (country in unique(brca_data$Country)){

            if (isFam){
                m <- filter(matched_plink_fam, SNP %in% subset(brca_data, Country == country)$Onc_ID)
            } else {
                m <- filter(matched_plink, SNP %in% subset(brca_data, Country == country)$Onc_ID)
            }

            print(dim(m))

            haplotypes <- findConsensus_plink(m)
            breaks <- findHaplotypeBreaks_plink(m, haplotypes)
            chr_pos <- mapNearestSNPs_plink(m, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)

            age_physicalDist <- mutationAge_physicalDist(chr_pos)
            age_geneticDist <- mutationAge_geneticDist(chr_pos, brca.mut.cM)


            age = rbind(age, c(group, age_physicalDist, age_geneticDist))
        }

        colnames(age) <- c("country", "age_physicalDist", "age_geneticDist")
    }

    if (nchar(filename1) > 0){
        colnames(age_table1) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename1, x = as.data.frame(age_table1), quote = F, row.names = F, sep = "\t")
    }
    if (nchar(filename2) > 0){
        colnames(age_table2) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename2, x = as.data.frame(age_table2), quote = F, row.names = F, sep = "\t")
    }
    if (nchar(filename3) > 0){
        colnames(age_table3) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename3, x = as.data.frame(age_table3), quote = F, row.names = F, sep = "\t")
    }
    if (nchar(filename4_1) > 0){
        colnames(age_table4_1) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename4_1, x = as.data.frame(age_table4_1), quote = F, row.names = F, sep = "\t")
    }
    if (nchar(filename4_2) > 0){
        colnames(age_table4_2) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename4_2, x = as.data.frame(age_table4_2), quote = F, row.names = F, sep = "\t")
    }
    if (nchar(filename5) > 0){
        colnames(age_table5) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename5, x = as.data.frame(age_table5), quote = F, row.names = F, sep = "\t")
    }
    if (nchar(filename6) > 0){
        colnames(age_table6) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename6, x = as.data.frame(age_table6), quote = F, row.names = F, sep = "\t")
    }
    if (nchar(filename6_cor) > 0){
        colnames(age_table6_cor) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename6_cor, x = as.data.frame(age_table6_cor), quote = F, row.names = F, sep = "\t")
    }
    if (nchar(filename8) > 0){
        colnames(age_table8) <- c("Ancestral haplotype", paste("Group", 1:num_groups, "age"))
        write.table(file = filename8, x = as.data.frame(age_table8), quote = F, row.names = F, sep = "\t")
    }
    #print(as.data.frame(age_table1))
    #print(as.data.frame(age_table2))

    age_table = setNames(as.data.frame(get(paste0("age_table", method))), c("Ancestral haplotype", paste("Group", 1:num_groups, "age")))
    return(age_table)
}
#age = computeMutationAge(gene, mut, fam=F, filename = "")

compareAgeMethod_generateFiles <- function(gen_start = 10, gen_stop = 510, gen_step = 20, sims = 10,
                                           samples = 100, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
                                           manual=FALSE, ancestral = "branchBoundIndep", minSampleBB=5,
                                           cs_correction="None", epsilon = 0.01, sim_type = "alleleFreqs"){
    # j = 1
    # ages1 <<- c()
    # ages1_mean <<- c()
    # ages2 <<- c()
    # ages3 <<- c()
    # ages3_mean <<- c()
    # ages5 <<- c()
    # ages5_new <<- c()
    # ages6 <<- c()
    # ages6_cor <<- c()
    # ages6_cs.cor <<- c()
    # ages6_cor_cs.cor <<- c()
    # ages7 <<- c()
    # ages7_group <<- c()
    
    if (manual){
        gen_start = 10
        gen_stop = 2510#510#10010
        gen_step = 50#10#20
        sims = 50#20#10
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
        sim_type = "knownBreaks"
        epsilon = 0.01
        minSampleBB=5
    }
    
    # hapmap2 <- read.table("/Users/lars/Downloads/annotHapMap2U.txt", header = T)
    # hapmap2.brca1 <- hapmap2 %>% filter(Chrom == "chr17")
    
    if (Sys.info()["sysname"] == "Linux"){
        popFreqs <<- read.table2("../pop_freqs_genotypes_brca1_ordered.txt", header = T)
        hapFreqs <<- read.table2("../../Scripts/simulate-population/famHaplotypes_simulationBasis_brca1_geno.txt", header = T)
    } else {
        popFreqs <<- read.table2("cache/pop_freqs_genotypes_brca1_ordered.txt", header = T)
        hapFreqs <<- read.table2("Scripts/simulate-population/famHaplotypes_simulationBasis_brca1_geno.txt", header = T)
        alleleFreqs <<- read.table2("cache/pop_freqs_minorAllele_brca1_ordered.txt", header = T)
        # ancestral_freqs <- c()
        # j = 1
        # for (i in 1:length(haplotypes2)){
        #     ancestral_freqs[j] = popFreqs[i, as.character(haplotypes2[i])]
        #     j = j+1
        # }
        # median_allele_freq = median(ancestral_freqs)
    }
        
    for (gen in seq(gen_start,gen_stop,gen_step)){
    #for (gen in seq(1010,1200,gen_step)){
        if (gen == gen_start){
            j = 1
            ages1 <<- c()
            ages1_mean <<- c()
            ages2 <<- c()
            ages3 <<- c()
            ages3_mean <<- c()
            ages5 <<- c()
            ages5_new <<- c()
            ages6 <<- c()
            ages6_cor <<- c()
            ages6_cs.cor <<- c()
            ages6_cor_cs.cor <<- c()
            ages7 <<- c()
        }
    
    # for (samples in seq(10,100,10)){
    #     gen = gen_start
        if (Sys.info()["sysname"] == "Linux"){
            pre_path = ""
            filename = paste0("simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen,"-seed_",seed)
        } else {
            pre_path = paste0("classify-simulated-population-age/", in_gene,"_simulated-",sim_type,"-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/")
            #filename = paste0("classify-simulated-population-age/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen,"-seed_",seed)
            #filename = paste0("classify-simulated-population-age/popFreqs_",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen,"-seed_",seed)
            filename = paste0(pre_path, "simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen,"-seed_",seed)
            #filename = paste0("classify-simulated-population-age/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",10,"_simulations-generations_",10,"_",10010,"_step",50,"-seed_",seed,"/simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",10,"_simulations-generations_",gen,"-seed_",seed)
            #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-20_samples-20_simulations-generations_10_100_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-20_samples-20_simulations-generations_",gen,"-seed_42")
            #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_",gen,"-seed_42")
            #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_10_500_step10-seed_24/simulated_population/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_",gen,"-seed_24")
            #filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-10_simulations-generations_10_10010_step50-seed_42/simulated_population/BRCA1_simulated-starGenealogy-100_samples-10_simulations-generations_",gen,"-seed_42")
        }
        brca1_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
        brca1_simulated_geno_plink <<- read.table2(file = paste0(filename, "-geno.txt"), header = T)
        
        for (i in 1:sims){
        #for (i in 1:1){
            mut = paste0("mutSim_BRCA1_",gen,"_",i); gene = "BRCA1_simulated"
            mut <<- mut; gene <<- gene
            readIndvAndFamData(gene, mut)
            dst <- firstBreakDist(fam=F)
            hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = 1, doHapPlot = F)
            brca.mut.cM <<- brca_cM_middle
            if (ancestral == "simulatedFounder"){
                #haplotypes2 <<- rep(0, ncol(matched_plink)-1)
                haplotypes <- read.table(paste0(filename, "-founder.txt"), header = T)
                #haplotypes2 <<- as.integer(haplotypes[haplotypes$SNP == mut,-1])
                haplotypes2 <<- as.integer(haplotypes[i,-1])
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
                breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
                chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
                #ancestral_method = "simulatedFounder"
            } else if (ancestral == "knownBreaks"){
                haplotypes2 <<- rep(0, ncol(matched_plink)-1)
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
                breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
                chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
                #ancestral_method = "knownBreaks"
            } else if (ancestral == "mostFreqBase"){
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
                haplotypes2 <<- findConsensus_plink(matched_plink2, cutoff = 0.5)
                breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
                chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
                #ancestral_method = "mostFreqBase"
            } else { # star_pop_freqs
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
                if (ancestral == "branchBoundIndep"){
                    res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = T, min_samples = minSampleBB)
                } else if (ancestral == "branchBound"){ #branchBound
                    res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = F, min_samples = minSampleBB)
                } else if (ancestral == "branchBoundIndep_alleleFreq"){
                    res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, 
                                                     indep_sides = T, min_samples = minSampleBB, alleleFreqs = alleleFreqs)
                } else if (ancestral == "branchBound_alleleFreq"){
                    res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, 
                                                     indep_sides = F, min_samples = minSampleBB, alleleFreqs = alleleFreqs)
                } else {
                    stop("Ancestral method does not exist")
                }
                haplotypes2 <<- res[[1]]
                breaks <<- res[[2]]
                chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
            }
            #ancestral_probs <<- sapply(1:length(haplotypes2), function(i) popFreqs[i, as.character(haplotypes2[i])])
            if (cs_correction=="gandolfo"){
                chr_pos <- adjustHapLengthsGandolfo(chr_pos, chr_coords, popFreqs, haplotypes2, epsilon = epsilon)
            } else if (cs_correction=="gandolfoAverage"){
                chr_pos <- adjustHapLengthsGandolfoAverage(chr_pos, chr_coords, popFreqs, haplotypes2, epsilon = epsilon)
            } else if (cs_correction=="adjustLengths"){
                chr_pos <- adjustHapLengths(chr_pos, chr_coords, popFreqs, matched_plink2, epsilon = epsilon)
            } else if (cs_correction=="adjustHaploFreqs"){
                chr_pos <- adjustHaploFreqs(chr_pos, matched_plink, chr_coords, epsilon)
            }
            
            ages1[j] <- mutationAge_geneticDist(chr_pos, brca.mut.cM)
            # ages1_mean[j] <- mutationAge_method1_mean(chr_pos, brca.mut.cM)

            # ages2[j] <- mutationAge_physicalDist(chr_pos)
            
            # ages3[j] <- mutationAge_method3_median(chr_pos)
            ages3_mean[j] <- mutationAge_method3(chr_pos)

            ages5[j] <- method5_fast(chr_pos)
                        
            # age5 = mutationAge_method5(chr_pos, brca.mut.cM, start=max(round(gen-gen*3/4),1), stop=round(gen+gen*3/4), step=300)
            # age5 = mutationAge_method5(chr_pos, brca.mut.cM, start=max(age5-300,1), stop=age5+300, step=100)
            # age5 = mutationAge_method5(chr_pos, brca.mut.cM, start=max(age5-60,1), stop=age5+60, step=10)
            # ages5[j] <- mutationAge_method5(chr_pos, brca.mut.cM, start=max(age5-5,1), stop=age5+5, step=1)
            
            # age5 = mutationAge_newMethod5(chr_pos, brca.mut.cM, start=max(round(gen-gen*3/4),1), stop=round(gen+gen*3/4), step=300)
            # age5 = mutationAge_newMethod5(chr_pos, brca.mut.cM, start=max(age5-300,1), stop=age5+300, step=100)
            # age5 = mutationAge_newMethod5(chr_pos, brca.mut.cM, start=max(age5-60,1), stop=age5+60, step=10)
            # ages5_new[j] <- mutationAge_newMethod5(chr_pos, brca.mut.cM, start=max(age5-5,1), stop=age5+5, step=1)

            age6 = mutationAge_method6(chr_pos, brca.mut.cM, chr_coords = chr_coords, popFreqs = popFreqs, haplotypes2 = haplotypes2)
            ages6[j] <- age6[1]
            ages6_cor[j] <- age6[4]
            
            age6 = mutationAge_method6(chr_pos, brca.mut.cM, chance.sharing.correction = T, epsilon = epsilon,
                                       chr_coords = chr_coords, popFreqs = popFreqs, haplotypes2 = haplotypes2)
            ages6_cs.cor[j] <- age6[1]
            ages6_cor_cs.cor[j] <- age6[4]
            
            # age7 = mutationAge_method7(chr_pos, brca.mut.cM, start=max(round(gen-gen*3/4),1), stop=round(gen+gen*3/4), step=300)
            # age7 = mutationAge_method7(chr_pos, brca.mut.cM, start=max(age5-300,1), stop=age5+300, step=100)
            # age7 = mutationAge_method7(chr_pos, brca.mut.cM, start=max(age5-60,1), stop=age5+60, step=10)
            # ages7[j] <- mutationAge_method7(chr_pos, brca.mut.cM, start=max(age5-5,1), stop=age5+5, step=1)
            
            # age7_group = mutationAge_method7_group(chr_pos, brca.mut.cM, start=max(round(gen-gen*3/4),1), stop=round(gen+gen*3/4), step=500)
            # age7_group = mutationAge_method7_group(chr_pos, brca.mut.cM, start=max(age5-300,1), stop=age5+300, step=100)
            # age7_group = mutationAge_method7_group(chr_pos, brca.mut.cM, start=max(age5-60,1), stop=age5+60, step=10)
            # ages7_group[j] <- mutationAge_method7_group(chr_pos, brca.mut.cM, start=max(age5-5,1), stop=age5+5, step=1)
            
            j = j+1
            
            #print(sum(matched_plink[,-1])/2)
            
            #r_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
            # for(i in generations){computeGenerationLikelihood(chr_pos,i)}
            # chr_pos
        }
        #stop()
    }
    #print(ages1)
    #print(ages5)
    #print(mean(ages1))
    #print(mean(ages5))
    #print(mean(ages6))
    #boxplot(ages1, ages5)
    
    #samples = rep(seq(10,100,10), each = sims)
    data = data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method1", samples=samples, predicted_age=ages1)
    #data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method1_mean", samples=samples, predicted_age=ages2))
    #data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method2", samples=samples, predicted_age=ages2))
    #data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method3", samples=samples, predicted_age=ages3))
    data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method3_mean", samples=samples, predicted_age=ages3_mean))
    #data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method5", samples=samples, predicted_age=ages5))
    #data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method5_new", samples=samples, predicted_age=ages5_new))
    data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method6", samples=samples, predicted_age=ages6))
    data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method6_cor", samples=samples, predicted_age=ages6_cor))
    data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method6_cs.cor", samples=samples, predicted_age=ages6_cs.cor))
    data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method6_cor_cs.cor", samples=samples, predicted_age=ages6_cor_cs.cor))
    #data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method7", samples=samples, predicted_age=ages7))
    #data = rbind(data, data.frame(actual_age=rep(seq(gen_start,gen_stop,gen_step), each = sims), method="method7_group", samples=samples, predicted_age=ages7_group))
    
    # write.table(data, file = paste0("compare-age-methods/",in_gene,"_age-comparisons_method1_1mean_2_3_3mean_5_6_6cor_6cscor_7-start_",
    #                                 gen_start,"-stop_",gen_stop,"-samples_",samples,"-genealogy_",genealogy,".txt"), 
    #             quote = F, row.names = F, col.names = T)
    
    if (cs_correction=="gandolfo"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_gandolfo_", epsilon)
    } else if (cs_correction=="adjustLengths"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_adjustLengths_", epsilon)
    } else if (cs_correction=="gandolfoAverage"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_gandolfoAverage_", epsilon)
    } else {
        custom = paste0("-minSample_",minSampleBB)
    }
    
    #dir.create(path = paste0("compare-age-methods/", ancestral, "-max5or25percent-cs_0.01"), showWarnings = F)
    dir.create(path = paste0(pre_path, "compare-age-methods/", ancestral, custom), showWarnings = F)
    write.table(data, file = paste0(pre_path, "compare-age-methods/",ancestral, custom,"/",in_gene,"_age-adjustedLengths-comparisons_method1_3mean_5_6_6cor_6cscor",
                                    "-start_",gen_start,"-stop_",gen_stop,"-samples_",samples,"-genealogy_",genealogy,".txt"), 
                quote = F, row.names = F, col.names = T)
    
    #compareAgeMethods_plots(filepath = paste0(pre_path, "compare-age-methods/",ancestral, custom), gen_start, gen_stop, gen_step, ancestral, sims, samples)
    
    #write.table(data, file = paste0("compare-age-methods/age-comparisons_method1_5_6-10_3510_step50_seed42.txt"), quote = F, row.names = F, col.names = T)
    
    ### RUNNING THIS METHOD ON SERVER
    # source("Scripts/Haplotype-projekt-ver2.R")
    # source("Scripts/Clustering.R")
    # source("Scripts/Mutation-Age.R")
    # source("Scripts/Gandolfo_Speed_Mutation_Age_Estimation.R")
    # readData()
    # testSimAge()
    
}

compareSampleSize_AgeMethods <- function(gen = 50){
    gen = 50
    # filepath = paste0("classify-simulated-population-age/BRCA1_simulated-test_sample_size-knownBreaks-starGenealogy-50_simulations-generations_",gen,"_",gen,"_step",gen,"-seed_42/compare-age-methods/simulatedFounder-minSample_0.25-epsilon_0.01/")
    # filepath = paste0("classify-simulated-population-age/BRCA1_simulated-test_sample_size-famHaplotypes-starGenealogy-50_simulations-generations_50_50_step50-seed_42/compare-age-methods/branchBoundIndep-minSample_0.25-epsilon_0.01/")
    filepath = paste0("classify-simulated-population-age/BRCA1_simulated-test_sample_size-famHaplotypes-starGenealogy-50_simulations-generations_50_50_step50-seed_42/compare-age-methods/mostFreqBase-minSample_0.25-epsilon_0.01/")
    data <- do.call(rbind, lapply(list.files(path = filepath, full.names = T, pattern = ".txt"), read.table, header = T))
    data = data %>% arrange(samples, method) # Should be 7000x4 in dim
    # Setup 
    methods = c("method1","method5","method9")
    methods = c("method1","method3_mean","method5","method8","method9")
    # methods = c("method1","method3_mean","method5","method6","method7")
    #methods = c("method1","method3_mean","method5","method6","method6_cor")
    #methods = c("method1","method1_mean","method2","method3","method3_mean","method5","method6","method6_cor","method6_cs.cor", "method6_cor_cs.cor","method7")
    data = filter(data, method %in% methods)
    
    data$method = as.character(data$method)
    data$method[data$method=="method1"] <- "method 1"
    data$method[data$method=="method3_mean"] <- "method 2"
    data$method[data$method=="method9"] <- "method 2-4"
    data$method[data$method=="method9"] <- "method 3"
    data$method[data$method=="method8"] <- "method 4"
    data$method[data$method=="method5"] <- "method 5"
    
    methods = c("method 1","method 2-4","method 5")
    methods = c("method 1","method 2","method 3","method 4","method 5")
    
    data = filter(data, samples %in% c(seq(10,100,10), seq(200,1000,100)))
    
    # For manuscript
    # data = data %>% arrange(method)
    # methods = c("method1","method2","method3","method4")
    # data$method = rep(methods, each=nrow(data)/length(methods))
    
    cols = colors_plot(color_names = methods, palette = "ggplot")
    # cols = colors_plot(color_names = methods, palette = "Set1")
    # cols = colors_plot(color_names = methods, palette = "sasha2")
    width = 21; height = 12
    
    gen_start = gen
    gen_stop = gen
    sims = 50
    #samples = 100
    
    for (m in methods){
        m1 = filter(data, method==m)
        p = m1 %>% ggplot(aes(x = factor(samples), y = predicted_age, fill = method)) +
            geom_boxplot() +
            scale_fill_manual(values = cols[m]) +
            #geom_point() + 
            geom_abline(intercept = gen_start, slope = 0, col = "darkblue") +
            scale_y_continuous(breaks = seq(gen_start-10*gen_start,gen_stop+10*gen_stop,10)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            xlab("sample size") + ylab("predicted age")
        #geom_smooth(method = "lm", se = F)
        print(p)
        ggsave(p, filename = paste0(filepath, "/boxplot_", m, "-generation_",gen_start, "-simulations_",sims, "-samplesSize-start_",min(data$samples),"-stop_",max(data$samples),".png"), width = width, height = height)
    }
    
    # combined
    p = ggplot(data, aes(x=factor(samples), y = predicted_age, fill = method)) + 
        geom_boxplot() +
        scale_fill_manual(values = cols) +
        geom_abline(intercept = gen_start, slope = 0, col = "black") + #, linetype=2) +
        scale_y_continuous(breaks = seq(gen_start-10*gen_start,gen_stop+10*gen_stop,10)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme_bw() +
        theme(text = element_text(size=14)) +
        theme(legend.title = element_text(face = "bold")) +
        xlab("Sample size") + ylab("Estimated age") + labs(fill="Method")
    print(p)
    # ggsave(p, filename = paste0(filepath, "/boxplot_all-methods-generation_",gen_start, "-simulations_",sims, "-samplesSize-start_",min(data$samples),"-stop_",max(data$samples),".png"), 
    #        width = width, height = height, scale = 0.7)
    ggsave(p, filename = paste0(filepath, "/boxplot_all-methods-generation_",gen_start, "-simulations_",sims, "-samplesSize-start_",min(data$samples),"-stop_",max(data$samples),".pdf"), 
           width = width, height = height, scale = 0.7)
    
    # Confidence intervals
    data$conf_interval = data$CI_max - data$CI_min
    p = data %>% #filter(method != "method 1") %>%
        ggplot(aes(x=factor(samples), y = conf_interval, fill = method)) + 
        geom_boxplot() +
        scale_fill_manual(values = cols) +
        # geom_abline(intercept = gen_start, slope = 0, col = "black") +
        scale_y_continuous(breaks = seq(0,gen_stop+20*gen_stop,20)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme_bw() +
        theme(text = element_text(size=14)) +
        theme(legend.title = element_text(face = "bold")) +
        xlab("Sample size") + ylab("Interval lengths") + labs(fill="Method") +
        ylim(c(0,125))
    print(p)
    ggsave(p, filename = paste0(filepath, "/confidence_intervals-boxplot_all-methods-generation_",gen_start, "-simulations_",sims, "-samplesSize-start_",min(data$samples),"-stop_",max(data$samples),".pdf"), 
           width = width, height = height, scale = 0.7)
}

compareAgeMethods_plots <- function(gen_start, gen_stop, gen_step, ancestral_method, sims, samples, sim_type, custom_string,
                                    pre_path=NULL, seed = 42, genealogy="starGenealogy", manual=F){
    
    breaks_step = gen_step
    
    if (manual){
        gen_start = 10
        gen_stop = 510 #1210 #510 #10010 #3010
        gen_step = 10 #10 #50 #20
        samples = 100
        seed = 42
        sims = 50 #10 #20
        breaks_step = 10 #20 #500 # 100 # gen_step
        genealogy = "starGenealogy"
        
            
        # Read interval data
        #filepath = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/compare-age-methods/")
        #filepath = paste0("classify-simulated-population-age/popFreqs_BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/compare-age-methods/")
        # filepath = paste0("classify-simulated-population-age/popFreqs_BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/compare-age-methods/mostFreqBase")
        # filepath = paste0("classify-simulated-population-age/popFreqs_BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/compare-age-methods/BranchBound")
        # filepath = paste0("classify-simulated-population-age/popFreqs_BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/compare-age-methods/BranchBoundIndep")
        # filepath = paste0("classify-simulated-population-age/popFreqs_BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/compare-age-methods/simulatedFounder")
        ancestral_method = "mostFreqBase"
        ancestral_method = "branchBound"
        ancestral_method = "branchBoundIndep"
        ancestral_method = "simulatedFounder"
        custom_string = "-minSample_2-epsilon_0.01"
        custom_string = "-max5"
        custom_string = "-minSample_5-cs_correction_gandolfo_0.001"
        custom_string = "-minSample_5-cs_correction_gandolfo_0.01"
        custom_string = "-minSample_5-cs_correction_adjustLengths_0.01"
        custom_string = "-minSample_5-cs_correction_adjustLengths_0.001"
        custom_string = "-minSample_5-epsilon_0.01"
        custom_string = "-minSample_0.25-epsilon_0.01"
        custom_string = "-minSample_auto-epsilon_0.01"
        sim_type = "alleleFreqs"
        sim_type = "popFreqs"
        sim_type = "knownBreaks"
        sim_type = "haplotypes"
        sim_type = "famHaplotypes"
        sim_type = "famHaplotypesNoHet"
    }
    if (is.null(pre_path)){
        pre_path = paste0("classify-simulated-population-age/BRCA1_simulated-",sim_type,"-",genealogy,"-",samples,"_samples-",sims,
                          "_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed)
    }
    filepath = paste0(pre_path, "/compare-age-methods/", ancestral_method, custom_string)
    data <- do.call(rbind, lapply(list.files(path = filepath, full.names = T, pattern = c("method1_3mean_5_6_6cor_6cscor_8_9", ".txt")), read.table, header = T))
    data = data %>% arrange(samples, method)
    dim(data)
    
    # Read saved data
    #data <- read.table("compare-age-methods/age-comparisons_method1-5-6_10-500-step10_seed42.txt", header = T)
    #data <- read.table("compare-age-methods/age-comparisons_method1-5-6_10-500-step10_seed24.txt", header = T)
    #data <- read.table("compare-age-methods/age-comparisons_method1-5-6_10-3010-step50_seed42.txt", header = T)
    #data <- read.table("compare-age-methods/age-comparisons_method1_5_6-10_3510_step50_seed42.txt", header = T)
    
    # Setup 
    methods = c("method1","method3_mean","method5","method6","method6_cor","method6_cs.cor","method6_cor_cs.cor", "method8", "method9")
    #methods = c("method1","method1_mean","method2","method3","method3_mean","method5","method6","method6_cor","method6_cs.cor","method6_cor_cs.cor","method7")
    #methods = c("method1","method1_mean","method2","method3","method3_mean","method5","method6","method6_cor","method7")
    #methods = c("method1","method3_mean","method5","method6","method6_cs.cor")
    # if (grepl("cs_correction", custom_string)){
    #     #methods = c("method3_mean","method5","method6","method6_cor")
    #     methods = c("method1","method3_mean","method5","method6")#,"method8")
    # } else {
    #     #methods = c("method1","method3_mean","method5","method6","method6_cor","method6_cs.cor","method6_cor_cs.cor", "method8", "method9")
    #     methods = c("method1","method3_mean","method5","method6","method6_cor","method6_cs.cor","method6_cor_cs.cor", "method8", "method9")
    # }
    #methods = c("method6","method6_cor","method6_cs.cor","method6_cor_cs.cor")
    #methods = c("method1","method3_mean","method5","method5_new","method6")
    #methods = c("method1", "method3_mean","method5","method5_new","method6")
    data = filter(data, method %in% methods)
    
    
    if (manual){
        # methods = c("method1","method3_mean","method5", "method8", "method9")
        methods = c("method1","method5","method9")
        data = filter(data, method %in% methods)
        data$method = as.character(data$method)
        data$method[data$method=="method1"] <- "method 1"
        data$method[data$method=="method3_mean"] <- "method 2"
        data$method[data$method=="method9"] <- "method 2-4"
        # data$method[data$method=="method9"] <- "method 3"
        data$method[data$method=="method8"] <- "method 4"
        data$method[data$method=="method5"] <- "method 5"
        # data$method <- replace(data$method, data$method %in% c("method3_mean", "method9", "method8"), values = c("method2", "method3", "method4"))
        
        # gen_step = 50; gen_stop = 1010
        gen_step = 50; gen_stop = 510
        # gen_step = 20; gen_stop = 210
        breaks_step = gen_step
        data = filter(data, actual_age %in% seq(10,gen_stop,gen_step))
        # data = filter(data, actual_age %in% c(seq(10,100,10),seq(110,gen_stop,gen_step)))
        
        # methods = c("method 1","method 2","method 3", "method 4", "method 5")
        methods = c("method 1","method 2-4", "method 5")
    }
    
    # For munascript
    #gen_step = 40
    #data = filter(data, method %in% methods, actual_age %in% seq(10,3010,gen_step))
    #data = filter(data, method %in% methods, actual_age <= 6560)
    # data = data %>% arrange(method)
    # methods = c("method1","method2","method3","method4")
    # data$method = rep(methods, each=nrow(data)/length(methods))

    cols = colors_plot(color_names = methods, palette = "ggplot")
    #width = 21; height = 12
    width = 16; height = 9
    #width = 30; height = 18
    
    # data = filter(data, predicted_age < 200)
    
    # stats
    stats <- data %>% group_by(actual_age, method) %>% summarise(mean_ages = mean(predicted_age))
    print(stats)
    
    #gen_start = 120
    #gen_stop = 6560
    # Scatter plot
    p = stats %>% ggplot(aes(x = factor(actual_age), y = mean_ages, col = method)) +
        geom_point() +
        scale_color_manual(values = cols) +
        # geom_abline(intercept = gen_start-gen_step, slope = gen_step, col = "darkblue") +
        geom_abline(intercept = gen_start-gen_step, slope = gen_step, col = "darkblue") +
        #scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,gen_step)) +
        scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
        scale_x_discrete(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        xlab("actual age") + ylab("mean age")
    print(p)
    #ggsave(p, filename = paste0("compare-age-methods/",ancestral_method,"-mean_age_all_methods-start_",gen_start, "-stop_",gen_stop,"-step_",gen_step,"-simulations_",sims, "-samples_",samples,".png"), width = width, height = height)
    # ggsave(p, filename = paste0(filepath, "/",ancestral_method,"-mean_age_all_methods-start_",gen_start, "-stop_",gen_stop,"-step_",gen_step,"-simulations_",sims, "-samples_",samples,".png"), width = width, height = height)
    ggsave(p, filename = paste0(filepath, "/",ancestral_method,"-mean_age_all_methods-start_",gen_start, "-stop_",gen_stop,"-step_",gen_step,"-simulations_",sims, "-samples_",samples,".pdf"), width = width, height = height)
    
    # Stats - median
    stats_median <- data %>% group_by(actual_age, method) %>% summarise(mean_ages = median(predicted_age))
    print(stats_median)
    
    p = stats_median %>% ggplot(aes(x = factor(actual_age), y = mean_ages, col = method)) +
        geom_point() +
        scale_color_manual(values = cols) +
        geom_abline(intercept = gen_start-gen_step, slope = gen_step, col = "darkblue") +
        #scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,gen_step)) +
        scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
        scale_x_discrete(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        xlab("actual age") + ylab("mean age")
    print(p)
    
    # Box plot
    p = ggplot(data, aes(x=factor(actual_age), y = predicted_age, fill = method)) +
    # gen_step=50
    # data2 = data %>% group_by(actual_age, method) %>% summarise_all(mean)
    # p = data %>% filter(actual_age %in% c(10,seq(50,1000,50))) %>%
        # ggplot(aes(x=factor(actual_age), y = predicted_age, fill = method)) +
        geom_boxplot() +
        scale_fill_manual(values = cols) +
        geom_abline(intercept = gen_start-gen_step, slope = gen_step, col = "darkblue") +
        # geom_abline(intercept = gen_start, slope = 0, col = "darkblue") +
        # geom_errorbar(aes(ymax=data2$CI_max, ymin=data2$CI_min, x=data$actual_age)) +
        # geom_point(aes(x=data$method, y=data$predicted_age)) +
        #scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,gen_step)) +
        scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
        scale_x_discrete(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(text = element_text(size=14)) +
        theme(legend.title = element_text(face = "bold")) +
        xlab("Actual age") + ylab("Estimated age") + labs(fill="Method")
    print(p)
    # ggsave(p, filename = paste0(filepath, "/",ancestral_method,"-boxplot_all-methods-start_",gen_start, "-stop_",gen_stop,"-step_",gen_step,"-simulations_",sims, "-samples_",samples,".png"), width = width, height = height)
    ggsave(p, filename = paste0(filepath, "/",ancestral_method,"-boxplot_all-methods-start_",gen_start, "-stop_",gen_stop,"-step_",gen_step,"-simulations_",sims, "-samples_",samples,".pdf"), 
           width = width, height = height, scale=0.8)
    
    
    # Confidence interval comparisons
    data$conf_interval = data$CI_max-data$CI_min
    p = ggplot(data, aes(x=factor(actual_age), y = conf_interval, fill = method)) +
        # gen_step=50
        # data2 = data %>% group_by(actual_age, method) %>% summarise_all(mean)
        # p = data %>% filter(actual_age %in% c(10,seq(50,1000,50))) %>%
        # ggplot(aes(x=factor(actual_age), y = predicted_age, fill = method)) +
        geom_boxplot() +
        scale_fill_manual(values = cols) +
        # geom_abline(intercept = gen_start-gen_step, slope = gen_step, col = "darkblue") +
        # geom_abline(intercept = gen_start, slope = 0, col = "darkblue") +
        # geom_errorbar(aes(ymax=data2$CI_max, ymin=data2$CI_min, x=data$actual_age)) +
        # geom_point(aes(x=data$method, y=data$predicted_age)) +
        #scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,gen_step)) +
        scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
        scale_x_discrete(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        theme(text = element_text(size=14)) +
        theme(legend.title = element_text(face = "bold")) +
        xlab("Actual age") + ylab("Interval length (in generations)") + labs(fill="Method")
    print(p)
    ggsave(p, filename = paste0(filepath, "/",ancestral_method,"-confidence_interval_lengths-boxplot-all_methods-start_",gen_start, "-stop_",gen_stop,"-step_",gen_step,"-simulations_",sims, "-samples_",samples,".pdf"), 
           width = width, height = height, scale=0.8)
    
    return()
    
    # Per sample scatter plot
    for (m in methods){
        m1 = filter(data, method==m)
        p = m1 %>% ggplot(aes(x = factor(actual_age), y = predicted_age, col = method)) +
            geom_point() + 
            scale_color_manual(values = cols[m]) +
            geom_abline(intercept = gen_start-gen_step, slope = gen_step, col = "darkblue") +
            #scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,gen_step)) +
            scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
            scale_x_discrete(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            xlab("actual age") + ylab("predicted age")
        print(p)
        ggsave(p, filename = paste0(filepath, "/",ancestral_method,"scatterplot_", m, "-start_",gen_start, "-stop_",gen_stop,"-step_",gen_step,"-simulations_",sims, "-samples_",samples,".png"), width = width, height = height)
    }
    
    # Per sample boxplot
    for (m in methods){
        m1 = filter(data, method==m)
        p = m1 %>% ggplot(aes(x = factor(actual_age), y = predicted_age, fill = method)) +
            geom_boxplot() +
            scale_fill_manual(values = cols[m]) +
            #geom_point() + 
            geom_abline(intercept = gen_start-gen_step, slope = gen_step, col = "darkblue") +
            #scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,gen_step)) +
            scale_y_continuous(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
            scale_x_discrete(breaks = seq(gen_start-20*gen_step,3*gen_stop,breaks_step)) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            xlab("actual age") + ylab("predicted age")
        #geom_smooth(method = "lm", se = F)
        print(p)
        ggsave(p, filename = paste0(filepath, "/",ancestral_method,"boxplot_", m, "-start_",gen_start, "-stop_",gen_stop,"-step_",gen_step,"-simulations_",sims, "-samples_",samples,".png"), width = width, height = height)
    }
}

comparing_founders <- function(gene = "BRCA1"){
    gene <<- gene
    if (gene == "BRCA1"){
        brca1_pheno_merged <- read.csv("input/139_ouh_june_2017/B1_Onco_phenotype_distribution_311215.csv", stringsAsFactors = F)
        brca1_pheno_merged$cluster_groups = 1
        brca_muts = count(brca1_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS
        bak_brca_pheno = brca1_pheno_merged
    } else {
        brca2_pheno_merged <- read.csv("input/139_ouh_june_2017/B2_Onco_phenotype_distribution_180816.csv", stringsAsFactors = F)
        brca2_pheno_merged$cluster_groups = 1
        brca_muts = count(brca2_pheno_merged, Mut1HGVS) %>% filter(n>100) %>% arrange(n) %>% .$Mut1HGVS
        bak_brca_pheno = brca2_pheno_merged   
    }
    
    dir.create(paste0("comparing_founders/", gene), showWarnings = F)
    removeSNPsInGene <<- TRUE
    isFam = TRUE
    fam = if (isFam) "fam" else "individual"
    k=2
    for (i in 1:(length(brca_muts)-1)){
        for (j in (i+1):length(brca_muts)){
            mut = as.character(brca_muts[i])
            mut2 = as.character(brca_muts[j])
            #mut = "c.427G>T"
            #mut2 = "c.3331_3334delCAAG"
            name_comb = paste0(mut, "_vs_", mut2)
            name_comb <- gsub(">", "_", name_comb)
            name_comb <<- gsub("\\?", "_", name_comb)
            
            dir.create(paste0("comparing_founders/", gene, "/", name_comb), showWarnings = F)
            
            if (gene == "BRCA1"){
                # Assume all carriers with same mut is from same country
                brca1_pheno_merged[brca1_pheno_merged$Mut1HGVS == mut, "Country"] <- mut
                brca1_pheno_merged[brca1_pheno_merged$Mut1HGVS == mut2, "Country"] <- mut2
                brca1_pheno_merged[brca1_pheno_merged$Mut1HGVS == mut2, "cluster_groups"] <- 2
                # Assume mut_i and mut_j has same mutation - to verify clustering
                brca1_pheno_merged[brca1_pheno_merged$Mut1HGVS == mut2, "Mut1HGVS"] <- mut
                # Make it global
                brca1_pheno_merged <<- brca1_pheno_merged
            } else {
                # Assume all carriers with same mut is from same country
                brca2_pheno_merged[brca2_pheno_merged$Mut1HGVS == mut, "Country"] <- mut
                brca2_pheno_merged[brca2_pheno_merged$Mut1HGVS == mut2, "Country"] <- mut2
                brca2_pheno_merged[brca2_pheno_merged$Mut1HGVS == mut2, "cluster_groups"] <- 2
                # Assume mut_i and mut_j has same mutation - to verify clustering
                brca2_pheno_merged[brca2_pheno_merged$Mut1HGVS == mut2, "Mut1HGVS"] <- mut
                # Make it global
                brca2_pheno_merged <<- brca2_pheno_merged
            }
            
            
            readIndvAndFamData(gene, mut)
            # dst <- firstBreakDist(fam=isFam)
            # #dst <- computeCustomDist(distMethod = firstUpperBreakDist, fam = T)
            # #dst <- computeCustomDist(distMethod = firstUnderBreakDist, fam = T)
            # hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = k, doHapPlot = F)
            # plot_dendrogram(hc, 2)
            # 
            # png(filename = paste0("comparing_founders/", gene, "/", name_comb, "/", name_comb, "-ward.D-",fam,"-clusters_",k,".png"), width = 1600, height = 900, res = 100)
            # plot_dendrogram(hc, 2)
            # dev.off()
            
            computeMutationAge(gene, mut, isFam=isFam, 
                               filename3 = paste0("comparing_founders/", gene, "/", name_comb, "/", name_comb, "-age-method_3-ward.D-", fam, "-clusters_",k,".txt"),
                               filename5 = paste0("comparing_founders/", gene, "/", name_comb, "/", name_comb, "-age-method_5-ward.D-", fam, "-clusters_",k,".txt"),
                               filename6 = paste0("comparing_founders/", gene, "/", name_comb, "/", name_comb, "-age-method_6-ward.D-", fam, "-clusters_",k,".txt"),
                               filename6_cor = paste0("comparing_founders/", gene, "/", name_comb, "/", name_comb, "-age-method_6cor-ward.D-", fam, "-clusters_",k,".txt"),
                               cluster = T, CI = F)
            
            
            if (gene == "BRCA1"){
                brca1_pheno_merged = bak_brca_pheno
            } else {
                brca2_pheno_merged = bak_brca_pheno
            }
            
        }
        
    }
}

# d <<- data.frame(A=c(1,2), B=c(3,4))
# for (i in 1:2){
#     d[,"A"] <- i
# }

comparing_founder_ages <- function(){
    gene = "BRCA1"
    method = "method6"
    # for (gene in c("BRCA1", "BRCA2")){
    #     for (method in c("method3","method5","method6")){
    for (gene in c("BRCA2")){
        for (method in c("method3","method6")){
            print(gene)
            print(method)
            #filepath = paste0("comparing_founders/BRCA1_method6_wardD_fam/")
            #filepath = paste0("comparing_founders/BRCA1_method5_wardD_fam/")
            #filepath = paste0("comparing_founders/BRCA1_method3_wardD_fam/")
            #filepath = paste0("comparing_founders/BRCA2_method3_wardD_fam/")
            #filepath = paste0("comparing_founders/BRCA2_method6_wardD_fam/")
            filepath = paste0("comparing_founders/",gene,"_",method,"_fam/")
            print(filepath)
            data <- do.call(rbind, lapply(list.files(path = filepath, full.names = T, pattern = "clusters_2.txt"), read.table, header = T, sep = "\t"))
            print(head(data))
            
            # matched assumed ancestral haplotype and group age
            inner_ages = c(filter(data, Consensus == "Group 1")$Group.1.age, filter(data, Consensus == "Group 2")$Group.2.age)
            # Group ages with other ancestral haplotype
            outer_ages = c(filter(data, Consensus == "Group 1")$Group.2.age, filter(data, Consensus == "Group 2")$Group.1.age)
            
            df <- data.frame(ages = c(inner_ages,outer_ages), category = c(rep("inner", length(inner_ages)), rep("outer", length(outer_ages))))
            #df <- data.frame(ages = inner_ages, category = rep("inner", length(inner_ages)))
            df$ages[df$ages > 5000] = 5000
            #df %>% arrange(desc(ages))
            
            p <- ggplot(df, aes(x=ages, fill = category)) + geom_density(alpha=0.5)
            print(p)
            ggsave(filename = paste0("comparing_founders/DensityPlot-", gene, "-", method, "-fam-preGrouped.png"), plot = p)
            
            # data2 <- do.call(rbind, lapply(list.files(path = filepath, full.names = T, pattern = "clusters_1.txt"), read.table, header = T, sep = "\t"))
            # print(head(data2))
            # combined_ages = data2$Group.1.age
            # 
            # df <- data.frame(ages = c(inner_ages,combined_ages), category = c(rep("inner", length(inner_ages)), rep("combined", length(combined_ages))))
            # df$ages[df$ages > 5000] = 5000
            # 
            # p <- ggplot(df, aes(x=ages, fill = category)) + geom_density(alpha=0.5)
            # print(p)
            # ggsave(filename = paste0("comparing_founders/DensityPlot-", gene, "-", method, "-wardD-fam-inner_vs_combined.png"), plot = p)
            
        }
    }
    
}

    
compareAgeMethod_generateFiles_parallel <- function(gen_start = 10, gen_stop = 510, gen_step = 20, sims = 10,
                                           samples = 100, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
                                           manual=FALSE, ancestral = "branchBoundIndep", minSampleBB=5,
                                           cs_correction="None", epsilon = 0.01, sim_type = "alleleFreqs",
                                           cores = NULL, pre_path=NULL, cmp_sample_size=F){
    
    custom = custom = paste0("-minSample_",minSampleBB,"-epsilon_", epsilon)
    if (cs_correction=="gandolfo"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_gandolfo_", epsilon)
    } else if (cs_correction=="gandolfoAverage"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_gandolfoAverage_", epsilon)
    } else if (cs_correction=="adjustLengths"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_adjustLengths_", epsilon)
    }
    if (is.null(pre_path)){
        pre_path = paste0("classify-simulated-population-age/", in_gene,"_simulated-",sim_type,"-",genealogy,"-",samples,"_samples-",sims,
                          "_simulations-generations_",gen_start,"_",gen_stop,"_step",gen_step,"-seed_",seed,"/")
    }
    filename = paste0(pre_path, "compare-age-methods/",ancestral, custom,"/",in_gene,"_age-comparisons_method1_3mean_5_6_6cor_6cscor_8_9",
                      "-start_",gen_start,"-stop_",gen_stop,"-samples_",samples,"-genealogy_",genealogy,".txt")
    
    # If ages already computed, only create plots
    if (file.exists(filename)){
        if (!cmp_sample_size){
            compareAgeMethods_plots(gen_start = gen_start, gen_stop = gen_stop, gen_step = gen_step, ancestral_method = ancestral, 
                                    sims = sims, samples = samples, sim_type = sim_type, custom_string = custom, seed = seed, 
                                    pre_path = pre_path)
        }
        return()
    }
    
    library(foreach)
    library(doMC)
    
    if (manual){
        gen_start = 100#10
        gen_stop = 101#2510#10010
        gen_step = 10#20#10#20
        sims = 10#20#10
        samples = 100
        in_gene = "BRCA1"
        genealogy="starGenealogy"
        genealogy="correlatedGenealogy"
        seed = 42
        ancestral = "branchBoundIndep"
        ancestral = "simulatedFounder"
        sim_type = "famHaplotypes"
        sim_type = "knownBreaks"
        #cs_correction="gandolfo"
        #cs_correction="adjustHaploFreqs"
        cs_correction="adjustLengths"
        cs_correction="none"
        epsilon = 0.01
        minSampleBB=5
    }
    
    if (Sys.info()["sysname"] == "Linux"){
        if (file.exists("cache/pop_freqs_genotypes_brca1_ordered.txt")){
            popFreqs <<- read.table2("cache/pop_freqs_genotypes_brca1_ordered.txt", header = T)
        } else if (file.exists("classify-simulated-population-age/pop_freqs_genotypes_brca1_ordered.txt")){
            popFreqs <<- read.table2("classify-simulated-population-age/pop_freqs_genotypes_brca1_ordered.txt", header = T)
        } else if (file.exists("../pop_freqs_genotypes_brca1_ordered.txt")){
            popFreqs <<- read.table2("../pop_freqs_genotypes_brca1_ordered.txt", header = T)
        } else {
            stop("pop_freqs_genotypes_brca1_ordered.txt does not exists!")
        }
    } else {
        popFreqs <<- read.table2("cache/pop_freqs_genotypes_brca1_ordered.txt", header = T)
    }
    
    if (!cmp_sample_size){
        if (is.null(cores)){
            doMC::registerDoMC(cores = detectCores())
        } else if (is.numeric(cores)){
            doMC::registerDoMC(cores = cores)
        } else {
            doMC::registerDoMC(cores = 6)
        }
    }
    
    # for (gen in seq(gen_start,gen_stop,gen_step)){
    data = foreach(gen=seq(gen_start,gen_stop,gen_step), .combine = rbind) %dopar% {
    # data = foreach(gen=seq(gen_start,30,gen_step), .combine = rbind) %dopar% {
    # for (gen in seq(gen_start,gen_stop,gen_step)){
        print(paste("Gen:", gen))
        
        j = 1
        
        # sims = 50
        filename = paste0(pre_path, "simulated_population/",in_gene,"_simulated-",genealogy,"-",samples,"_samples-",sims,"_simulations-generations_",gen,"-seed_",seed)
        brca1_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
        brca1_simulated_geno_plink <<- read.table2(file = paste0(filename, "-geno.txt"), header = T)
        
        # sims = 2
        sims2 = sims
        ages = list()
        CI_min = list()
        CI_max = list()
        # parameters.tw = list()
        for (i in 1:sims){
        # for (i in 1:2){
            print(paste0("Samples_", samples, "_sim_", i))
            mut = paste0("mutSim_BRCA1_",gen,"_",i); gene = "BRCA1_simulated"
            mut <<- mut; gene <<- gene
            readIndvAndFamData(gene, mut)
            # Adjust number of simulations if some failed
            if (nrow(matched_plink) < 2){
                sims2 = sims2-1
                next()
            }
            dst <- firstBreakDist(fam=F)
            #dst <- computeCustomDist(distMethod = firstBreakDist_cM, fam=F)
            hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = 1, doHapPlot = F)
            #plot_dendrogram(hc, k=1)
            brca.mut.cM <<- brca_cM_middle
            if (ancestral == "simulatedFounder"){
                haplotypes <- read.table(paste0(filename, "-founder.txt"), header = T)
                #haplotypes2 <<- as.integer(haplotypes[haplotypes$SNP == mut,-1])
                haplotypes2 <<- as.integer(haplotypes[i,-1])
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
                breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
                chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
            } else if (ancestral == "knownBreaks"){
                haplotypes2 <<- rep(0, ncol(matched_plink)-1)
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
                breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
                chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
            } else if (ancestral == "mostFreqBase"){
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
                haplotypes2 <<- findConsensus_plink(matched_plink2, cutoff = 0.5)
                breaks <- findHaplotypeBreaks_plink(matched_plink2, haplotypes2)
                chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
            } else {
                matched_plink2 <<- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID)
                if (minSampleBB == "auto") minSampleBB <- floor((gen+56.465)/9.975)/100 # famHaplotypes
                if (minSampleBB == "auto2") minSampleBB <- floor((gen+31.52)/16.43)/100 # famHaplotypesNoHet
                if (minSampleBB == "auto3") minSampleBB <- floor((gen+32.383)/4.014)/100 # haplotypes
                if (ancestral == "branchBoundIndep"){
                    res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = T, min_samples = minSampleBB)
                } else if (ancestral == "branchBound"){ #branchBound
                    res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = F, min_samples = minSampleBB)
                } else if (ancestral == "branchBoundIndep_alleleFreq"){
                    res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = T, min_samples = minSampleBB, alleleFreqs = alleleFreqs)
                } else if (ancestral == "branchBound_alleleFreq"){
                    res = branch_and_bound_ancestral(matched_plink, matched_plink2, cutoff = 0.5, indep_sides = F, min_samples = minSampleBB, alleleFreqs = alleleFreqs)
                } else {
                    stop("Ancestral method does not exist")
                }
                haplotypes2 <<- res[[1]]
                breaks <<- res[[2]]
                chr_pos <- mapNearestSNPs_plink(matched_plink2, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
            }
            if (cs_correction=="gandolfo"){
                chr_pos <- adjustHapLengthsGandolfo(chr_pos, chr_coords, popFreqs, haplotypes2, epsilon = epsilon)
            } else if (cs_correction=="gandolfoAverage"){
                chr_pos <- adjustHapLengthsGandolfoAverage(chr_pos, chr_coords, popFreqs, haplotypes2, epsilon = epsilon)
            } else if (cs_correction=="adjustLengths"){
                chr_pos <- adjustHapLengths(chr_pos, chr_coords, popFreqs, matched_plink2, epsilon = epsilon)
            } else if (cs_correction=="adjustHaploFreqs"){
                chr_pos <- adjustHaploFreqs(chr_pos, matched_plink, chr_coords, epsilon)
            }
            
            #age1 <- mutationAge_geneticDist(chr_pos, brca.mut.cM)
            age1 <- round(mutationAge_method1_median_with_bootstrap(chr_pos), digits = 2)

            #age3_mean <- mutationAge_method3(chr_pos)
            age3_mean <- round(mutationAge_method3_with_bootstrap(chr_pos), digits = 2)
            
            # age5 = mutationAge_method5(chr_pos, brca.mut.cM, start=max(round(gen-gen*3/4),1), stop=round(gen+gen*3/4), step=300)
            # age5 = mutationAge_method5(chr_pos, brca.mut.cM, start=max(age5-300,1), stop=age5+300, step=100)
            # age5 = mutationAge_method5(chr_pos, brca.mut.cM, start=max(age5-60,1), stop=age5+60, step=10)
            # age5 = mutationAge_method5(chr_pos, brca.mut.cM, start=max(age5-5,1), stop=age5+5, step=1)
            
            #age5 = method5_fast(chr_pos)
            age5 <- round(method5_fast_with_bootstrap(chr_pos), digits = 2)
            
            age6 = mutationAge_method6(chr_pos, brca.mut.cM, chr_coords = chr_coords, popFreqs = popFreqs, haplotypes2 = haplotypes2)
            # ages6[j] <- age6[1]
            # ages6_cor[j] <- age6[4]
            
            age6_cor = mutationAge_method6(chr_pos, brca.mut.cM, chance.sharing.correction = T, epsilon = epsilon,
                                           chr_coords = chr_coords, popFreqs = popFreqs, haplotypes2 = haplotypes2)
            # ages6_cs.cor[j] <- age6[1]
            # ages6_cor_cs.cor[j] <- age6[4]
            
            age8 = mutationAge_method8_tweedie(chr_pos, brca.mut.cM)
            # parameters.tw[[j]] = paste0(age8[4:6], collapse = " ")
            
            age9 = mutationAge_method9_rawGamma(chr_pos, brca.mut.cM)
            
            ages[[j]] = c(age1[1], age3_mean[1], age5[1], age6[1], age6[4], age6_cor[1], age6_cor[4], age8[1], age9[1])
            CI_min[[j]] =  c(age1[2], age3_mean[2], age5[2], age6[2], age6[5], age6_cor[2], age6_cor[5], age8[2], age9[2])
            CI_max[[j]] =  c(age1[3], age3_mean[3], age5[3], age6[3], age6[6], age6_cor[3], age6_cor[6], age8[3], age9[3])
            
            j = j+1
        }
        
        data = data.frame(actual_age=rep(gen, each = sims2),
                          method=c("method1", "method3_mean", "method5", "method6", "method6_cor", "method6_cs.cor", "method6_cor_cs.cor", "method8", "method9"),
                          sim=rep(1:sims2, each = 9),
                          samples=samples,
                          predicted_age=unlist(ages),
                          CI_min=unlist(CI_min),
                          CI_max=unlist(CI_max)
                          ) %>% arrange(method, sim)
        # data$tweedie.parameters = c(rep("", 8*sims2), unlist(parameters.tw))
        print(data)
        # return(data)
    }
    
    if (cs_correction=="gandolfo"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_gandolfo_", epsilon)
    } else if (cs_correction=="adjustLengths"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_adjustLengths_", epsilon)
    } else if (cs_correction=="gandolfoAverage"){
        custom = paste0("-minSample_",minSampleBB,"-cs_correction_gandolfoAverage_", epsilon)
    } else {
        custom = paste0("-minSample_",minSampleBB,"-epsilon_", epsilon)
    }
    
    #dir.create(path = paste0(pre_path, "compare-age-methods/"), showWarnings = F)
    dir.create(path = paste0(pre_path, "compare-age-methods/", ancestral, custom), showWarnings = F, recursive = T)
    write.table(data, 
                file = paste0(pre_path, "compare-age-methods/",ancestral, custom,"/",in_gene,
                                    "_age-comparisons_method1_3mean_5_6_6cor_6cscor_8_9",
                                    "-start_",gen_start,"-stop_",gen_stop,"-samples_",samples,
                                    "-genealogy_",genealogy,".txt"), 
                quote = F, row.names = F, col.names = T)
    
    # Tweedie parameters
    # write.table(x = setNames(as.data.frame(do.call(rbind,strsplit(data[,8][data[,8]!=""], split = " "))), c("mu","phi","p")),
    #             file = paste0(pre_path, "compare-age-methods/",ancestral, custom,"/",in_gene,
    #                           "_age-tweedie_parameters",
    #                           "-start_",gen_start,"-stop_",gen_stop,"-samples_",samples,
    #                           "-genealogy_",genealogy,".txt"),
    #             quote = F, row.names = F, col.names = T)
    
    if (!cmp_sample_size){
        compareAgeMethods_plots(gen_start = gen_start, gen_stop = gen_stop, gen_step = gen_step, ancestral_method = ancestral,
                                sims = sims, samples = samples, sim_type = sim_type, custom_string = custom, seed = seed,
                                pre_path = pre_path)
    }
}

run_compareAgeMethod_generateFiles_parallel <- function(g_start=NULL, g_stop=NULL, g_step=NULL, sim=NULL, minSamplesBB=0.25,
                                                        s_type = "knownBreaks", ancestral_method="simulatedFounder", 
                                                        eps = 0.01, cs_correction = "none", samples = 100){
    if (Sys.info()["sysname"] == "Linux"){
        setwd("/work/sduvarcall/haplotype-project")
        cores=26
        # cores=51
    } else {
        setwd("~/Documents/PhD/haplotype-project")
        # cores=6
        cores=1
    }
    source("Scripts/Haplotype-projekt-ver2.R")
    source("Scripts/Clustering.R")
    source("Scripts/Mutation-Age.R")
    source("Scripts/Mutation-age-tripleA.R")
    source("Scripts/Gandolfo_Speed_Mutation_Age_Estimation.R")
    source("Scripts/maxLikelihood.R")
    ## Read Data
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    morgan.brca1 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    #morgan.brca2 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
    alleleFreqs_brca1 <<- as.numeric(read.table2("cache/pop_freqs_minorAllele_brca1_ordered.txt", header = T))
    
    if (is.null(g_start)){
        #g_start = 10; g_stop = 2510; g_step = 20; sim = 10
        #g_start = 10; g_stop = 2510; g_step = 50; sim = 10
        #g_start = 10; g_stop = 2010; g_step = 50; sim = 10
        g_start = 10; g_stop = 1210; g_step = 10; sim = 50
    }
    
    system.time(compareAgeMethod_generateFiles_parallel(gen_start = g_start, gen_stop = g_stop, gen_step = g_step, sims = sim,
                                                        samples = samples, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
                                                        manual=FALSE, ancestral = ancestral_method, minSampleBB=minSamplesBB,
                                                        cs_correction=cs_correction, epsilon = eps, sim_type = s_type,
                                                        cores = cores))
    
    #for (s_type in c("famHaplotypes", "alleleFreqs")){
    #for (s_type in c("alleleFreqs")){
    # for (s_type in c("famHaplotypes")){
    # for (s_type in c("knownBreaks")){
    # for (s_type in c("knownBreaks")){
    #     #for (ancestral_method in c("simulatedFounder", "branchBoundIndep")){
    #     # for (ancestral_method in c("branchBoundIndep")){
    #     # for (ancestral_method in c("simulatedFounder")){
    #     for (ancestral_method in ancestral){
    #         for (n_samples in minSamplesBB){
    #             #for (cs_cor in c("adjustLengths", "gandolfo", "None")){
    #             for (cs_cor in c("none")){
    #             # for (cs_cor in c("adjustLengths")){
    #                 # for (e in c(0.01, 0.001,0.0001)){
    #                 for (e in c(0.01)){
    #                 #for (e in c(0.01)){
    #                     system.time(compareAgeMethod_generateFiles_parallel(gen_start = g_start, gen_stop = g_stop, gen_step = g_step, sims = sim,
    #                                                                         samples = 100, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
    #                                                                         manual=FALSE, ancestral = ancestral_method, minSampleBB=n_samples,
    #                                                                         cs_correction=cs_cor, epsilon = e, sim_type = s_type,
    #                                                                         cores = cores))
    #                 }
    #             }
    #         }
    #     }
    # }
    print("Done")
    
    
    # system.time(compareAgeMethod_generateFiles_parallel(gen_start = 10, gen_stop = 510, gen_step = 20, sims = 10,
    #                                                     samples = 100, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
    #                                                     manual=FALSE, ancestral = "branchBoundIndep", minSampleBB=5,
    #                                                     cs_correction="gandolfo", epsilon = 0.001, sim_type = "famHaplotypes",
    #                                                     cores = 6))
    # 
    # system.time(compareAgeMethod_generateFiles_parallel(gen_start = 10, gen_stop = 2010, gen_step = 50, sims = 10,
    #                                                     samples = 100, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
    #                                                     manual=FALSE, ancestral = "simulatedFounder", minSampleBB=5,
    #                                                     cs_correction="None", epsilon = 0.01, sim_type = "alleleFreqs",
    #                                                     cores = 6))
    # 
    # system.time(compareAgeMethod_generateFiles_parallel(gen_start = 10, gen_stop = 2010, gen_step = 50, sims = 10,
    #                                                     samples = 100, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
    #                                                     manual=FALSE, ancestral = "simulatedFounder", minSampleBB=5,
    #                                                     cs_correction="gandolfo", epsilon = 0.01, sim_type = "alleleFreqs",
    #                                                     cores = 6))
    # 
    # system.time(compareAgeMethod_generateFiles_parallel(gen_start = 10, gen_stop = 2010, gen_step = 50, sims = 10,
    #                                                     samples = 100, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
    #                                                     manual=FALSE, ancestral = "simulatedFounder", minSampleBB=5,
    #                                                     cs_correction="adjustLengths", epsilon = 0.01, sim_type = "alleleFreqs",
    #                                                     cores = 6))
}

run_compareSampleSize_generateFiles_parallel <- function(sample_sizes=NULL, gen=NULL, sim=NULL, minSamplesBB=0.25,
                                                        s_type = "knownBreaks", ancestral_method="simulatedFounder", 
                                                        eps = 0.01, cs_correction = "none"){
    if (Sys.info()["sysname"] == "Linux"){
        setwd("/work/sduvarcall/haplotype-project")
        cores=25
        # cores=51
    } else {
        setwd("~/Documents/PhD/haplotype-project")
        # cores=6
        cores=1
    }
    source("Scripts/Haplotype-projekt-ver2.R")
    source("Scripts/Clustering.R")
    source("Scripts/Mutation-Age.R")
    source("Scripts/Mutation-age-tripleA.R")
    source("Scripts/Gandolfo_Speed_Mutation_Age_Estimation.R")
    source("Scripts/maxLikelihood.R")
    ## Read Data
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    morgan.brca1 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    #morgan.brca2 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
    alleleFreqs_brca1 <<- as.numeric(read.table2("cache/pop_freqs_minorAllele_brca1_ordered.txt", header = T))
    
    doMC::registerDoMC(cores = cores)
    # for (samples in sample_sizes){
    foreach(samples=sample_sizes) %dopar% {
        pre_path = paste0("classify-simulated-population-age/BRCA1_simulated-test_sample_size-",s_type,"-starGenealogy-",sim,"_simulations-generations_",gen,"_",gen,"_step",gen,"-seed_42/")
        system.time(compareAgeMethod_generateFiles_parallel(gen_start = gen, gen_stop = gen, gen_step = gen, sims = sim,
                                                        samples = samples, seed = 42, in_gene = "BRCA1", genealogy="starGenealogy", 
                                                        manual=FALSE, ancestral = ancestral_method, minSampleBB=minSamplesBB,
                                                        cs_correction=cs_correction, epsilon = eps, sim_type = s_type,
                                                        cores = 1, pre_path = pre_path, cmp_sample_size = T))
    }
    
    print("Done")
}

test_parameters <- function(){
    if (Sys.info()["sysname"] == "Linux"){
        setwd("/work/sduvarcall/haplotype-project")
        cores=24
    } else {
        setwd("~/Documents/PhD/haplotype-project")
        cores=6
    }
    source("Scripts/Haplotype-projekt-ver2.R")
    source("Scripts/Clustering.R")
    source("Scripts/Mutation-Age.R")
    source("Scripts/Mutation-age-tripleA.R")
    source("Scripts/maxLikelihood.R")
    source("Scripts/Gandolfo_Speed_Mutation_Age_Estimation.R")
    ## Read Data
    chr_coords_all <<- read.table("input/139_ouh_june_2017/qry_139_ouh_geno_snps_with_position-fixed-with_centiMorgan_coords_from_linear_interpolation.txt", header = T, sep = "\t")
    morgan.brca1 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    #morgan.brca2 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
    alleleFreqs_brca1 <<- as.numeric(read.table2("cache/pop_freqs_minorAllele_brca1_ordered.txt", header = T))
    
    samples=100; genealogy="starGenealogy"
    #g_start = 10; g_stop = 2510; g_step = 20; sim = 10
    #g_start = 10; g_stop = 2510; g_step = 50; sim = 50
    #g_start = 10; g_stop = 3010; g_step = 50; sim = 20
    #g_start = 10; g_stop = 510; g_step = 10; sim = 20
    g_start = 10; g_stop = 1210; g_step = 10; sim = 50
    #g_start = 10; g_stop = 510; g_step = 50; sim = 10
    #g_start = 50; g_stop = 51; g_step = 10; sim = 20; genealogy="correlatedGenealogy"; samples=200
    # g_start = 100; g_stop = 101; g_step = 10; sim = 50; genealogy="correlatedGenealogy"; samples=100
    
    
    #for (s_type in c("famHaplotypes", "alleleFreqs")){
    #for (s_type in c("alleleFreqs")){
    for (s_type in c("knownBreaks")){
    #for (s_type in c("famHaplotypes")){
        #for (ancestral_method in c("branchBoundIndep")){
        for (ancestral_method in c("simulatedFounder")){
            #for (n_samples in c(2, 5, 10, 15, 20, 25)){
            for (n_samples in c(5)){
                #for (cs_cor in c("adjustLengths", "gandolfo", "gandolfoAverage", "none")){
                for (cs_cor in c("none")){
                    for (e in c(0.01)){
                    #for (e in c(0.01, 0.001, 0.0001, 0.00001)){
                    #for (e in c(0.001, 0.0001, 0.00001)){
                        system.time(compareAgeMethod_generateFiles_parallel(gen_start = g_start, gen_stop = g_stop, gen_step = g_step, sims = sim,
                                                                            samples = samples, seed = 42, in_gene = "BRCA1", genealogy=genealogy, 
                                                                            manual=FALSE, ancestral = ancestral_method, minSampleBB=n_samples,
                                                                            cs_correction=cs_cor, epsilon = e, sim_type = s_type,
                                                                            cores = cores))
                    }
                }
            }
        }
    }
    print("Done")
}

compare_gamma_and_tweedie_CI <- function(){
    gen_start = 10
    gen_step = 50
    k = 600#400
    
    data = read.table("classify-simulated-population-age/BRCA1_simulated-knownBreaks-correlatedGenealogy-100_samples-50_simulations-generations_100_101_step10-seed_42/compare-age-methods/simulatedFounder-minSample_5-epsilon_0.01/BRCA1_age-adjustedLengths-comparisons_method1_3mean_5_6_6cor_6cscor-start_100-stop_101-samples_100-genealogy_correlatedGenealogy.txt", header = T)

    methods = c("3_mean",5,6,8,9)
    df = filter(data, method %in% paste0("method",methods), actual_age < k) %>% arrange(method, actual_age)
    count(df, method, actual_age) %>% summarise(sum(n))
    
    # If data from methods above
    # df = data
    
    df$method = as.character(df$method)
    df$method[df$method=="method3_mean"] = "method3"
    # df = filter(df, predicted_age < 400)
    df2 = df %>% group_by(actual_age,method) %>% summarise_all(median)
    # df2 = df %>% group_by(actual_age,method,samples) %>% summarise(predicted_age=median(predicted_age), 
    #                                                                CI_min=min(CI_min), CI_max=max(CI_max))
    df3 = df %>% group_by(actual_age,method) %>% summarise_all(mean)
    
    
    # df2 %>% #mutate(predicted_age = predicted_age-actual_age) %>%
    # df3 %>% filter(method %in% paste0("method ", c(3,4))) %>%
    df2 %>% filter(method %in% paste0("method ", c(2,3,4))) %>%
        ggplot(aes(x = actual_age, y = predicted_age, ymin = CI_min, ymax = CI_max, fill=method, colour = method)) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        #geom_smooth() +
        geom_ribbon(alpha=0.1) +
        geom_abline(intercept = 0, slope = 1, col = "black") +
        #geom_point() +
        geom_line() +
        theme_bw() +
        xlab("Predicted age") +
        ylab("Actual age")
    

    df2 %>% 
        ggplot(aes(x=method, y=predicted_age, ymin = CI_min, ymax = CI_max, fill=method, colour = method)) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        #geom_boxplot() +
        #geom_ribbon(alpha=0.1) +
        #geom_abline(intercept = 0, slope = 1, col = "black") +
        geom_point() +
        geom_point(aes(x=df2$method, y=df2$CI_min), alpha=0.5) +
        geom_point(aes(x=df2$method, y=df2$CI_max), alpha=0.5) +
        #geom_line() +
        theme_bw() +
        xlab("Predicted age") +
        ylab("Actual age")
    
    #df %>% filter(predicted_age < 400) %>%
    #ggplot(aes(x=method, y=predicted_age, fill = method)) +
    p = ggplot() +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        #geom_boxplot() +
        #geom_ribbon(alpha=0.1) +
        #geom_histogram(bins=50) +
        geom_boxplot(data = df, mapping = aes(x=method, y=predicted_age, color = method), outlier.size = 1) +
        geom_abline(intercept = 100, slope = 0, col = "gray", linetype = 2) +
        #stat_summary(fun.data = mean, geom="errorbar") +
        geom_errorbar(aes(ymax=df3$CI_max, ymin=df3$CI_min, x=df3$method, y=df3$predicted_age)) +
        geom_point(aes(x=df3$method, y=df3$predicted_age)) +
        # geom_errorbar(aes(ymax=df2$CI_max, ymin=df2$CI_min, x=df2$method, y=df2$predicted_age), color="gray") +
        # geom_point(aes(x=df2$method, y=df2$predicted_age), color="gray") +
        # geom_line() +
        theme_bw() +
        #guides(fill=F) +
        xlab("Method") +
        ylab("Predicted age")
    print(p)
    ggsave(filename = "classify-simulated-population-age/BRCA1_simulated-knownBreaks-correlatedGenealogy-100_samples-50_simulations-generations_100_101_step10-seed_42/compare-age-methods/boxplot_with_confidence_intervals.png", plot = p)

    ################################################
    ### deprecated
    ################################################
    data = read.table("classify-simulated-population-age/BRCA1_simulated-knownBreaks-starGenealogy-100_samples-20_simulations-generations_10_510_step10-seed_42/compare-age-methods/simulatedFounder-minSample_2-epsilon_0.01/BRCA1_age-adjustedLengths-comparisons_method1_3mean_5_6_6cor_6cscor-start_10-stop_510-samples_100-genealogy_starGenealogy.txt", header = T)
    data = read.table("classify-simulated-population-age/BRCA1_simulated-knownBreaks-starGenealogy-100_samples-20_simulations-generations_10_3010_step50-seed_42/compare-age-methods/simulatedFounder-minSample_2-epsilon_0.01/BRCA1_age-adjustedLengths-comparisons_method1_3mean_5_6_6cor_6cscor-start_10-stop_3010-samples_100-genealogy_starGenealogy.txt", header = T)
    data = read.table("classify-simulated-population-age/BRCA1_simulated-knownBreaks-starGenealogy-100_samples-50_simulations-generations_10_2510_step50-seed_42/compare-age-methods/simulatedFounder-minSample_2-epsilon_0.01/BRCA1_age-adjustedLengths-comparisons_method1_3mean_5_6_6cor_6cscor-start_10-stop_2510-samples_100-genealogy_starGenealogy.txt", header = T)
    data = read.table("classify-simulated-population-age/BRCA1_simulated-knownBreaks-starGenealogy-100_samples-50_simulations-generations_10_1210_step10-seed_42/compare-age-methods/simulatedFounder-minSample_2-epsilon_0.01/BRCA1_age-adjustedLengths-comparisons_method1_3mean_5_6_6cor_6cscor-start_10-stop_1210-samples_100-genealogy_starGenealogy.txt", header = T)
    methods = c(6,8,9)
    df = filter(data, method %in% paste0("method",methods), actual_age < k) %>% arrange(method, actual_age)
    df = cbind(df,
               # CI_min = filter(data, grepl("CI_min",method)) %>% arrange(method, actual_age) %>% .$predicted_age,
               # CI_max = filter(data, grepl("CI_max",method)) %>% arrange(method, actual_age) %>% .$predicted_age)
               CI_min = filter(data, grepl("CI_min",method), actual_age < k) %>% arrange(method, actual_age) %>% .$predicted_age,
               CI_max = filter(data, grepl("CI_max",method), actual_age < k) %>% arrange(method, actual_age) %>% .$predicted_age
               )
    df = filter(df, method %in% paste0("method",c(8,9)))
    df2 = df %>% group_by(actual_age,method) %>% summarise_all(median)
    
    
    df2 %>% #mutate(predicted_age = predicted_age-actual_age) %>%
        ggplot(aes(x = actual_age, y = predicted_age, ymin = CI_min, ymax = CI_max, fill=method, colour = method)) +
        scale_color_brewer(palette = "Set1") +
        scale_fill_brewer(palette = "Set1") +
        #geom_smooth() +
        geom_ribbon(alpha=0.1) +
        geom_abline(intercept = 0, slope = 1, col = "black") +
        #geom_point() +
        geom_line() +
        theme_bw() +
        xlab("Predicted age") +
        ylab("Actual age")
    
            
    # data %>% filter(!grepl("cor", method)) %>% #group_by(actual_age,method) %>% summarise_all(mean) %>%#filter(method == "method8") %>%
    #     ggplot(aes(x = actual_age, y = predicted_age, color = method)) +
    #     geom_abline(intercept = 0, slope = 1, col = "darkblue") +
    #     #geom_point() +
    #     geom_smooth() +
    #     theme_bw()
}



