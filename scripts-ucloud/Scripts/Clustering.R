library(amap)
#library(ggdendro)
library(dendextend)
library(RColorBrewer)
library(proxy)
library(factoextra)
library(NbClust)
library(dbscan)
library(fpc)

#setwd("~/Documents/PhD/haplotype-project")
#source("Scripts/Haplotype-projekt-ver2.R")
#source("Maximum-Likelihood.R")

#mut = "c.5266dupC"; gene = "BRCA1"
#mut = "c.1016delA"; gene = "BRCA1"
#mut = "c.5123C>A"; gene = "BRCA1"
#mut = "c.7617+1G>A"; gene = "BRCA2"
#mut = "c.6275_6276delTT"; gene = "BRCA2"
#mut = "c.3331_3334delCAAG"; gene = "BRCA1"

###############################################################################
### Read Data
###############################################################################

readFamData <- function(){
    # Load family geno data for brca1
    if (file.exists("cache/famHaplotypes-BRCA1-geno.RData") || file.exists("cache/famHaplotypes-BRCA1-geno_plink_format.RData")){
        #load("cache/famHaplotypes-BRCA1-geno.RData", envir = .GlobalEnv)
        load("cache/famHaplotypes-BRCA1-geno_plink_format.RData", envir = .GlobalEnv)
    }else{
        famHaplotypes_brca1_geno <<- read.table(file = "input/139_ouh_june_2017/famHaplotypes-BRCA1-geno.txt", sep = "\t", header = T)
        save(famHaplotypes_brca1_geno, file = "cache/famHaplotypes-BRCA1-geno.RData")
        famHaplotypes_brca1_geno_plink <<- read.table(file = "input/139_ouh_june_2017/famHaplotypes-BRCA1-geno_fam_plink_format.txt", sep = "\t", header = T)
        save(famHaplotypes_brca1_geno_plink, file = "cache/famHaplotypes-BRCA1-geno_plink_format.RData")
    }

    # To be done - single-fam-members for BRCA1

    # Load family geno data for brca2
    if (file.exists("cache/famHaplotypes-BRCA2-geno.RData") || file.exists("cache/famHaplotypes-BRCA2-geno_plink_format.RData")){
        #load("cache/famHaplotypes-BRCA2-geno.RData", envir = .GlobalEnv)
        load("cache/famHaplotypes-BRCA2-geno_plink_format.RData", envir = .GlobalEnv)
    }else{
        famHaplotypes_brca2_geno <<- read.table(file = "input/139_ouh_june_2017/famHaplotypes-BRCA2-geno.txt", sep = "\t", header = T)
        save(famHaplotypes_brca2_geno, file = "cache/famHaplotypes-BRCA2-geno.RData")
        famHaplotypes_brca2_geno_plink <<- read.table(file = "input/139_ouh_june_2017/famHaplotypes-BRCA2-geno_fam_plink_format.txt", sep = "\t", header = T)
        save(famHaplotypes_brca2_geno_plink, file = "cache/famHaplotypes-BRCA2-geno_plink_format.RData")
    }

    # get geno information for single-fam-members
    # singleHaplotypes_pheno <- brca_data %>% group_by(FamCode) %>% filter(n()==1)
    # singleHaplotypes_geno <- extractSamples(singleHaplotypes_pheno, brca2_geno_merged, chr_coords)
    # singleHaplotypes_geno_plink <- extractSamples(singleHaplotypes_pheno, brca2_geno_plink, chr_coords)
}

extractFamPlink <- function(req_fam_members=2){
    if (gene == "BRCA1"){
        #famHaplotypes_geno <<- famHaplotypes_brca1_geno
        famHaplotypes_geno_plink <<- famHaplotypes_brca1_geno_plink
    } else if (gene == "BRCA2") {
        #famHaplotypes_geno <<- famHaplotypes_brca2_geno
        famHaplotypes_geno_plink <<- famHaplotypes_brca2_geno_plink
    } else if (gene == "BRCA2_toySample"){
        famHaplotypes_geno_plink <<- brca2_toy_fam_geno
    } else {
        # Fix for simulated data with no family info
        famHaplotypes_geno_plink <<- data.frame(SNP=0,SNP1=0)
        matched_plink_fam <<- data.frame(SNP=0,SNP1=0)
        #famHaplotypes_geno_plink <<- brca_geno_plink
        #matched_plink_fam <<- brca_geno_plink
        brca_data.fam <<- brca_data
        return()
    }

    # Check number of families with more than one fam member - Seems to be a little buggy for some reason!
    numFam = length(famHaplotypes_geno_plink[na.omit(fmatch(famHaplotypes_geno_plink$FamCode, brca_data$FamCode)),]$FamCode)
    print(numFam)

    # Extract pheno and geno data from families with more than one member
    MultiplefamMembers_pheno <- brca_data %>% group_by(FamCode) %>% filter(n()>=req_fam_members) %>% select(FamCode, Onc_ID)
    print(nrow(MultiplefamMembers_pheno))
    #MultiplefamMembers_pheno <- brca_data %>% count(FamCode) %>% filter(n>1) %>% .$FamCode
    #famHaplotypes_geno <- famHaplotypes_geno[fmatch(unique(MultiplefamMembers_pheno$FamCode), famHaplotypes_geno$FamCode), c(1, 2, na.omit(fmatch(chr_coords$SNP, colnames(famHaplotypes_geno))))]
    #famHaplotypes_geno_plink <- famHaplotypes_geno_plink[fmatch(unique(MultiplefamMembers_pheno$FamCode), famHaplotypes_geno_plink$FamCode), c(1, 2, na.omit(fmatch(chr_coords$SNP, colnames(famHaplotypes_geno_plink))))]
    #matched_fam <- famHaplotypes_geno[na.omit(fmatch(MultiplefamMembers_pheno$Onc_ID, famHaplotypes_geno$SNP)), c(1, 2, na.omit(fmatch(chr_coords$SNP, colnames(famHaplotypes_geno))))]
    matched_plink_fam <- famHaplotypes_geno_plink[na.omit(fmatch(MultiplefamMembers_pheno$Onc_ID, famHaplotypes_geno_plink$SNP)), c(1, 2, na.omit(fmatch(chr_coords$SNP, colnames(famHaplotypes_geno_plink))))]

    # Extract pheno and geno data from families with one member
    # famCode <- filter(brca_data, Onc_ID %in% singleHaplotypes_geno$SNP) %>% select(FamCode)
    # singleHaplotypes_geno <- cbind(famCode, singleHaplotypes_geno)
    # singleHaplotypes_geno_plink <- cbind(famCode, singleHaplotypes_geno_plink)

    # Combine single and multi families
    # matched_fam <<- rbind(famHaplotypes_geno, singleHaplotypes_geno)
    # matched_plink_fam <<- rbind(famHaplotypes_geno_plink, singleHaplotypes_geno_plink)

    #matched_fam <<- famHaplotypes_geno
    #matched_plink_fam <<- famHaplotypes_geno_plink

    #rownames(matched_fam) <<- matched_fam$FamCode
    #rownames(matched_plink_fam) <<- matched_plink_fam$FamCode
    #rownames(matched_fam) <- matched_fam$SNP
    #rownames(matched_plink_fam) <- matched_plink_fam$SNP

    #matched_fam <<- matched_fam
    matched_plink_fam <<- matched_plink_fam

    brca_data.fam <<- filter(brca_data, Onc_ID %in% matched_plink_fam$SNP)
    #print(dim(matched))
    #print(dim(brca_data))
    #return(matched_plink_fam)
}

readIndvAndFamData <- function(gene, mut){
    ## Read input data into global variables
    # readData()
    ## Override matched object with fam information - if preferred
    # readFamData()

    ## Prepare data frames
    prepare_dataframe(mut, gene)

    extractFamPlink()
}

###############################################################################
### Hierachical clustering
###############################################################################
hiearchical_clustering <- function(dst, clustMethod="complete", k=3, doHapPlot=F, printInfo=F){
    hc <- hclust(dst, clustMethod)
    # hc <- hclust(dst, "ward.D2")
    # hc <- hclust(dst, "ward.D")
    # hc <- hclust(dst, "complete")
    # hc <- hclust(dst, "ave")

    cluster_groups <<- cutree(hc, k=k)
    #names(cluster_groups) <- snp_names

    #brca_data.filtered = brca_data
    #brca_data.filtered <- filter(brca_data, Country == "SPAIN" | Country == "DENMARK")
    brca_data.filtered["cluster_groups"] <<- cluster_groups
    brca_data <<- brca_data.filtered

    # Print group information
    if (printInfo){
        for (i in unique(cluster_groups)){
            #print(subset(brca_data.filtered, cluster_groups == i) %>% select(FamCode, Onc_ID, Country, cluster_groups))
            print(as.data.frame(subset(brca_data.filtered, cluster_groups == i) %>% count(Country)))
        }
    }

    if (doHapPlot){
        for (i in unique(cluster_groups)){
            if (nrow(subset(brca_data.filtered, cluster_groups == i)) > 5){
                haplotype_mutation(mut, prepare = F, group = i)
                consensus <- paste0("cluster_",i,"_of_",k)
                ggsave(paste0("plots/", brca_name, "-", mut, "-", consensus, "-", num_samples, "_samples.png"),
                       scale = 1, width = 16, height = 10)
            }
        }
    }
    return(hc)
}

plot_dendrogram <- function(hc, k = 3, leaf = F, sub_k = F, plot_type="country"){
    # Draw dendrogram
    if (leaf && sub_k){
        dend <- get_subdendrograms(as.dendrogram(hc), sub_k)[[leaf]]
    } else {
        dend <- as.dendrogram(hc)
    }
    FamCode <- brca_data$FamCode[order.dendrogram(dend)]
    #Country <- brca_data.filtered$Country[order.dendrogram(dend)]
    #Country <- brca_data.filtered$Country
    #cols <- 1:length(unique(Country))
    # if (length(unique(Country)) < 10){
    #     #cols <- brewer.pal(max(length(unique(Country)), 3), "Set1")[1:length(unique(Country))]
    #     cols <- brewer.pal(max(length(unique(Country)), 3), "Dark2")[1:length(unique(Country))]
    # } else {
    #     #cols <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(Country)))
    #     cols <- colorRampPalette(brewer.pal(9, "Dark2"))(length(unique(Country)))
    # }
    # names(cols) <- unique(Country)
    # cols <<- cols
    #cols <<- colors_plot(pheno_data = brca_data.filtered, palette = "Dark2")
    
    if (plot_type == "country"){
        
        cols <- colors_plot(pheno_data = brca_data[order.dendrogram(dend), ], palette = "sasha")
        Country_ordered <- brca_data$Country[order.dendrogram(dend)]
        label_cols <- cols[as.character(Country_ordered)]
        
        #title = "Country"
        title = NULL
        
        # set("branches_col", cols) %>%
        # dend %>% set("branches_k_color", k = k) %>%
        #     set("labels", FamCode) %>% set("labels_colors", country_cols) %>%
        #     set("labels_cex", 1) %>% set("branches_lwd", 3) %>% plot
        #dend %>% rect.dendrogram(k=k, border = 8, lty = 5, lwd = 2) # Draw squares around clusters
        #legend('topright', names(cols), col=cols, fill = cols, cex = 1.2, border = F, y.intersp = 1, box.lwd = F, bty = "o", inset = c(0.01,0)) # lty=1 for lines in legend
        #legend(x=1600, y=120, names(cols), col=cols, fill = cols, cex = 1.2, border = F, y.intersp = 1, box.lwd = F, bty = "o", inset = c(0.01,0)) # lty=1 for lines in legend
        
        #plot(hc, hang = -1, cex=0.6, labels = famCode_labels)
        #rect.hclust(hc, k=k, border="darkblue")
    } else if (plot_type == "color_labels_as_branches"){
        cols <- colors_plot(color_names = 1:max(brca_data$cluster_groups))
        label_cols <- cols[brca_data$cluster_groups[order.dendrogram(dend)]]
        
        title = "Group"
        
    } else if (plot_type == "ethnicity"){
        ethnicity = as.character(brca_data$ethnicity)
        ethnicity[ethnicity==""] = "NA"
        cols <- colors_plot(color_names = unique(ethnicity[order.dendrogram(dend)]))
        label_cols <- cols[as.character(ethnicity[order.dendrogram(dend)])]
        title="Ethnicity"
    }
    branch_colors = colors_plot(color_names = 1:k)
    #dend <- color_branches(dend, clusters = cluster_groups[order.dendrogram(dend)], col = unique(branch_colors[cluster_groups[order.dendrogram(dend)]]))
    dend <- color_branches(dend, clusters = cluster_groups[order.dendrogram(dend)], groupLabels = unique(cluster_groups[order.dendrogram(dend)]))
    dend %>% #set("branches_k_color", k = k, col = branch_colors) %>%
        set("labels", FamCode) %>% set("labels_colors", label_cols) %>%
        set("labels_cex", 1) %>% set("branches_lwd", 3) %>% plot
    legend('topright', names(cols), col=cols, fill = cols, cex = 1.2, border = F, y.intersp = 1, box.lwd = F, bty = "o", inset = c(0.01,0), title = title)

    # label_cols <- colors_plot(palette = "sasha")[as.character(distinct(chr_pos, sample_id, Country)$Country)]
    # country_freq <- names(sort(table(names(label_cols)), decreasing = T))
    # val_col <- label_cols[country_freq]
    # legend("bottom", names(val_col), col = val_col, fill = val_col)
}

###############################################################################
### Euclidean distance
###############################################################################
euclideanDist <- function(fam=F, snp_rm = 1000){
    #snp_rm = 0
    if (fam){
        dst <- dist(matched_plink_fam[, (2+snp_rm):(ncol(matched_plink)-snp_rm)])
        brca_data.filtered <<- brca_data.fam
    } else {
        dst <- dist(matched_plink[, (2+snp_rm):(ncol(matched_plink)-snp_rm)])
        brca_data.filtered <<- brca_data.single
    }
    # hc <- hclust(dst, "ward.D2")
    # hc <- hclust(dst, "ward.D")
    # hc <- hclust(dst, "complete")
    # hc <- hclust(dst, "mcquitty")
    # hc <- hclust(dst, "ave")
    return(dst)
}

###############################################################################
### Euclidean distance - modified according to BRCA mut distance
###############################################################################
euclideanDist_BRCAdistWeighted <- function(fam=F, snp_rm = 1000){
    euclidean_brca_dist <- function(x,y){
        #x = matched_plink[1,2:ncol(matched_plink)]
        #y = matched_plink[2,2:ncol(matched_plink)]
        #print((x-y)/brca_dist)
        diff = (x-y)^2
        #diff[diff==0 & brca_dist < -1] = -1
        #diff[diff==0 & brca_dist > -1] = -1/brca_dist[diff==0 & brca_dist > -1]
        diff[diff==0] = -1/brca_dist[diff==0]
        dist <- sqrt(max(sum(diff), 0))
        #print(dist)
        return(dist)
        #return( sqrt(sum( (x-y)^2 )) )
    }
    chr_pos.dist <- mapSNPsToCoords(chr_coords, matched_plink)
    # Distance in Kb = 1000, Mb = 1000000
    #brca_dist = (abs(chr_pos[chr_pos$SNP %in% names(matched_plink),]$brca_distance) / 1000000)
    brca_dist <<- (abs(chr_pos.dist$brca_distance) / 1000)[(1+snp_rm):(nrow(chr_pos.dist)-snp_rm)]
    if (fam){
        dst = proxy::dist(matched_plink_fam[, (3+snp_rm):(ncol(matched_plink_fam)-snp_rm)], method = euclidean_brca_dist)
        brca_data.filtered <<- brca_data.fam
    } else {
        dst = proxy::dist(matched_plink[, (2+snp_rm):(ncol(matched_plink)-snp_rm)], method = euclidean_brca_dist)
        brca_data.filtered <<- brca_data.single
    }

    # hc <- hclust(dst, "ward.D2")
    # hc <- hclust(dst, "ward.D")
    # hc <- hclust(dst, "complete")
    # hc <- hclust(dst, "ave")

    #hiearchical_clustering(hc, k=4, F)
    #hiearchical_clustering(hc, k=3, T)
    return(dst)
}

###############################################################################
### amap distances (parallel computed) - e.g. Centered Pearson, euclidean
###############################################################################
amapDist <- function(distMethod="euclidean", fam=F, snp_rm = 1000){
    if (fam){
        #dst = Dist(matched_plink_fam[,(3+snp_rm):(ncol(matched_plink_fam)-snp_rm)], method = "correlation")
        dst = Dist(matched_plink_fam[,(3+snp_rm):(ncol(matched_plink_fam)-snp_rm)], method = distMethod)
        brca_data.filtered <<- brca_data.fam
    } else{
        #dst = Dist(matched_plink[,(2+snp_rm):(ncol(matched_plink)-snp_rm)], method = "correlation")
        dst = Dist(matched_plink[,(2+snp_rm):(ncol(matched_plink)-snp_rm)], method = distMethod)
        brca_data.filtered <<- brca_data.single
    }

    # hc <- hclust(dst, "ward.D2")
    # hc <- hclust(dst, "ward.D")
    # hc <- hclust(dst, "complete")
    # hc <- hclust(dst, "ave")
    #
    # hiearchical_clustering(hc, k=2, F)
    # hiearchical_clustering(hc, k=3, T)
    return(dst)
}

###############################################################################
### Define similarity based on length of similar elements around brca gene
###############################################################################
firstUpperBreakDist <- function(x,y){
    #x = matched_plink[1,-1]
    #y = matched_plink[2,-1]
    # Find homozygote break
    condition = (x-y == -2 | x-y == 2)
    # Split vector according to positive or negative distance to brca gene
    cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
    # get length (in Kb) to brca gene for first homozygous break on either side
    min_index=min(which(cond_pos==T), length(cond_pos))
    #print(which(cond_pos==T))
    pos.len = chr_pos.dist[chr_pos.split+min_index,"brca_distance"] / 1000000
    # Convert from cM to morgan
    #pos.len = chr_pos.dist[chr_pos.split+min_index,6] * 100
    # Compute similarity based on length with no homozygous break
    similarity = 1/(pos.len)
    return(similarity)
}

firstUpperBreakDist_cM <- function(x,y){
    #x = matched_plink[1,-1]
    #y = matched_plink[2,-1]
    # Find homozygote break
    condition = (x-y == -2 | x-y == 2)
    # Split vector according to positive or negative distance to brca gene
    cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
    # get length (in Kb) to brca gene for first homozygous break on either side
    min_index=min(which(cond_pos==T), length(cond_pos))
    #print(which(cond_pos==T))
    pos.len = chr_pos.dist[chr_pos.split+min_index,"brca_distance_cM"] * 100
    # Convert from cM to morgan
    #pos.len = chr_pos.dist[chr_pos.split+min_index,6] * 100
    # Compute similarity based on length with no homozygous break
    similarity = 1/(pos.len)
    return(similarity)
}

firstUnderBreakDist <- function(x,y){
    #x = matched_plink[1,2:ncol(matched_plink)]
    #y = matched_plink[2,2:ncol(matched_plink)]
    # Find homozygote break
    condition = (x-y == -2 | x-y == 2)
    # Split vector according to positive or negative distance to brca gene
    cond_neg = condition[1:chr_pos.split.left]
    # get length (in Kb) to brca gene for first homozygous break on either side
    max_index=max(which(cond_neg==T), 1)
    neg.len = abs(chr_pos.dist[max_index,"brca_distance"]) / 1000000
    # Convert from cM to morgan
    #neg.len = abs(chr_pos.dist[max_index,6]) * 100
    # Compute similarity based on length with no homozygous break
    similarity = 1/(neg.len)
    return(similarity)
}

firstUnderBreakDist_cM <- function(x,y){
    #x = matched_plink[1,2:ncol(matched_plink)]
    #y = matched_plink[2,2:ncol(matched_plink)]
    # Find homozygote break
    condition = (x-y == -2 | x-y == 2)
    # Split vector according to positive or negative distance to brca gene
    cond_neg = condition[1:chr_pos.split.left]
    # get length (in Kb) to brca gene for first homozygous break on either side
    max_index=max(which(cond_neg==T), 1)
    neg.len = abs(chr_pos.dist[max_index,"brca_distance_cM"]) * 100
    # Convert from cM to morgan
    #neg.len = abs(chr_pos.dist[max_index,6]) * 100
    # Compute similarity based on length with no homozygous break
    similarity = 1/(neg.len)
    return(similarity)
}

firstBreakDist_cM <- function(x,y){
    # Find homozygote break
    condition = (x-y == -2 | x-y == 2)
    # Split vector according to positive or negative distance to brca gene
    cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
    cond_neg = condition[1:chr_pos.split.left]
    # get length (in Kb) to brca gene for first homozygous break on either side
    min_index=min(which(cond_pos==T), length(cond_pos))
    #print(which(cond_pos==T))
    #pos.len = chr_pos.dist[chr_pos.split+min_index,"brca_distance"] / 1000000
    pos.len = chr_pos.dist[chr_pos.split+min_index,"brca_distance_cM"] #* 100
    max_index=max(which(cond_neg==T), 1)
    #print(which(cond_neg==T))
    #neg.len = abs(chr_pos.dist[max_index,"brca_distance"]) / 1000000
    neg.len = abs(chr_pos.dist[max_index,"brca_distance_cM"]) #* 100
    # Compute similarity based on length with no homozygous break
    similarity = 1/(pos.len+neg.len)
    return(similarity)
}

firstBreakEuclideanDist <- function(x,y){
    # Find homozygote break
    condition = (x-y == -2 | x-y == 2)
    # Split vector according to positive or negative distance to brca gene
    cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
    cond_neg = condition[1:chr_pos.split.left]
    # get length (in Kb) to brca gene for first homozygous break on either side
    min_index=min(which(cond_pos==T), length(cond_pos))
    #print(which(cond_pos==T))
    pos.len = chr_pos.dist[chr_pos.split+min_index,"brca_distance"] / 1000000
    max_index=max(which(cond_neg==T), 1)
    #print(which(cond_neg==T))
    neg.len = abs(chr_pos.dist[max_index,"brca_distance"]) / 1000000
    # Compute similarity based on length with no homozygous break
    similarity = 1/(sqrt(pos.len^2+neg.len^2))
    return(similarity)
}

firstBreakEuclideanDist_cM <- function(x,y){
    # Find homozygote break
    condition = (x-y == -2 | x-y == 2)
    # Split vector according to positive or negative distance to brca gene
    cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
    cond_neg = condition[1:chr_pos.split.left]
    # get length (in Kb) to brca gene for first homozygous break on either side
    min_index=min(which(cond_pos==T), length(cond_pos))
    #print(which(cond_pos==T))
    pos.len = chr_pos.dist[chr_pos.split+min_index,"brca_distance_cM"] * 100
    max_index=max(which(cond_neg==T), 1)
    #print(which(cond_neg==T))
    neg.len = abs(chr_pos.dist[max_index,"brca_distance_cM"]) * 100
    # Compute similarity based on length with no homozygous break
    similarity = 1/(sqrt(pos.len^2+neg.len^2))
    return(similarity)
}


# The following four distance measures are given in: https://www.biorxiv.org/content/biorxiv/early/2014/11/28/004184.full.pdf
matchstates <- function(x,y){
    intersect_locus <- 2-abs(x-y)
    nucleotides_locus <- 4
    dist = 1 - 2 * (intersect_locus / nucleotides_locus)
    # combine
    combined_dist = sum(dist) / length(x)
    return(combined_dist)
}

genpofad <- function(x,y){
    intersect_locus <- 2-abs(x-y)
    nucleotides_locus <- 2 # max for one individual
    dist = 1 - intersect_locus / nucleotides_locus
    # combine
    combined_dist = sum(dist) / length(x)
    return(combined_dist)
}

mrca <- function(x,y){
    intersect_locus <- 2-abs(x-y)
    dist = rep(0, length(intersect_locus))
    dist[intersect_locus==0] = 1
    # combine
    combined_dist = sum(dist) / length(x)
    return(combined_dist)
}

nei <- function(x,y){
    freq_x = rep(0, length(x))
    freq_x[x==0 | x == 2] = 1
    freq_x[x==1] = 0.5

    freq_y = rep(0, length(y))
    freq_y[y==0 | y == 2] = 1
    freq_y[y==1] = 0.5

    dist = 1 - sqrt(freq_x * freq_y)
    # combine
    combined_dist = sum(dist) / length(x)
    return(combined_dist)
}

computeCustomDist <- function(distMethod = firstUnderBreakDist, fam=F, filename="", project_path=""){

    individualDist <- function(filename="", project_path=""){
        #snp_names <<- matched_plink$SNP
        #famCode_labels <<- filter(brca_data, Onc_ID %in% matched_plink$SNP) %>% select(FamCode) %>% .[,1]
        #rownames(matched_plink) <<- matched_plink$SNP
        brca_data.filtered <<- brca_data.single
        matched <<- matched_single
        matched_plink <<- matched_plink_single

        # Define snp distances to BRCA mutation
        chr_pos.dist <<- mapSNPsToCoords(chr_coords, matched_plink)
        # Define position splitting snps on either side of BRCA mutation
        chr_pos.split <<- nrow(filter(chr_pos.dist, brca_distance < 0))
        chr_pos.split.left <<- nrow(filter(chr_pos.dist, brca_distance <= 0))

        if (file.exists(filename)){
            dst <- readRDS(paste0(project_path, filename))
        } else {
            dst = proxy::dist(matched_plink[,2:ncol(matched_plink)], method=distMethod)
            if (nchar(filename) > 0){
                saveRDS(dst, file = paste0(project_path, filename))
            }
        }

        return(dst)
    }

    famDist <- function(filename="", project_path=""){
        #snp_names <<- matched_plink_fam$SNP
        #famCode_labels <<- matched_plink_fam$FamCode
        brca_data.filtered <<- brca_data.fam
        # Remove famCode column
        #matched <<- matched_fam[,2:ncol(matched_fam)]
        matched_plink <<- matched_plink_fam[,2:ncol(matched_plink_fam)]

        # Define snp distances to BRCA mutation
        chr_pos.dist <<- mapSNPsToCoords(chr_coords, matched_plink_fam)
        # Removes SNPs in brca gene, as for now we're not sure which side they belong to
        if (removeSNPsInGene){
            chr_pos.dist <<- filter(chr_pos.dist, position_b37 < brca_start | position_b37 > brca_stop)
        }
        # Define position splitting snps on either side of BRCA mutation
        chr_pos.split <<- nrow(filter(chr_pos.dist, brca_distance < 0))

        if (file.exists(filename)){
            dst <- readRDS(paste0(project_path, filename))
        } else {
            dst = proxy::dist(matched_plink_fam[,3:ncol(matched_plink_fam)], method=distMethod)
            if (nchar(filename) > 0){
                saveRDS(dst, file = paste0(project_path, filename))
            }
        }

        return(dst)
    }

    if(fam==T){
        dst <- famDist(filename, project_path)
    } else { # fam == F
        dst <- individualDist(filename, project_path)
    }
    return(dst)
}



###############################################################################
### Define similarity based on length of similar elements around brca
###############################################################################
mapSNPsToCoords <- function(chr_coords, matched){
    chr_pos = filter(chr_coords, SNP %in% names(matched)) %>%
        arrange(position_b37) %>%
        mutate(brca_distance = as.integer(as.character(position_b37)) - brca_middle) %>%
        mutate(brca_distance_cM = cM - brca_cM_middle)
    return(chr_pos)
}

firstBreakDist <- function(matched_custom=NULL, fam=F, filename="", project_path=""){

    distanceSimilarity <- function(x,y){
        #x = matched_plink[1,2:ncol(matched_plink)]
        #y = matched_plink[2,2:ncol(matched_plink)]
        # Find homozygote break
        #condition = (x-y == -2 | x-y == 2)
        condition = (x-y == -2 | x-y == 2 | x == -1 | y == -1)
        # Split vector according to positive or negative distance to brca gene
        cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
        cond_neg = condition[1:chr_pos.split.left]
        # get length (in Kb) to brca gene for first homozygous break on either side
        min_index=min(which(cond_pos==T), length(cond_pos))
        #print(which(cond_pos==T))
        pos.len = chr_pos.dist[chr_pos.split+min_index,"brca_distance"] / 1000000
        #pos.len = chr_pos.dist[chr_pos.split+min_index,"brca_distance_cM"] * 100
        max_index=max(which(cond_neg==T), 1)
        #print(which(cond_neg==T))
        neg.len = abs(chr_pos.dist[max_index,"brca_distance"]) / 1000000
        #neg.len = abs(chr_pos.dist[max_index,"brca_distance_cM"]) * 100
        # Compute similarity based on length with no homozygous break
        similarity = 1/(pos.len+neg.len)
        return(similarity)
    }

    individualDist <- function(filename="", project_path=""){
        #snp_names <<- matched_plink$SNP
        #famCode_labels <<- filter(brca_data, Onc_ID %in% matched_plink$SNP) %>% select(FamCode) %>% .[,1]
        #rownames(matched_plink) <<- matched_plink$SNP
        brca_data.filtered <<- brca_data.single
        matched <<- matched_single
        matched_plink <<- matched_plink_single

        # Define snp distances to BRCA mutation
        chr_pos.dist <<- mapSNPsToCoords(chr_coords, matched_plink)
        # Removes SNPs in brca gene, as for now we're not sure which side they belong to
        if (removeSNPsInGene){
            chr_pos.dist <<- filter(chr_pos.dist, position_b37 < brca_start | position_b37 > brca_stop)
        }
        # Define position splitting snps on either side of BRCA mutation
        chr_pos.split.left <<- nrow(filter(chr_pos.dist, brca_distance <= 0))
        chr_pos.split <<- nrow(filter(chr_pos.dist, brca_distance < 0))

        if (file.exists(paste0(project_path, filename))){
            dst <- readRDS(paste0(project_path, filename))
        } else {
            dst = proxy::dist(matched_plink[,2:ncol(matched_plink)], method=distanceSimilarity)
            if (nchar(filename) > 0){
                saveRDS(dst, file = paste0(project_path, filename))
            }
        }

        return(dst)
    }

    famDist <- function(filename="", project_path=""){
        #snp_names <<- matched_plink_fam$SNP
        #famCode_labels <<- matched_plink_fam$FamCode
        brca_data.filtered <<- brca_data.fam
        # Remove famCode column
        #matched <<- matched_fam[,2:ncol(matched_fam)]
        matched_plink <<- matched_plink_fam[,2:ncol(matched_plink_fam)]

        # Define snp distances to BRCA mutation
        chr_pos.dist <<- mapSNPsToCoords(chr_coords, matched_plink_fam)
        # Removes SNPs in brca gene, as for now we're not sure which side they belong to
        if (removeSNPsInGene){
            chr_pos.dist <<- filter(chr_pos.dist, position_b37 < brca_start | position_b37 > brca_stop)
        }
        # Define position splitting snps on either side of BRCA mutation
        chr_pos.split.left <<- nrow(filter(chr_pos.dist, brca_distance <= 0))
        chr_pos.split <<- nrow(filter(chr_pos.dist, brca_distance < 0))

        print(filename)
        print(paste0(project_path, filename))
        
        if (file.exists(paste0(project_path, filename))){
            dst <- readRDS(paste0(project_path, filename))
        } else {
            dst = proxy::dist(matched_plink_fam[,3:ncol(matched_plink_fam)], method=distanceSimilarity)
            if (nchar(filename) > 0){
                saveRDS(dst, file = paste0(project_path, filename))
            }
        }

        return(dst)
    }

    customDist <- function(matched_custom, filename=""){
        #brca_data.filtered <<- brca_data
        chr_pos.dist <<- mapSNPsToCoords(chr_coords, matched_custom)
        # Removes SNPs in brca gene, as for now we're not sure which side they belong to
        if (removeSNPsInGene){
            chr_pos.dist <<- filter(chr_pos.dist, position_b37 < brca_start | position_b37 > brca_stop)
        }
        chr_pos.split.left <<- nrow(filter(chr_pos.dist, brca_distance <= 0))
        chr_pos.split <<- nrow(filter(chr_pos.dist, brca_distance < 0))
        if (file.exists(filename)){
            dst <- readRDS(filename)
        } else {
            dst = proxy::dist(matched_custom[,2:ncol(matched_custom)], method=distanceSimilarity)
            # length(dst)
            # if (fam){
            #     dst = proxy::dist(matched_custom[,3:ncol(matched_custom)], method=distanceSimilarity)
            # } else {
            #     dst = proxy::dist(matched_custom[,2:ncol(matched_custom)], method=distanceSimilarity)
            # }
            if (nchar(filename) > 0){
                saveRDS(dst, file = filename)
            }
        }
        return(dst)
    }

    if (!is.null(matched_custom)){
        dst <- customDist(matched_custom, filename)
    } else if(fam==T){
        dst <- famDist(filename, project_path)
    } else { # fam == F
        dst <- individualDist(filename, project_path)
    }

    # hc <- hclust(dst, "ward.D2")
    # hc <- hclust(dst, "ward.D")
    # hc <- hclust(dst, "complete")

    # hiearchical_clustering(hc, k=2, F)
    # hiearchical_clustering(hc, k=3, F)
    # hiearchical_clustering(hc, k=4, F)
    # hiearchical_clustering(hc, k=5, F)
    # hiearchical_clustering(hc, k=6, F)
    # hiearchical_clustering(hc, k=7, F)
    # hiearchical_clustering(hc, k=4, T)
    return(dst)
}

runHClust <- function(mut = "c.7617+1G>A", gene = "BRCA2"){
    mut = "c.7617+1G>A"; gene = "BRCA2"
    #mut = "c.5946delT"; gene = "BRCA2"
    #mut = "c.3264dupT"; gene = "BRCA2"
    #mut = "c.8537_8538delAG"; gene = "BRCA2"
    #mut = "c.5722_5723delCT"; gene = "BRCA2"
    # mut = "c.-200-?_80+?del"; gene = "BRCA1"
    mut = "c.68_69delAG"; gene = "BRCA1"
    #mut = "c.4327C>T"; gene = "BRCA1"
    #mut = "c.211A>G"; gene = "BRCA1"
    # mut = "c.2681_2682delAA"; gene = "BRCA1"
    #mut = "mutSim_BRCA1_100_1"; gene = "BRCA1_simulated"
    #mut = "mutSim_BRCA2_100_1"; gene = "BRCA2_simulated"
    mut <<- mut; gene <<- gene
    readIndvAndFamData(gene, mut)
    #dst <- euclideanDist_BRCAdistWeighted()
    # dst <- firstBreakDist(fam=F)#, filename)
    # dst <- computeCustomDist(firstUpperBreakDist, fam=F)
    # dst <- computeCustomDist(firstUnderBreakDist, fam=F)
    # dst <- computeCustomDist(matchstates, fam=F)
    dst <- firstBreakDist(fam=T)
    k = 2
    #hc <- hiearchical_clustering(dst, clustMethod = "complete", k = k, doHapPlot = F)
    # hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = k, doHapPlot = F) # WRONG
    hc <- hiearchical_clustering(dst, clustMethod = "ward.D2", k = k, doHapPlot = F)
    #hc <- hiearchical_clustering(dst, clustMethod = "complete", k = 2, doHapPlot = T)
    plot_dendrogram(hc, k = k)
    #plot_dendrogram(hc, k = k, printLeaf = 1)
    group = 2
    chr_pos = haplotype_mutation(mut, prepare = F, group = group, cutoff = 0.5, ancestral_method = "branchBoundIndep", min_samples = 0.25)
    chr_pos = haplotype_mutation(mut, prepare = F, group = group, cutoff = 0.5, ancestral_method = "branchBoundIndep", min_samples = 5)
    chr_pos = haplotype_mutation(mut, prepare = F, group = group, cutoff = 0.5, ancestral_method = "branchBound")
    chr_pos = haplotype_mutation(mut, prepare = F, group = group, cutoff = 0.5, ancestral_method = "mostFreqBase")
    plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = group)
    #plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = group, col="position_b37")
    #plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = group, col="centimorgan")
    plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
    p = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group=group)[[1]]
    # p = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group=group)
    # q = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group=group)
    # p2 = ggarrange(p[[1]], q[[1]], ncol = 2, nrow = 1)
    # ggsave(filename = "BRCA2-c.7617+1G_A-SPAIN-ancestral.png", plot = p2, scale = 1, width = 16, height = 9)
    system.time(computeMutationAge(gene = gene, mut = mut, isFam = T, method = 5))
}

runHClust_simulated <- function(mut = "mutSim_BRCA1_1", gene = "BRCA1_simulated"){
    simulations = 20; generations = 100
    brca1_simulated_pheno_merged <<- read.table(file = paste0("Scripts/simulate-population/BRCA1_simulated_starGenealogy_200-samples_",simulations,"-muts_",generations,"-generations_pheno.txt"), header = T)
    brca1_simulated_geno_plink <<- read.table2(file = paste0("Scripts/simulate-population/BRCA1_simulated_starGenealogy_200-samples_",simulations,"-muts_",generations,"-generations_geno.txt"), header = T)

    mut = "mutSim_BRCA1_6"; gene = "BRCA1_simulated"
    mut <<- mut; gene <<- gene
    readIndvAndFamData(gene, mut)
    dst <- firstBreakDist(fam=F)#, filename)
    #dst <- firstBreakDist(fam=T)
    k = 1
    #hc <- hiearchical_clustering(dst, clustMethod = "complete", k = k, doHapPlot = F)
    hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = k, doHapPlot = F)
    plot_dendrogram(hc, k = k)
    group = 1
    chr_pos = haplotype_mutation(mut, prepare = F, group = group, cutoff = 0.5)
    plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = group)
    #plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
    p = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group=group)
    p[[1]]
    p[[4]]
}

runHClust_mark <- function(){
    loadProject("LPR5")
    mut = project_mut_info$Mut1HGVS
    gene
    prepare_dataframe(mut, gene, gene_info_obj = project_mut_info)
    dim(matched_plink)
    dim(open_projects$LRP5$project_geno)
    dim(brca_data)
    dim(chr_coords)
    dst <- firstBreakDist(fam=F)
    k = 3
    hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = k, doHapPlot = F)
    plot_dendrogram(hc, k = k)
    # for (group in 1:k){
    #     chr_pos = haplotype_mutation(mut, prepare = F, group = group, cutoff = 0.5, ancestral_method = "branchBoundIndep")
    #     p = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group=group)[[1]]
    # }
    chr_pos = haplotype_mutation(mut, prepare = F, group = 1, cutoff = 0.5, ancestral_method = "mostFreqBase")
    chr_pos = haplotype_mutation(mut, prepare = F, group = 3, cutoff = 0.5, ancestral_method = "branchBoundIndep")
    plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = group)
    p = add_second_legend(plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group=group)[[1]], chr_pos)
    p = add_second_legend(p, chr_pos)
    
    # Country grouping
    chr_pos = haplotype_mutation(mut, prepare = F, country = "case", cutoff = 0.5, ancestral_method = "branchBoundIndep")
    chr_pos = haplotype_mutation(mut, prepare = F, country = "case", cutoff = 0.5, ancestral_method = "mostFreqBase")
    plot_breakpoints_country(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, otherCountry = F)
    p = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, otherCountry = F)[[1]]
    
}

# compare_dendrograms()
compare_dendrograms <- function(dst, dst2, clustMethod1, clustMethod2){
    # Compute distance matrix
    #dst <- firstBreakDist(fam=F)
    # Compute 2 hierarchical clusterings
    #hc1 <- hclust(dst, method = "ward.D")
    #hc2 <- hclust(dst, method = "ward.D2")
    hc1 <- hclust(dst, method = clustMethod1)
    hc2 <- hclust(dst2, method = clustMethod2)

    # Create two dendrograms
    cols <<- colors_plot(pheno_data = brca_data.filtered, palette = "Dark2")

    dend1 <- as.dendrogram(hc1)
    # FamCode <- brca_data.filtered$FamCode[order.dendrogram(dend1)]
    # Country_ordered <- brca_data.filtered$Country[order.dendrogram(dend1)]
    # country_cols <- cols[as.character(Country_ordered)]
    # dend1 <- dend1 %>% set("branches_k_color", k = k) %>%
    #     set("labels", FamCode) %>% set("labels_colors", country_cols) %>%
    #     set("labels_cex", 1) %>% set("branches_lwd", 3)

    dend2 <- as.dendrogram(hc2)
    # FamCode <- brca_data.filtered$FamCode[order.dendrogram(dend2)]
    # Country_ordered <- brca_data.filtered$Country[order.dendrogram(dend2)]
    # country_cols <- cols[as.character(Country_ordered)]
    # dend2 <- dend2 %>% set("branches_k_color", k = k) %>%
    #     set("labels", FamCode) %>% set("labels_colors", country_cols) %>%
    #     set("labels_cex", 1) %>% set("branches_lwd", 3)

    # Create a list to hold dendrograms
    dend_list <- dendlist(dend1, dend2)

    tanglegram(dend1, dend2, sort = T,
       highlight_distinct_edges = TRUE, # Turn-off dashed lines
       common_subtrees_color_lines = TRUE, # Turn-off line colors
       common_subtrees_color_branches = TRUE, # Color common branches
       main = paste("entanglement =", round(entanglement(dend_list), 2)))
}

kmeans_clust <- function(){
    km = kmeans(dst, 2, 100)
    brca_data.filtered$cluster_groups = km$cluster
    # Print group information
    for (i in unique(brca_data.filtered$cluster_groups)){
        print(subset(brca_data.filtered, cluster_groups == i) %>% select(FamCode, Onc_ID, Country, cluster_groups))
    }
    brca_data <<- brca_data.filtered
    for (i in unique(brca_data.filtered$cluster_groups)){
        if (nrow(subset(brca_data.filtered, cluster_groups == i)) > 5){
            chr_pos = haplotype_mutation(mut, prepare = F, group = i)
            #plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = i)
            plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, i)
        }
    }
}

numCluster_estimation <- function(max.clust = 12,
                                  distMethod = "firstBreakDist",
                                  clustMethod = "ward.D",
                                  fam = "fam",
                                  dist = dst,
                                  req_members = 3,
                                  index = c("mcclain", "cindex", "silhouette", "dunn")){

    print(index)
    print(length(dist))

    filename = paste0(project_path, "/ClusterValidationValues/", gene, "-", mut_name, "-ClusterValidationPlot-", length(index),
                      "_scores-distMethod-", distMethod, "-clustMethod-", clustMethod, "-", fam, ".txt")
    # Check if data already cached, if so, return it.
    if (file.exists(filename)){
        num_clust <- read.table(filename, header = T)[,2]
        numCluster <<- table(num_clust)
        numCluster_optimal = max(as.integer(names(which.max(numCluster))), 1)
        return(numCluster_optimal)
    }

    # dist=dst
    # No need for data matrix
    num_clust = rep(0, length(index))
    numCluster_plots = list()
    for (i in 1:length(index)){

        #matched_pli = matched_plink

        print(index[i])
        # if (index[i] == "silhouette"){
        #     v <- rep(0, max.clust)
        #     for(j in 2:max.clust){
        #         hc <- hiearchical_clustering(dist, clustMethod = clustMethod, k = j, doHapPlot = F)
        #         v[j] <- mean(cluster::silhouette(cluster_groups, dist, FUN = mean)[,3])
        #     }
        #     num_clust[i] = which.max(v)
        #
        #     # Set up data frame for plotting
        #     df = data.frame(clusters = 1:max.clust, y = v)
        # } else
        if (index[i] %in% c("mcclain", "cindex", "silhouette", "dunn")) {
            nb=NbClust(data = NULL, distance = NULL, diss = dist, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
            #nb=NbClust(matched_plink[,(2+1310):(ncol(matched_plink)-792)], distance = NULL, diss = dst, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
            print(nb)
            num_clust[i] = nb$Best.nc[1]

            # Set up data frame for plotting
            df = data.frame(clusters = 1:max.clust, y = c(0, nb$All.index))
        } else {
            #nb=NbClust(matched_plink[,2:ncol(matched_plink)], distance = NULL, diss = dist, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
            if (gene == "BRCA2"){
                # Starting around 5550000 bases from brca gene
                nb=NbClust(matched_plink[,(2+1310):(ncol(matched_plink)-792)], distance = NULL, diss = dist, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
                #nb=NbClust(matched_plink[,2295:2305], distance = NULL, diss = dst, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
            } else { # BRCA1
                # Starting around 5550000 bases from brca gene
                nb=NbClust(matched_plink[,(2+1170):(ncol(matched_plink)-1990)], distance = NULL, diss = dist, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
            }
            print(nb)

            # might be wrong to just find max value - should probably be looking for knee graphically
            if (index[i] %in% c("dindex", "hubert")){
                num_clust[i] = as.integer(names(which.max(nb$All.index)))
            } else {
                num_clust[i] = nb$Best.nc[1]
            }

            # Set up data frame for plotting
            df = data.frame(clusters = 1:max.clust, y = c(0, nb$All.index))
        }

        plotname = paste0(project_path, "/ClusterValidationPlots/", gene, "-", mut_name, "-ClusterValidationPlot-", index[i],
                          "_score-distMethod-", distMethod, "-clustMethod-", clustMethod, "-", fam, ".png")

        if (!file.exists(plotname)){
            linecolor = "steelblue"
            p <- ggpubr::ggline(df, x = "clusters", y = "y", group = 1,
                                color = linecolor, ylab = paste("Average", index[i], "score"), xlab = "Number of clusters k",
                                main = paste("Optimal number of clusters -", index[i], "score")) +
                scale_x_continuous(breaks = c(1:max.clust)) +
                geom_vline(xintercept = num_clust[i], linetype = 2, color = linecolor)
                #geom_vline(xintercept = which.max(v), linetype = 2, color = linecolor)
            #print(p)
            cur_dev <- dev.cur()
            ggsave(filename = plotname, plot = p, scale = 1, width = 17.5, height = 5.714)
            dev.set(cur_dev)
            #numCluster_plots[[i]] = p
        }
    }
    # cache data, for faster retrieval
    write.table(x = cbind(index, num_clust), file = filename, quote = F, sep = "\t", row.names = F)

    print(table(num_clust))
    numCluster <<- table(num_clust)
    numCluster_optimal = max(as.integer(names(which.max(table(num_clust)))), 1)
    #numCluster_plots <<- numCluster_plots

    # filter(brca_data.single, Onc_ID %in% matched_plink$SNP)$Country
    # filter(brca_data.single, Onc_ID %in% matched_plink$SNP)$Country[nb$Best.partition==1] %>% unique()
    # filter(brca_data.single, Onc_ID %in% matched_plink$SNP)$Country[nb$Best.partition==2] %>% unique()
    # filter(brca_data.single, Onc_ID %in% matched_plink$SNP)$Country[nb$Best.partition==3] %>% unique()
    # filter(brca_data.single, Onc_ID %in% matched_plink$SNP)$Country[nb$Best.partition==4] %>% unique()
    #return(list(numCluster_optimal, numCluster_plots))
    return(numCluster_optimal)
}

numCluster_estimation_extended <- function(max.clust = 4, clustMethod = "complete"){
    mut <<- "c.7617+1G>A"; gene <<- "BRCA2"
    readIndvAndFamData(gene, mut)
    #dst <- firstBreakDist(fam=F)
    dst <- firstBreakDist(fam=T)

    # With data matrix
    # index = "all" # Does not work, since one of the following does not work
    index = c("kl", "ch", "hartigan", "ball", "ptbiserial", "gap", "gamma", "gplus", "tau", "sdindex", "sdbw", "db", "pseudot2", "hubert", "dindex")
    #index = c("kl")

    num_clust = rep(0, length(index))
    for (i in 1:length(index)){

        print(index[i])
        if (gene == "BRCA2"){
            # Starting around 5550000 bases from brca gene
            nb=NbClust(matched_plink[,(2+1310):(ncol(matched_plink)-792)], distance = NULL, diss = dst, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
            #nb=NbClust(matched_plink[,2295:2305], distance = NULL, diss = dst, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
        } else { # BRCA1
            # Starting around 5550000 bases from brca gene
            nb=NbClust(matched_plink[,(2+1170):(ncol(matched_plink)-1990)], distance = NULL, diss = dst, min.nc=2, max.nc=max.clust, method = clustMethod, index = index[i])
        }
        print(nb)

        if (index[i] %in% c("dindex", "hubert")){
            num_clust[i] = as.integer(names(which.max(nb$All.index)))
        } else {
            num_clust[i] = nb$Best.nc[1]
        }

        linecolor = "steelblue"
        df = data.frame(clusters = 1:max.clust, y = c(0, nb$All.index))
        p <- ggpubr::ggline(df, x = "clusters", y = "y", group = 1,
                            color = linecolor, ylab = paste("Average", index[i], "score"), xlab = "Number of clusters k",
                            main = "Optimal number of clusters") + scale_x_continuous(breaks = c(1:max.clust)) +
             geom_vline(xintercept = num_clust[i], linetype = 2, color = linecolor)
        print(p)
    }

    print(table(num_clust))
    # df = data.frame(clusters = 1:max.clust, y = num_clust)
    # p <- ggpubr::ggline(df, x = "clusters", y = "y", group = 1,
    #                     color = linecolor, ylab = ylab, xlab = "Number of clusters k",
    #                     main = "Optimal number of clusters")

    # Metrics
    # Success: "kl", "ch", "hartigan", "ball", "ptbiserial", "gap", "gamma", "gplus", "tau", "sdindex", "sdbw", "db", "duda", "pseudot2"
    # Failed: "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "hubert", "dindex", "frey", "beale"
    # Maybe: "ratkowsky" (kinda works, but suggests unlimited clusters - further testing needed),
}

gap_statistic <- function(){
    #set.seed(123)
    fviz_nbclust(matched_plink, hcut, method = "gap_stat", diss = dst, nboot = 100) + labs(subtitle = "Gap statistic method")
}

pvclust_compute <- function(){
    customDist <- function(matched_plink){
        distanceSimilarity <- function(x,y){
            #x = matched_plink[1,2:ncol(matched_plink)]
            #y = matched_plink[2,2:ncol(matched_plink)]
            # Find homozygote break
            condition = (x-y == -2 | x-y == 2)
            # Split vector according to positive or negative distance to brca gene
            cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
            cond_neg = condition[1:chr_pos.split]
            # get length (in Kb) to brca gene for first homozygous break on either side
            min_index=min(which(cond_pos==T), length(cond_pos))
            #print(which(cond_pos==T))
            pos.len = chr_pos.dist[chr_pos.split+min_index,4] / 1000
            max_index=max(which(cond_neg==T), 1)
            #print(which(cond_neg==T))
            neg.len = abs(chr_pos.dist[max_index,4]) / 1000
            # Compute similarity based on length with no homozygous break
            similarity = 1/(pos.len+neg.len)
            return(similarity)
        }

        dst = proxy::dist(t(matched_plink[,-1]), method=distanceSimilarity)
        return(dst)
    }
    chr_pos.dist <<- mapSNPsToCoords(chr_coords, matched_plink)
    chr_pos.split <<- nrow(filter(chr_pos.dist, brca_distance < 0))

    library(pvclust)
    library(parallel)
    # Non parallel
    pv = pvclust(matched_plink[,-1], method.hclust = "complete", method.dist = customDist, nboot = 10)

    # Parallel
    cl <- makeCluster(detectCores()-1)
    clusterExport(cl, list("chr_pos.dist", "chr_pos.split", "matched_plink"))
    pv = pvclust(t(matched_plink[,-1]), method.hclust = "complete", method.dist = customDist, nboot = 100, parallel = cl)
    stopCluster(cl)

    # Print AU and BP scores
    print(pv)
    plot(pv, hang = -1)
    pvrect(pv, alpha = 0.95)
}

Recluster_subgroups <- function(group, k){
    ### NOTE: This works. But essentially it's the same cluster as the big one, just "zoomed in"
    group = 2
    # Backup original brca_data
    brca_data_bak <- brca_data
    brca_data.filtered <<- subset(brca_data, cluster_groups == group)
    matched_plink2 = filter(matched_plink, SNP %in% brca_data.filtered$Onc_ID)
    dst <- firstBreakDist(matched_custom = matched_plink2)
    k = 6
    #hc <- hiearchical_clustering(dst, clustMethod = "complete", k = k, doHapPlot = F)
    hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = k, doHapPlot = F)

    # Restore brca_data object
    brca_data <<- brca_data_bak

    plot_dendrogram(hc, k = k)

    # brca_data <<- brca_data.filtered
    # subgroup = 1
    # chr_pos = haplotype_mutation(mut, prepare = F, group = subgroup, cutoff = 0.5)
    # plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = subgroup)
    # plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, subgroup)
}

#test_dbscan(8,0.2)
#test_dbscan(10,40)
test_dbscan <- function(minPts, eps){
    set.seed(123)
    minPts <<- minPts; eps <<- eps
    db <- fpc::dbscan(data = dst, eps = eps, MinPts = minPts, method = "dist")

    kNNdistplot(x = dst, k = minPts)

    dbscan_plot(db)

    # brca_data.filtered$cluster_groups = db$cluster
    # # Print group information
    # for (i in unique(brca_data.filtered$cluster_groups)){
    #     print(subset(brca_data.filtered, cluster_groups == i) %>% select(FamCode, Onc_ID, Country, cluster_groups))
    # }
    # brca_data <<- brca_data.filtered
    # for (i in unique(brca_data.filtered$cluster_groups)){
    #     if (nrow(subset(brca_data.filtered, cluster_groups == i)) > 5){
    #         chr_pos <<- haplotype_mutation(mut, prepare = F, group = i)
    #         print(plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = i))
    #         plot(plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, i)[[4]])
    #     }
    # }

    return(db)
}

dbscan_plot <- function(db){
    p = fviz_cluster(db, data = matched_plink[,1312:(ncol(matched_plink)-792)], stand = FALSE,
                     #fviz_cluster(db, data = matched_plink[,1172:(ncol(matched_plink)-1990)], stand = FALSE,
                     ellipse = FALSE, show.clust.cent = FALSE, pointsize = 2.5,
                     geom = "point", ggtheme = theme_classic())
    #print(p)
    return(p)
}


hdbscan_plot <- function (hdb_object = hdb, data = dst)
{
    db = dbscan::dbscan(x = data, eps = 0.2, minPts = 5)
    db$cluster = hdb_object$cluster
    dbscan_plot(db)
}

debug_fam_vs_individual <- function(){
    # Figuring stuff out
    # Group 1
    filter(combined, FamCode == "CIMF027658")[2690:2710]
    filter(matched_fam, SNP %in% c("149341"))[2690:2710]
    filter(matched, SNP %in% c("149341", "149346", "149351", "149388", "149398", "149399"))[,2690:2710]

    # Group 2
    filter(combined, FamCode == "CIMF027815")[1:10]
    filter(matched_fam, SNP %in% c("149622"))[2690:2710]
    filter(matched, SNP %in% c("149622", "149653", "149655", "149658", "149685"))[,2690:2710]
    #filter(brca_data, FamCode == "CIMF027803") # 3
    #filter(brca_data, FamCode == "CIMF027815") # 5
    filter(brca_data, Onc_ID == "146809") %>% select(FamCode)
    filter(brca_data, FamCode == "CIMF018512")
    filter(matched, SNP %in% c("146809", "146810", "146811", "146812"))[(chr_pos.split+370):(chr_pos.split+374)]
    filter(matched, SNP %in% c("146809"))[(chr_pos.split+60):(chr_pos.split+70)]
    filter(matched_plink, SNP %in% c("146809"))[(chr_pos.split+60):(chr_pos.split+70)]
    filter(matched_plink_fam, SNP %in% c("146809"))[(chr_pos.split+370):(chr_pos.split+374)]

    filter(brca_data, Onc_ID == "146334") %>% select(FamCode)
    filter(brca_data, FamCode == "CIMF002010")
    filter(matched, SNP %in% c("146334","146346"))[(chr_pos.split+370):(chr_pos.split+374)]
    filter(matched_plink, SNP %in% c("146334"))[(chr_pos.split+60):(chr_pos.split+70)]
    filter(matched_plink_fam, SNP %in% c("146334"))[(chr_pos.split+370):(chr_pos.split+374)]



    filter(brca_data, Onc_ID == "259467")
    filter(brca_data, Onc_ID == "146337")
    filter(matched_plink_fam, SNP == "259467")[1:10]
    filter(matched_plink_fam, SNP == "146337")[1:10]

    # x=filter(matched_plink_fam, SNP %in% c("149622"))[2:4600]
    # y=filter(matched_plink_fam, SNP %in% c("149341"))[2:4600]
    # x=filter(matched_plink, SNP %in% c("149622"))[2:4600]
    # y=filter(matched_plink, SNP %in% c("149341"))[2:4600]

    # x=filter(matched_plink_fam, SNP %in% c("149566"))[2:4600]
    # y=filter(matched_plink_fam, SNP %in% c("305697"))[2:4600]
    # x=filter(matched_plink, SNP %in% c("149566"))[2:4600]
    # y=filter(matched_plink, SNP %in% c("305697"))[2:4600]

    x=filter(matched_plink_fam, SNP %in% c("146334"))[2:4600]
    y=filter(matched_plink_fam, SNP %in% c("146809"))[2:4600]
    x=filter(matched_plink, SNP %in% c("146334"))[2:4600]
    y=filter(matched_plink, SNP %in% c("146809"))[2:4600]


    condition = (x-y == -2 | x-y == 2)
    # Split vector according to positive or negative distance to brca gene
    cond_pos = condition[(chr_pos.split+1):nrow(chr_pos.dist)]
    cond_neg = condition[1:chr_pos.split]
    # get length (in Kb) to brca gene for first homozygous break on either side
    min_index=min(which(cond_pos==T), length(cond_pos))
    #print(which(cond_pos==T))
    pos.len = chr_pos.dist[chr_pos.split+min_index,4] / 1000000
    max_index=max(which(cond_neg==T), 1)
    #print(which(cond_neg==T))
    neg.len = abs(chr_pos.dist[max_index,4]) / 1000000
    print(pos.len)
    print(neg.len)
    1/(pos.len+neg.len)

    test_matched_plink <- filter(matched_plink, SNP %in% matched_plink_fam$SNP) %>% arrange(SNP)
    chr_pos.dist <- mapSNPsToCoords(chr_coords, test_matched_plink)
    chr_pos.split = nrow(filter(chr_pos.dist, brca_distance < 0))
    dst1 = proxy::dist(test_matched_plink[,2:ncol(test_matched_plink)], method=distanceSimilarity)

    test_matched_plink_fam <- arrange(matched_plink_fam, SNP)
    chr_pos.dist <- mapSNPsToCoords(chr_coords, test_matched_plink_fam)
    chr_pos.split = nrow(filter(chr_pos.dist, brca_distance < 0))
    dst2 = proxy::dist(test_matched_plink_fam[,2:ncol(test_matched_plink_fam)], method=distanceSimilarity)

    which((dst1==dst2)==F)[30:50]
    dst1[which((dst1==dst2)==F)][30:50]
    dst2[which((dst1==dst2)==F)][30:50]
}

plot_dendrogram_alternatives <- function(){
    ### ggplot2 and ggdendro for dendrogram
    hcd <- dendro_data(as.dendrogram(hc))
    FamCode <- brca_data.filtered[fmatch(label(hcd)$label, brca_data.filtered$Onc_ID),"FamCode"]
    Country <- brca_data.filtered[fmatch(label(hcd)$label, brca_data.filtered$Onc_ID),"Country"]
    seq2 <- filter(segment(hcd), x > 60)
    p <- ggplot(segment(hcd)) +
        geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
        geom_text(data=label(hcd), aes(label=FamCode, x=x, y=0, colour=Country, angle=90, vjust = 0.5, hjust = 1.1), size=3) +
        ylim(-2, max(segment(hcd)$yend)) +
        scale_color_brewer(palette="Set1")
    print(p)


    # hcd %>% set("labels", FamCode) %>% set("labels_col", c("blue")) %>% plot

    # load code of A2R function
    # source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
    # # Changing labels
    # hcd <- as.dendrogram(hc)
    # famLab <- function(n){
    #     if (is.leaf(n)) {
    #
    #     }
    # }
    # FamCode <- brca_data[fmatch(label(hcd)$label, brca_data$Onc_ID),"FamCode"]
    # label(hcd)$label <- FamCode
    # # colored dendrogram
    # op = par(bg = "#EFEFEF")
    # A2Rplot(hc, k = 2, boxes = F, col.up = "gray50", col.down = c("#FF6B6B", "#4ECDC4", "#556270", "green"))


    # hcd <- as.dendrogram(hc)
    # labelColors = c("#CDB380", "#036564", "#EB6841", "#EDC951")
    # clusMember = cutree(hc, 2)
    # # function to get color labels
    # colLab <- function(n) {
    #     print(n)
    #     if (is.leaf(n)) {
    #         a <- attributes(n)
    #         labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    #         attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    #     }
    #     n
    # }
    # clusDendro = dendrapply(hcd, colLab)
    # # make plot
    # plot(clusDendro, main = "Cool Dendrogram")
}

test_simulated_data <- function(){
    mut = "mutSim_BRCA1_100_1"; gene = "BRCA1_simulated"
    gen = strsplit(mut, "_")[[1]][3]
    #filename = paste0("classify-simulated-population-age/training-BRCA1_simulated-starGenealogy-100_samples-1_simulations-generations_100_100_step10-seed_42/simulated_population/training-BRCA1_simulated-starGenealogy-100_samples-1_simulations-generations_100-seed_42")
    filename = paste0("classify-simulated-population-age/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_100_100_step10-seed_9012/simulated_population/BRCA1_simulated-starGenealogy-100_samples-20_simulations-generations_",gen,"-seed_9012")
    brca1_simulated_pheno_merged <<- read.table(file = paste0(filename, "-pheno.txt"), header = T)
    brca1_simulated_geno_plink <<- read.table(file = paste0(filename, "-geno.txt"), header = T)
    mut <<- mut; gene <<- gene
    readIndvAndFamData(gene, mut)
    dst <- firstBreakDist(fam=F)#, filename)
    #dst <- firstBreakDist(fam=T)
    k = 1
    #hc <- hiearchical_clustering(dst, clustMethod = "complete", k = k, doHapPlot = F)
    hc <- hiearchical_clustering(dst, clustMethod = "ward.D", k = k, doHapPlot = F)
    plot_dendrogram(hc, k = k)
    group = 1
    chr_pos = haplotype_mutation(mut, prepare = F, group = group, cutoff = 0.5)
    plot_breakpoints_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group = group, col = "genetic_dist")
    #plot_nearest_brca_break_group(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group)
    chr_pos2 = mapNearestSNPs_plink(matched_plink, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    p = plot_haplotype(chr_pos, brca_chr, brca_name, brca_start, brca_stop, mut, group=group, col = "genetic_dist")
    p[[4]]
    
    chr_coords_mod = chr_coords %>% mutate(cM_mod = (cM - brca_cM_middle)/100)
    num_breaks <- c()
    # From -0.17 to 0.15 with step 0.01, only representing end of interval
    intervals <- seq(-0.16, 0.15, 0.01)
    index_start = 0
    for (i in intervals){
        index_end = which.min(abs(chr_coords_mod$cM_mod - i))
        print(paste0(index_end, " ", chr_coords_mod$cM_mod[index_end]))
        num_breaks <- c(num_breaks, sum(breaks[, index_start:index_end]))
        index_start = index_end
    }
    plot(intervals, num_breaks)
}
