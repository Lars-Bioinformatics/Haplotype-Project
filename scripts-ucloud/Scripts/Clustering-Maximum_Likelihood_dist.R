library(parallel)

#setwd("~/Documents/PhD/haplotype-project")
print(getwd())
source("Scripts/Haplotype-projekt-ver2.R")
source("Scripts/Maximum-Likelihood.R")

#mut = "c.1016delA"; gene = "BRCA1"
mut = "c.7617+1G>A"; gene = "BRCA2"

###############################################################################
### Maximum likelihood distance
###############################################################################
## Get data needed for maximum likelihood
readData()
prepare_dataframe(mut, gene)
chr_pos <<- mapSNPsToCoords(chr_coords, matched)
pop_freqs <<- computePopFreqs()

#matched <- matched %>% mutate_all(as.character) %>% mutate_if("AA","A")
#matched <- sapply(matched, as.character)

#brca2_geno_singleHomo <- read.table("input/139_ouh_june_2017/139_ouh_brca2_onco_geno_single_homozygote.txt", header = T, stringsAsFactors = T)
#matched <<- extractSamples(brca_data, brca2_geno_singleHomo, chr_coords)

maxLikeDist <- function(h,a){
    # h: ancestral ref
    # a: haplotype sequence under study
    # G: # generations to test
    h <- getHaplotypeVector(h)
    a <- getHaplotypeVector(a)
    #h <<- h
    #a <<- a
    #G <<- 1
    factor = 1.15
    
    #print(a)
    #print(h)
    
    h_m <- h[1:(length(h)-1)]
    h_n <- h[2:length(h)]
    a_m <- a[1:(length(a)-1)]
    a_n <- a[2:length(a)]
    
    k = (prob2_vectorized(h_m, h_n, a_m, a_n) / prob1_vectorized(h_m, a_m)) * factor
    L1 = prob1_vectorized(h_m[1], a_m[1]) * prod(k)
    
    #print(k)
    
    k = (prob2_vectorized(a_m, a_n, h_m, h_n) / prob1_vectorized(a_m, h_m)) * factor
    L2 = prob1_vectorized(a_m[1], h_m[1]) * prod(k)
    
    k = (prob2_vectorized(h_m, h_n, h_m, h_n) / prob1_vectorized(h_m, h_m)) * factor
    L3 = prob1_vectorized(h_m[1], h_m[1]) * prod(k)

    k = (prob2_vectorized(a_m, a_n, a_m, a_n) / prob1_vectorized(a_m, a_m)) * factor
    L4 = prob1_vectorized(a_m[1], a_m[1]) * prod(k)
    
    ### Distance - own idea
    D = 1/mean(c(L1, L2))
    #print(c(L1,L2,D))
    
    ### Distance - David
    #D = -(L1+L2-L3-L4)
    #print(c(L1,L2,L3,L4,D))
    
    return(D)
}

getHaplotypeVector <- function(vec){
    vec <- sapply(vec, as.character)
    vec[vec=="AA"] = "A"
    vec[vec=="CC"] = "C"
    vec[vec=="GG"] = "G"
    vec[vec=="TT"] = "T"
    vec[vec=="DD"] = "D"
    vec[vec=="II"] = "I"
    return(vec)
}

testDist <- function(x,y){
    return(sqrt(abs(sum(x-y))))
}

#mat = matrix(1:9, nrow=3, ncol=3)
#proxy::dist(mat, method = testDist)

customDistance <- function(){
    start.time = Sys.time()
    snp_rm = 1000
    G <<- 1
    #dst <<- proxy::dist(matched[, (2+snp_rm):(ncol(matched)-snp_rm)], method = maxLikeDist)
    dst <<- proxy::dist(matched[1:2, (2+snp_rm):(ncol(matched)-snp_rm)], method = maxLikeDist)
    #dst <<- proxy::dist(matched[1:10, 2:ncol(matched)], method = maxLikeDist)
    #dst <<- proxy::dist(matched[1:20, 2:ncol(matched)], method = maxLikeDist)
    print(dst)
    #save(dst, file = paste0("dist_matrix/", brca_name, "-", mut, "-", "distance_matrix_DavDist.RRData"))
    save(dst, file = paste0("dist_matrix/", brca_name, "-", mut, "-", "distance_matrix.RRData"))
    print(Sys.time() - start.time)
}
customDistance()

customDistance_manual <- function(){
    start.time = Sys.time()
    snp_rm = 1000
    dist_mat = matrix(0, nrow=nrow(matched), ncol=nrow(matched))
    G <<- 1
    n = 2
    #n = nrow(dist_mat)
    #for (i in 1:(nrow(dist_mat)-1)){
    for (i in 1:(n-1)){
        if (i %% 5 == 0){
            print(paste0("i: ", i))
            print(Sys.time() - start.time)
        }
        #h <<- getHaplotypeVector(matched[i,2:ncol(matched)])
        h <<- getHaplotypeVector(matched[i,(2+snp_rm):(ncol(matched)-snp_rm)])
        #for (j in (i+1):(nrow(dist_mat))){
        for (j in (i+1):n){
            #a <<- getHaplotypeVector(matched[j,2:ncol(matched)])
            a <<- getHaplotypeVector(matched[j,(2+snp_rm):(ncol(matched)-snp_rm)])
            #print(paste(length(a), length(h)))
            #print(paste(i,j))
            #dist_mat[j,i] <- 1/mean(c(maxLikeDist(h, a), maxLikeDist(a, h)))
            #dist_mat[j,i] <- maxLikeDist(h, a)
            dist_mat[j,i] <- maxLikeDist(a, h)
            #print(c(i,j))
        }
    }
    print(as.dist(dist_mat[1:n,1:n]))
    dist_mat <- as.dist(dist_mat)
    save(dist_mat, file = paste0("dist_matrix/", brca_name, "-", mut, "-", "distance_matrix.RData"))
    print(Sys.time() - start.time)
}
#customDistance_manual()

customDistance_parManual <- function(){
    start.time = Sys.time()
    snp_rm = 1000
    G <<- 1
    #n = 5
    n = nrow(matched)
    dist_mat = matrix(0, nrow=n, ncol=n)
    cl <- makeCluster(detectCores())
    #cl = makeCluster(3)
    clusterExport(cl, list("maxLikeDist", "getHaplotypeVector", "prob2_vectorized", "matched",
                           "prob1_vectorized", "compareStrings_vectorized", "G", "getPopFreqs_vectorized",
                           "theta_vectorized", "theta_mn_vectorized", "chr_pos", "pop_freqs"))
    for (i in 1:(n-1)){
        if (i %% 5 == 0){
            print(paste0("i: ", i))
            print(Sys.time() - start.time)
        }
        h <<- getHaplotypeVector(matched[i,(2+snp_rm):(ncol(matched)-snp_rm)])
        clusterExport(cl, list("h"))
        # Apply rowise
        print(length(cl))
        print(n-i)
        if (length(cl) > (n-i)){
            #vec = apply(matched[(i+1):n, (2+snp_rm):(ncol(matched)-snp_rm)], 1, function(x) maxLikeDist(x, h))
            vec = parRapply(cl[1:(n-i)], matched[(i+1):n, (2+snp_rm):(ncol(matched)-snp_rm)], function(x) maxLikeDist(x, h))
        }else{
            vec = parRapply(cl, matched[(i+1):n, (2+snp_rm):(ncol(matched)-snp_rm)], function(x) maxLikeDist(x, h))
        }
        print(vec)
        dist_mat[,i] = c(rep(0,n-length(vec)),vec)
    }
    stopCluster(cl)
    dist_mat <- as.dist(dist_mat)
    print(dist_mat)
    #save(dist_mat, file = paste0("dist_matrix/", brca_name, "-", mut, "-", "distance_matrix_parallel.RData"))
    save(dist_mat, file = paste0("dist_matrix/", brca_name, "-", mut, "-", "distance_matrix_parallel_DavDist.RData"))
    print(Sys.time() - start.time)
    return(dist_mat)
}
dist_mat <- customDistance_parManual()

hc <- hclust(dist_mat, method = "ward.D")
hc <- hclust(dist_mat, method = "complete")
#hc <- hclust(dst, method = "ward.D")
#hc <- hclust(dst, method = "complete")
plot(hc)
k = 2
cluster_groups <<- cutree(hc, k=k)
names(cluster_groups) <- matched$SNP
rect.hclust(hc, k=k, border="red")

brca_data["cluster_groups"] <- cluster_groups
brca_data.filtered <- filter(brca_data, Country == "SPAIN" | Country == "DENMARK")
# Print group information
for (i in unique(cluster_groups)){
    print(subset(brca_data.filtered, cluster_groups == i) %>% select(Onc_ID, Country, cluster_groups))
}

test_small <- function(){
    danes <- brca_data %>% filter(Country == "DENMARK") %>% sample_n(5)
    spanes <- brca_data %>% filter(Country == "SPAIN") %>% sample_n(5)
    small_brca_data <- bind_rows(danes, spanes)
    matched <<- extractSamples(small_brca_data, haplotypes_brca, chr_coords)
    chr_pos <<- mapSNPsToCoords(chr_coords, matched)
    pop_freqs <<- computePopFreqs()
    customDistance()
    
    #hc <- hclust(dst, method = "complete")
    hc <- hclust(dst, method = "ward.D")
    plot(hc, hang = -1)
    k = 3
    cluster_groups <<- cutree(hc, k=k)
    names(cluster_groups) <- matched$SNP
    rect.hclust(hc, k=k, border="red")
    
    small_brca_data["cluster_groups"] <- cluster_groups
    # Print group information
    for (i in unique(cluster_groups)){
        print(subset(small_brca_data, cluster_groups == i) %>% select(Onc_ID, Country, cluster_groups))
    }
}
#test_small()
