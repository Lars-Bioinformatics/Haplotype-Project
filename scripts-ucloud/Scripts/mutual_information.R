library(FNN)
library(muti)
library(infotheo)

dim(matched_plink)
x = as.numeric(matched_plink[1,2100:2300])
y = as.numeric(matched_plink[2,2100:2300])

popFreqs = read.table(file = "cache/pop_freqs_genotypes_brca2_ordered.txt", header = T)
names(popFreqs) = c("0","1","2")
head(popFreqs)
p = popFreqs[2101:2301,]
dim(popFreqs)

x2=c()
y2=c()
for (i in 1:nrow(p)){
    x2[i] = p[i,as.character(x[i])]
    y2[i] = p[i,as.character(y[i])]
}
x2
y2

#x[x==0] = 1
#x[x==2] = 2
#y[y==0] = 1
#y[y==2] = 2

is.vector(x2)

mutinfo(x2, y2, k = 20, direct=T)
muti(x2, y2, normal = T)


set.seed(123)
TT <- 30
x1 <- rnorm(TT)
y1 <- x1 + rnorm(TT)
x1
y1


mutinfo(x1, y1, k = 10, direct=T)






sim_mut_info <- function(){
    chr_coords2 = filter(chr_coords, (cM - brca_cM_middle)>0)
    dim(chr_coords2)
    col_index = match(chr_coords2$SNP, names(matched_plink))
    row_index = match(subset(brca_data, cluster_groups == 1)$Onc_ID, matched_plink$SNP)
    m = matched_plink[,col_index]
    dim(m)
    # popFreqsBrca2 <- read.table2("cache/pop_freqs_genotypes_brca2_ordered.txt", header = T)[index,]
    # dim(popFreqsBrca2)
    #for (j in seq(1,(nrow(chr_coords2)-1),100)){
    for (j in 1:1){
        bits = c()
        snp1 = rep(0, nrow(m))
        snp1[row_index] = 2
        for (i in (j+1):(nrow(chr_coords2)-1)){
        #for (i in 2:10){
            #snp1_aa = if (popFreqsBrca2[j,"0"]==0) 0 else popFreqsBrca2[j,"0"] * log2(popFreqsBrca2[j,"0"])
            #snp1_ab = if (popFreqsBrca2[j,"1"]==0) 0 else popFreqsBrca2[j,"1"] * log2(popFreqsBrca2[j,"1"])
            #snp1_bb = if (popFreqsBrca2[j,"2"]==0) 0 else popFreqsBrca2[j,"2"] * log2(popFreqsBrca2[j,"2"])
            #snp1 = -(snp1_aa + snp1_ab + snp1_bb)
            # snp1_freq = plyr::count(snp1)$freq/length(snp1)
            # snp1 = -sum(snp1_freq * log2(snp1_freq))
            # snp2_aa = if (popFreqsBrca2[i,"0"]==0) 0 else popFreqsBrca2[i,"0"] * log2(popFreqsBrca2[i,"0"])
            # snp2_ab = if (popFreqsBrca2[i,"1"]==0) 0 else popFreqsBrca2[i,"1"] * log2(popFreqsBrca2[i,"1"])
            # snp2_bb = if (popFreqsBrca2[i,"2"]==0) 0 else popFreqsBrca2[i,"2"] * log2(popFreqsBrca2[i,"2"])
            # snp2 = -(snp2_aa + snp2_ab + snp2_bb)
            # 
            # #print(paste(snp1, snp2, snp1*snp2))
            # #print(snp1*snp2)
            # 
            # c=plyr::count(m[,c(1,i)])
            # freq=c$freq/nrow(m)
            # joint = -sum(freq*log2(freq))
            # mi = sum(snp1, snp2) - joint
            #print(mi)
            #print(paste(snp1+snp2, mi))
            
            snp2 = m[,i]
            mi = mutinfo(snp1,snp2,k=100)
            
            val = if (mi<0) 0 else mi
            bits = c(bits,val)
        }
        plot(x=1:length(bits), y = bits)
    }
    
    # Example
    A = c(0,1,0,2,0,1,1,0,1,0,1,0,2,1,0,0,1,0,1,0)
    B = c(2,2,2,0,1,1,2,2,0,2,1,1,0,1,2,1,2,2,1,2)
    freqA=plyr::count(A)$freq/length(A)
    freqB=plyr::count(B)$freq/length(B)
    snpA = -sum(freqA * log2(freqA))
    snpB = -sum(freqB * log2(freqB))
    #prob_indep=sum(snpA,snpB)
    
    df=data.frame(A,B)
    freqAB=plyr::count(df)$freq/nrow(df)
    snpAB = -sum(freqAB * log2(freqAB))
    
    # H(X) + H(Y) - H(X,Y)
    mi = snpA+snpB-snpAB
    mi
    
    mutinfo(A,B,k=8)
    muti(A,B,normal = T, lag=seq(1,10))
}


test_mutInfo <- function(){
    chr_coords2 = filter(chr_coords, (cM - brca_cM_middle)>0)
    dim(chr_coords2)
    #index = match(chr_coords2$SNP, names(matched_plink))
    #m = matched_plink[,index]
    all_data = rbind(brca1_geno_plink, brca2_geno_plink)
    index = match(chr_coords2$SNP, names(all_data))
    m = all_data[,c(1,index)]
    dim(m)
    snp_m = rep(0, nrow(m))
    snp_m[which(m$SNP %in% brca_data$Onc_ID)] = 2
    snp1 = m[,900]
    
    freqA=plyr::count(snp_m)$freq/length(snp_m)
    freqB=plyr::count(snp1)$freq/length(snp1)
    snpA = -sum(freqA * log2(freqA))
    snpB = -sum(freqB * log2(freqB))
    
    df=data.frame(snp_m,snp1)
    freqAB=plyr::count(df)$freq/nrow(df)
    snpAB = -sum(freqAB * log2(freqAB))
    
    # H(X) + H(Y) - H(X,Y)
    mi = snpA+snpB-snpAB
    mi
    
    mutinfo(snp_m,snp1,k=100)
    muti(snp_m,snp1, normal = T)
    mutinformation(snp_m, snp1)
    
    # LD
    all_data = rbind(famHaplotypes_brca1_geno_plink, famHaplotypes_brca2_geno_plink)
    index = match(chr_coords2$SNP, names(all_data))
    m = all_data[,index]
    dim(m)
    snp1 = 1
    snp2 = 2
    countA = plyr::count(m[,snp1])
    freqA=setNames(countA$freq/length(m[,snp1]), countA$x)
    countB = plyr::count(m[,snp2])
    freqB=setNames(countB$freq/length(m[,snp2]), countB$x)
    #countAB = plyr::count(m[, c(snp1,snp2)])
    #names(countAB) = c("snp1", "snp2", "freq")
    #LD = mutate(countAB, LD=countAB$freq/nrow(m) - freqA[as.character(snp1)]*freqB[as.character(snp2)])
    
    matched_plink2 <- filter(matched_plink, SNP %in% subset(brca_data, cluster_groups == 1)$Onc_ID) %>% select(as.character(chr_coords2$SNP))
    dim(matched_plink2)
    
    for (snp2 in 2:nrow(chr_coords2)){
        countA = plyr::count(m[,snp1])
        freqA=setNames(countA$freq/length(m[,snp1]), countA$x)
        countB = plyr::count(m[,snp2])
        freqB=setNames(countB$freq/length(m[,snp2]), countB$x)
        countAB = plyr::count(matched_plink2[, c(snp1,snp2)])
        names(countAB) = c("snp1", "snp2", "freq")
        LD = mutate(countAB, LD=countAB$freq/nrow(matched_plink2) - freqA[as.character(snp1)]*freqB[as.character(snp2)])
        LD = mutate(LD, LD_norm = ifelse(test = LD<0, 
                                         yes = max(-freqA[as.character(snp1)]*freqB[as.character(snp2)], -(1-freqA[as.character(snp1)])*(1-freqB[as.character(snp2)])),
                                         no = min(freqA[as.character(snp1)]*(1-freqB[as.character(snp2)]), (1-freqA[as.character(snp1)])*freqB[as.character(snp2)])
                                         )
                    )
        print(snp2)
        print(LD)
    }
}
