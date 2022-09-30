
#theta = 0.002

getData <- function(){
    mut = "c.7617+1G>A"; gene = "BRCA2"
    readIndvAndFamData(gene, mut)
    matched_plink_fam = matched_plink_fam[,-1]
    matched_plink_fam = filter(matched_plink_fam, SNP %in% subset(brca_data, Country == "DENMARK")$Onc_ID)
    haplotypes <<- findConsensus_plink(matched_plink_fam)
    breaks <<- findHaplotypeBreaks_plink(matched_plink_fam, haplotypes)
    chr_pos <<- mapNearestSNPs_plink(matched_plink_fam, breaks, chr_coords, brca_middle, brca_start, brca_stop, brca_chr, brca_name)
    # TESTING
    #chr_pos <<- subset(chr_pos, Country=="DENMARK")
    #chr_pos <<- subset(chr_pos, Country=="SPAIN")
    #chr_pos <<- chr_pos %>% mutate(positive_pos_adj=positive_pos-brca_middle) %>% mutate(negative_pos_adj=abs(negative_pos-brca_middle))
    morgan.brca1 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_Phase3/genetic_map_chr17_combined_b37.txt", header = T) # phase3 map
    morgan.brca2 <<- read.table("input/139_ouh_june_2017/genetic_maps/1000GP_Phase3/genetic_map_chr13_combined_b37.txt", header = T) # phase3 map
}


S2 <- function(dst, n){
    # Get distance from snp to mutation in Mb
    #dst_new = abs(dst - brca_middle) / 1000000
    #theta = 1/2 * (1 - exp(-2*dst_new))
    #return((1-theta)^n)

    #brca1.mut.cM = morgan.brca1[which.min(abs(morgan.brca1$position-brca_middle)),3]
    #brca2.mut.cM = morgan.brca2[which.min(abs(morgan.brca2$position-brca_middle)),3]
    theta = abs(dst - brca.mut.cM)/100
    # print(dst)
    # print(theta)
    # print(n)
    #print((1-theta)^n)
    #print(brca2.mut.cM)
    return((1-theta)^n)

}

f2 <- function(dst, n, right_side=T){
    i = 1
    index = match(dst, chr_coords$cM)
    repeat {
        dst2 = if (right_side) chr_coords$cM[index-1] else dst2 = chr_coords$cM[index+1]
        if ((dst-dst2) != 0) break
        i = i + 1
    }
    #print(S(dst2, n) - S(dst, n))
    return(S(dst2, n) - S(dst, n)) # Correct!
}

group1_prob <- function(group, n, right_side=T){
    ## L_a(n) = f(k,n)^y -> last marker, y individuals in group1, n generations
    ## L_b(n) = y*S(k,n)*f(k,n)^(y-1)
    ## Return L_a(n) + L_b(n)
    y = nrow(group)
    if (right_side){
        #k = group[1,2]
        k = group[1,4]
    } else {
        #k = group[1,3]
        k = group[1,5]
    }
    L_a = f(k,n,right_side)^y
    L_b = y*S(k,n)*f(k,n,right_side)^(y-1)
    return(L_a + L_b)
}

group2_prob <- function(group, n, right_side=T){
    if (right_side){
        #index = 2
        index = 4
    } else {
        #index = 3
        index = 5
    }
    L_i = prod(f(group[,index], n, right_side))
    return(L_i)
}

runLikelihood <- function(){
    longest_hap = chr_pos %>% count(positive_pos) %>% filter(n>1) %>% summarise(max=max(positive_pos))
    group1_right = chr_pos %>% filter(positive_pos >= longest_hap$max)
    group1_right$positive_pos = longest_hap$max
    group2_right = chr_pos %>% filter(positive_pos < longest_hap$max)

    longest_hap = chr_pos %>% count(negative_pos) %>% filter(n>1) %>% summarise(min=min(negative_pos))
    group1_left = chr_pos %>% filter(negative_pos <= longest_hap$min)
    group1_left$negative_pos = longest_hap$min
    group2_left = chr_pos %>% filter(negative_pos > longest_hap$min)

    maxLikelihood = 0
    generations = 0
    for (n in seq(1, 101, 10)){
        group1_prob_rightSide = group1_prob(group1_right, n, right_side = T)
        group2_prob_rightSide = group2_prob(group2_right, n, right_side = T)
        L_rightSide = group1_prob_rightSide * group2_prob_rightSide

        group1_prob_leftSide = group1_prob(group1_left, n, right_side = F)
        group2_prob_leftSide = group2_prob(group2_left, n, right_side = F)
        L_leftSide = group1_prob_leftSide * group2_prob_leftSide

        L = L_rightSide * L_leftSide
        if (L > maxLikelihood){
            maxLikelihood = L
            generations = n
        }
        print(paste("n:", n, "L:", L))
    }
    print(paste("generations:", generations, "Likelihood:", maxLikelihood))
}
#runLikelihood()

group1_prob2 <- function(k, n, y, right_side=T){
    L_a = f(k,n,right_side)^y
    L_b = y*S(k,n)*f(k,n,right_side)^(y-1)
    return(L_a + L_b)
}

group2_prob2 <- function(k, n, right_side=T){
    L_i = f(k, n, right_side)
    return(L_i)
}

computeMethod7GroupLikelihood <- function(group1_right, group2_right, group1_left, group2_left, gen){
    # Left side
    print(paste("generation:", gen))
    y = nrow(group1_left)
    group1_prob_left = sapply(group1_left$negative_cM, group1_prob2, n=gen, y=y, right_side=F)
    group2_prob_left = sapply(group2_left$negative_cM, group2_prob2, n=gen, right_side=F)
    prob_left_norm = sapply(c(group1_prob_left, group2_prob_left), getNumAbove)
    
    likelihood_left = productNorm(prob_left_norm[2,])
    zeros_left = sum(prob_left_norm[1,]) + likelihood_left[1]
    likelihood_left = likelihood_left[2]
    
    # Right side
    y = nrow(group1_right)
    group1_prob_right = sapply(group1_right$positive_cM, group1_prob2, n=gen, y=y, right_side=F)
    group2_prob_right = sapply(group2_right$positive_cM, group2_prob2, n=gen, right_side=F)
    prob_right_norm = sapply(c(group1_prob_right, group2_prob_right), getNumAbove)
    
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


#100/(2*median(abs(c(chr_pos$positive_pos, chr_pos$negative_pos)-brca_middle)) / 1000000)
#100/(median(abs(c(chr_pos$positive_pos, chr_pos$negative_pos)-brca_middle)) / 1000000)
