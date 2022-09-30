#setwd("~/Documents/PhD/haplotype-project")
source("Scripts/Haplotype-projekt-ver2.R")

#mut = "c.1016delA"; gene = "BRCA1"
mut = "c.7617+1G>A"; gene = "BRCA2"

#getData()
getData <- function(){
    readData()
    
    #readFamData()
    #readIndvAndFamData(gene, mut)
    #hc <- hiearchical_clustering(firstBreakDist(fam=T), clustMethod = "complete", k = 2, doHapPlot = F)
    #famHaplotypes <<- filter(famHaplotypes_brca2_geno, FamCode %in% (filter(brca_data.filtered, cluster_groups==1)$FamCode))
    #famHaplotypes <<- famHaplotypes[,3:ncol(famHaplotypes)]
    
    prepare_dataframe(mut, gene)
    famHaplotypes <<- findMutFamHaplotypes(mut, gene) # Family haplotypes
    dan_fam <- filter(brca_data, Country == "DENMARK") %>% group_by(FamCode) %>% filter(n()>1) %>% distinct(FamCode)
    famHaplotypes <<- famHaplotypes %>% filter(FamCode %in% dan_fam$FamCode)
    
    chr_pos <<- mapSNPsToCoords(chr_coords, famHaplotypes)
    h <<- findConsensus(famHaplotypes) # ancestral haplotype
    h <<- c(h[fmatch(chr_pos$SNP, names(h))])
    #h <<- h[2700:4000]
    h <- h[1312:3806]
    pop_freqs <<- computePopFreqs()
}

#likelihood <- runMaximumLikelihood_vectorized()
runMaximumLikelihood_vectorized <- function(){
    start.time = Sys.time()
    
    likelihood <<- list()
    #i = 1
    #g <- seq(5, 200, 10)
    for (G in seq(5, 200, 10)) {
    #for (G in c(1,5,15)) {
        start.time2 = Sys.time()
        print(G)
        L <- Likelihood_vectorized(G)
        like <- data.frame(G, L)
        #print(like)
        #like <- c(G,L)
        names(like) = c("Generation", "Likelihood")
        likelihood[[G]] <- like
        #i = i+1
        print(Sys.time() - start.time2)
        #like <- c(G,L)
        #like
    }
    likelihood <- bind_rows(likelihood)
    
    print(likelihood)
    end.time = Sys.time() - start.time
    print(end.time)
    return(likelihood)
}

#runMaximumLikelihood()
runMaximumLikelihood <- function(){
    start.time = Sys.time()
    
    #G <- seq(5, 200, 10)
    G = c(1,5,15)
    L <- Likelihood(G)
    likelihood <- data.frame(G,L)
    
    #likelihood <<- list()
    #i = 1
    #g <- seq(5, 200, 10)
    #likelihood <- foreach(G=g) %dopar% {
    # for (G in seq(5, 200, 10)) {
    #     start.time2 = Sys.time()
    #     print(G)
    #     #Likelihood(G = 30)
    #     L <- Likelihood(G)
    #     #L = 0.02
    #     like <- data.frame(G, L)
    #     #like <- c(G,L)
    #     names(like) = c("Generation", "Likelihood")
    #     likelihood[[G]] <- like
    #     #i = i+1
    #     print(Sys.time() - start.time2)
    #     #like <- c(G,L)
    #     #like
    # }
    #likelihood <- bind_rows(likelihood)
    
    print(likelihood)
    end.time = Sys.time() - start.time
    print(end.time)
}

# for (f in seq(1.5,1.7,0.01)){
#     print(paste("factor:", f))
#     #Likelihood_vectorized(G = 5, factor = f)
#     Likelihood_vectorized(G = 10, factor = f)
# }

# testLikelihood()
testLikelihood <- function(){
    for (g in seq(20,100,20)){
        print(paste("Testing Generation:", g))
        #Likelihood_vectorized(G = g, factor = 1.313)
        #Likelihood_vectorized(G = g, factor = 1.285)
        Likelihood_vectorized(G = g, factor = 1.32)
    }
}

Likelihood_vectorized <- function(G, factor = 1.32){
    G <<- G
    
    h_m <- h[1:(length(h)-1)]
    h_n <- h[2:length(h)]
    
    L = 1
    #for(i in 1:1){
    #for(i in 10:11){
    for(i in 1:nrow(famHaplotypes)){
        if (L == 0 || L == Inf) next
        #a <<- as.character(unlist(famHaplotypes[i,]))
        #names(a) <- names(famHaplotypes)
        
        a <- unlist(famHaplotypes[i,])
        a <- c(a[fmatch(chr_pos$SNP, names(a))])
        #a <- a[2700:4000]
        a <- a[1312:3806]
        
        a_m <- a[1:(length(a)-1)]
        a_n <- a[2:length(a)]
        
        l = prob1_vectorized(h_m[1], a_m[1])
        k = (prob2_vectorized(h_m, h_n, a_m, a_n) / prob1_vectorized(h_m, a_m)) * factor
        
        if ((l * prod(k)) > 1){
            print("LOLLOLLOLLOLLOLLOLLOLLOLLOL")
            print(l * prod(k))
            print(i)
            print("LOLLOLLOLLOLLOLLOLLOLLOLLOL")
        }
        
        L = L * abs(log(l * prod(k)))
        print(L)
    }
    print(paste("L:", L))
    return(L)
}

# testPairwiseLikelihood()
testPairwiseLikelihood <- function(){
    for (g in seq(10,300,20)){
        print(paste("Testing Generation:", g))
        #Likelihood_pairwise(G = g, factor = 1.7)
        #Likelihood_pairwise(G = g, factor = 1.5)
        #Likelihood_pairwise(G = g, factor = 1.335)
        Likelihood_pairwise(G = g, factor = 1.2)
        #Likelihood_pairwise(G = g, factor = 1)
    }
}

Likelihood_pairwise <- function(G, factor = 1.32){
    G <<- G
    
    L = 1
    for (run in 1:2){
        #for(j in 1:1){
        #for(j in 1:3){
        for (j in 1:nrow(famHaplotypes)){
            
            L_temp = 1
            
            h <- unlist(famHaplotypes[j,])
            h <- c(h[fmatch(chr_pos$SNP, names(h))])
            #h <- h[1312:3806]
            if (run==1){
                h <- h[1312:2694]
            } else {
                h <- h[2695:3806]
            }
            h_m <- h[1:(length(h)-1)]
            h_n <- h[2:length(h)]
            
            #for(i in 2:2){
            #for(i in 1:3){
            for(i in 1:nrow(famHaplotypes)){
                #if (i == j) next
                if (L_temp == 0 || L_temp == Inf) next
                #a <<- as.character(unlist(famHaplotypes[i,]))
                #names(a) <- names(famHaplotypes)
                
                # print(paste("i:", i, "j:", j))
                # print(L_temp)
                # print(L)
                
                a <- unlist(famHaplotypes[i,])
                a <- c(a[fmatch(chr_pos$SNP, names(a))])
                #a <- a[2700:4000]
                #a <- a[1312:3806]
                if (run==1){
                    a <- a[1312:2694]
                } else {
                    a <- a[2695:3806]
                }
                
                a_m <- a[1:(length(a)-1)]
                a_n <- a[2:length(a)]
                
                l = prob1_vectorized(h_m[1], a_m[1])
                k = (prob2_vectorized(h_m, h_n, a_m, a_n) / prob1_vectorized(h_m, a_m)) * factor
                
                if ((l * prod(k)) > 1){
                    print("LOLLOLLOLLOLLOLLOLLOLLOLLOL")
                    print(l * prod(k))
                    print(i)
                    print("LOLLOLLOLLOLLOLLOLLOLLOLLOL")
                }
                
                L_temp = L_temp * abs(log(l * prod(k)))
                # print(L_temp)
            }
            if (L_temp == Inf){
                print("OBS: is Inf")
            } else {
                L = L * log(L_temp)
            }
            print(paste("L:", L))
        }
    }
    print(L)
    return(L)
}

Likelihood <- function(G){
    G <<- G
    L = 1
    if (brca_name == "BRCA1"){
        factor = 1.35
    } else {
        factor = 1.35
    }
    #for(i in 1:nrow(famHaplotypes)){
    for(i in 1:3){
    #for(i in 1:1){
        a <<- unlist(famHaplotypes[i,])
        a <<- c(a[fmatch(chr_pos$SNP, names(a))])
        a <<- a[2700:4000]

        l = prob1(1)
        kage = numeric(length(a)-1)
        for(m in 1:(length(a)-1)){
        #for(m in 4102:4102){
        #for(m in 4000:(length(a)-1)){
            #m <<- m
            #print(paste0("m: ", m))
            #print(names(h[m]))
            #print(paste0("likelihood: ", l))
            # if (is.na(l)){
            #     stop("likelihood is NA!")
            # }
            k = (prob2(m, m+1) / prob1(m)) * factor
            #print(paste0("Next val: ", k))
            l = l * k
            #kage[m] = k
        }
        
        L = L * l
        print(paste0("i: ", i))
        #print(l)
        print(L)
    }
    print(L)
    return(L)
}

prob1_vectorized <- function(h_m, a_m){
    prob1_vector <- numeric(length(a_m))
    index = compareStrings_vectorized(a_m, h_m)
    if (length(a_m[index]) > 0){
        prob1_vector[index] = theta_vectorized(h_m[index])*getPopFreqs_vectorized(h_m[index]) + (1 - theta_vectorized(h_m[index]))
    }
    index = !compareStrings_vectorized(a_m, h_m)
    if (length(a_m[index]) > 0){
        prob1_vector[index] = theta_vectorized(h_m[index])*getPopFreqs_vectorized(a_m[index])
    }
    return(prob1_vector)
}

prob1 <- function(m){
    ifelse(a[m] == h[m],
           return( theta(m)*getPopFreqs(h[m]) + (1 - theta(m)) ),
           return( theta(m)*getPopFreqs(a[m]) )
    )
}

prob2_vectorized <- function(h_m, h_n, a_m, a_n){
    
    prob2_vector <- numeric(length(a_m))
    
    index = compareStrings_vectorized(a_m, h_m) & compareStrings_vectorized(a_n, h_n)
    if (length(a_m[index])>0){
        prob2_vector[index] = 
            theta_vectorized(h_m[index]) * getPopFreqs_vectorized(h_m[index]) * getPopFreqs_vectorized(h_n[index]) 
            + (1-theta_vectorized(h_m[index])) * (1 - theta_mn_vectorized(h_m[index], h_n[index]) + 
            theta_vectorized(h_n[index]) * getPopFreqs_vectorized(h_n[index]))
    }
    
    index = compareStrings_vectorized(a_m, h_m) & !compareStrings_vectorized(a_n, h_n)
    if (length(a_m[index])>0){
        prob2_vector[index] = 
            theta_vectorized(h_m[index]) * getPopFreqs_vectorized(h_m[index]) * getPopFreqs_vectorized(a_n[index]) + 
            (1-theta_vectorized(h_m[index])) * (theta_mn_vectorized(h_m[index],h_n[index]) * getPopFreqs_vectorized(a_n[index]))
    }
    
    index = !compareStrings_vectorized(a_m, h_m) & compareStrings_vectorized(a_n, h_n)
    if (length(a_m[index])>0){
        prob2_vector[index] =
            theta_vectorized(h_m[index]) * getPopFreqs_vectorized(a_m[index]) * getPopFreqs_vectorized(h_n[index])
    }
    
    index = !compareStrings_vectorized(a_m, h_m) & !compareStrings_vectorized(a_n, h_n)
    if (length(a_m[index])>0){
        prob2_vector[index] =
            theta_vectorized(h_m[index]) * getPopFreqs_vectorized(a_m[index]) * getPopFreqs_vectorized(a_n[index])
    }
    
    return(prob2_vector)
}

# Assuming my_m (snp mutation freq) = 0
prob2 <- function(m, n){
    if(compareStrings(a[m], h[m]) & compareStrings(a[n], h[n])){
        #print("prob2: 1")
        return( theta(m) * getPopFreqs(h[m]) * getPopFreqs(h[n]) +
            (1-theta(m)) * (1 - theta_mn(m,n) + theta(n) * getPopFreqs(h[n])) )

    }else if(compareStrings(a[m], h[m]) & !compareStrings(a[n], h[n])){
        #print("prob2: 2")
        return( theta(m) * getPopFreqs(h[m]) * getPopFreqs(a[n]) + 
            (1-theta(m)) * (theta_mn(m,n) * getPopFreqs(a[n])) )

    }else if(!compareStrings(a[m], h[m]) & compareStrings(a[n], h[n])){

        return( theta(m) * getPopFreqs(a[m]) * getPopFreqs(h[n]) )

    }else if(!compareStrings(a[m], h[m]) & !compareStrings(a[n], h[n])){

        return( theta(m) * getPopFreqs(a[m]) * getPopFreqs(a[n]) )
        
    }
}

theta_vectorized <- function(h_m){
    #X_m = abs(filter(chr_pos, chr_pos$SNP %in% names(h_m))$brca_distance) / 1000000
    #X_m = abs(chr_pos[chr_pos$SNP %in% names(h_m),]$brca_distance) / 1000000
    X_m = abs(chr_pos[chr_pos$SNP %in% names(h_m),]$brca_distance) / 1000 # distance in Kb (kilobases)
    return(1-exp(-X_m*G/1000))
}

# Probability of recombination between locus m at distance X kb from BRCA mut
theta <- function(m){
    X_m = abs(filter(chr_pos, SNP == names(h[m]))$brca_distance) / 1000000
    # print(1-exp(-G*X_m))
    return(1-exp(-G*X_m))
}

theta_mn_vectorized <- function(h_m, h_n){
    #X_m = abs(filter(chr_pos, chr_pos$SNP %in% names(h_m))$brca_distance) / 1000000
    #X_m = abs(chr_pos[chr_pos$SNP %in% names(h_m),]$brca_distance) / 1000000
    X_m = abs(chr_pos[chr_pos$SNP %in% names(h_m),]$brca_distance) / 1000 # distance in Kb (kilobases)
    
    #X_n = abs(filter(chr_pos, chr_pos$SNP %in% names(h_n))$brca_distance) / 1000000
    #X_n = abs(chr_pos[chr_pos$SNP %in% names(h_n),]$brca_distance) / 1000000
    X_n = abs(chr_pos[chr_pos$SNP %in% names(h_n),]$brca_distance) / 1000 # distance in Kb (kilobases)
    return( 1-exp(-G * abs(X_n - X_m) / 1000) )
}

theta_mn <- function(m, n){
    X_m = abs(filter(chr_pos, SNP == names(h[m]))$brca_distance) / 1000000
    X_n = abs(filter(chr_pos, SNP == names(h[n]))$brca_distance) / 1000000
    return( 1-exp(-G * abs(X_n - X_m)) )
}

mapSNPsToCoords <- function(chr_coords, matched){
    chr_pos = filter(chr_coords, SNP %in% names(matched)) %>% 
        arrange(position_b37) %>%
        mutate(brca_distance = as.integer(as.character(position_b37)) - brca_middle)
    return(chr_pos)
}

compareStrings_vectorized <- function(a, b){
    bools <- logical(length(a))
    bools[a == b] = TRUE
    bools[nchar(a) > 1 | nchar(b) > 1] = TRUE
    return(bools)
}

compareStrings <- function(a, b){
    # ifelse(test, yes, no) - vectorized function
    ifelse(nchar(a) == 1 & nchar(b) == 1,
           return(a == b), 
           #return((strsplit(b, split = "")[[1]] %in% strsplit(a, split = "")[[1]])[1])
           return(sum(strsplit(b, split = "")[[1]] %in% strsplit(a, split = "")[[1]]))
           )
}

computePopFreqs <- function(){
    # Load cached pop frequencies, if available
    if(file.exists("cache/pop_freqs.RData")){
        load("cache/pop_freqs.RData")
        return(pop_freqs)
    }
    
    ### generate population frequencies ###
    brca_geno_reference = bind_rows(brca1_geno_merged, brca2_geno_merged)
    
    start.time = Sys.time()
    pop_freqs = c()
    for (pos in 2:ncol(brca_geno_reference)){
    #for (pos in 2:7){
        # Combine a SNP column to one string and split string into one-character vector, 
        # then count the number of occurences of each element
        counts <- table(unlist(strsplit(as.character(brca_geno_reference[, pos]), split = "")))

        if (names(counts[1]) == "-"){
            counts <- counts[2:length(counts)]
        }
        
        freq = counts/sum(counts)
        
        # #print(paste("pos:", pos))
        # #print(counts)
        # #print(freq)
        pop_freqs = c(pop_freqs, list(freq))
    }
    names(pop_freqs) = colnames(brca_geno_reference[2:ncol(brca_geno_reference)])
    print(Sys.time()-start.time)
    head(pop_freqs)
    save(pop_freqs, file = "cache/pop_freqs.RData")
    return(pop_freqs)
}

getPopFreqs_vectorized <- function(b){
    a = numeric(length(b))
    
    # Not exactly sure what to do about missing data, or
    # locations with multiple possible alleles
    a[nchar(b) > 1 | b == "-"] = 0.5 #= 1 #= 0.5
    
    # a[nchar(b) > 1] = 1
    # 
    # index <- b[b == "-"]
    # if (length(index) > 0){
    #     names(index) <- names(b)[b=="-"]
    #     int <- numeric(length(index))
    #     for (i in 1:length(index)){
    #         int[i] = max(pop_freqs[[names(index[i])]])
    #     }
    #     a[b == "-"] = int
    # }
    
    #index <- which(nchar(b)==1 & b != "-")
    
    
    index <- b[nchar(b)==1 & b != "-"]
    if (length(index) > 0){
        names(index) <- names(b)[nchar(b)==1 & b != "-"]
        int <- numeric(length(index))
        #print("lol")
        #print(index)
        #print(b)
        #print(names(index))
        for (i in 1:length(index)){
            #print(i)
            #print(names(index[i]))
            int[i] = pop_freqs[[names(index[i])]][index[i]]
        }
        a[nchar(b) == 1 & b != "-"] = int
    }
    
    return(a)
}

getPopFreqs <- function(b){
    if (nchar(b) > 1){
        return(0.5)
        #return(1)
    } else if (b == "-"){
        # If missing value, assume most frequent (could be biased)
        return(max(pop_freqs[[names(b)]]))
    } else {
        return(pop_freqs[[names(b)]][b])
    }
}

convertPopFreqsToTable <- function(){
    df = as.data.frame(matrix(nrow=length(pop_freqs), ncol = 8, data = 0))
    names(df) = c("SNP", "A", "C", "G", "T", "D", "I", "-")
    df[,1] = names(pop_freqs)
    for (i in 1:length(pop_freqs)){
        for (j in 1:length(pop_freqs[[i]])){
            print(pop_freqs[[i]][j])
            index = which(names(pop_freqs[[i]][j])==names(df))
            df[i,index] = pop_freqs[[i]][j]
        }
    }
    df[1:10,]
    df <- df[,1:7]
}

#getPopulationFreqs(brca1_geno_merged)
#getData()

#print(pop_freqs[[names(h[11])]])

#Likelihood()

# ifelse(nchar(a) == 1 & nchar(b) == 1, 
#        print(a == b), 
#        print((strsplit(b, split = "")[[1]] %in% strsplit(a, split = "")[[1]])[1])
# )
