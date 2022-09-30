marker_before_break = function(chr_pos){
    get_x2 <- function(row,n,right_side,chr_pos){
        #print(row)
        #j = which(sid == chr_pos$sample_id)
        j = row
        i = 1
        x = if (right_side) chr_pos[j, "positive_cM"] else chr_pos[j, "negative_cM"]
        
        if (x == max(chr_coords$cM) || x == min(chr_coords$cM)){
            #return(x-0.00001)
            return(x)
        }
        
        index = match(x, chr_coords$cM)
        repeat{
            index2 = if (right_side) index-i else index+i
            # Check that we don't cross the mutation and measure on wrong side
            if (chr_coords$cM[index2] > brca.mut.cM && !right_side){
                x2 = max(chr_coords$cM[chr_coords$cM<brca_cM_middle])
                return(x2)
            }
            if (chr_coords$cM[index2] < brca.mut.cM && right_side){
                x2 = min(chr_coords$cM[chr_coords$cM>brca_cM_middle])
                return(x2)
            }
            
            # Find marker closest to breakpoint with a cM distance > 0
            cmp = matched_plink2[j, index2] - haplotypes2[index2]
            if (haplotypes2[index2] != 1 && cmp == 0){
                x2 = chr_coords$cM[index2]
                # Check for genetic distance above 0
                if (x2 - x != 0) break
            }
            i = i + 1
        }
        
        return(x2)
    }
    chr_pos$negative_cM_preMarker = sapply(1:nrow(chr_pos), get_x2, right_side=F, chr_pos = chr_pos)
    chr_pos$positive_cM_preMarker = sapply(1:nrow(chr_pos), get_x2, right_side=T, chr_pos = chr_pos)
    # If break at closest marker
    chr_pos$negative_cM[chr_pos$negative_cM==max(chr_coords$cM[chr_coords$cM < brca.mut.cM])] = max(chr_coords$cM[chr_coords$cM < max(chr_coords$cM[chr_coords$cM < brca.mut.cM])-0.00001])
    chr_pos$positive_cM[chr_pos$positive_cM==min(chr_coords$cM[chr_coords$cM > brca.mut.cM])] = min(chr_coords$cM[chr_coords$cM > min(chr_coords$cM[chr_coords$cM > brca.mut.cM])+0.00001])
    adjusted_lengths = data.frame(break_cM = c(chr_pos$positive_cM-brca.mut.cM, abs(chr_pos$negative_cM-brca.mut.cM)),
                                  before_break_cM = c(chr_pos$positive_cM_preMarker-brca.mut.cM, abs(chr_pos$negative_cM_preMarker-brca.mut.cM)))
}

maxLikelihoodGeneration <- function(n, y){
    prob = sum(log(ifelse((1-(1-abs(y$break_cM-y$before_break_cM)/100)^n)==0,
                          (1-abs(y$before_break_cM)/100)^n, 
                          (1-abs(y$before_break_cM)/100)^n*(1-(1-abs(y$break_cM-y$before_break_cM)/100)^n))))
    #prob = sum(log((1-abs(y$before_break_cM)/100)^n * (1-(1-abs(y$break_cM-y$before_break_cM)/100)^n)))
    # p1 = sum(log((1-y$before_break_cM/100)^n))
    # p2 = (1-(1-abs(y$break_cM-y$before_break_cM)/100)^n)
    # prob = sum(p1, sum(log(p2[p2!=0])))
    # print(n)
    # print(p1)
    return(prob)
}

method5_fast <- function(chr_pos,i=NULL){
    #print(i)
    if (is.null(i)){
        y = marker_before_break(chr_pos)
    } else {
        y = chr_pos[i,]
    }
    #y[y$break_cM==y$before_break_cM,]$break_cM = y[y$break_cM==y$before_break_cM,]$break_cM+100
    l = suppressWarnings(optimize(f = maxLikelihoodGeneration, maximum = T, interval = c(0,10000), y = y))
    #print(paste("Method 5 age:", l$maximum))
    return(l$maximum)
}

test_maxLikelihood <- function(){
    y = marker_before_break(chr_pos)
    #optim(par = c(100), fn = maxLikelihoodGeneration, y = y, method = "Brent", lower = 0, upper = 1000)
    optimize(f = maxLikelihoodGeneration, maximum = T, interval = c(0,10000), y = y)
    l=sapply(1:1000, maxLikelihoodGeneration, y)
    which.max(l)
    l[which.max(l)]
}

method5_fast_with_bootstrap <- function(chr_pos){
    set.seed(626)
    #cl = makeCluster(detectCores())
    #clusterExport(cl, c("marker_before_break", "chr_coords", "brca.mut.cM"))
    y = marker_before_break(chr_pos)
    #bootcorr <- boot(data = y, statistic = method5_fast, R=100, cl = cl, parallel = "multicore", ncpus = detectCores())
    bootcorr <- boot(data = y, statistic = method5_fast, R=100)
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
