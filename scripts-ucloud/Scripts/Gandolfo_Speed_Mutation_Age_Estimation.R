# l.lengths = abs(chr_pos$negative_cM - brca_cM_middle)
# r.lengths = abs(chr_pos$positive_cM - brca_cM_middle)
# confidence.coefficient = 0.95
# chance.sharing.correction = F

run_gandolfo <- function(chr_pos, brca_cM_middle, confidence.coefficient = 0.95, chance.sharing.correction = F){
    gandolfo_age(l.lengths = chr_pos$negative_cM,
                 r.lengths = chr_pos$positive_cM,
                 confidence.coefficient = 0.95,
                 chance.sharing.correction = T,
                 chr_coords = chr_coords,
                 popFreqs = popFreqs,
                 haplotypes2 = haplotypes2)
}

gandolfo_age <- function(l.lengths = chr_pos$negative_cM,
                         r.lengths = chr_pos$positive_cM,
                         confidence.coefficient = 0.95,
                         chance.sharing.correction = F,
                         e = 0.01,
                         popFreqs, haplotypes2, chr_coords) {

    #Chance sharing correction
    if (chance.sharing.correction == TRUE){
        # my input
        #median.allele.frequency = 0.5
        ancestral_freqs = sapply(1:nrow(chr_coords), function(i) popFreqs[i, as.character(haplotypes2[i])])
        median.allele.frequency = median(na.omit(ancestral_freqs))
        median.allele.frequency
        
        length.of.chromosome = max(chr_coords[,4])-min(chr_coords[,4])
        #length.of.chromosome = max(r.lengths)-min(l.lengths)
        markers.on.chromosome = length(chr_coords[,4])
        #markers.on.chromosome = nrow(filter(chr_coords, cM <= max(r.lengths), cM >= min(l.lengths)))
        
        # original
        #e = 0.01
        p = (median.allele.frequency)^2 + (1-median.allele.frequency)^2
        phi = (length.of.chromosome/100)/markers.on.chromosome
        loci = log(e)/log(p)
        cs.correction = loci*phi
    }
    if (chance.sharing.correction == FALSE){cs.correction = 0}
    
    # Adjust arm lengths relative to mutation
    l.lengths = abs(l.lengths - brca_cM_middle)
    r.lengths = abs(r.lengths - brca_cM_middle)
    
    #Age estimation and confidence intervals
    cc = confidence.coefficient
    l.lengths = (1/100)*l.lengths 
    r.lengths = (1/100)*r.lengths
    n = length(l.lengths)
    
    #Assuming an 'independent' genealogy
    if (n < 10){i.cs.correction = 0}
    if (n >= 10){i.cs.correction = cs.correction}
    
    length.correction = (sum(l.lengths) + sum(r.lengths) - 2*(n-1)*i.cs.correction)/(2*n)
    
    sum.lengths = sum(l.lengths) + sum(r.lengths) + 2*length.correction - 2*(n-1)*i.cs.correction
    b.c = (2*n-1)/(2*n)
    i.tau.hat <- (b.c*2*n)/sum.lengths
    g_l <- qgamma(shape=2*n,scale=1/(2*n*b.c),((1-cc)/2))
    g_u <- qgamma(shape=2*n,scale=1/(2*n*b.c),(cc+(1-cc)/2))		
    i.l = g_l*i.tau.hat
    i.u = g_u*i.tau.hat
    
    #Assuming a 'correlated' genealogy
    length.correction = (sum(l.lengths) + sum(r.lengths) - 2*(n-1)*cs.correction)/(2*n)
    
    longest.l.lengths = match(sort(l.lengths,decreasing=TRUE)[1], l.lengths)
    l.lengths[longest.l.lengths] <- l.lengths[longest.l.lengths] + length.correction + cs.correction
    longest.r.lengths = match(sort(r.lengths,decreasing=TRUE)[1], r.lengths)
    r.lengths[longest.r.lengths] <- r.lengths[longest.r.lengths] + length.correction + cs.correction
    
    lengths = l.lengths + r.lengths
    lengths = lengths - 2*cs.correction
    var_lengths = if (length(lengths)>1) var(lengths) else 0
    rho.hat = (n*(mean(lengths))^2 - var_lengths*(1+2*n))/(n*(mean(lengths))^2 + var_lengths*(n-1))
    #rho.hat = (n*(mean(lengths))^2 - var(lengths)*(1+2*n))/(n*(mean(lengths))^2 + var(lengths)*(n-1))
    n.star = n/(1+(n-1)*rho.hat)
    if (n.star > n) {n.star = n}
    if (n.star < -n) {n.star = -n}		
    b.c = (2*n.star-1)/(2*n.star)
    c.tau.hat = (b.c*2*n)/sum(lengths)
    if (rho.hat < -2/(n-1)){n.star = n/(1+(n-1)*abs(rho.hat))}
    if (-2/(n-1) <= rho.hat & rho.hat < -1/(n-1)){n.star = n}
    g_l = qgamma(shape=2*n.star,scale=1/(2*n.star*b.c),(1-cc)/2)
    g_u = qgamma(shape=2*n.star,scale=1/(2*n.star*b.c),cc+(1-cc)/2)
    c.l = g_l*c.tau.hat
    c.u = g_u*c.tau.hat
    
    #Print results
    print(paste("Assuming an 'independent' genealogy: age estimate =", 
                round(i.tau.hat, digits=1) ,
                "generations, with confidence interval", 
                paste("(",round(i.l, digits=1),",",
                      round(i.u, digits=1),")",sep="")), 
          quote=FALSE);
    
    print(paste("Assuming a 'correlated' genealogy: age estimate =", 
                round(c.tau.hat, digits=1) ,
                "generations, with confidence interval", 
                paste("(",round(c.l, digits=1),",",
                      round(c.u, digits=1),")",sep="")), 
          quote=FALSE)
    # Return (independent_age, indp_conf_lower, indp_conf_upper, correlated_age, corr_conf_lower, corr_conf_upper)
    return(c(round(i.tau.hat, digits=1), round(i.l, digits=1), round(i.u, digits=1), round(c.tau.hat, digits=1), round(c.l, digits=1), round(c.u, digits=1)))
}
