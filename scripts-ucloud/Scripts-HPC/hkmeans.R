hkmeans_mod <- function (x, k, hc.metric = "euclidean", hc.method = "ward.D2", 
          iter.max = 10, km.algorithm = "Hartigan-Wong") 
{
    #res.hc <- stats::hclust(stats::dist(x, method = "euclidean"), method = hc.method)
    res.hc <- hclust(firstBreakDist(fam=F), method = hc.method)
    grp <- stats::cutree(res.hc, k = k)
    clus.centers <- stats::aggregate(x, list(grp), mean)[, -1]
    res.km <- kmeans(x, centers = clus.centers, iter.max = iter.max, 
                     algorithm = km.algorithm)
    class(res.km) <- "hkmeans"
    res.km$data <- x
    res.km$hclust <- res.hc
    res.km
}
