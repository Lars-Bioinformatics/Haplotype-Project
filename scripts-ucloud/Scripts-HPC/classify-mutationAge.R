# library(e1071)
# library(randomForest)
# library(Metrics)

abacus = T
if (abacus){
    dists <- read.table(file = "training-distances-generations_10-100.txt", header = F)
    res <- read.table(file = "res.txt", header = F)
} else {
    dists <- read.table(file = "classify-simulated-population-age/training-distances-generations_10-100.txt", header = F)
    res <- read.table(file = "classify-simulated-population-age/res.txt", header = F)
}

run_loocv <- function(parallel=T){
    if (parallel){
        library(parallel)
        cl <- makeCluster(detectCores(), outfile="progress.txt")
        #cl <- makeCluster(2, outfile="~/log.txt")
        clusterExport(cl, varlist = c("dists", "res"), envir = .GlobalEnv)
        results <- as.data.frame(t(parSapply(cl = cl, X = 1:24, FUN = loocv)))
        #results <- as.data.frame(t(parSapply(cl = cl, X = 1:nrow(dists), FUN = loocv)))
        stopCluster(cl)
    } else {
        #results <- as.data.frame(t(sapply(1:nrow(dists), loocv)))
        results <- as.data.frame(t(sapply(1:32, loocv)))
    }

    names(results) <- c("sample", "predicted_age", "actual_age")
    print(results)
    return(results)
}


### Leave one out cross validation ###
loocv <- function(sample_id){
    start_time <- Sys.time()
    set.seed(42)
    library(e1071)
    library(randomForest)
    library(Metrics)

    #results <- {}

    # print(sample_id)
    # sample_id = 1

    all <- 1:nrow(dists)
    nvar <- ncol(dists)
    #for (d in 1:nrow(dists)){
    for (d in sample_id:sample_id){

        externalSample <- d
        inpairs <- all[-d]

        train.mat <- dists[inpairs, ]
        train.res <- res[inpairs,]

        ### Random Forest
        gna <- c(1:nvar)

        mt <- sqrt(nvar)

        train <- data.frame(class = train.res, p=train.mat)
        model <- randomForest(class ~ ., data = train, mtry=mt, importance = T, na.action = na.omit)#, ntrees = 5000)

        imp11<-cbind(importance(model),c(1:nvar))
        imp<-imp11[rev(order(imp11[,2])),]
        rf_rank<-as.numeric(imp[,3])

        test.res <- train.res

        nsel <- length(rf_rank)

        RMSEsvm <- {}
        for (ct in 2:nsel){
        #for (ct in 2:100){
            ans <- {}
            gn3 <- rf_rank[c(1:ct)]

            all2 <- 1:99
            for (a in 1:(nrow(dists)-1)){
            #for (a in 1:1){
                train1 <- all2[-a]
                test1 <- all2[a]

                classtrain <- train.res[train1]
                train2 <- train.mat[train1, gn3]

                classtest <- test.res[test1]
                test2 <- train.mat[test1, gn3]

                train <- data.frame(class = classtrain, p=train2)
                model2 <- svm(class ~ ., data = train)

                test <- data.frame(p=test2)
                ans[a] <- predict(model2, test)

            }

            #MEANsvm <- mean(ans-train.res)
            RMSEsvm <- c(RMSEsvm, rmse(ans, train.res))

            # if (d %% 10 == 0){
            #     print(paste("S-SVM SAM SAMPLE-OUT: ",d, "GENE",ct, "Out of ",nsel))
            # }
            # plot(1:ct,RMSEsvm[1:ct],type="b",ylim=c(40,100),main=paste("S-SVM SAM SAMPLE-OUT: ",d, "GENE",ct, "Out of ",nsel))
            # lines(1:ct,sen2[1:ct],col="red")
            # lines(1:ct,spe2[1:ct],col="blue")


        }
        #print(RMSEsvm)
        # We pick num features according to which minimizes the value most (+1 as we start with two features)
        n_features <- which.min(RMSEsvm)+1
        gn4 <- rf_rank[c(1:n_features)]

        classtrain <- train.res
        train3 <- train.mat[, gn4]

        classtest <- res[externalSample,]
        test3 <- dists[externalSample, gn4]

        train <- data.frame(class = classtrain, train3)
        model3 <- svm(class ~ ., data = train)

        fx <- predict(model3, test3)
        #print(fx)

        #results[externalSample] <- fx
        result_vector <- c(d, fx, classtest)
        #names(result_vector) <- c("sample", "predicted_age", "actual_age")
        print(result_vector)

    }

    #write(sample_id, file = "~/logfile.txt", append = T)

    print(paste("Sample_id:", sample_id, "- Time:", Sys.time()-start_time))
    return(result_vector)
}

# Select variables to look at
#dists = dists[, 1:1000]


p.val = {}
#for (i in 1:2){
for (i in 1:ncol(dists)){
    model <- lm(res[,1] ~ dists[,i])
    obj <- summary(model)
    p.val[i] <- obj$coefficients[2,4]
}
length(p.val)
p.adj <- p.adjust(p.val, n = length(p.val), method = "bonferroni")
selected <- which(p.adj < 0.0000000001)
print(length(which(p.adj < 0.0000000001)))

dists <- dists[, selected]

# system.time(run_loocv(parallel = T))
# system.time(run_loocv(parallel = F))

results <- run_loocv(parallel = T)
write.table(results, file = "prediction_results_loocv_10-100_generations.txt", row.names = F, quote = F, sep = "\t")
