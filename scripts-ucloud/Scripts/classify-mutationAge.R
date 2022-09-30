library(e1071)
library(randomForest)
library(Metrics)
library(tsne)
library(data.table)
library(caret)

abacus = F
if (abacus){
    setwd("/work/sduvarcall/haplotype-project")
    # First dataset
    #training_dists <- read.table(file = "training-distances-generations_10-100.txt", header = F)
    #res <- read.table(file = "res.txt", header = F)

    # Second dataset
    training_dists <- read.table(file = "training-distances-generations_10-500.txt", header = F)
    test_dists <- read.table("test-distances-generations_10-500.txt", header = F)
    res <- read.table(file = "results_10-500.txt", header = F)
} else {
    
    # Similarity Matrix
    # 10-100-step1 - 100 samples
    sims = 10
    samples = 100
    gen_start_tr = 10
    gen_stop_tr = 100
    gen_step_tr = 1
    gen_start_te = gen_start_tr
    gen_stop_te = gen_stop_tr
    gen_step_te = gen_step_tr
    seed_training = 42
    seed_test = 1024
    
    # 10-100-step1 - 20 samples
    sims = 20
    samples = 20
    gen_start_tr = 10
    gen_stop_tr = 100
    gen_step_tr = 10
    gen_start_te = gen_start_tr
    gen_stop_te = gen_stop_tr
    gen_step_te = gen_step_tr
    seed_training = 42
    seed_test = 24
    
    # training: 10-520-step10, test: 10-500-step10
    sims = 20
    samples = 100
    gen_start_tr = 10
    gen_stop_tr = 520
    gen_step_tr = 10
    gen_start_te = gen_start_tr
    gen_stop_te = 500
    gen_step_te = gen_step_tr
    seed_training = 31
    seed_test = 24
    
    # training: 10-550-step10, test: 10-500-step10
    sims = 20
    samples = 100
    gen_start_tr = 10
    gen_stop_tr = 550
    gen_step_tr = 10
    gen_start_te = gen_start_tr
    gen_stop_te = 500
    gen_step_te = gen_step_tr
    seed_training = 9012
    seed_test = 24
    
    # 10-500-step10 - 100 samples, 20 sims
    sims = 20
    samples = 100
    gen_start_tr = 10
    gen_stop_tr = 500
    gen_step_tr = 10
    gen_start_te = gen_start_tr
    gen_stop_te = gen_stop_tr
    gen_step_te = gen_step_tr
    seed_training = 42
    seed_test = 24
    
    # training: 10-1000-step10, test: 10-500-step10
    sims = 20
    samples = 100
    gen_start_tr = 10
    gen_stop_tr = 1000
    gen_step_tr = 10
    gen_start_te = gen_start_tr
    gen_stop_te = 500
    gen_step_te = gen_step_tr
    seed_training = 42
    seed_test = 24
    
    training_dists <- as.data.frame(data.table::fread(paste0("classify-simulated-population-age/",
                                  "BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start_tr,"_",gen_stop_tr,"_step",gen_step_tr,"-seed_",seed_training,"/",
                                  "classification_data/",
                                  "Distances-BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start_tr,"_",gen_stop_tr,"_step",gen_step_tr,"-seed_",seed_training,".txt"), 
                                  header = F))
    test_dists <- as.data.frame(data.table::fread(paste0("classify-simulated-population-age/",
                                    "BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start_te,"_",gen_stop_te,"_step",gen_step_te,"-seed_",seed_test,"/",
                                    "classification_data/",
                                    "Distances-BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start_te,"_",gen_stop_te,"_step",gen_step_te,"-seed_",seed_test,".txt"),
                             header = F))
    res <- read.table(file = paste0("classify-simulated-population-age/",
                                    "BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start_tr,"_",gen_stop_tr,"_step",gen_step_tr,"-seed_",seed_training,"/",
                                    "classification_data/results-generations_",gen_start_tr,"_",gen_stop_tr,"_step",gen_step_tr,".txt"), header = F)[,1]
    res2 <- read.table(file = paste0("classify-simulated-population-age/",
                                    "BRCA1_simulated-starGenealogy-",samples,"_samples-",sims,"_simulations-generations_",gen_start_te,"_",gen_stop_te,"_step",gen_step_te,"-seed_",seed_test,"/",
                                    "classification_data/results-generations_",gen_start_te,"_",gen_stop_te,"_step",gen_step_te,".txt"), header = F)[,1]

}

run_loocv <- function(parallel=T){
    # Select variables to look at
    #training_dists = training_dists[, 1:1000]
    p.val = {}
    #for (i in 1:2){
    for (i in 1:ncol(training_dists)){
        model <- lm(res[,1] ~ training_dists[,i])
        obj <- summary(model)
        p.val[i] <- obj$coefficients[2,4]
    }
    length(p.val)
    #p.adj <- p.adjust(p.val, n = length(p.val), method = "bonferroni") # bonferroni-holm (BH) better version of bonferroni!
    p.adj <- p.adjust(p.val, n = length(p.val), method = "fdr") # fdr = BH
    selected <-  which(p.adj < 0.000000000000001)
    #selected <-  which(p.adj < 0.000000002)
    print(length(which(p.adj < 0.000000000000001)))
    
    training_dists <- training_dists[, selected]
    dim(training_dists)
    
    if (parallel){
        library(parallel)
        cl <- makeCluster(detectCores(), outfile="progress.txt")
        #cl <- makeCluster(2, outfile="~/log.txt")
        clusterExport(cl, varlist = c("training_dists", "res"), envir = .GlobalEnv)
        #results <- as.data.frame(t(parSapply(cl = cl, X = 1:8, FUN = loocv)))
        results <- as.data.frame(t(parSapply(cl = cl, X = 1:nrow(training_dists), FUN = loocv)))
        stopCluster(cl)
    } else {
        #results <- as.data.frame(t(sapply(1:nrow(training_dists), loocv)))
        results <- as.data.frame(t(sapply(1:32, loocv)))
    }

    names(results) <- c("sample", "predicted_age", "actual_age")
    print(results)
    write.table(results, file = "prediction_results_loocv_10-100_generations.txt", row.names = F, quote = F, sep = "\t")
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

    all <- 1:nrow(training_dists)
    nvar <- ncol(training_dists)
    #for (d in 1:nrow(training_dists)){
    for (d in sample_id:sample_id){

        externalSample <- d
        inpairs <- all[-d]

        train.mat <- training_dists[inpairs, ]
        train.res <- res[inpairs,]

        ### Random Forest
        #gna <- c(1:nvar)

        squared_nvar <- sqrt(nvar)

        train <- data.frame(class = train.res, p=train.mat)
        model <- randomForest(class ~ ., data = train, mtry=squared_nvar, importance = T, na.action = na.omit)#, ntrees = 5000)

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
            for (a in 1:(nrow(training_dists)-1)){
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
        test3 <- training_dists[externalSample, gn4]

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


results <- run_loocv(parallel = T)


result <- {}
#test_dists2 <- training_dists[, selected]
test_dists2 <- training_dists
all <- 1:nrow(test_dists2)
d = 50
for (d in 1:nrow(test_dists2)){
    externalSample <- d
    inpairs <- all[-d]

    #test_dists2 <- training_dists[,sample(1:ncol(training_dists), 1000)]

    train.mat <- test_dists2[inpairs, ]
    train.res <- res[inpairs,]
    train <- data.frame(class = train.res, p=train.mat)
    #mod <- lm(class ~ ., data = train)
    mod <- svm(class ~ ., data = train, kernel = "radial")
    #mod <- svm(class ~ ., data = train, kernel = "linear")
    #mod <- svm(class ~ ., data = train, kernel = "polynomial")
    #mod <- svm(class ~ ., data = train, kernel = "sigmoid")
    result[d] <- predict(mod, newdata = data.frame(p=test_dists2[d,]))
    print(result[d])
    print(res[d])
}
plot(res,result)
cor.test(result,res)
summary(mod)


# Rearrange columns
sums = apply(training_dists, 2, sum)
names(sums) = 1:length(sums)
sums2 = sort(sums)
training_dists2 = training_dists[,as.integer(names(sums2))]

sums = apply(test_dists, 2, sum)
names(sums) = 1:length(sums)
sums2 = sort(sums)
test_dists2 = test_dists[,as.integer(names(sums2))]

training_dists2 = training_dists
test_dists2 = test_dists
# training_dists2 = test_dists
# test_dists2 = training_dists

set.seed(31)

# Real data
dst <- firstBreakDist(fam=T, matched_custom = matched_plink2)
test_dists2 = data.frame(t(as.vector(dst)))
names(test_dists2) = paste0("V", 1:length(test_dists2))
train = data.frame(class = res, p = training_dists2)#[,sample(1:ncol(training_dists), length(test_dists2))])
names(train) = c("class", paste0("p.V", 1:length(test_dists2)))

# Testing intervals
res.bak = res; res2.bak = res2
# 10-100
training_dists2 = training_dists[1:200,]
test_dists2 = test_dists[1:200,]
res = res.bak; res2 = res2.bak
res = res[1:200]
res2 = res2[1:200]

# 110-200
training_dists2 = training_dists[201:400,]
test_dists2 = test_dists[201:400,]
res = res.bak; res2 = res2.bak
res = res[201:400]
res2 = res2[201:400]

# 210-300
training_dists2 = training_dists[401:600,]
test_dists2 = test_dists[401:600,]
res = res.bak; res2 = res2.bak
res = res[401:600]
res2 = res2[401:600]

# 310-400
training_dists2 = training_dists[601:800,]
test_dists2 = test_dists[601:800,]
res = res.bak; res2 = res2.bak
res = res[601:800]
res2 = res2[601:800]

# 410-500
training_dists2 = training_dists[801:1000,]
#test_dists2 = test_dists[801:1000,]
test_dists2 = training_dists[801:1000,]
res = res.bak; res2 = res2.bak
res = res[801:1000]
res2 = res2[801:1000]



# Training data
#train = data.frame(class = res[,1], p = training_dists[,sample(1:ncol(training_dists), 190)])
train = data.frame(class = res, p = training_dists2)
#train = data.frame(class = res, p = training_dists)

# Model
model <- svm(class ~ ., data = train, kernel = "radial")
#model <- svm(class ~ ., data = train, kernel = "p", degree = 2)
#model <- lm(class  ~ ., data = train)
model <- randomForest(class ~ ., data=train, mtry=round(sqrt(ncol(train)-1)), importance=TRUE, na.action=na.omit)
#model <- randomForest(class ~ ., data=train, mtry=99, importance=TRUE, na.action=na.omit)
#model <- pcr(class ~ ., data=train, validation = "CV")

# Predictions
#predictions <- predict(model, newdata = data.frame(p=test_dists))
#predictions <- predict(model, newdata = data.frame(p=test_dists2), ncomp = 10) # pcr
predictions <- predict(model, newdata = data.frame(p=test_dists2))
outcome <- data.frame(predicted=predictions, actual=res2)

# Plot outcome
plot(outcome[,2], outcome[,1])
abline(a = 1, b = 1, col = "red")

# Correlation test
cor.test(outcome[,2], outcome[,1])
# Average difference
rmse(outcome[,2], outcome[,1])

# Boxplot
boxplot(outcome[,2], outcome[,1])


# xgboost
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(42)
model_xgb <- caret::train(class ~ .,
                          data = train,
                          method = "xgbTree",
                          preProcess = c("scale", "center"),
                          trControl = trainControl(method = "repeatedcv", 
                                                   number = 5, 
                                                   repeats = 3, 
                                                   savePredictions = TRUE, 
                                                   verboseIter = FALSE))
model_xgb

set.seed(42)
model_rf <- caret::train(class ~ .,
                         data = train,
                         method = "rf",
                         preProcess = c("scale", "center"),
                         trControl = trainControl(method = "repeatedcv", 
                                                  number = 5, 
                                                  repeats = 3, 
                                                  savePredictions = TRUE, 
                                                  verboseIter = FALSE))
model_rf

set.seed(42)
model_rf <- caret::train(class ~ .,
                         data = train,
                         method = "rf",
                         tunegrid = expand.grid(.mtry=99))

set.seed(42)
model <- caret::train(class ~ .,
                         data = train,
                         method = "rf",
                         preProcess = c("scale", "center"),
                         trControl = trainControl(method = "repeatedcv", 
                                                  number = 5, 
                                                  repeats = 3, 
                                                  savePredictions = TRUE, 
                                                  verboseIter = FALSE))