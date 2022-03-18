# Remove highly correlated proteins, in terms of the cosine similarity of their
#drug binding profile, using hierarchical clustering
#
#This function removes highly correlated proteins from the drug-protein
#affinity matrix, based on hierarchical clustering, and only keeps on protein of
#the clusters of highly correlated proteins. Rows should contain drugs and columns should contain targets.
removeCorrelatedTargets <- function(x, cutoff = 0.8, cluster_method = "ward.D2") {
    # calculate distance
    # cosine similarity maybe preferred for sparse matrix
    cosineSimi <- function(x){
        x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
    }
    distMat <- stats::as.dist(1-cosineSimi(t(x)))
    
    #hierarchical clustering
    hc <- stats::hclust(distMat, method = cluster_method)
    clusters <- stats::cutree(hc, h = 1-cutoff)
    x.re <- x[,!duplicated(clusters)]

    #record the removed features
    mapList <- lapply(colnames(x.re), function(i) {
        members <- names(clusters[clusters == clusters[i]])
        members[members != i]
    })
    names(mapList) <- colnames(x.re)

    return(list(reduced = x.re,
                mapReduce = mapList))
}


# Function for processing glm results (called within main LASSO regression function)
processGlm <- function(results, X, y, lambda = "lambda.min") {
    modelList <- list()
    lambdaList <- rep(NA,length(results))
    varExplain.all <- rep(NA,length(results))
    varExplain.cv <- rep(NA,length(results))
    coefTab <- data.frame(a = rep(colnames(X), times = ncol(y)),
                          b = rep(colnames(y), each = ncol(X)))
    
    for (i in seq(length(results))) {
        res <- results[[i]]
        lambdaList[i] <- res[[lambda]]
        #process coefficient model
        coefModel <- glmnet::coef.glmnet(res, s = lambda) #remove intercept row
        coefModel <- as.vector(Reduce(cbind, coefModel)[-1,])
        coefTab[[paste0("r",i)]] <- coefModel
        #calculate variance explained
        y.pred <- predict(res, s = lambda, newx = X)
        varExp <- stats::cor(as.vector(y),as.vector(y.pred))^2
        varExplain.all[i] <- varExp
    }
    
    #generate result matrix, the sign of coefficient is reversed in order to let 
    #higher coefficient indicates higher target important. 
    resMat <- -as.matrix(coefTab[,!colnames(coefTab) %in% c("a","b")])
    coefTab[["freq"]] <- rowSums(!resMat == 0)/ncol(resMat)
    coefTab[["med"]] <- rowMedians(resMat)
    
    freqMat <- coefMat <- matrix(data=NA, nrow = ncol(X), ncol = ncol(y),
                                 dimnames = list(colnames(X),colnames(y)))
    
    for (eachCol in colnames(freqMat)) {
        freqMat[,eachCol] <- coefTab[coefTab$b %in% eachCol,]$freq
        coefMat[,eachCol] <- coefTab[coefTab$b %in% eachCol,]$med
    }
    
    list(coefMat = coefMat, freqMat = freqMat,
         lambdaList = lambdaList, varExplain.all = varExplain.all,
         inputX = X, inputY = y)
}