# Remove highly correlated proteins, in terms of the cosine similarity of their drug binding profile, using hierarchical clustering
#
#This function removes highly correlated proteins from the drug-protein affinity matrix, based on hierarchical clustering, and only keeps on protein of the clusters of highly correlated proteins.
removeCorrelatedTargets <- function(x, cutoff = 0.8, distance = "cosine", cluster_method = "ward.D2") {
  # calculate distance matrix
  if (distance == "binary") {
    #maybe also useful is the input is a sparse matrix
    distMat <- stats::dist(t(x), method = "binary")
  } else if (distance == "pearson") {
    #otherwise, using pearson correlation
    distMat <- stats::as.dist(1-cor(x))
  } else if (distance == "euclidean") {
    distMat <- stats::dist(t(x), method = "euclidean")
  } else if (distance == "cosine") {
    # cosine similarity maybe preferred for sparse matrix
    cosineSimi <- function(x){
      x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
    }
    distMat <- stats::as.dist(1-cosineSimi(t(x)))
  } else if (distance == "canberra") {
    distMat <- stats::as.dist(as.matrix(dist(t(x), method = "canberra"))/nrow(x))
  }

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
  #add identifiers to coefMat, change it to coefTab
  coefTab <- tibble(a = rep(colnames(X), times = ncol(y)),
                    b = rep(colnames(y), each = ncol(X)))

  for (i in seq(length(results))) {
    res <- results[[i]]
    lambdaList[i] <- res[[lambda]]

    #process coefficient model
    coefModel <- coef(res, s = lambda) #remove intercept row
    coefModel <- as.vector(Reduce(cbind, coefModel)[-1,])
    coefTab[[paste0("r",i)]] <- coefModel

    #calculate variance explained
    y.pred <- predict(res, s = lambda, newx = X)
    varExp <- stats::cor(as.vector(y),as.vector(y.pred))^2
    varExplain.all[i] <- varExp
  }

  #generate result matrix
  resMat <- -as.matrix(dplyr::select(coefTab, -a , -b))
  coefTab <- dplyr::mutate(coefTab,
                           freq = rowSums(!resMat == 0)/ncol(resMat),
                           med = rowMedians(resMat),
                           sd = apply(resMat, 1, sd))

  #selection frequency matrix
  freqMat <- dplyr::select(coefTab, a, b, freq) %>%
    tidyr::pivot_wider(names_from = b, values_from = freq) %>%
    data.frame() %>% remove_rownames() %>% column_to_rownames("a") %>% as.matrix()

  #mean coefficient matrix
  coefMat <- dplyr::select(coefTab, a, b, med) %>%
    tidyr::pivot_wider(names_from = b, values_from = med) %>%
    data.frame() %>% remove_rownames() %>% column_to_rownames("a") %>% as.matrix()

  list(coefMat = coefMat, freqMat = freqMat,
       lambdaList = lambdaList, varExplain.all = varExplain.all,
       inputX = X, inputY = y)
}

