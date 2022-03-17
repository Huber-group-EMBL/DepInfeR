#' Function for pre-processing drug-protein affinity dataset
#'
#' This function is used to preprocess the drug-protein affinity dataset
#' including the following steps:
#' - log-transform kd values (KdAsInput = TRUE)
#' - arctan-transform log(kd) values (KdAsInput = TRUE)
#' - check target similarity and remove highly correlated proteins
#' (removeCorrelated = TRUE)
#' - specify targets that should be kept in the matrix (keepTargets = NULL)
#'
#' All steps within this function are optional depending on input data.
#' The transformation steps should be performed
#' if the affinity matrix consists of kd values.
#' If there are highly correlated features within the affinity matrix,
#' they can be removed using the provided function.
#'
#' @param targetsMat Drug-protein affinity matrix with kd values (or optionally other
#' affinity measurement values at roughly normal distribution). Rows should contain drugs and columns should contain targets.
#' @param KdAsInput A boolean value indicating whether the drug-protein
#' affinity matrix contains kd values which should be log- and arctan-transformed.
#' The default value is TRUE.
#' @param removeCorrelated A boolean value indicating whether highly
#' correlated proteins should be summarized in target groups. The default value is TRUE.
#' @param keepTargets  A character variable that specifies important proteins
#' that should be retained in the matrix.
#' @param cutoff A similarity cutoff value for clustering proteins into one target group.
#' @import matrixStats rlist tibble
#' @export
#' @return A list of two element: 1)\code{targetMatrix} Pre-processed drug-protein
#' affinity matrix; 2)\code{targetCluster}, a list that contains the targets
#' show high correlations with each other.
#'
#' @examples
#' data(targetMatrix)
#' processTarget(targetsMat = targetMatrix, KdAsInput = TRUE , removeCorrelated = TRUE)
#'

processTarget <- function(targetsMat, KdAsInput = TRUE, removeCorrelated = TRUE,
                          keepTargets = NULL, cutoff=0.8) {
    if (KdAsInput) {
        targetsMat <- -log10(targetsMat) #log transform kd values
        targetsMat[is.na(targetsMat)] <- -10 #fill NA values with a very small pKd value (-10)
        arcTrans <- function(x,b = 2, g = 1) {     #define arctan function
            y <- (atan((x + b) * g) + pi/2)/pi
        }

        targetsMat <- arcTrans(targetsMat, b = 2, g = 3) #apply arctan transformation
    }

    if (removeCorrelated) {
        cosineSimi <- function(x) {
            x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))
        }
        simiMat <- cosineSimi(t(targetsMat)) #save similarity matrix
        #remove highly correlated features
        #allow manually specify important proteins to keep
        #sort target matrix per target importance to keep most important proteins in the dataset
        targetOrder <- colnames(targetsMat)[order(colSums(-targetsMat))]

        if (!is.null(keepTargets)) { #manually specify some important proteins
            targetOrder <- c(keepTargets, targetOrder[!targetOrder %in% keepTargets])
        }

        targetsMat <- targetsMat[,targetOrder]
        res <- removeCorrelatedTargets(targetsMat, cutoff = cutoff,
                                       distance = "cosine",
                                       cluster_method = "ward.D2")
        resTarMat <- res$reduced
        mapReduce_kd <- res$mapReduce
        ProcessTargetResults <- list(targetMatrix = resTarMat,
                                     targetCluster = mapReduce_kd)

    } else {
        ProcessTargetResults <- list(targetMatrix = targetsMat,
                                     targetCluster = NULL)
    }
    return(ProcessTargetResults)
}



#' Main function to run LASSO regression
#'
#' @param TargetMatrix Pre-processed drug-protein affinity matrix. Rows should contain drugs and columns should contain targets.
#' @param ResponseMatrix Pre-processed drug-response matrix. Rows should contain drugs and columns should contain samples.
#' @param cores A numeric variable specifying the number of cores. Multi-core parallelization may only work for Mac OS and Linux.
#' @param repeats A numeric variable specifying the number of regression repeats.
#' @return Pre-processed drug-protein affinity matrix
#' @export
#' @import glmnet rlang stats BiocParallel
#'
#' @examples
#' data(responseInput) #load drug response matrix
#' data(targetInput) #load drug-target affinity matrix
#' runLASSOregression(TargetMatrix = targetInput, ResponseMatrix = responseInput, repeats = 5)
#' # Please refer to the package vignette for more detailed information about this function.
#'
#'
runLASSOregression <- function(TargetMatrix, ResponseMatrix, cores = 1, repeats = 100) {
  
  #function for multi-target LASSO with repeated cross-validation
  runGlm.multi.para <- function(i, X, y, folds=3, lambda = "lambda.min",
                                standardize = FALSE) {
    res <- cv.glmnet(X, y, family = "mgaussian",
                     nfolds = folds, alpha = 1,
                     standardize = standardize)
    res
  }
  
  #script for LASSO regression on saved input matrices
  multicoreParam <- MulticoreParam(workers = cores)
  
  allResults <- bplapply(seq(repeats), runGlm.multi.para, 
                         TargetMatrix, ResponseMatrix, 
                         BPPARAM = multicoreParam)
  
  #Run function for processing glm results
  processGlm(allResults, TargetMatrix, ResponseMatrix)
}
