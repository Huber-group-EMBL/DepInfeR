# Documentation of datasets''


#' drug_response_GDSC
#'
#' This is the processed Genomics of Drug Sensitivity in Cancer (GDSC) drug sensitivity dataset.
#' The raw dataset was downloaded from \url{https://www.cancerrxgene.org/downloads/bulk_download}. 
#' The post-processing steps can be found at: \url{https://www.huber.embl.de/users/jlu/depInfeR/process_GDSC.html}.
#'
#' @docType data
#' @usage data(drug_response_GDSC)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(drug_response_GDSC)
"drug_response_GDSC"


#' mutation_GDSC
#'
#' This cancer type and genomic background annotation for cancer cell lines, 
#' use for the analysis of the GDSC dataset in the package vignette.
#' The raw dataset was downloaded from \url{https://www.cancerrxgene.org/downloads/bulk_download}. 
#' The post-processing steps can be found at: \url{https://www.huber.embl.de/users/jlu/depInfeR/process_GDSC.html}.
#'
#' @docType data
#' @usage data(mutation_GDSC)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(mutation_GDSC)
"mutation_GDSC"


#' targetsGDSC
#'
#' This drug-protein affinity profiling data for the analysis of the GDSC dataset - 
#' a subset of the data provided by Klaeger et al. 2017.
#' The raw data can be found in the supplementary file of the paper (Table_S1 & Table_S2): \url{https://science.sciencemag.org/content/358/6367/eaan4368/tab-figures-data}.
#' The post-processing steps can be found at: \url{https://www.huber.embl.de/users/jlu/depInfeR/process_kinobeads.html}.
#' 
#' @docType data
#' @usage data(targetsGDSC)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(targetsGDSC)
"targetsGDSC"

#' targetMatrix
#'
#' A toy data set that contains drug-target affinity matrix for examples and test of processTarget function. 
#' Rows contain drugs and columns contain targets.
#'
#' @docType data
#' @usage data(targetMatrix)
#' @format an object of matrix
#' @examples
#' data(targetMatrix)
"targetMatrix"



#' responseInput
#'
#' A toy data set that contains processed drug response matrix for examples and test of runLASSORegression function. 
#' Rows contain drugs and columns contain samples.
#'
#' @docType data
#' @usage data(responseInput)
#' @format an object of matrix
#' @examples
#' data(responseInput)
"responseInput"

#' targetInput
#'
#' A toy data set that contains processed drug-target affinity matrix for examples and test of runLASSOregression function. 
#' Rows contain drugs and columns contain targets.
#'
#' @docType data
#' @usage data(targetInput)
#' @format an object of matrix
#' @examples
#' data(targetInput)
"targetInput"


