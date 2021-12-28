# Documentation of datasets''


#' drug_response_GDSC
#'
#' This is the Genomics of Drug Sensitivity in Cancer (GDSC) drug sensitivity dataset, downloaded from https://www.cancerrxgene.org/.
#'
#' @docType data
#' @usage data(drug_response_GDSC)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(drug_response_GDSC)
"drug_response_GDSC"


#' targetsGDSC
#'
#' This drug-protein affinity profiling data for the analysis of the GDSC dataset - a subset of the data provided by Klaeger et al. 2017.
#'
#' @docType data
#' @usage data(targetsGDSC)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(targetsGDSC)
"targetsGDSC"


#' mutation_GDSC
#'
#' This cancer type and genomic background annotation for cancer cell lines, use for the analysis of the GDSC dataset in the package vignetee.
#'
#' @docType data
#' @usage data(mutation_GDSC)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(mutation_GDSC)
"mutation_GDSC"


#' responseInput
#'
#' Processed drug response matrix for examples and test of runLASSOregression
#'
#' @docType data
#' @usage data(responseInput)
#' @format an object of matrix
#' @examples
#' data(responseInput)
"responseInput"

#' targetInput
#'
#' Processed drug-target affinity matrix for examples and test of runLASSOregression
#'
#' @docType data
#' @usage data(targetInput)
#' @format an object of matrix
#' @examples
#' data(targetInput)
"targetInput"


#' targetMatrix
#'
#' Drug-target affinity matrix
#'
#' @docType data
#' @usage data(targetMatrix)
#' @format an object of matrix
#' @examples
#' data(targetMatrix)
"targetMatrix"

