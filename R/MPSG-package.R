#' MPSG: Meta-Analysis Model for Pleiotropy Selection with Group Structure
#'
#' Implementation of the method described in Sugier et al. (2025),
#' Biostatistics, 26(1), kxaf037, doi:10.1093/biostatistics/kxaf037.
#'
#' @keywords internal
#' @useDynLib MPSG, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
"_PACKAGE"
