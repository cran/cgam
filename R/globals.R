# This goes at the top of testpar-fit.R
#' @keywords internal
#' @name cgam-globalVariables
utils::globalVariables(c("object", "sse1", "etahat", "ahatc", "face", "qv", "edf", "edfs", "gcv", "gcvs", 
                         "sighat", "edfu_use", "covmat", "covmat_inv", "phi1", "mat1"))
`%||%` <- function(a, b) if (!is.null(a)) a else b
