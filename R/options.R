.onLoad <- function(libname, pkgname) {
  op <- options()
  op.cgam <- list(
    cgam.parallel = FALSE,
    cgam.cores = max(1, parallel::detectCores(logical = FALSE) - 1)
  )
  toset <- !(names(op.cgam) %in% names(op))
  if (any(toset)) options(op.cgam[toset])
  invisible()
}

.get_cgam_option <- function(name, default = NULL) {
  getOption(paste0("cgam.", name), default)
}
