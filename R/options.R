.onLoad <- function(libname, pkgname) {
  op <- options()
  ncores <- parallel::detectCores(logical = FALSE)
  
  op.cgam <- list(
    cgam.parallel = FALSE,
    cgam.cores = if (is.na(ncores)) 1L else 2L #max(1L, ncores - 1L)
  )

  toset <- !(names(op.cgam) %in% names(op))
  if (any(toset)) options(op.cgam[toset])
  invisible()
}

.get_cgam_option <- function(name, default = NULL) {
  getOption(paste0("cgam.", name), default)
}
