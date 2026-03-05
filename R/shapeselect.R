#-----------------------------------------------------------------------------------------
#' Control parameters for \code{shapeselect()}
#'
#' Creates a list of parameters used by \code{shapeselect()} to control
#' basis construction, information-criterion calibration, simulation settings,
#' parallel computing etc.
#'
#' @param knots_lst Optional list specifying knots for spline basis functions,
#'   one element per predictor. Use `NULL` to let \code{shapeselect()} (and its
#'   helper routines) choose defaults.
#' @param boundary.knots Optional list specifiying knots for numeric predictors.
#' @param space Character string indicating the knots being equally spaced for spline 
#'   bases in underlying fitting routines (e.g., `"E"`). 
#' @param nsim Integer; number of simulations used  when estimating 
#'  mean degrees of freedom terms used in variable and shape selection criteria CBIC.
#' @param nfolds Integer; number of folds used for cross-validation when
#'   \code{type = "classo"} (shrinkage selection). It is not in use for now.
#' @param eps Numeric tolerance used for detecting near-zero signal.
#' @param pnt Logical; whether to use penalized spline estimator.
#' @param multicore Logical/integer; parallel setting passed to \code{shapeselect()} via the
#'   `parallel` control entry. Defaults to \code{getOption("cgam.parallel")}.
#' @param a Optional numeric; calibration constant used in the modified CBIC penalty term.
#'   If `NULL`, \code{shapeselect()} chooses a value based on `n`, `p`, and the implied
#'   total number of basis functions.
#' @param maxedf Logical; if `TRUE`, enforce an effective degrees of freedom cap in internal
#'   forward-step routines (as implemented).
#' @param add_edfu Numeric; optional additive adjustment to the unconstrained EDF used in the
#'   information criterion calculations.
#' @param add_kts Integer; optional adjustment to the number of knots used in spline bases.
#'   If `NULL`, \code{shapeselect()} chooses a default depending on `pnt`.
#' @param shrinketa Logical; if `TRUE` and \code{type = "classo"}, shrink fitted component
#'   contributions by the selected shrinkage parameters. It is not in use for now.
#' @param catnrep Logical; if `TRUE`, print the iteration counter during the selection path.
#' @param standardize Logical; if `TRUE`, standardize numeric predictors by their standard deviations.
#' @param center Logical; if `TRUE`, center numeric predictors when scaling them.
#' @param verbose Logical; if `TRUE`, print the model-path tracking matrix as the algorithm runs.
#' @param set4 Logical; optional flag passed through to internal routines used to restrict or
#'   modify the candidate shape set to be increasing, decreasing, convex and concave; if false, then
#'   combinations of monotonicity and convexity will be added
#'
#' @return A named list suitable for the \code{control} argument of \code{shapeselect()}.
#'
#' @examples
#' ctrl <- shapeselect_control(nsim = 100, nfolds = 10, verbose = TRUE)
#' ctrl
#'
#' @export
#' 
shapeselect_control <- function(knots_lst = NULL, boundary.knots = NULL, space = "E", nsim = 200, nfolds = 5, 
                                eps = 1e-8, pnt = TRUE, multicore = getOption("cgam.parallel"), 
                                a = NULL, maxedf = FALSE, add_edfu = 0, add_kts = NULL,
                                shrinketa = FALSE, catnrep = FALSE, standardize = TRUE,
                                center = TRUE, verbose = FALSE, set4 = TRUE)
{
  rslt <- list(knots_lst = knots_lst, boundary.knots = boundary.knots,
               space = space, nsim = nsim, nfolds = nfolds, 
               eps = eps, pnt = pnt, parallel = multicore,
               a = a, maxedf = maxedf, add_edfu = add_edfu, add_kts = add_kts,
               shrinketa = shrinketa, catnrep = catnrep, standardize = standardize,
               center = center, verbose = verbose, set4 = set4)
  return (rslt)
}

#' Shape selection for additive components with shape or order constraints
#'
#' Fits an additive model by selecting variables with shape or order constraints
#' for each predictor using Cone Information Bayesian Criterion (CBIC).
#' The function can perform forward selection or stepwise (forward/backward) updates
#' over a set of candidate shapes appropriate to each predictor type.
#'
#' @param xmat A data frame of predictors. Columns may be numeric or factors.
#'   Column names are used in output when available. If x is numeric, then the shape options are 
#'   increasing, decreasing, convex, concave or combination of monotonicity and convexity; if x is a 
#'   factor, then it can have monotonic or convex ordering options; if x is a factor, the option can also be "in" or "out-of" the model, or "in", "out-of" the model, or
#'   following a "tree" ordering, depending on the `types` parameter.
#' @param y A response vector of length `nrow(xmat)`. It can be either numeric or binary.
#' @param types Optional character vector of length `ncol(xmat)` specifying the predictor type
#'   for each column of `xmat`. This is used to determine which shape set is eligible for each
#'   predictor. Values used for this parameter include `"num"`/`"numeric"` (continuous),
#'   `"ord"`/`"ordinal"` (ordinal/categorical), `"tree"` (categorical),
#'  and `"nom"`/`"nominal"` (categorical). Default is `NULL`, but the user must provide a value for each column of xmat;
#'  otherwise the function will throw an error message.
#' @param direction Direction of the search. `"forward"` performs forward selection only.
#'   `"both"` allows backward revisions (stepwise) after additions.
#' @param type Search strategy. `"stepwise"` adds predictors while allowing backward revisions.
#'   `"classo"` then applies shrinkage/selection (e.g., via cross-validated lasso). The available option 
#'   is `"stepwise"` for now.
#' @param infocrit Information criterion used to decide adding or removing a predictor, and shape updates. Currently supports
#'   `"CBIC"` and `"CAIC"`.
#' @param family A GLM family. May be a family function or a character string naming one. Available families are `gaussian` and `binomial`.
#' @param control A list of control parameters. Any entries not supplied fall back to
#'   `shapeselect_control()`. Additional named arguments in `...` override both defaults and
#'   `control`.
#' @param ... Additional control parameters (passed as overrides). These are merged into `control`
#'   after defaults, so they take highest priority.
#'
#' @return An object of class `c("shapeselect","cgam")`, a list containing (among others):
#' \describe{
#'   \item{mod_tried}{A matrix recording the selection path (chosen shapes, criterion values, adding or removing a predictor).}
#'   \item{shps}{Numeric shape codes for each predictor (0 indicates excluded/flat).}
#'   \item{shapes}{Character labels for the selected shapes.}
#'   \item{etahat}{Fitted linear predictor values.}
#'   \item{muhat}{Fitted mean response on the data scale (`family$linkinv(etahat)`).}
#'   \item{etacomps}{Matrix of fitted component contributions by predictor (rows).}
#'   \item{cic_vec}{Criterion values along the path.}
#'   \item{cvec}{Shrinkage/selection vector (used when `type = "classo"`).}
#'   \item{ddkeep_x_lst, ddkeep_z_lst}{Stored design/basis components for numeric and categorical predictors used for prediction and plotting.}
#' }
#'
#' @details
#' The algorithm begins with an intercept-only fit, then evaluates each inactive predictor under
#' a set of allowable shapes (depending on `types`). At each step it selects the predictor/shape
#' combination that improves the chosen information criterion. When `direction = "both"`, it
#' checks non-chosen shapes for already-selected predictors via backward steps.
#'
#' Shape and order choices in this implementation include 
#' 1) smooth effects: monotone increasing/decreasing, convex/concave when `set4 = TRUE`, and combinations `set4 = FALSE`, if the corresponding value 
#' in `types` is `"num"`/`"numeric"` (x is continuous).
#' 2) non-smooth effects: monotone increasing/decreasing and convex/concave (when `set4 = TRUE`), or their combinations (when `set4 = FALSE`), 
#'    if the corresponding value in `types` is `"ord"`/`"ordinal"` (i.e., x is ordinal, categorical, or numeric but treated as 
#'    discrete with an inherent ordering).
#' in `types` is `"ord"`/`"ordinal"` (x is ordinal/categorical/having numeric values but treated as discrete values with ordering).
#' 3) in or out-of the model, if the corresponding value in `types` is `"nom"`/`"nominal"` (x is categorical)
#' 4) in, out-of the model or following a tree ordering, if the corresponding value in `types` is `"tree"` (i.e., x is categorical or numeric but treated as 
#'    discrete with a tree ordering).
#' @examples
#' library(MASS)
#' data("Rubber", package = "MASS")
#' y <- Rubber[, 1]
#' xmat <- Rubber[, -1]
#' ans <- shapeselect(xmat, y, types = c("num", "num"), pnt = FALSE)
#' 
#' print (ans$shapes)
#' 
#' plot(ans)
#' @export

#-----------------------------------------------------------------------------------------
#coneA is used; one is included in bigmat
#main routine: shapeselect
#xmat must be a data frame
#tree: the first unique value of x is placebo
#add tree ordering when the type of x is 'tree': shps: 0 (out), 1 (tree), 2 (no order)
#-----------------------------------------------------------------------------------------
shapeselect <- function(xmat, y, types = NULL, direction = c("both", "forward"),
                        type = "stepwise", infocrit = c("CBIC", "CAIC"), 
                        family = gaussian, control = list(),...)
{	
  cl <- match.call()
  type <- match.arg(type)
  infocrit <- match.arg(infocrit)
  direction <- match.arg(direction)
  
  if (is.null(types)) stop("'types' must be provided and have length ncol(xmat).")
  
  defaults <- shapeselect_control()
  extras <- list(...)
  control <- modifyList(defaults, control)
  control <- modifyList(control, extras)
  
  #unpack the control list
  knots_lst <- control$knots_lst
  boundary.knots <- control$boundary.knots
  space <- control$space
  nsim <- control$nsim
  nfolds <- control$nfolds
  eps <- control$eps
  pnt <- control$pnt
  parallel <- control$parallel
  a <- control$a
  maxedf <- control$maxedf
  add_edfu <- control$add_edfu
  add_kts <- control$add_kts
  shrinketa <- control$shrinketa
  catnrep <- control$catnrep
  standardize <- control$standardize
  center <- control$center
  verbose <- control$verbose
  set4 <- control$set4
  
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family))
    stop("'family' not recognized!")
  
  #new
  cicfamily <- CicFamily(family)
  llh.fun <- cicfamily$llh.fun
  
  fml <- family$family
  if(fml == "gaussian"){
    wt.iter <- FALSE
  } else {wt.iter <- TRUE}
  xnms <- colnames(xmat)
  if(is.null(xnms)){xnms <- paste('x', 1:NCOL(xmat), sep='')}
  #new: use types to check if x is for spline or ordinal basis or just nominal
  
  #use 18 for linear shape
  #use 17 for flat
  #shps <- 9:18
  shps <- 9:17
  #shps_char <- c('s.incr','s.decr','s.conv','s.conc','s.incr.conv','s.incr.conc','s.decr.conv','s.decr.conc', 'flat', 'linear')
  shps_char <- c('s.incr','s.decr','s.conv','s.conc','s.incr.conv','s.incr.conc','s.decr.conv','s.decr.conc', 'flat')
  mod_0 <- data.frame(x = character(), 'shape/order' = numeric(), cic = numeric(), check.names = FALSE)
  n <- length(y)
  cic_vec <- NULL
  active <- NULL
  x_rm <- NULL #keep track of redundant z predictors
  im <- inactive <- 1:NCOL(xmat)
  #--------------------------------
  #step 0: project y onto 1
  #--------------------------------
  one <- rep(1, n)
  
  xmat0 <- xmat #used xmat0 later in plot
  #new
  #xmat[] <- lapply(xmat, function(col) if(is.numeric(col)) scale(col, center = FALSE) else col)
  p <- NCOL(xmat)
  n <- NROW(xmat)
  
  sdx <- rep(1, p)
  mx <- rep(0, p)
  #new: new rule for a
  mult <- ifelse(pnt, 2, 1)
  m_total <- 0
  
  for(ix in 1:p){
    if(is.numeric(xmat[,ix]) & types[ix] %in% c("num", "numeric")){
      #if(is.numeric(xmat[,ix])){
      #xmat[,i] <- (xmat[,i] - min(xmat[,i])) / (max(xmat[,i]) - min(xmat[,i]))
      
      kts0 <- NULL
      if(is.list(boundary.knots)){
        kts0 <- boundary.knots[[ix]]
        if (!is.numeric(kts0) | any(!is.finite(kts0)) | length(kts0) < 2){
          warning(
            "Invalid 'boundary.knots': must be a numeric vector of at least two finite values. ",
            "Using range(x) instead.",
            call. = FALSE
          )
          bounndary.knots[ix] <- list(NULL)
          kts0 <- NULL
        }
      }
      
      #new: check scale first
      if(standardize){
        xmat[,ix] <- scale(xmat[,ix], center=center)
        #sdx <- c(sdx, attr(xmat[,ix], "scaled:scale"))
        sdx[ix] <- attr(xmat[,ix], "scaled:scale")
        if(center){
          mx[ix] <- attr(xmat[, ix], 'scaled:center')
        }
        
        #new: scale user-defined boundary knots
        if(length(kts0) >= 2){
          if(center){
            kts0 <- kts0 - mx[ix]
          }
          kts0 <- kts0 / sdx[ix]
          
          boundary.knots[[ix]] <- kts0
        }
      }
      #new:
      n1 <- length(unique(xmat[,ix]))
      m_total <- m_total + mult * (trunc(n1^(1/7)) + 6)
    }
    if(inherits(xmat[,ix], "factor") & types[ix] %in% "tree"){
      xmat[,ix] <- as.numeric(as.character(xmat[,ix]))
    }
    if(inherits(xmat[,ix], "factor") & types[ix] %in% "ord"){
      xmat[,ix] <- as.numeric(as.character(xmat[,ix]))
    }
    attr(xmat[, ix], "type") <- types[ix]
  }
  
  # for(ix in 1:p){
  #   if(is.numeric(xmat[,ix]) & types[ix] %in% c("num", "numeric")){
  #     #if(is.numeric(xmat[,ix])){
  #     #xmat[,i] <- (xmat[,i] - min(xmat[,i])) / (max(xmat[,i]) - min(xmat[,i]))
  #     xmat[,ix] <- scale(xmat[,ix], center=center)
  #     #sdx <- c(sdx, attr(xmat[,ix], "scaled:scale"))
  #     sdx[ix] <- attr(xmat[,ix], "scaled:scale")
  #     if(center){
  #       mx[ix] <- attr(xmat[, ix], 'scaled:center')
  #     }
  #     #new:
  #     n1 <- length(unique(xmat[,ix]))
  #     m_total <- m_total + mult * (trunc(n1^(1/7)) + 6)
  #   }
  #   if(inherits(xmat[,ix], "factor") & types[ix] %in% "tree"){
  #     xmat[,ix] <- as.numeric(as.character(xmat[,ix]))
  #   }
  #   if(inherits(xmat[,ix], "factor") & types[ix] %in% "ord"){
  #     xmat[,ix] <- as.numeric(as.character(xmat[,ix]))
  #   }
  #   attr(xmat[, ix], "type") <- types[ix]
  # }
  
  #new rule for a
  if(is.null(a)){
    a <- 1
    # if(p >= n){
    #   a <- 2
    # }
    # #new: decide value of a according to the number of spline bases functions
    # if(m_total > n & n > p){
    #   a <- 1
    # }
    # if(n > m_total){
    #   a <- 0.5
    # }
  }
  
  #new rule for add_kts
  if(is.null(add_kts)){
    if(pnt){
      add_kts <- 0
    }else{
      add_kts <- -2
    }
  }
  
  #np will change later when including convex/concave x
  np <- 1
  zmat <- NULL
  etahat <- muhat <- mean(y)
  rss <- sum((y - etahat)^2)
  #need a larger value to start cic/bic?
  if(infocrit %in% c("CAIC", "CIC")) {
    if(!wt.iter){
      #test!
      cic <- log(rss) +  2*1/n
    } else {
      #cat('fml:', fml, '\n')
      llh <- llh.fun(y=y, muhat=rep(1/2, n), etahat=rep(0, n), n=n, fml=fml)
      cic <- llh + 2*1/n
    }
  }
  
  if(infocrit == "CBIC") {
    if(!wt.iter){
      #cic <- log(rss) + log(n) * 1 / n
      #new
      #cic <- log(rss) + a * log(n) * 1 / n^(6/7)
      degree <- 6/7#1
      cic <- log(rss) + a * log(n) * 1 / n^(degree)
    } else {
      m0 <- stats::glm(y ~ 1, family=family)
      llh <- llh.fun(y=y, muhat=rep(1/2, n), etahat=rep(0, n), n=n, fml=fml)
      #cic <- llh + log(n) * 1 / n
      #new
      #cic <- llh + a * log(n) * 1 / n^(6/7)
      degree <- 6/7#1
      cic <- llh + a * log(n) * 1 / n^(degree)
    }
  }
  
  cic0 <- oldcic <- cic
  cic_vec <- c(cic_vec, cic)
  mod_0[1, 1] <- "1 vector"
  mod_0[1, 2] <- "constant"
  mod_0[1, 3] <- cic 
  y0 <- y
  
  if(!wt.iter){
    mu <- mean(y)
  }
  
  #----------------------------------------------------------------------
  #step 1: regress y vs x1, y vs x2, ...., y vs xp
  #----------------------------------------------------------------------
  x_chosen <- shps_chosen <- NULL
  #keep track of how many spline bases have been chosen
  capm <- 0
  ddkeep_z <- ddkeep_z_all <- ddkeep_x_all <- ddkeep_x <- NULL
  #new: keep track of penalty matrix
  dmat_keep <- NULL 
  ps_keep <- 0
  capl <- NCOL(xmat) #capl includes z
  #8 shapes for continuous x; 2 shapes for z
  #include linear option
  #1 is not in the projection space for residual
  cic_st <- oldcic
  rslt_step1 <- forward_cic(xmat, y, cic_st, infocrit, type, knots_lst, boundary.knots, 
                            nsim, ddkeep_z_all, wt.iter, family, parallel, space, pnt, a,
                            maxedf, add_edfu, add_kts, set4)
  #rslt_step1 <- forward_cic(xmat, y, cic_st, infocrit, type, knots_lst, nsim, 
  #                          ddkeep_z_all, wt.iter, family, parallel, space, pnt, a,
  #                          maxedf, add_edfu, add_kts, set4)
  cic_vals <- rslt_step1$cic_vals
  #new: avoid wrong choice when e.g. conv and incr.conv has the diff of 1e-16
  cic_vals <- round(cic_vals, 8)
  #need to add a warning: is the smallest cic_val is for the flat shape for all predictors, then stop the algorithm
  ddkeep_x <- rslt_step1$ddkeep_x
  ddkeep_z <- rslt_step1$ddkeep_z
  #etahatkeep is cumulative estimate; it's not smooth component chosen at each step
  etahatkeep <- rslt_step1$etahatkeep
  coefskeep <- rslt_step1$coefskeep
  mskeep <- rslt_step1$mskeep
  oldetahat <- etahatkeep
  oldcoefs <- coefskeep
  ktskeep <- rslt_step1$ktskeep
  cic_st_next_step <- rslt_step1$cic_st_next_step
  dmat_keep <- rslt_step1$dmat_keep
  ps_keep <- rslt_step1$ps_keep
  
  #NA's will be for categorical predictors
  cic <- min(cic_vals, na.rm = TRUE)
  #cic_vals = cic_vals[-9, ,drop = FALSE] #no need to remove -9; include flat
  if(cic_st_next_step < oldcic){
    x_index <- which(cic_vals == cic, arr.ind = TRUE)
  } else {
    #flat for all predictors
    x_index <- which(cic_vals == round(oldcic, 8), arr.ind = TRUE)
  }
  x_index <- rbind(x_index)
  
  #if(NROW(x_index) > 1 & any(x_index[, 1] == 9)){ 
  if(NROW(x_index) > 1 & any(x_index[, 1] %in% c(5, 9))){ 
    #if(NROW(x_index) > 1 & any(x_index[, 1] == nrow(cic_vals))){ 
    #stop('check!')
    #more than 1 x's best shape is flat; should we add them all at once?
    x_index <- x_index[1, ]
  } else if(NROW(x_index) > 1 & length(unique(x_index[, 2])) == 1){
    #more than 1 shape of 1 x has the same cic_val; choose the one with smallest null dfmean
    #dfmean_cbic_vals <- rslt_step1$dfmean_cbic_vals
    col_ps <- unique(x_index[, 2])
    dfmean_x_chosen <- dfmean_cbic_vals[, col_ps]
    row_ps <- x_index[, 1]
    x_index <- x_index[which.min(dfmean_x_chosen[row_ps]), ]
  }
  
  x_chosen <- c(x_chosen, inactive[x_index[2]])
  #is_tree <- inherits(xmat_inactive[, x_index[2]], "factor") & attr(xmat_inactive[, x_index[2]], "type") %in% c("tree")
  if(!inherits(xmat[, x_index[2]], "factor")){
    if(attr(xmat[, x_index[2]], "type") %in% c("ord", "ordinal")){
      shps <- c(1:8, 0)
      if(wt.iter & infocrit == "CBIC"){
        shps <- c(1:4, 0)
      }
    } else if (attr(xmat[, x_index[2]], "type") %in% c("tree")) {
      shps <- c(1:2, 0)
    } else {
      shps <- 9:17
      if(wt.iter & infocrit == "CBIC"){
        shps <- c(9:12, 0)
      }
    }
    shp_ch_i <- shps[x_index[1]]
    shps_chosen <- c(shps_chosen, shps[x_index[1]])
  } else {
    if (cic == oldcic) {shp_ch_i <- 0} 
    if (cic < oldcic) {shp_ch_i <- 1}
    shps_chosen <- c(shps_chosen, shp_ch_i)
  }
  
  #record chosen model's cic, x, shape
  cic_vec <- c(cic_vec, cic)
  colnames(cic_vals) <- colnames(xmat)
  rownames(cic_vals) <- shps_char
  #----------------------------------------------------------------------------------------------------
  #step 2 to step k: regress y ~ x1 vs x2, ...., y ~ x1 vs xp
  #stopping criteria: p = 0; oldcic < cic
  #if classo = T, then stop when p = 0? go forward and backward to find the smallest CIC
  #----------------------------------------------------------------------------------------------------
  #mod_record_all <- list()
  nrep <- 0
  varlst <- varlst_z <- NULL #to keep track of delta's
  #varlst <- NULL
  #varlst_z <- 1 #start with one vector
  #new keep track of predictors tested 
  #1 if x is in the model; 0 otherwise
  totx <- NCOL(xmat)
  mod_tried <- matrix(0, nrow = 1000, ncol = (totx + 3))
  #colnames(mod_tried) = c(xnms, "# of x's", "cic", "new add")
  xnms_with_id <- paste0(xnms, '|', 1:totx)
  colnames(mod_tried) <- c(xnms_with_id, "# of x's", infocrit, "new add")#, "logrss", "edf0", "edf_max")
  
  #used to create bigmat in the order cgam uses
  is_z <- sapply(xmat, function(e) is.factor(e))
  obs <- 1:ncol(xmat)
  zid <- obs[is_z]
  xid <- obs[!is_z]
  #new: 
  mskeep_lst <- ddkeep_x_lst <- ddkeep_z_lst <- vector("list", length = totx)
  dmat_keep_lst <- vector("list", length = totx)
  ps_keep_vec <- rep(0, length = totx)
  if(is.null(knots_lst)){
    knots_lst <- vector("list", length = totx)
  }
  check1 <- cic_st_next_step < oldcic & length(inactive) > 0
  
  rm_id <- NULL
  x_rm <- NULL
  ansi_back <- NULL
  #new:
  coefskeep_x_lst <- vector("list", length = totx)
  coefskeep_z_lst <- vector("list", length = totx)
  varlst_all <- NULL
  elapsed_times <- NULL
  ncols <- NULL
  sttime <- proc.time()[3]
  while(check1){
    new_add <- inactive[x_index[2]]
    x_rm_nrep <- inactive[rm_id]
    x_rm <- c(x_rm, x_rm_nrep)
    inactive <- inactive[-c(x_index[2], rm_id)]
    active <- im[!im %in% inactive]
    
    nrep <- nrep + 1
    if(catnrep){
      cat('nrep: ', nrep, '\n')
    }
    end_time <- proc.time()[3]
    elapsed_time <- end_time - sttime
    #cat('elapsed time: ', elapsed_time, '\n')
    sttime <- end_time
    elapsed_times <- c(elapsed_times, elapsed_time)
    #new:
    if(nrep > 1) {
      mod_tried[nrep, ] <- mod_tried[nrep-1, ]
    }
    new_add_num <- id_ch <- new_add 
    #use she chosen shape for id_ch instead? 
    mod_tried[nrep, id_ch] <- rev(shps_chosen)[1]#1
    mod_tried[nrep, totx + 1] <- length(active)
    #------------
    #test more!
    #------------
    #mod_tried[nrep, totx + 2] <- cic
    mod_tried[nrep, totx + 2] = cic_st_next_step
    
    mod_tried[nrep, totx + 3] <- id_ch
    nac <- mod_tried[nrep, totx + 1] 
    
    #test more
    if(!is.null(rm_id)){
      nrep <- nrep + 1
      mod_tried[nrep, ] <- mod_tried[nrep-1, ]
      mod_tried[nrep, totx + 3] <- x_rm_nrep
    }
    
    if(verbose){
      formatted_df <- mod_tried[1:nrep, , drop = FALSE]
      cat(capture.output(print(formatted_df)), sep = "\n")
    }
    
    #etahat_all = etahat_all + etahatkeep
    ddkeep_z_all <- cbind(ddkeep_z_all, ddkeep_z)
    ddkeep_x_all <- cbind(ddkeep_x_all, ddkeep_x)#may not follow the order of x, may be 1,3,2.
    
    np <- NCOL(ddkeep_z_all) #no intercept in it
    capm <- NCOL(ddkeep_x_all)
    
    if(!is.null(ddkeep_x)) {
      ddkeep_x_lst[[new_add_num]] <- ddkeep_x
      mskeep_lst[[new_add_num]] <- mskeep
      if(!is.null(ktskeep)){
        knots_lst[[new_add_num]] <- ktskeep
      }else{
        knots_lst[new_add_num] <- list(NULL)
      }
      #new:
      dmat_keep_lst[[new_add_num]] <- dmat_keep
      ps_keep_vec[new_add_num] <- ps_keep
      #varlst not record z
      varlst <- c(varlst, rep(new_add_num, NCOL(ddkeep_x)))
      #new
      #varlst_all <- c(varlst_all, rep(new_add_num, NCOL(ddkeep_x)))
      #coefskeep_x_lst[[new_add_num]] <- (coefskeep[-1])[which(varlst_all %in% new_add_num)]
    }
    
    if(!is.null(ddkeep_z)) {
      ddkeep_z_lst[[new_add_num]] <- ddkeep_z
      mskeep_lst[[new_add_num]] <- mskeep
      varlst_z <- c(varlst_z, rep(new_add_num, NCOL(ddkeep_z)))
    }
    
    #new: go backwards to see if other shapes of previously added predictors will decrease cic
    if(direction == "both"){
      if(nac >= 2 & (!shp_ch_i %in% c(0, 17))){
        chs <- which(mod_tried[nrep-1, 1:totx] != 0)
        #in backward step, use smaller oldcic?
        oldcic <- cic_st_next_step
        #new: used to determin degree of n in CBIC
        shps_all <- mod_tried[nrep, 1:totx]
        #cat('start backward step: ', '\n')
        ansi_back <- backward_cic(xmat, y, varlst, varlst_z, ddkeep_x_all, ddkeep_x_lst, knots_lst,
                                  boundary.knots, 
                                  ddkeep_z_all, ddkeep_z_lst, mskeep_lst, mod_tried, nrep, totx,
                                  oldcic, infocrit, dfmean_old, dfmeanmax_old,
                                  penalty_old, oldcoefs, oldetahat, 
                                  active, type, nsim, wt.iter, family,
                                  parallel, space, pnt, dmat_keep_lst, ps_keep_vec,
                                  a, add_edfu, add_kts, shps_all)
        
        # ansi_back <- backward_cic(xmat, y, varlst, varlst_z, ddkeep_x_all, ddkeep_x_lst, knots_lst,
        #                           ddkeep_z_all, ddkeep_z_lst, mskeep_lst, mod_tried, nrep, totx,
        #                           oldcic, infocrit, dfmean_old, dfmeanmax_old,
        #                           penalty_old, oldcoefs, oldetahat, 
        #                           active, type, nsim, wt.iter, family,
        #                           parallel, space, pnt, dmat_keep_lst, ps_keep_vec,
        #                           a, add_edfu, add_kts, shps_all)
        #cat('finished backward step: ', '\n')
        if(ansi_back$changed){
          back_inactive <- ansi_back$back_inactive
          if(length(back_inactive) > 0) {
            inactive <- sort(c(inactive, back_inactive))
            active <- im[!im %in% inactive]
          }
          
          #create new ddkeep_x_all, ddkeep_z_all
          ddkeep_x_lst <- ansi_back$ddkeep_x_lst
          #new:
          coefskeep_x_lst <- ansi_back$coefskeep_x_lst
          coefskeep_z_lst <- ansi_back$coefskeep_z_lst
          #print (coefskeep_z_lst)
          #new:
          dmat_keep_lst <- ansi_back$dmat_keep_lst
          ps_keep_vec <- ansi_back$ps_keep_vec
          ddkeep_z_lst <- ansi_back$ddkeep_z_lst
          #print (ddkeep_z_lst)
          knots_lst <- ansi_back$knots_lst
          mod_tried <- ansi_back$mod_tried
          mskeep_lst <- ansi_back$mskeep_lst
          cic <- ansi_back$cic_st
          cic_st_next_step <- ansi_back$cic_st_next_step
          #xtxkeep <- ansi_back$xtx
          #cat('backward step xtx:', dim(xtxkeep), '\n')
          #face <- ansi_back$face
          
          etahatkeep <- ansi_back$etahatkeep
          coefskeep <- ansi_back$coefskeep
          oldcoefs <- coefskeep
          oldetahat <- etahatkeep
          
          seq_added <- mod_tried[1:nrep, totx + 3]
          
          # ddkeep_x_lst_ord will only store dd's for chosen x
          ddkeep_x_lst_ord <- ddkeep_x_lst[seq_added]
          ddkeep_z_lst_ord <- ddkeep_z_lst[seq_added]
          
          #keep the adding order of x, may not be 1,2,3; may be 1,3,2
          ddkeep_x_all <- do.call(base::cbind, ddkeep_x_lst_ord) 
          ddkeep_z_all <- do.call(base::cbind, ddkeep_z_lst_ord)
          
          #new:
          #dmat_keep_all <- as.matrix(bdiag(Filter(Negate(is.null), dmat_keep_lst)))
          
          #need to create varlst and varlst_z again
          seq_added <- mod_tried[1:nrep, totx + 3]
          
          obs <- 1:totx
          obs_ord <- obs[seq_added]
          
          if(!is.null(ddkeep_x_all)){
            varlst_lst <- lapply(1:length(obs_ord), function(i, xlst) {if(NCOL(xlst[[i]]) >= 1) {rep(obs_ord[i], NCOL(xlst[[i]]))}}, xlst = ddkeep_x_lst_ord)
            #print (varlst_lst)
            varlst <- do.call(base::c, varlst_lst)
          } else {varlst <- NULL}
          
          if(!is.null(ddkeep_z_all)){
            varlst_z_lst <- lapply(1:length(obs_ord), function(i, zlst) {if(NCOL(zlst[[i]]) >= 1) {rep(obs_ord[i], NCOL(zlst[[i]]))}}, zlst = ddkeep_z_lst_ord)
            varlst_z <- do.call(base::c, varlst_z_lst)
          } else {varlst_z <- NULL}
        }
      }
    }   
    
    xmat_inactive <- xmat[, -active, drop=FALSE]
    if(length(inactive) == 0) {break}
    
    #oldcic <- cic 
    #cic_st <- cic
    oldcic <- cic_st_next_step
    cic_st <- cic_st_next_step
    
    knots_lst_inactive <- NULL
    boundary.knots_inactive <- NULL
    
    if(!is.null(knots_lst)){
      knots_lst_inactive <- knots_lst[-active]
    }
    
    if(!is.null(boundary.knots)){
      boundary.knots_inactive <- boundary.knots[-active]
    }
    
    #cic_st should be the same as oldcic? 
    #cic_st is a little larger than cic_st_fixed_shps; it is cic val for 17; used in classo to make every x to be chosen??
    #new: use shapes chosen so far to determine the denominator of CBIC
    shps_all <- mod_tried[nrep, 1:totx]
    shps_non_zero <- shps_all[!shps_all %in% c(0, 17)]
    #new! dmat_keep_all should be of the same order as ddkeep_x_all
    #cannot use this: may add some x twice if the x was out of the model in the previous step: 2 5 3 1 6 2
    seq_added <- mod_tried[1:nrep, totx + 3]
    #test!
    seq_added <- seq_added[!duplicated(seq_added, fromLast = TRUE)]
    ansi <- forward_cic_step2(xmat_inactive, y, cic_st, ddkeep_z_all, ddkeep_x_all, knots_lst_inactive, 
                              boundary.knots_inactive,
                              oldcic, oldetahat, oldcoefs, infocrit, dfmean_old, dfmeanmax_old, penalty_old,
                              active, inactive, type, nsim, wt.iter, family, parallel, space, pnt, 
                              dmat_keep_lst, ps_keep_vec, a, maxedf, add_edfu, add_kts, shps_non_zero, seq_added, set4)
    
    
    #new! dmat_keep_all should be of the same order as ddkeep_x_all
    #seq_added <- mod_tried[1:nrep, totx + 3]
    #ansi <- forward_cic_step2(xmat_inactive, y, cic_st, ddkeep_z_all, ddkeep_x_all, knots_lst_inactive, 
    #                          oldcic, oldetahat, oldcoefs, infocrit, dfmean_old, dfmeanmax_old, penalty_old,
    #                          active, inactive, type, nsim, wt.iter, family, parallel, space, pnt, 
    #                          dmat_keep_lst, ps_keep_vec, a, maxedf, add_edfu, add_kts, shps_non_zero, seq_added, set4)
    
    #cat('finished 2nd step', ncol(ddkeep_x_all), '\n')
    ncols <- c(ncols, ncol(ddkeep_x_all))
    cic_vals <- ansi$cic_vals
    #new: avoid tie of e.g. conv and s.incr.conv
    cic_vals <- round(cic_vals, 8)
    cic <- min(cic_vals, na.rm = TRUE)
    
    ddkeep_x <- ansi$ddkeep_x
    ddkeep_z <- ansi$ddkeep_z
    etahatkeep <- ansi$etahatkeep 
    coefskeep <- ansi$coefskeep
    oldetahat <- etahatkeep
    oldcoefs <- coefskeep
    rm_id <- ansi$rm_id
    mskeep <- ansi$mskeep
    #penalty_old <- ansi$penalty 
    ktskeep <- ansi$ktskeep
    cic_st_next_step <- ansi$cic_st_next_step #used edf0 not edf_m
    dmat_keep <- ansi$dmat_keep
    ps_keep <- ansi$ps_keep
    #loop through predictors in the inactive set
    if(cic_st_next_step < oldcic) {
      x_index <- which(cic_vals == min(cic_vals, na.rm = TRUE), arr.ind = TRUE)
    } else {
      x_index <- which(cic_vals == round(oldcic, 8), arr.ind = TRUE)
    }
    x_index <- rbind(x_index)
    
    #test more: tie
    #flat should be the 5th if y is binomial and CBIC
    #if(NROW(x_index) > 1 & any(x_index[, 1] == 9)){ 
    if(NROW(x_index) > 1 & any(x_index[, 1] %in% c(5, 9))){ 
      #if(NROW(x_index) > 1 & any(x_index[, 1] == nrow(cic_vals))){ 
      #stop('check!')
      #more than 1 x's best shape is flat; should we add them all at once?
      x_index <- x_index[NROW(x_index), ]
    } else if(NROW(x_index) > 1 & length(unique(x_index[, 2])) == 1){
      #more than 1 shape of 1 x has the same cic_val; choose the one with smallest null dfmean
      dfmean_cbic_vals <- ansi$dfmean_cbic_vals
      if(!anyNA(dfmean_cbic_vals)){
        col_ps <- unique(x_index[, 2])
        dfmean_x_chosen <- dfmean_cbic_vals[, col_ps]
        row_ps <- x_index[, 1]
        x_index <- x_index[which.min(dfmean_x_chosen[row_ps]), ]
      } else {
        #temp; no simulation for dfmean_cbic; choose any of the tie
        x_index <- x_index[NROW(x_index), ]
      }
    }
    
    x_index <- rbind(x_index)
    x_chosen <- c(x_chosen, inactive[x_index[2]])
    # print (x_chosen)
    # print (head(xmat_inactive))
    # print (x_index)
    #is_tree <- inherits(xmat_inactive[, x_index[2]], "factor") & attr(xmat_inactive[, x_index[2]], "type") %in% c("tree")
    if(!inherits(xmat_inactive[, x_index[2]], "factor")){
      if(attr(xmat_inactive[, x_index[2]], "type") %in% c("ord", "ordinal")){
        #cat('check shapes', '\n')
        shps <- c(1:8, 0)
        if(wt.iter & infocrit == "CBIC"){
          shps <- c(1:4, 0)
        }
      } else if(attr(xmat_inactive[, x_index[2]], "type") %in% c("tree")) {
        shps <- c(1:2, 0)
      }else {
        shps <- 9:17
        if(wt.iter & infocrit == "CBIC"){
          shps <- c(9:12, 0)
        }
      }
      shp_ch_i <- shps[x_index[1]]
      shps_chosen <- c(shps_chosen, shps[x_index[1]])
    } else {
      shp_ch_i <- ifelse(is.null(ddkeep_z), 0, 1)
      shps_chosen <- c(shps_chosen, shp_ch_i)
    }
    #record chosen model's cic, x, shape
    cic_vec <- c(cic_vec, cic)
    colnames(cic_vals) <- colnames(xmat_inactive)
    rownames(cic_vals) <- shps_char
    
    #new
    if(type == "classo") {
      check1 <- length(inactive) > 0
    }
    if(type == "stepwise") {
      check1 <- cic_st_next_step < oldcic & length(inactive) > 0 #works for diabetes together with cic_st = cic
      #check1 = cic < oldcic & length(inactive) > 0
    }
  }
  #cat('final step xtx:', dim(xtxkeep), '\n')
  #new
  ch1 <- is.null(ansi_back)
  ch2 <- FALSE
  if(!is.null(ansi_back)){
    if(!ansi_back$changed) {
      ch2 <- TRUE
    }
  }
  if(ch1 | ch2){
    if(length(varlst) > 0){
      ux_added <- unique(varlst)
      for(x_added in ux_added){
        if(np == 0){
          coefskeep_x_lst[[x_added]] <- (coefskeep[-1])[which(varlst %in% x_added)]
        }
        if(np > 0){
          coefskeep_x_lst[[x_added]] <- (coefskeep[-c(1:(1+np))])[which(varlst %in% x_added)]
        }
      }
    }
    #new
    if(length(varlst_z) > 0){
      uz_added <- unique(varlst_z)
      for(z_added in uz_added){
        coefskeep_z_lst[[z_added]] <- (coefskeep[2:(1+np)])[which(varlst_z %in% z_added)]
      }
    }
  }
  
  #the order of gmat in ansi_back is not following the order of seq_added: need to get coefskeep here
  #for etacomps 
  #new:
  dv <- 0
  if(pnt & nrep >= 1){
    dmat <- Matrix::bdiag(Filter(Negate(is.null), dmat_keep_lst))
    dmat <- as.matrix(dmat)
    if(prod(dim(dmat)) != 0){
      np_ix <- NCOL(ddkeep_z_all) + 1
      if(np_ix > 0){
        dm_zero <- matrix(0, nrow = nrow(dmat), ncol = np_ix)
        dmat <- cbind(dm_zero, dmat)
      }
      dv <- crossprod(dmat)
    }
  }
  
  #1 df for the intercept 
  #dfmeankeep <- dfmean_old #+ 1
  #penalty <- penalty_old
  if(nrep == 0) {
    #mod_tried <- mod_0
    mod_tried <- mod_tried[1, ,drop = FALSE]
    rownames(mod_tried) <- paste("step", 0, sep = "")
    mod_tried[, totx+2] <- mod_0$cic
    yhat <- mu
  } else if (nrep >= 1) {
    mod_tried <- mod_tried[seq(nrep), ,drop = FALSE]	#
    
    #test!
    #nrep <- min(totx, nrep)
    mod_tried <- mod_tried[!duplicated(mod_tried[, totx + 3]), ,drop=FALSE]
    nrep <- nrow(mod_tried)
    
    rownames(mod_tried) <- paste("step", seq(nrep), sep = "")
    seq_added <- mod_tried[, totx + 3]
    #in case some z is redundant; rm_id is not NULL
    yhat <- etahatkeep 
  }
  
  #yhat = etahatkeep + mu
  if(!wt.iter){
    r2 <- 1 - sum((y0 - yhat)^2) / sum((y0 - mu)^2)
  } else {
    r2 <- NULL #change to be pseudo r2 later
  }
  #----------------------------------------------------------------------------------------------------
  #etacomps store fitted components NOT orthogonal to 1 vector or [1,x] if convex/concave
  #----------------------------------------------------------------------------------------------------
  zcoefs <- xcoefs <- NULL
  coefskeep_withone <- coefskeep #not ordered
  coefskeep <- coefskeep[-1] #remove intercept coef
  #new xcoefs and zcoefs should appear in the order of seq_add; otherwise etacomps is wrong!
  if(!is.null(varlst_z)){
    nz <- NCOL(ddkeep_z_all)
    zcoefs <- coefskeep[1:nz]
    xcoefs <- coefskeep[-c(1:nz)]
    #test!
    #zcoefs <- zcoefs[order(varlst_z)]
    #xcoefs <- xcoefs[order(varlst)]
  } else {
    xcoefs <- coefskeep
    #test!
    #xcoefs <- xcoefs[order(varlst)]
  }
  
  rm_rep_x <- NULL
  if(nrep >= 1) {
    seq_added <- mod_tried[, totx + 3]
    shps <- (mod_tried[nrep, 1:totx]) #not ordered
    shps_added <- shps[seq_added]
    #print (shps_added)
    xnms_added <- xnms[seq_added]
    xnms_cont <- xnms[!shps %in% 1]
    xmat_cont <- xmat[, !shps %in% 1, drop = FALSE]
    shps_cont <- shps[!shps %in% 1]
  } else if (nrep == 0) {
    xnms_cont <- xnms
    xmat_cont <- xmat
    #new:
    shps <- rep(0, totx)
    shps_cont <- shps
    seq_added <- NULL 
    #print (shps_added)
  }
  
  #print (length(mskeep_lst))
  #if a chosen shp is 17/flat, then let etacomp = 0; it will be removed at fmat
  etacomps <- matrix(0, nrow = totx, ncol = n)
  coefskeep_lst <- vector("list", length = totx)
  #coefskeep_x_lst <- vector("list", length = totx)
  #coefskeep_z_lst <- vector("list", length = totx)
  if(nrep >= 1) {
    totx_ch <- length(seq_added)
    for(i in 1:totx_ch){
      #print (i)
      ps <- seq_added[i]
      xnmi <- xnms_added[i]
      #print (ps)
      xi <- xmat[, ps]
      shpi <- shps_added[i]
      #if(shpi == 17 | ps %in% x_rm){#flat
      if(shpi == 17 | shpi == 0){#flat
        next
      } else {
        #one <- cbind(rep(1, n))
        #xi = as.numeric(xi)
        #if(shpi == 18) {
        #  etacomps[ps, ] <- ddkeep_z_all[, varlst_z %in% ps, drop = FALSE] %*% zcoefs[which(varlst_z %in% ps)]
        #  coefskeep_lst[[ps]] <- zcoefs[which(varlst_z %in% ps)]
        #} else 
        
        if ((shpi == 1 & inherits(xi, "factor")) | (shpi == 2 & attr(xi, "type") == "tree")) {#need to use 1 for incr later
          #etacomps[ps, ] <- (ddkeep_z_all[, which(varlst_z %in% ps), drop = FALSE]) %*% zcoefs[which(varlst_z %in% ps)]
          
          etacomps[ps, ] <- ddkeep_z_lst[[ps]] %*% coefskeep_z_lst[[ps]]
          #coefskeep_lst[[ps]] <- zcoefs[which(varlst_z %in% ps)]
          
          #coefskeep_x_lst
          #coefskeep_z_lst[[ps]] <- zcoefs[which(varlst_z %in% ps)]
          #capk <- capk + length(zcoefs[which(varlst_z %in% ps)])
          #zmat <- cbind(zmat, ddkeep_z_all[, which(varlst_z %in% ps), drop = FALSE])
        } else if (shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5) {
          #etacomps[ps, ] <- (ddkeep_x_all[, which(varlst %in% ps), drop = FALSE]) %*% xcoefs[which(varlst %in% ps)] 
          #etacomps[ps, ] <- etacomps[ps, ] + (xi-mean(xi)) * zcoefs[which(varlst_z %in% ps)]
          
          etacomps[ps, ] <- ddkeep_x_lst[[ps]] %*% coefskeep_x_lst[[ps]]
          #etacomps[ps, ] <- etacomps[ps, ] + ddkeep_z_lst[[ps]] %*% coefskeep_z_lst[[ps]]
          etacomps[ps, ] <- etacomps[ps, ] + ddkeep_z_lst[[ps]] * coefskeep_z_lst[[ps]]
          
          #coefskeep_lst[[ps]] <- c(xcoefs[which(varlst %in% ps)], zcoefs[which(varlst_z %in% ps)])
          #coefskeep_x_lst[[ps]] <- xcoefs[which(varlst %in% ps)]
          #coefskeep_z_lst[[ps]] <- zcoefs[which(varlst_z %in% ps)]
          #print (ps)
          #print (coefskeep_z_lst)
          #zmat_x <- cbind(zmat_x, (xi-mean(xi)))
        } else {
          #etacomps[ps, ] <- (ddkeep_x_all[, which(varlst %in% ps), drop = FALSE]) %*% xcoefs[which(varlst %in% ps)]
          
          etacomps[ps, ] <- ddkeep_x_lst[[ps]] %*% coefskeep_x_lst[[ps]]
          
          #coefskeep_lst[[ps]] <- xcoefs[which(varlst %in% ps)]
          #coefskeep_x_lst[[ps]] <- xcoefs[which(varlst %in% ps)]
        }
      }
    }
  }
  #k-fold CV to find shrinkage for each component
  if(nrep >= 1 & type == "classo"){
    fmat <- t(etacomps)
    colnames(fmat) <- xnms
    if(anyNA(fmat)) {
      print (length(coefskeep))
      print (NCOL(ddkeep_x_all))
      print (NCOL(ddkeep_z_all))
      print (head(fmat))
      stop('check fmat!')
    }
    
    one <- rep(1, n)
    normx <- sqrt(drop(one %*% (fmat^2)))
    nosignal <- normx / sqrt(n) < eps
  }
  
  nng.mod <- NULL
  cv_lasso <- NULL
  mod_final <- mod_0
  if(nrep == 0) {
    cvec <- rep(0, totx)
    #mod_final = mod_0
    cv_rslt <- list(cvec = cvec)
    mod_tried <- mod_0
    shps_added <- NULL 
    #colnames(mod_tried)[totx + 4] = "shrink parameter"
  } else if(nrep >= 1) {
    if(type == "classo") {
      cvec0 <- rep(0, NCOL(fmat))
      if(all(nosignal)) {
        cvec <- cvec0
        cv_rslt <- list(cvec = cvec)
      } else {
        #cv_rslt <- cv.classo(fmat, y0, nfolds = nfolds, type = "glmnet", family = family)
        #cvec <- cv_rslt$cvec
        
        if(any(nosignal)){
          nosignal_id <- which(nosignal)
          fmat <- fmat[, -nosignal_id, drop = FALSE]
        }
        cv_rslt <- cv.classo(fmat, y0, nfolds = nfolds, type = "glmnet", family = family)
        cvec <- cv_rslt$cvec
        if(any(nosignal)){
          cvec0[-nosignal_id] <- cvec
          cvec <- cvec0
        }
      }
      cvec00 <- rep(0, NROW(mod_tried))
      mod_tried <- cbind(mod_tried, cvec)
      colnames(mod_tried)[NCOL(mod_tried)] <- "shrink parameter"
      #new: change etacomps and etahatkeep here! each etacomps should be multiplied by its cvec
      #need more test!
      #use for predict.classo
      totx_ch <- length(seq_added)
      for(i in 1:totx_ch){
        ps <- seq_added[i]
        shpi <- shps_added[i]
        if(shpi == 17 | shpi == 0){#flat
          next
        } else {
          if(cvec[ps] == 0){
            if((shpi == 1 & inherits(xmat[,ps], "factor")) | (shpi == 2 & attr(xmat[,ps], "type") == "tree") | shpi %in% c(11, 12) | shpi %in% c(3, 4)){
              coefskeep_z_lst[ps] <- list(NULL)
            }
            if(shpi >= 9 & shpi <= 16 | shpi >= 1 & shpi <= 8){
              coefskeep_x_lst[ps] <- list(NULL)
            }
            etacomps[ps, ] <- 0
            etahatkeep <- colSums(etacomps) + coefskeep_withone[1]
            yhat <- etahatkeep
          }
          
          if(cvec[ps] > 0 & shrinketa){
            etacomps[ps, ] <- cvec[ps] * etacomps[ps, ]
            etahatkeep <- colSums(etacomps) + coefskeep_withone[1]
            yhat <- etahatkeep
          }
        }
      }
    } else {
      seq_added <- mod_tried[1:nrep, totx + 3]
      x_kept <- seq_added
      cvec <- rep(0, totx)
      if(any(shps_added %in% c(0, 17))) {
        rm_ps <- which(shps_added %in% c(0, 17))
        x_kept <- x_kept[-rm_ps]
      }
      cvec[x_kept] <- 1
      cv_rslt <- NULL
    }
    if(type == "classo") {
      shps[shps == 17] <- 0
      shps_added_char <- sapply(shps_added, ShapeToChar_shapeselect)
      cvec_added <- mod_tried[, NCOL(mod_tried)]
      
      if(any(cvec_added < eps)) {
        ps_out <- which(cvec_added < eps)
        shps_added_char[ps_out] <- 'out-of-mod'
        #new: test more
        shps[ps_out] <- 0
      }
    } else {
      shps[shps == 17] <- 0
      shps_added_char <- sapply(shps_added, ShapeToChar_shapeselect)
    }
  }
  
  #mod_final <- mod_tried[nrow(mod_tried), ,drop = FALSE]
  #new: create some terms matching terms in cgam, to use vcov
  #if(!wt.iter){
  if(nrep >= 1){
    ddx <- do.call(base::cbind, ddkeep_x_lst)
    zkeep_id <- intersect(obs[sapply(ddkeep_z_lst, function(e) !is.null(e))], zid)
    zmat <- do.call(base::cbind, ddkeep_z_lst[zkeep_id]) # no one, exclude conv/conc, used for plot.classo
    capk <- NCOL(zmat)
    xckeep_id <- setdiff(obs[sapply(ddkeep_z_lst, function(e) !is.null(e))], zid)
    xcmat <- do.call(base::cbind, ddkeep_z_lst[xckeep_id]) #conv/conc/linear in the linear space# test linear
    zmat <- cbind(zmat, xcmat)
  } else {
    ddx <- NULL
    zmat <- NULL
    capk <- 0
  }
  
  #match cgam
  etahat <- yhat #added ybar if gaussian
  muhat <- family$linkinv(etahat)
  wt <- rep(1, n) #temp, replace it with what cgam has later
  if(nrep >= 1){
    bigmat <- t(cbind(one, zmat, ddx))
    d0 <- NCOL(zmat) + 1 #add constant back
    edf <- sum(abs(coefskeep) > 1e-12) #include 1
  } else {
    bigmat <- NULL
    d0 <- 1
    edf <- 0
  }
  
  #new: change shps to be words
  tags <- rep('x', p)
  tags[which(types %in% c('nom', 'nominal'))] <- 'z'
  tags[which(types %in% 'tree')] <- 'tree'
  #shapes in character
  shapes <- mapply(ShapeToChar_shapeselect, shps, tags)
  rslt <- list(mod_tried = mod_tried, etahatkeep = etahatkeep, ddkeep_x_all = ddkeep_x_all, 
               knots_lst =  knots_lst, ddkeep_z_all = ddkeep_z_all, 
               mskeep_lst = mskeep_lst, varlst = varlst, varlst_z = varlst_z, coefskeep = coefskeep, 
               r2 = r2, cic_vec = cic_vec, d0 = d0, wt = wt, y = y0, edf = edf,
               etahat = etahat, muhat = muhat, etacomps = etacomps, coefskeep_lst = coefskeep_lst,
               cvec = cvec, cv_rslt = cv_rslt, family = family, dv = dv,
               shps_added = shps_added, bigmat = bigmat, capk = capk, ddkeep_x_lst = ddkeep_x_lst,
               ddkeep_z_lst = ddkeep_z_lst,
               zcoefs = zcoefs, xcoefs = xcoefs, nfolds = nfolds, xmat = xmat0, xmat_cont = xmat_cont, 
               shps = shps, shapes = shapes, shps_cont = shps_cont, xnms = xnms, xnms_cont = xnms_cont, rm_id=rm_id,
               ps_keep_vec = ps_keep_vec, type = type, space = space, pnt = pnt, 
               sdx = sdx, mx = mx, coefskeep_withone = coefskeep_withone, seq_added = seq_added, dmat_keep = dmat_keep,
               coefskeep_x_lst = coefskeep_x_lst, coefskeep_z_lst = coefskeep_z_lst, 
               types = types, center = center, standardize=standardize, boundary.knots=boundary.knots,
               elapsed_times = elapsed_times, a = a) #, ncols=ncols, a = a, m_total=m_total)
  class(rslt) <- c("shapeselect", "cgam")
  return (rslt)
}

#------------------------------------------
#create penalty matrix
#difference from quadratic or cubic splines
#for i-spl and c-spl
#2nd deri of c-spl is 1st deri of i-spl
#------------------------------------------
makedmat_1D_ic <- function(m, knots, shpi = 9) {
  q <- 1
  dmat <- matrix(0, nrow = (m - q), ncol = m)
  nkts <- length(knots)
  
  if(shpi %in% c(9, 10)){
    vals0 <- iSpline(knots, knots = knots[-c(1, nkts)], degree = 1, derivs = 1L, scale=FALSE)
    vals0 <- diag(vals0)
    if(shpi == 9){
      vals <- vals0 
    }
    if(shpi == 10){
      vals <- -vals0
    }
  }
  
  if(shpi %in% c(11, 12)){
    vals0 <- cSpline(knots, knots = knots[-c(1, nkts)], degree = 1, derivs = 2L, scale=FALSE)
    vals0 <- diag(vals0)
    if(shpi == 11){
      vals <- vals0 
    }
    if(shpi == 12){
      vals <- -vals0
    }
  }
  
  if(shpi == 13) {
    vals0 <- cSpline(knots, knots = knots[-c(1, nkts)], degree = 1, derivs = 2L, scale=FALSE)
    vals0 <- diag(vals0)
    vals <- c(vals0, 0) #the last column in delta is x
  }
  
  if(shpi == 14) {
    vals0 <- cSpline(knots, knots = knots[-c(1, nkts)], degree = 1, derivs = 2L, scale=FALSE)
    vals0 <- -diag(vals0)
    vals <- c(vals0, 0)     
  }
  
  if(shpi == 15) {
    vals0 <- cSpline(knots, knots = knots[-c(1, nkts)], degree = 1, derivs = 2L, scale=FALSE)
    vals0 <- diag(vals0)
    vals <- c(vals0, 0) 
  }
  
  if(shpi == 16) {
    vals0 <- cSpline(knots, knots = knots[-c(1, nkts)], degree = 1, derivs = 2L, scale=FALSE)
    vals0 <- -diag(vals0)
    vals <- c(vals0, 0) 
  }
  
  #vals0 <- diag(vals0)
  
  # if (shpi == 9 | shpi == 11) {
  #   vals <- vals0
  # } else if (shpi == 10 | shpi == 12) {
  #   vals <- -vals0
  # } else if (shpi == 13 | shpi == 14) {
  #   vals <- c(vals0, 0)
  # } else if (shpi == 15 | shpi == 16) {
  #   vals <- c(-vals0, 0)
  # }
  
  for (i in 2:m) {
    dmat[i-1, i-1] <- -1 * vals[i-1]; dmat[i-1, i] <- 1 * vals[i]
  }
  return (dmat)
}

#--------------------------------------------------------------------
#forward step used in classo or stepwise
#change something like cic_vals to be with 10 rows: 1 for tree ordering
#--------------------------------------------------------------------
forward_cic <- function(xmat, y, cic_st, infocrit = "CAIC", type = "stepwise", 
                        knots_lst = NULL, boundary.knots = NULL, 
                        nsim=50, ddkeep_z_all=NULL, wt.iter=FALSE, 
                        family = gaussian(), parallel=FALSE, space="E", pnt=FALSE,
                        a=1, maxedf = FALSE, add_edfu=0, add_kts=0, set4=FALSE)
{
  cicfamily <- CicFamily(family)
  llh.fun <- cicfamily$llh.fun
  #new: use log link in gamma
  linkfun <- cicfamily$linkfun
  etahat.fun <- cicfamily$etahat.fun
  gr.fun <- cicfamily$gr.fun
  wt.fun <- cicfamily$wt.fun
  zvec.fun <- cicfamily$zvec.fun
  muhat.fun <- cicfamily$muhat.fun
  ysim.fun <- cicfamily$ysim.fun
  deriv.fun <- cicfamily$deriv.fun
  dev.fun <- cicfamily$dev.fun
  capl <- NCOL(xmat)  
  n <- NROW(xmat)
  #new: add tree like z: 0 out, 1 tree, 2 unordered
  #cic_vals <- matrix(NA, nrow = 10, ncol = capl)
  #dfmean_cbic_vals <- matrix(NA, nrow = 10, ncol = capl)
  #no 18
  cic_vals <- matrix(NA, nrow = 9, ncol = capl)
  dfmean_cbic_vals <- matrix(NA, nrow = 9, ncol = capl)
  one <- cbind(rep(1, n))
  ddkeep_x <- NULL
  ddkeep_z <- NULL
  dmat_keep <- NULL
  ps_keep <- 0
  dfmean <- dfmeankeep <- dfmeanmaxkeep <- 0 #shpi = 17
  penalty <- penaltykeep <- 0
  ktskeep <- NULL
  mskeep <- NULL
  etahatkeep <- NULL
  coefskeep <- NULL
  oldcic <- cic_st
  cic_st_next_step <- oldcic
  #dfmean_keep_cbic <- 1e+3 #?
  #new:
  #dfmean_caic_vec <- NULL
  dm <- ps <- NULL
  xtxkeep <- face <- facekeep <- NULL
  for(ix in seq_len(capl))
  {
    #print (ix)
    x <- xmat[, ix]
    #if(inherits(x, "factor") & attr(x, 'type') != 'tree'){
    if(inherits(x, "factor") & attr(x, 'type') == 'nom'){
      dd <- model.matrix(~ x)[, -1, drop=FALSE]
      ms <- apply(dd, 2, mean)
      #ms <- one %*% rbind(ms)
      
      dd <- apply(dd, 2, function(e) e - mean(e))
      #compute cic
      #xm <- dd
      xm <- cbind(one, dd) #add intercept
      #new: update np
      np_ix <- NCOL(xm)
      #ybar is removed from y
      #add ybar back 
      zvec <- y
      if(!wt.iter){
        bhat <- solve(crossprod(xm), t(xm)) %*% zvec
        etahat <- xm %*% bhat
      } else {
        #m0 <- glm(y ~ factor(x), family=family$family)
        m0 <- stats::glm(y ~ xm - 1, family = family) #xm includes intercept, needs to remove 1 from glm
        #bhat <- coef(m0)[-1]
        bhat <- coef(m0)
        etahat <- m0$linear.predictors
        muhat <- family$linkinv(etahat)
      }
      
      rss <- sum((y - etahat)^2)
      edf0 <- 0 
      
      dfmean <- length(bhat) + edf0
      
      if(wt.iter){
        llh <- llh.fun(y, muhat, etahat, phihat=NULL, n, weights=NULL, fml = family$family)
        #test!
      } else {llh <- log(rss)}
      
      if (infocrit == "CBIC") {
        #cic_val <- llh + log(n) * (dfmean) / n
        #new
        #cic_val <- llh + a * log(n) * (dfmean) / n^(6/7)
        degree <- 6/7#1
        #degree <- 1
        cic_val <- llh + a * log(n) * (dfmean) / n^degree
      }
      
      if (infocrit == "CAIC") {
        cic_val <- llh + 2 * (dfmean) / n
      }
      
      #if this x reduces cic, then let its etahat to be kept; update cic_st to be cic
      #the new dd will be added to null space; np_add the the # of new columns added
      cic_vals[1, ix] <- oldcic #start with the cic from the yhat = ybat
      if(cic_val < cic_st){
        etahatkeep <- etahat
        coefskeep <- bhat
        ddkeep_z <- dd
        ddkeep_x <- NULL
        #dmat_keep <- NULL
        mskeep <- ms
        cic_st <- cic_val
        
        cic_st_next_step_0 <- llh + 2 * (dfmean) / n
        
        if(cic_st_next_step_0 < oldcic) {
          cic_st_next_step <- cic_st_next_step_0
        }
        
        np_add <- NCOL(dd)
        capm_add <- 0
        
        cic_vals[1, ix] <- cic_val
      }
      #cic_vals[1, ix] <- cic_val
    } else {
      if(attr(x, 'type') %in% c('num', 'numeric')){
        if(set4){ #new
          shps <- 9:12
        }else{
          shps <- 9:16
        }
      }
      if(attr(x, 'type') %in% c('ord', 'ordinal')){
        if(set4){#new
          shps <- 1:4
        }else{
          shps <- 1:8
        }
      }
      if(attr(x, 'type') %in% c('tree')){ #1: tree order #2: unordered
        shps <- 1:2
      }
      #new: if we use cbic and wt.iter then we only have 4 options
      if(wt.iter & infocrit == "CBIC"){
        if(attr(x, 'type') %in% c('num', 'numeric')){
          shps <- 9:12
        }
        if(attr(x, 'type') %in% c('ord', 'ordinal')){
          shps <- 1:4
        }
        if(attr(x, 'type') %in% c('tree')){
          shps <- 1:2
        }
      }
      #shps <- c(9:16, 18) #include 17 later
      #shps <- 9:18 #17 is flat
      #shps <- c(9:16)
      nobs <- length(shps)
      #loop through shapes and compute cic
      gmat_lst <- vector("list", length = nobs)
      np_ix_lst <- vector("list", length = nobs)
      #mu0_lst <- vector("list", length = nobs)
      ms_lst <- vector("list", length = nobs)
      dd_lst <- vector("list", length = nobs)
      #new: penalty matrix
      dmat_lst <- vector("list", length = nobs)
      dm_lst <- vector("list", length = nobs)
      #new:
      knots_lst <- vector("list", length = nobs)
      ps_vec <- rep(0, nobs)
      cic_vals_max <- cic_vals_edf0 <- NULL
      etahat_lst <- vector("list", length = nobs)
      bhat_lst <- vector("list", length = nobs)
      xtx_lst <- vector("list", length = nobs)
      face_lst <- vector("list", length = nobs)
      
      #record gmat, np_ix, etc
      #no 17
      for(i in 1:nobs)
      {
        #shpi <- c(9:16)[i]
        #shpi <- c(9:16, 18)[i]
        shpi <- shps[i]
        #print (c(ix, shpi))
        if(shpi == 18) {
          #dd <- cbind(x)
          ms <- mean(x)
          #check more!
          dd <- cbind(x - ms)
          gmat <- cbind(one, dd)
          np_ix <- 2
          #test more
          dm <- NULL
          ps <- 0
          ps_vec[nobs] <- ps
          #dd <- x
        } else {
          if(shpi <= 8){
            
            if(attr(x, 'type') %in% c('ord', 'ordinal')){
              dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = 0, space = space)
              ms <- dd_ans$ms
              dd <- t(dd_ans$amat)
              ps <- 0
              dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
              
              if(shpi > 2 & shpi < 5){
                gmat <- cbind(one, x - mean(x), dd)
                np_ix <- 2
              } else {
                gmat <- cbind(one, dd)
                np_ix <- 1
              }
            }
            
            #new
            if(attr(x, 'type') %in% c('tree')){ 
              #ms <- NULL
              dd <- t(tree.fun(x))
              ms <- apply(dd, 2, mean)
              dd <- apply(dd, 2, function(e) e - mean(e))
              gmat <- cbind(one, dd)
              if(shpi == 1){
                #dd <- t(tree.fun(x))
                dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
                ps <- 0
                #ms <- apply(dd, 2, mean)
                #dd <- apply(dd, 2, function(e) e - mean(e))
                #gmat <- cbind(one, dd)
                np_ix <- 1
              }
              
              if(shpi == 2){
                #dd <- model.matrix(~ x)[, -1, drop=FALSE]
                #dd <- t(tree.fun(x))
                #ms <- apply(dd, 2, mean)
                #dd <- apply(dd, 2, function(e) e - mean(e))
                #gmat <- cbind(one, dd)
                dm <- NULL
                ps <- 0
                np_ix <- 1 + ncol(dd)
                #ms <- apply(dd, 2, mean)
                #dd <- apply(dd, 2, function(e) e - mean(e))
              }
            }
            
            #ps <- 0
            #dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
            #dm_lst[[ix]] <- dm #no penalty
          }
          
          if(shpi > 8){
            mult <- ifelse(pnt, 2, 1)
            spl_degree_i <- ifelse(shpi %in% c(9, 10), 2, 3) #use 2 for i-spl, 3 for c-spl
            n1 <- length(unique(x))
            #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/5)) + 5, trunc(n1^(1/7)) + 5)
            #new: see if we can fix choosing 13 - 16 too much
            #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/5)) + 5, trunc(n1^(1/7)) + 5)
            nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/7)) + 6 + add_kts, trunc(n1^(1/9)) + 6 + add_kts)
            edfu <- ifelse(spl_degree_i == 2, nkts/mult, nkts/mult + 1) + add_edfu
            
            kts <- seq.int(min(x), max(x), length.out = nkts)
            
            #new: boundary knots
            kts0 <- boundary.knots[[ix]]
            if(length(kts0) >= 2){
              x_rng <- range(x, na.rm = TRUE)
              kts0_rng <- range(kts0, na.rm = TRUE)
              
              tol <- 1e-10
              reset_lwr <- kts0_rng[1] > x_rng[1] + tol
              reset_upp <- kts0_rng[2] < x_rng[2] - tol
              
              if(reset_lwr){
                kts_lwr <- x_rng[1]
              } else {
                kts_lwr <- kts0_rng[1]
              }
              if(reset_upp){
                kts_upp <- x_rng[2]
              } else {
                kts_upp <- kts0_rng[2]
              }
              
              kts[1] <- kts_lwr
              kts[nkts] <- kts_upp
            }
            
            dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = kts, space = space)
            #dd is orth to 1
            ms <- dd_ans$ms
            dd <- t(dd_ans$amat)
            #kts <- dd_ans$knots
            
            if(shpi > 10 & shpi < 13){
              gmat <- cbind(one, x - mean(x), dd)
              np_ix <- 2
            } else {
              gmat <- cbind(one, dd)
              np_ix <- 1
            }
            #print (kts)
            #new:
            if(pnt){
              nc <- NCOL(dd)
              dm <- makedmat_1D_ic(m = nc, knots = kts, shpi = shpi) #different from makedmat_1D for testpar
              qvi <- crossprod(dd)
              dvi <- crossprod(dm)
              #ps <- 1e-3
              #9, 10, 11, 12, 13, 14, 15, 16
              if(shpi %in% c(9, 11, 13)){
                log_root_res <- uniroot(f = function(log_pen,...) .search_ps_shapeselect(exp(log_pen),...),
                                        interval = c(-8, 8), 
                                        qv0 = qvi, dv = dvi, dd = dd, edfu = edfu,
                                        extendInt = "yes",
                                        tol = .Machine$double.eps^0.25)
                # convert back to the original scale
                #temp
                if(abs(log_root_res$root) > 8) {
                  ps <- floor(exp(8))
                } else {
                  ps <- exp(log_root_res$root)
                }
                #print (ps)
                #print (log_root_res$root)
                #temp!
                #if(abs(ps) < 1e-5) {ps <- 1e-5}
                ps_vec[which(shps %in% shpi)] <- ps
              } else if (shpi == 10){
                ps <- ps_vec[1]
                ps_vec[which(shps %in% shpi)] <- ps
              } else if (shpi == 12) {
                ps <- ps_vec[3]
                ps_vec[which(shps %in% shpi)] <- ps
              } else {
                ps <- ps_vec[5]
                ps_vec[which(shps %in% shpi)] <- ps
              }
              dm <- sqrt(ps) * dm
            }
          }
        }
        
        dd_lst[[i]] <- dd
        
        if(!is.null(ms)){
          ms_lst[[i]] <- ms
        }
        
        gmat_lst[[i]] <- gmat
        np_ix_lst[[i]] <- np_ix
        
        if(pnt & !is.null(dm)){
          dmat_lst[[i]] <- dm
          dm_lst[[i]] <- dm
          if(shpi > 8){
            knots_lst[[i]] <- kts
          }
        }
      }
      
      #compute dfmean_max, which is the average of all dfmean_max over nsim simulations
      dfmean_rslt <- compute_edf(family=family, gmat_lst, dmat_lst=dmat_lst,
                                 mu0_lst=NULL, np_ix_lst, weights=NULL, 
                                 nsim=nsim, wt.iter=wt.iter, parallel=parallel,
                                 pnt=pnt)
      dfmean_caic <- dfmean_rslt$dfmean_caic
      dfmean_cbic <- dfmean_rslt$dfmean_cbic
      #new
      #if(nobs > 4){
      #dfmean_cbic_vals[-9, ix] <- dfmean_cbic 
      #} 
      #if(nobs == 4){
      dfmean_cbic_vals[1:nobs, ix] <- dfmean_cbic 
      #}
      #new:
      #dfmean_caic_vec <- c(dfmean_caic_vec, dfmean_caic)
      
      #go through all shapes to compute rss and cic/bic
      #cic_vals[9, ix] <- oldcic #start with yhat = ybar
      for(shpi in shps){
        #ps_shpi <- which(c(9:16, 18) %in% shpi)
        #ps_shpi <- which(c(9:16) %in% shpi)
        ps_shpi <- which(shps %in% shpi)
        gmat <- gmat_lst[[ps_shpi]]
        ms <- ms_lst[[ps_shpi]]
        dd <- dd_lst[[ps_shpi]]
        np_ix <- np_ix_lst[[ps_shpi]]
        
        #new
        # if(shpi <= 8){
        #   pnt <- FALSE
        # }
        
        if(pnt){
          dm <- dmat_lst[[ps_shpi]]
          ps <- ps_vec[ps_shpi]
        }
        
        nc = ncol(gmat)
        amat = diag(nc)
        ps_zero = 1:np_ix
        amat[ps_zero, ps_zero] = 0
        bvec = cbind(rep(0, nc))
        xtx = crossprod(gmat)
        
        if(pnt & !is.null(dm)){
          dm_zero <- matrix(0, nrow = nrow(dm), ncol = np_ix)
          dmat <- cbind(dm_zero, dm)
          dv <- crossprod(dmat) #absort penalty already
          xtx <- xtx + dv
        }
        
        if (shpi == 18 | shpi == 2 & attr(x, "type") == 'tree') {
          zvec <- y
          if(!wt.iter){
            bhat <- solve(crossprod(gmat), t(gmat)) %*% zvec
            etahat <- gmat %*% bhat
          } else {
            #gmat has one
            m0 <- stats::glm(y ~ gmat - 1, family=family)
            bhat <- coef(m0)
            etahat <- m0$linear.predictors
          }
        } else {
          zvec <- y
          if(!wt.iter){
            if (!pnt | pnt & all(dm == 0)) {
              dsend = gmat[, -c(1:np_ix), drop = FALSE]
              zsend = gmat[, 1:np_ix, drop = FALSE]
              ans = coneB(zvec, dsend, zsend)
              etahat = ans$yhat
              bhat = ans$coefs
              face = ans$face
            } else {
              #ans <- coneB(zvec, dsend, zsend)
              #etahat <- ans$yhat
              #bhat <- ans$coefs
              xty = crossprod(gmat, zvec)
              # calculate eigenvalues
              ev <- eigen(xtx, only.values = TRUE)$values
              # check if the smallest eigenvalue is too close to zero or negative
              if (min(ev) < 1e-8) {
                diag(xtx) <- diag(xtx) + 1e-8
              }
              #ans = solve.QP(xtx, xty, t(amat), bvec)
              ans = solve.QP(xtx, xty, amat, bvec)
              bhat = ans$solution
              bhat = round(bhat, 10)
              etahat = gmat %*% bhat
              face = ans$iact
            }
          } else {
            gmat0 <- gmat
            weights <- rep(1, n)
            etahat <- etahat.fun(n, y, fml = family$family)
            gr <- gr.fun(y, etahat, weights, fml = family$family)
            wt <- wt.fun(y, etahat, n, weights, fml = family$family)
            cvec <- wt * etahat - gr
            zvec <- zvec.fun(cvec, wt, y, fml = family$family)
            muhat <- family$linkinv(etahat)
            diff <- 1
            niter <- 0
            sm <- 1e-5
            #sm <- 1e-7
            #mdiff <- TRUE
            if (family$family == "binomial") {
              mdiff <- abs(max(muhat) - 1) > sm
            } else {mdiff <- TRUE}
            while (diff > sm & mdiff & niter < 10) 
            {
              oldmu <- muhat
              niter <- niter + 1
              gr <- gr.fun(y, etahat, weights=NULL, fml = family$family)
              wt <- wt.fun(y, etahat, n, weights=NULL, fml = family$family)
              
              cvec <- wt * etahat - gr
              zvec <- zvec.fun(cvec, wt, y, fml = family$family)
              gmat <- sweep(gmat0, MARGIN = 1, sqrt(wt), FUN = '*')
              #test
              if(!pnt | pnt & all(dm == 0)){
                dsend <- gmat[, -c(1:np_ix), drop = FALSE]
                zsend <- gmat[, 1:np_ix, drop = FALSE]
                ansi <- coneB(zvec, dsend, zsend)
                bhat <- ansi$coefs
              }
              
              if(pnt & !is.null(dm) & !all(dm == 0)){
                xtx = crossprod(gmat)
                xtx = xtx + dv
                xty = crossprod(gmat, zvec)
                # calculate eigenvalues
                ev <- eigen(xtx, only.values = TRUE)$values
                # check if the smallest eigenvalue is too close to zero or negative
                if (min(ev) < 1e-8) {
                  diag(xtx) <- diag(xtx) + 1e-8
                }
                ansi = solve.QP(xtx, xty, amat, bvec)
                #ansi = solve.QP(xtx, xty, t(amat), bvec)
                bhat = ansi$solution
                face = ansi$iact
              }
              
              bhat = round(bhat, 10)
              etahat <- gmat0 %*% bhat
              muhat <- family$linkinv(etahat)
              mdiff <- abs(max(muhat) - 1)
              if (family$family == "binomial") {
                mdiff <- abs(max(muhat) - 1) > sm
              } else {mdiff <- TRUE}
              
              diff <- mean((muhat - oldmu)^2)
              #mdiff <- abs(max(muhat) - 1) > sm
            }
          }
        }
        
        rss <- sum((y - etahat)^2)
        
        #new
        etahat_lst[[ps_shpi]] <- etahat
        bhat_lst[[ps_shpi]] <- bhat
        
        xtx_lst[[ps_shpi]] <- xtx
        face_lst[[ps_shpi]] <- face
        
        if(wt.iter){
          llh <- llh.fun(y, muhat, etahat, phihat=NULL, n, weights=NULL, fml=family$family)
        } else {llh <- log(rss)}
        
        if (infocrit == "CBIC") {
          #cic_val <- llh + log(n) * (dfmean_caic) / n
          #cic_val_edf0 <- llh + log(n) * (dfmean_cbic[ps_shpi]) / n
          
          #new
          #cic_val <- llh + a * log(n) * (dfmean_caic) / n^(6/7)
          #cic_val_edf0 <- llh + a * log(n) * (dfmean_cbic[ps_shpi]) / n^(6/7)
          
          if(shpi %in% c(9, 10)){
            degree <- 6/7
          } else if (shpi > 10 & shpi < 17){
            degree <-  6/7#8/9 # 6/7#8/9
          } else if (shpi > 0 & shpi < 9) {
            degree <- 6/7#1 # 6/7#1
          } else {
            stop('check!')
          }
          cic_val <- llh + a * log(n) * (dfmean_caic) / n^(degree)
          cic_val_edf0 <- llh + a * log(n) * (dfmean_cbic[ps_shpi]) / n^(degree)
        } 
        
        if (infocrit == "CAIC") {
          cic_val <- llh +  2 * (dfmean_caic) / n 
          cic_val_edf0 <- llh +  2 * (dfmean_cbic[ps_shpi]) / n
        }
        
        #new:
        if(maxedf){
          cic_vals_max <- c(cic_vals_max, cic_val)
        } else {
          cic_vals_max <- c(cic_vals_max, cic_val_edf0)
        }
        
        cic_vals_edf0 <- c(cic_vals_edf0, cic_val_edf0)
        
        #test!
        if(shpi > 8){
          cic_vals[(shpi - 8), ix] <- cic_val
        }
        if(shpi <= 8){
          cic_vals[shpi, ix] <- cic_val
          #cic_vals[which(shps %in% shpi), ix] <- cic_val
        }
        #cic_vals[(shpi - 8), ix] <- cic_val_edf0
      }
      
      cic_vals[9, ix] <- oldcic #flat
      #cic_vals[(1+nobs), ix] <- oldcic #flat
      #new:
      if(any(cic_vals_max < cic_st)){
        #test!
        #cic_vals[-9, ix] <- cic_vals_edf0 #no flat
        
        if (wt.iter & infocrit == "CBIC") {
          cic_vals[1:length(cic_vals_edf0), ix] <- cic_vals_edf0
        } else {
          cic_vals[1:length(cic_vals_edf0), ix] <- cic_vals_edf0
          #cic_vals[-9, ix] <- cic_vals_edf0 #no flat, will not work when it is tree
          #cic_vals[-(1+nobs), ix] <- cic_vals_edf0 #no flat
        }
        
        ps_shpi <- which(cic_vals_edf0 == min(cic_vals_edf0))
        
        if(length(ps_shpi) > 1) {
          ps_shpi <- rev(ps_shpi)[1]
        }
        
        etahatkeep <- etahat_lst[[ps_shpi]] #etahat #will include 1 for binomial; but etacomps will not use this
        coefskeep <- bhat_lst[[ps_shpi]] # bhat
        mskeep <- ms_lst[[ps_shpi]] # ms
        ktskeep <- knots_lst[[ps_shpi]] # kts; will be NULL for ordinal
        dd <- dd_lst[[ps_shpi]] 
        dm <- dmat_lst[[ps_shpi]] 
        ps <- ps_vec[ps_shpi]
        
        shpi <- shps[ps_shpi]
        
        if(shpi != 18) {
          if(shpi == 2 & attr(x, 'type') == 'tree'){
            ddkeep_x <- NULL
            ddkeep_z <- dd
          } else {
            ddkeep_x <- dd
          }
          if(pnt){
            dmat_keep <- dm
            ps_keep <- ps
          }
        }
        
        if(shpi != 2 | attr(x, 'type') != 'tree'){
          ddkeep_z <- NULL
          np_add <- 0
          if (shpi > 10 & shpi < 13 | shpi == 18 | shpi > 2 & shpi < 5){
            ddkeep_z <- cbind(x - mean(x))
            np_add <- 1
          }
        }
        
        cic_st_next_step_0 <- cic_vals_edf0[ps_shpi]
        cic_st_next_step <- cic_st_next_step_0
        #test
        cic_st <- cic_st_next_step_0
        
        #xtxkeep <- xtx_lst[[ps_shpi]]
        #facekeep <- face_lst[[ps_shpi]]
      }
    }
  }
  
  rslt <- list(cic_vals = cic_vals, ddkeep_x = ddkeep_x, ddkeep_z = ddkeep_z, 
               dmat_keep = dmat_keep, ps_keep = ps_keep, mskeep = mskeep, ktskeep = ktskeep,
               etahatkeep = etahatkeep, coefskeep = coefskeep, dfmean_cbic_vals = dfmean_cbic_vals,
               cic_st = cic_st, cic_st_next_step = cic_st_next_step)
  return (rslt)
}


#--------------------------------------------------------------------------------------------------------------------------------
forward_cic_step2 <- function(xmat_inactive, y, cic_st, ddkeep_z_all = NULL, ddkeep_x_all = NULL, knots_lst_inactive = NULL,
                              boundary.knots_inactive = NULL, 
                              oldcic = NULL, oldetahat = NULL, oldcoefs = NULL, infocrit = "CAIC", dfmean_old = NULL,
                              dfmeanmax_old = NULL, penalty_old = NULL,
                              active, inactive, type = "stepwise", nsim=50, wt.iter=FALSE,
                              family = gaussian(), parallel=FALSE, space="E", pnt=FALSE, dmat_keep_lst=NULL,
                              ps_keep_vec=NULL, a=1, maxedf = FALSE,
                              add_edfu = 0, add_kts = 0, shps_non_zero, seq_added, set4 = FALSE) {
  cicfamily <- CicFamily(family)
  llh.fun <- cicfamily$llh.fun
  #new: use log link in gamma
  linkfun <- cicfamily$linkfun
  etahat.fun <- cicfamily$etahat.fun
  gr.fun <- cicfamily$gr.fun
  wt.fun <- cicfamily$wt.fun
  zvec.fun <- cicfamily$zvec.fun
  muhat.fun <- cicfamily$muhat.fun
  ysim.fun <- cicfamily$ysim.fun
  deriv.fun <- cicfamily$deriv.fun
  dev.fun <- cicfamily$dev.fun
  
  capl <- NCOL(xmat_inactive)
  n <- NROW(xmat_inactive)
  cic_st_0 <- cic_st #largetest cic_st before putting in any new x; cic val for 17 in classo
  cic_st_next_step <- oldcic #for stepwise, cic_st and oldcic have the same starting val; cic_st will be replaced when decreased
  #new: add tree
  #cic_vals <- matrix(NA, nrow = 10, ncol = capl)
  #dfmean_cbic_vals <- matrix(NA, nrow = 10, ncol = capl)
  #no 18
  cic_vals <- matrix(NA, nrow = 9, ncol = capl)
  dfmean_cbic_vals <- matrix(NA, nrow = 9, ncol = capl)
  etahatkeep <- oldetahat
  coefskeep <- oldcoefs
  ddkeep_z <- NULL
  ddkeep_x <- NULL
  ktskeep <- NULL
  mskeep <- NULL
  np_add <- 0
  rm_id <- NULL
  one <- cbind(rep(1, n))
  #new:
  dmat_keep <- dm <- dmat <- dmat_keep_all <- NULL
  
  ps_keep <- ps <- 0
  if(pnt){
    dmat_keep_all <- Matrix::bdiag(Filter(Negate(is.null), dmat_keep_lst[seq_added]))
    #dmat_keep_all <- Matrix::bdiag(dmat_keep_lst[seq_added])
    dmat_keep_all <- as.matrix(dmat_keep_all)
    if(prod(dim(dmat_keep_all)) == 0) {
      dmat_keep_all <- NULL
    }
  }
  for(ix in seq_len(capl)){
    #print (ix)
    x <- xmat_inactive[, ix]
    #for categorical variable: oldcic is for 0; cic_val is for 1. they need to be compared with each other
    if(inherits(x, "factor") & attr(x, 'type') == 'nom'){
      #dd is orth to 1
      dd <- model.matrix(~ x)[, -1, drop=FALSE]
      #new:
      #if(NCOL(ddkeep_z_all) == 0){
      ms <- apply(dd, 2, mean)
      dd <- apply(dd, 2, function(e) e - mean(e))
      #}
      #not sure how to do this in predict.classo....
      # if(NCOL(ddkeep_z_all) > 0){
      #   dz_temp <- cbind(1, ddkeep_z_all)
      #   pmz <- dz_temp %*% MASS::ginv(crossprod(dz_temp)) %*% t(dz_temp)
      #   ms <- pmz %*% dd
      #   #coefms <- solve(crossprod(dz_temp)) %*% t(dz_temp) %*% dd
      #   dd <- dd - ms
      #   nosignal <- apply(dd, 2, function(e) all(round(e, 8) < 1e-8))
      #   if(all(nosignal)){
      #     rm_id <- c(rm_id, ix)
      #     cic_vals[1, ix] <- oldcic
      #     next
      #   } else if(any(!nosignal)){
      #     kept_id <- which(!nosignal)
      #     dd <- dd[, kept_id, drop = FALSE]
      #     ms <- ms[, kept_id, drop = FALSE]
      #   }
      # }
      #new: some column in dd may be the same as some column in ddkeep_z_all
      #skip z if all columns in z are already in the space of ddkeep_z_all
      # if(NCOL(ddkeep_z_all) > 0){
      #   cat('ms changed!', '\n')
      #   #pmz = ddkeep_z_all %*% solve(crossprod(ddkeep_z_all), t(ddkeep_z_all))
      #   pmz <- ddkeep_z_all %*% MASS::ginv(crossprod(ddkeep_z_all)) %*% t(ddkeep_z_all)
      #   ms_add <- pmz %*% dd
      #   dd <- dd - ms_add
      #   #test!
      #   ms_addt <- t(ms_add)
      #   ms <- ms_addt + ms
      #   ms <- t(ms)
      #   sgn_dd <- apply(dd, 2, function(e) all(round(e, 8) < 1e-8))
      #   if(all(sgn_dd)) {
      #     rm_id <- c(rm_id, ix)
      #     #new: make sure cic_vals is not all NA
      #     cic_vals[1, ix] <- oldcic
      #     next
      #   } else if(any(!sgn_dd)){
      #     #cat('ms changed!', '\n')
      #     kept_id <- which(!sgn_dd)
      #     dd <- dd[, kept_id, drop = FALSE]
      #     #ms <- one %*% rbind(ms) + ms_add
      #     #test!
      #     ms <- ms[, kept_id, drop = FALSE]
      #   }
      # }
      xm <- cbind(one, ddkeep_z_all, dd, ddkeep_x_all)
      gmat <- xm
      #zvec <- y - etahat_all #residual vector
      zvec <- y
      capm <- NCOL(ddkeep_x_all)
      np_ix <- NCOL(ddkeep_z_all) + NCOL(dd) + 1
      dmat <- dmat_keep_all
      #vm <- gmat[, 1:np_ix]
      if(capm == 0){
        #no x has been chosen
        edf0 = 0
        dfmean = edf0 + ncol(xm)
        dfmean_caic <- dfmean
        dfmean_cbic <- dfmean
      } else {
        dfmean_rslt <- compute_edf(family = family, gmat_lst = list(gmat), dmat_lst=list(dmat),
                                   mu0_lst = NULL,
                                   np_ix_lst = list(np_ix), weights=NULL, nsim=nsim, wt.iter=wt.iter,
                                   parallel=parallel, pnt=pnt) #oldcic is without z (edf is the average of the chosen shape without max)
        dfmean_caic <- dfmean_rslt$dfmean_caic #only has 1 val; caic and cbic are the same
        dfmean_cbic <- dfmean_rslt$dfmean_cbic #only has 1 val
      }
      #new:
      #dfmean_caic_vec <- c(dfmean_caic_vec, dfmean_caic)
      
      if(capm == 0){
        if(!wt.iter){
          bhat = solve(crossprod(xm), t(xm)) %*% zvec
          etahat <- xm %*% bhat
        } else {
          m0 = stats::glm(y ~ xm - 1, family=family)
          bhat = coef(m0)
          etahat <- m0$linear.predictors
          muhat <- family$linkinv(etahat)
        }
      }
      
      if(capm > 0){ #no x
        #gmat <- xm
        nc = ncol(gmat)
        amat = diag(nc)
        ps_zero = 1:np_ix
        amat[ps_zero, ps_zero] = 0
        #amat[1:np_ix, 1:np_ix] = 0
        bvec = cbind(rep(0, nc))
        xtx = crossprod(gmat)
        
        #test! avoid linear, shp = 18
        if(pnt & !is.null(dmat)){
          dm_lst = vector("list", length = 2)
          dm0 = matrix(0, nrow = nrow(dmat), ncol = np_ix)
          dm_lst[[1]] = dm0
          dm_lst[[2]] = dmat
          dmat = as.matrix(bdiag(dm_lst))
          dv = crossprod(dmat)
          xtx = xtx + dv
        }
        
        if(!wt.iter){
          #test
          if(!pnt | pnt & all(dmat == 0)){
            dsend = gmat[, -c(1:np_ix), drop = FALSE]
            zsend = gmat[, 1:np_ix, drop = FALSE]
            ans = coneB(zvec, dsend, zsend)
            bhat = ans$coefs
            bhat = round(bhat, 10)
            etahat = gmat %*% bhat
            face = ans$face
          } else {
            xty = crossprod(gmat, zvec)
            # calculate eigenvalues
            ev <- eigen(xtx, only.values = TRUE)$values
            # check if the smallest eigenvalue is too close to zero or negative
            if (min(ev) < 1e-8) {
              diag(xtx) <- diag(xtx) + 1e-8
            }
            ans = solve.QP(xtx, xty, amat, bvec)
            #ans = solve.QP(xtx, xty, t(amat), bvec)
            bhat = ans$solution
            bhat = round(bhat, 10)
            #cat('finished bhat in check z:', '\n')
            etahat = gmat %*% bhat
            face = ans$iact
          }
        } else {
          gmat0 <- gmat
          weights <- rep(1, n)
          etahat <- etahat.fun(n, y, fml = family$family)
          gr <- gr.fun(y, etahat, weights, fml = family$family)
          wt <- wt.fun(y, etahat, n, weights, fml = family$family)
          cvec <- wt * etahat - gr
          zvec <- zvec.fun(cvec, wt, y, fml = family$family)
          muhat <- family$linkinv(etahat)
          sm <- 1e-5
          #sm <- 1e-7 #will cause solve error
          
          if (family$family == "binomial") {
            mdiff <- abs(max(muhat) - 1) > sm
          } else {mdiff <- TRUE}
          
          diff <- 1
          niter <- 0
          while (diff > sm & mdiff & niter < 10)
          {
            oldmu <- muhat
            niter <- niter + 1
            gr <- gr.fun(y, etahat, weights=NULL, fml = family$family)
            wt <- wt.fun(y, etahat, n, weights=NULL, fml = family$family)
            
            cvec <- wt * etahat - gr
            zvec <- zvec.fun(cvec, wt, y, fml = family$family)
            gmat <- sweep(gmat0, MARGIN = 1, sqrt(wt), FUN = '*')
            
            #test
            if(!pnt | pnt & all(dmat == 0)){
              dsend <- gmat[, -c(1:np_ix), drop = FALSE]
              zsend <- gmat[, 1:np_ix, drop = FALSE]
              ansi <- coneB(zvec, dsend, zsend)
              bhat <- ansi$coefs
            }
            
            if(pnt & !is.null(dmat) & !all(dmat == 0)){
              xtx = crossprod(gmat)
              xtx = xtx + dv
              xty = crossprod(gmat, zvec)
              # calculate eigenvalues
              ev <- eigen(xtx, only.values = TRUE)$values
              # check if the smallest eigenvalue is too close to zero or negative
              if (min(ev) < 1e-8) {
                diag(xtx) <- diag(xtx) + 1e-8
              }
              ansi = solve.QP(xtx, xty, amat, bvec)
              #ansi = solve.QP(xtx, xty, t(amat), bvec)
              bhat = ansi$solution
              face = ansi$iact
            }
            
            bhat = round(bhat, 10)
            etahat <- gmat0 %*% bhat
            muhat <- family$linkinv(etahat)
            mdiff <- abs(max(muhat) - 1)
            if (family$family == "binomial") {
              mdiff <- abs(max(muhat) - 1) > sm
            } else {mdiff <- TRUE}
            
            diff <- mean((muhat - oldmu)^2)
          }
        }
      }
      #check
      #vm <- cbind(one, ddkeep_z_all, dd)
      rss <- sum((y - etahat)^2)
      if(wt.iter){
        llh <- llh.fun(y, muhat, etahat, phihat=NULL, n, weights=NULL, fml=family$family)
        #test!
        #llh <- 2 * (llh * n / 2) #2 times neg-loglike
      } else {llh <- log(rss)}
      
      #include 1 when computing cic
      if (infocrit == "CBIC") {
        #cic_val <- llh + log(n) * (dfmean_cbic) / n
        #new:
        #cic_val <- llh + a * log(n) * (dfmean_cbic) / n^(6/7)
        
        all_shps <- c(shps_non_zero, 1)
        if(any(all_shps %in% c(9, 10))){ #smooth incr/decr
          degree <- 6/7
        } else if (any(all_shps > 10 & all_shps < 17)){ #smooth conv/conc
          degree <- 6/7#8/9 #
        } else if (any(all_shps > 0 & all_shps < 9)) { #ord or parametric
          degree <- 6/7#1 #
        } else {
          stop('check!')
        }
        
        cic_val <- llh + a * log(n) * (dfmean_cbic) / n^degree
      }
      if (infocrit == "CAIC") {
        cic_val <- llh + 2 * (dfmean_caic) / n
      }
      
      #start with shape 0; oldcic
      #cic_vals[1, ix] <- cic_st_0
      cic_vals[1, ix] <- oldcic
      if(cic_val < cic_st){
        etahatkeep <- etahat
        coefskeep <- bhat
        
        ddkeep_z <- dd
        mskeep <- ms
        ddkeep_x <- NULL
        dmat_keep <- NULL
        ps_keep <- 0
        ktskeep <- NULL
        cic_st <- cic_val
        cic_st_next_step_0 <- cic_val
        cic_st_next_step <- cic_st_next_step_0
        np_add <- NCOL(ddkeep_z)
        cic_vals[1, ix] <- cic_val
        
        #xtxkeep <- xtx
        #facekeep <- face
      }
      #cat(ix, 'z', cic_val, '\n')
    } else {
      #test more
      #cic_st <- oldcic
      #shps <- 9:18 #17 is flat; 17 is included
      if(attr(x, 'type') %in% c('num', 'numeric')){
        #new
        if(set4){
          shps <- c(9:12, 17)
        }else{
          shps <- 9:17
        }
        #shps <- 9:17
      }
      if(attr(x, 'type') %in% c('ord', 'ordinal')){
        #new
        if(set4){
          shps <- c(1:4, 0)
        }else{
          shps <- c(1:8, 0)
        }
        #shps <- c(1:8, 0)
      }
      if(attr(x, 'type') %in% c('tree')){
        shps <- c(1, 2, 0)
      }
      #new: if we use cbic and wt.iter then we only have 4 options
      if(wt.iter & infocrit == "CBIC"){
        if(attr(x, 'type') %in% c('num', 'numeric')){
          shps <- c(9:12, 17)
        }
        if(attr(x, 'type') %in% c('ord', 'ordinal')){
          shps <- c(1:4, 0)
        }
        if(attr(x, 'type') %in% c('tree')){
          shps <- c(1, 2, 0)
        }
      }
      nobs <- length(shps)
      gmat_lst <- vector("list", length = nobs)
      np_ix_lst <- vector("list", length = nobs)
      #mu0_lst <- vector("list", length = nobs)
      ms_lst <- vector("list", length = nobs)
      dd_lst <- vector("list", length = nobs)
      kts_lst <- vector("list", length = nobs)
      dmat_lst <- vector("list", length = nobs)
      dm_lst <- vector("list", length = nobs)
      ps_vec <- rep(0, length = nobs)
      
      #dmat <- dmat_keep_all
      #record gmat, np_ix, etc
      for(i in 1:nobs)
      {
        shpi <- shps[i]
        #print (c(ix, shpi))
        if (shpi == 17 | shpi == 0) {
          dd <- NULL
          ms <- NULL
          kts <- NULL
          gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all)
          np_ix <- NCOL(ddkeep_z_all) + 1
          dm <- NULL
          ps <- 0
          dmat <- dmat_keep_all
          cic_val_edf0 <- oldcic
          cic_val <- oldcic
        } else if (shpi == 18) {
          #dd <- cbind(x)
          dd <- cbind(x - mean(x))
          kts <- NULL
          ms <- mean(x)
          gmat <- cbind(one, ddkeep_z_all, dd, ddkeep_x_all)
          np_ix <- NCOL(ddkeep_z_all) + 2
          #dmat_lst[[i]] <- 0
          #not_skip_cic <- TRUE
          dm <- NULL
          ps <- 0
          dmat <- dmat_keep_all
        } else {
          
          if(shpi <= 8) {
            
            if(attr(x, 'type') %in% c('ord')){
              dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = 0, space = space)
              ms <- dd_ans$ms
              dd <- t(dd_ans$amat)
              if(shpi > 2 & shpi < 5){
                gmat <- cbind(one, ddkeep_z_all, x-mean(x), ddkeep_x_all, dd)
                np_ix <- NCOL(ddkeep_z_all) + 2
              } else {
                gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all, dd)
                np_ix <- NCOL(ddkeep_z_all) + 1
              }
              dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
            }
            
            if(attr(x, 'type') %in% c('tree')){
              #ms <- NULL
              dd <- t(tree.fun(x))
              ms <- apply(dd, 2, mean)
              dd <- apply(dd, 2, function(e) e - mean(e))
              if(shpi == 1){
                #dd <- t(tree.fun(x))
                #ms <- apply(dd, 2, mean)
                #dd <- apply(dd, 2, function(e) e - mean(e))
                gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all, dd)
                np_ix <- NCOL(ddkeep_z_all) + 1
                dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
              }
              if(shpi == 2){
                #dd <- model.matrix(~ x)[, -1, drop=FALSE]
                #dd <- t(tree.fun(x))
                #ms <- apply(dd, 2, mean)
                #dd <- apply(dd, 2, function(e) e - mean(e))
                gmat <- cbind(one, ddkeep_z_all, dd, ddkeep_x_all)
                np_ix <- NCOL(ddkeep_z_all) + 1 + NCOL(dd)
                dm <- NULL
              }
            }
            
            ps <- 0
            dmat <- dmat_keep_all
            if(!is.null(dmat)){
              dm_lst_tmp <- vector('list', length = 2)
              dm_lst_tmp[[1]] <- dmat
              dm_lst_tmp[[2]] <- dm
              dmat <- as.matrix(Matrix::bdiag(dm_lst_tmp))
            }else{
              dmat <- dm
            }
          }
          
          if(shpi > 8) {
            mult <- ifelse(pnt, 2, 1)
            spl_degree_i <- ifelse(shpi %in% c(9, 10), 2, 3) #use 2 for i-spl, 3 for c-spl
            n1 <- length(unique(x))
            #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/5)) + 5, trunc(n1^(1/7)) + 5)
            #new
            #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/5)) + 5, trunc(n1^(1/7)) + 5)
            nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/7)) + 6 + add_kts, trunc(n1^(1/9)) + 6 + add_kts)
            edfu <- ifelse(spl_degree_i == 2, nkts/mult, nkts/mult + 1) + add_edfu
            
            kts <- seq.int(min(x), max(x), length.out = nkts)
            #new: boundary knots
            kts0 <- boundary.knots_inactive[[ix]]
            if(length(kts0) >= 2){
              x_rng <- range(x, na.rm = TRUE)
              kts0_rng <- range(kts0, na.rm = TRUE)
              
              tol <- 1e-10
              reset_lwr <- kts0_rng[1] > x_rng[1] + tol
              reset_upp <- kts0_rng[2] < x_rng[2] - tol
              
              if(reset_lwr){
                kts_lwr <- x_rng[1]
              } else {
                kts_lwr <- kts0_rng[1]
              }
              if(reset_upp){
                kts_upp <- x_rng[2]
              } else {
                kts_upp <- kts0_rng[2]
              }
              
              kts[1] <- kts_lwr
              kts[nkts] <- kts_upp
            }
            
            dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = kts, space = space)
            #dd is orth to 1
            ms <- dd_ans$ms
            dd <- t(dd_ans$amat)
            
            #gmat <- cbind(zmat, ddkeep_z_all, ddkeep_x_all, dd)
            if(shpi > 10 & shpi < 13){
              gmat <- cbind(one, ddkeep_z_all, x-mean(x), ddkeep_x_all, dd)
              np_ix <- NCOL(ddkeep_z_all) + 2
            } else {
              gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all, dd)
              np_ix <- NCOL(ddkeep_z_all) + 1
            }
            
            if(pnt){
              nc <- NCOL(dd)
              dm <- makedmat_1D_ic(m = nc, knots = kts, shpi = shpi) #different from makedmat_1D for testpar
              qvi <- crossprod(dd)
              dvi <- crossprod(dm)
              #ps <- 1e-6
              #9, 10, 11, 12, 13, 14, 15, 16
              if(shpi %in% c(9, 11, 13)){
                log_root_res <- uniroot(f = function(log_pen,...) .search_ps_shapeselect(exp(log_pen),...),
                                        interval = c(-8, 8),
                                        qv0 = qvi, dv = dvi, dd = dd, edfu = edfu,
                                        extendInt = "yes",
                                        tol = .Machine$double.eps^0.25)
                # convert back to the original scale
                #temp
                if(abs(log_root_res$root) > 8) {
                  ps <- floor(exp(8))
                } else {
                  ps <- exp(log_root_res$root)
                }
                #print (ps)
                #ps <- exp(log_root_res$root)
                ps_vec[which(shps %in% shpi)] <- ps
              } else if (shpi == 10){
                ps <- ps_vec[1]
                ps_vec[which(shps %in% shpi)] <- ps
              } else if (shpi == 12) {
                ps <- ps_vec[3]
                ps_vec[which(shps %in% shpi)] <- ps
              } else {
                ps <- ps_vec[5]
                ps_vec[which(shps %in% shpi)] <- ps
              }
              dm <- sqrt(ps) * dm
              
              dmat <- dmat_keep_all
              if(!is.null(dmat)){
                dm_lst_tmp <- vector('list', length = 2)
                dm_lst_tmp[[1]] <- dmat
                dm_lst_tmp[[2]] <- dm
                dmat <- as.matrix(bdiag(dm_lst_tmp))
              }else{
                dmat <- dm
              }
            }
          }
        }
        
        if(shpi == 17 | shpi == 0) {
          dd_lst[i] <- list(NULL)
          ms_lst[i] <- list(NULL)
          kts_lst[i] <- list(NULL)
          #dmat_lst[[i]] <- dmat
        } else {
          dd_lst[[i]] <- dd
          if(!is.null(ms)){
            ms_lst[[i]] <- ms
          }
          #dmat_lst[[i]] <- dmat
          if(shpi == 18){
            kts_lst[i] <- list(NULL)
          } else {
            #will be NULL for ordinal
            if(shpi > 8){
              kts_lst[[i]] <- kts
            }
          }
        }
        
        if(!is.null(dmat)){
          dmat_lst[[i]] <- dmat #diff from forward_cic_step
        }
        
        if(pnt & !is.null(dm)){
          dm_lst[[i]] <- dm
        }
        gmat_lst[[i]] <- gmat
        np_ix_lst[[i]] <- np_ix
      }
      
      #new: compute llh first
      #compute dfmean_max, which is the average of all dfmean_max over nsim simulations
      #need to remove 17 when compute dfmean
      #ps_17 <- which(shps %in% c(17, 0))
      #dfmean_rslt <- compute_edf(family=family, gmat_lst=gmat_lst[-ps_17], dmat_lst=dmat_lst[-ps_17],
      #                           mu0_lst=NULL, np_ix_lst=np_ix_lst[-ps_17],
      #                           weights=NULL, nsim=nsim, wt.iter=wt.iter, parallel=parallel, pnt=pnt)
      #dfmean_caic <- dfmean_rslt$dfmean_caic
      #dfmean_cbic <- dfmean_rslt$dfmean_cbic
      #dfmean_cbic_vals[-9, ix] <- dfmean_cbic
      
      #new
      # if(nobs > 5){
      #  dfmean_cbic_vals[-9, ix] <- dfmean_cbic
      # }
      
      #if(nobs == 5){
      #no 17
      #dfmean_cbic_vals[1:length(dfmean_cbic), ix] <- dfmean_cbic
      #}
      
      #test!
      #dfmean_cbic_vals[1:length(dfmean_cbic), ix] <- dfmean_cbic
      
      #new: compute the cic lower bound for each x first
      # which uses the total number of predictors as dfmean
      #dfmean_caic_vec <- c(dfmean_caic_vec, dfmean_caic)
      #go through all shapes to compute rss and cic/bic
      cic_vals_max <- cic_vals_edf0 <- NULL # rep(cic_st_0, nobs)
      llh_vec <- NULL
      etahat_lst <- vector("list", length = nobs)
      bhat_lst <- vector("list", length = nobs)
      
      face_lst <- vector("list", length = nobs)
      xtx_lst <- vector("list", length = nobs)
      
      for(shpi in shps){
        #print (shpi)
        #ps_shpi <- which(c(9:18) %in% shpi)
        #ps_shpi <- which(c(9:17) %in% shpi)
        ps_shpi <- which(shps %in% shpi)
        gmat <- gmat_lst[[ps_shpi]]
        ms <- ms_lst[[ps_shpi]]
        dd <- dd_lst[[ps_shpi]]
        dm <- dm_lst[[ps_shpi]]
        np_ix <- np_ix_lst[[ps_shpi]]
        kts <- kts_lst[[ps_shpi]]
        
        if(pnt){
          dmat <- dmat_lst[[ps_shpi]]
          ps <- ps_vec[ps_shpi]
        }
        
        nc <- ncol(gmat)
        amat <- diag(nc)
        ps_zero <- 1:np_ix
        amat[ps_zero, ps_zero] <- 0
        bvec <- cbind(rep(0, nc))
        xtx <- crossprod(gmat)
        
        if(pnt & !is.null(dmat)){
          dm_zero <- matrix(0, nrow = nrow(dmat), ncol = np_ix)
          dmat <- cbind(dm_zero, dmat)
          dv <- crossprod(dmat) #absort penalty already
          xtx <- xtx + dv
        }
        
        if (shpi == 17 | shpi == 0) {
          #etahat <- rep(mean(y), n)
          etahat <- oldetahat
          #rss <- sum((y - etahat)^2)
          #cic_val <- oldcic
          cic_val <- cic_st_0
          cic_val_edf0 <- cic_st_0
          bhat <- oldcoefs
          #print (length(bhat))
          dd <- NULL
          ms <- NULL
          #dm <- NULL
        } else {
          if(np_ix > 0 & NCOL(gmat) > np_ix){
            zvec <- y
            if(!wt.iter){
              #test
              if(!pnt | pnt & all(dmat == 0)){
                dsend = gmat[, -c(1:np_ix), drop = FALSE]
                zsend = gmat[, 1:np_ix, drop = FALSE]
                ans = coneB(zvec, dsend, zsend)
                bhat = ans$coefs
                bhat = round(bhat, 10)
                etahat = gmat %*% bhat
                face = ans$face
              } else {
                xty = crossprod(gmat, zvec)
                # calculate eigenvalues
                ev <- eigen(xtx, only.values = TRUE)$values
                # check if the smallest eigenvalue is too close to zero or negative
                if (min(ev) < 1e-8) {
                  diag(xtx) <- diag(xtx) + 1e-8
                }
                #ans = solve.QP(xtx, xty, t(amat), bvec)
                ans = solve.QP(xtx, xty, amat, bvec)
                bhat = ans$solution
                bhat = round(bhat, 10)
                etahat = gmat %*% bhat
                face = ans$iact
              }
            } else {
              gmat0 <- gmat
              weights <- rep(1, n)
              etahat <- etahat.fun(n, y, fml = family$family)
              gr <- gr.fun(y, etahat, weights, fml = family$family)
              wt <- wt.fun(y, etahat, n, weights, fml = family$family)
              cvec <- wt * etahat - gr
              zvec <- zvec.fun(cvec, wt, y, fml = family$family)
              muhat <- family$linkinv(etahat)
              sm <- 1e-5
              #sm <- 1e-7
              if (family$family == "binomial") {
                mdiff <- abs(max(muhat) - 1) > sm
              } else {mdiff <- TRUE}
              diff <- 1
              niter <- 0
              
              while (diff > sm & mdiff & niter < 10)
              {
                oldmu <- muhat
                niter <- niter + 1
                gr <- gr.fun(y, etahat, weights=NULL, fml = family$family)
                wt <- wt.fun(y, etahat, n, weights=NULL, fml = family$family)
                
                cvec <- wt * etahat - gr
                zvec <- zvec.fun(cvec, wt, y, fml = family$family)
                gmat <- sweep(gmat0, MARGIN = 1, sqrt(wt), FUN = '*')
                
                #test
                if(!pnt | pnt & all(dmat == 0)){
                  dsend <- gmat[, -c(1:np_ix), drop = FALSE]
                  zsend <- gmat[, 1:np_ix, drop = FALSE]
                  ansi <- coneB(zvec, dsend, zsend)
                  bhat <- ansi$coefs
                }
                
                if(pnt & !is.null(dmat) & !all(dmat == 0)){
                  xtx = crossprod(gmat)
                  xtx = xtx + dv
                  xty = crossprod(gmat, zvec)
                  # calculate eigenvalues
                  ev <- eigen(xtx, only.values = TRUE)$values
                  # check if the smallest eigenvalue is too close to zero or negative
                  if (min(ev) < 1e-8) {
                    diag(xtx) <- diag(xtx) + 1e-8
                  }
                  ansi = solve.QP(xtx, xty, amat, bvec)
                  #ansi = solve.QP(xtx, xty, t(amat), bvec)
                  bhat = ansi$solution
                  face = ansi$iact
                }
                
                bhat = round(bhat, 10)
                etahat <- gmat0 %*% bhat
                muhat <- family$linkinv(etahat)
                mdiff <- abs(max(muhat) - 1)
                if (family$family == "binomial") {
                  mdiff <- abs(max(muhat) - 1) > sm
                } else {mdiff <- TRUE}
                
                diff <- mean((muhat - oldmu)^2)
              }
            }
          } else if (np_ix > 0 & NCOL(gmat) == np_ix) { #shpi = 18 and no constrained x has been chosen
            #zvec <- y
            #bhat <- solve(crossprod(gmat), t(gmat)) %*% zvec
            #etahat <- gmat %*% bhat
            zvec <- y
            if(!wt.iter){
              #bhat <- solve(crossprod(gmat), t(gmat)) %*% zvec
              bhat <- ginv(crossprod(gmat)) %*% t(gmat) %*% zvec
              etahat <- gmat %*% bhat
            } else {
              m0 <- stats::glm(y ~ gmat - 1, family = family)
              etahat <- m0$linear.predictors
              bhat <- coef(m0)
              muhat <- family$linkinv(etahat)
            }
            #etahat <- gmat %*% solve(crossprod(gmat), t(gmat)) %*% zvec
          }
          
          rss <- sum((y - etahat)^2)
          if(wt.iter){
            llh <- llh.fun(y, muhat, etahat, phihat=NULL, n, weights=NULL, fml = family$family)
            #test!
            #llh <- 2 * (llh * n / 2) #2 times neg-loglike
          } else {llh <- log(rss)}
          
          #new
          llh_vec <- c(llh_vec, llh)
          
          if(shpi != 17 & shpi != 0){ #remove this if?
            #ps_shpi_bic = which(c(9:16, 18) %in% shpi)
            if(shpi > 8){
              ps_shpi_bic = which(c(9:16) %in% shpi)
            }
            if(shpi <= 8){
              ps_shpi_bic = which(c(1:8) %in% shpi)
            }
            if (infocrit == "CBIC") {
              #cic_val <- llh + log(n) * (dfmean_caic) / n
              #cic_val_edf0 <- llh + log(n) * (dfmean_cbic[ps_shpi_bic]) / n
              
              #new:
              #cic_val <- llh + a * log(n) * (dfmean_caic) / n^(6/7)
              #cic_val_edf0 <- llh + a * log(n) * (dfmean_cbic[ps_shpi_bic]) / n^(6/7)
              
              all_shps <- c(shps_non_zero, shpi)
              if(any(all_shps %in% c(9, 10))){ #smooth incr/decr
                degree <- 6/7#
              } else if (any(all_shps > 10 & all_shps < 17)){ #smooth conv/conc
                degree <- 6/7#8/9
              } else if (any(all_shps > 0 & all_shps < 9)) { #ord or parametric
                degree <- 6/7#1
              } else {
                stop('check!')
              }
              
              #new
              dfmean_lwr <- length(active) + 1
              cic_val <- llh + a * log(n) * (dfmean_lwr) / n^(degree)
              cic_val_edf0 <- llh + a * log(n) * (dfmean_lwr) / n^(degree)
            }
            
            if (infocrit == "CAIC") {
              #new
              dfmean_lwr <- length(active) + 1
              cic_val <- llh + 2 * (dfmean_lwr) / n
              cic_val_edf0 <- llh + 2 * (dfmean_lwr) / n
            }
          }
        }
        
        #cat('include x', '\n')
        if(shpi > 8){
          #cic_vals[(shpi - 8), ix] <- cic_val
          #new: use the starting cic, bc we use lower bound for edf0 now, cic_vals could be much smaller than they should
          cic_vals[(shpi - 8), ix] <- cic_st_0
          
        }
        if(shpi <= 8) {
          if(shpi == 0) {
            #cic_vals[9, ix] <- cic_val
            #new
            cic_vals[9, ix] <- cic_st_0
            #cic_vals[nobs, ix] <- cic_val
          } else {
            #cic_vals[ps_shpi, ix] <- cic_val
            #cic_vals[shpi, ix] <- cic_val
            #new
            cic_vals[shpi, ix] <- cic_st_0
          }
        }
        #print (ps_shpi)
        etahat_lst[[ps_shpi]] <- etahat
        bhat_lst[[ps_shpi]] <- bhat
        #face_lst[[ps_shpi]] <- face
        #xtx_lst[[ps_shpi]] <- xtx
        
        #cic_vals_max <- c(cic_vals_max, cic_val)
        if(maxedf){ #not use anymore
          cic_vals_max <- c(cic_vals_max, cic_val)
        } else {
          cic_vals_max <- c(cic_vals_max, cic_val_edf0)
        }
        cic_vals_edf0 <- c(cic_vals_edf0, cic_val_edf0)
        #cic_vals[(shpi - 8), ix] <- cic_val_edf0
      }
      
      #new:
      ps_to_sims <- NULL
      dfmean_cbic_vals[1:(length(shps)-1), ix] <- length(active) + 1
      shps_to_sims <- shps[-length(shps)]
      
      if(any(cic_vals_max < cic_st)){
        
        ps_shps_tosim <- which(cic_vals_max < cic_st)
        shps_to_sims <- shps_to_sims[ps_shps_tosim]
        ps_17 <- which(shps %in% c(17, 0))
        
        #no need to worry about ps_17, because it will be the last one?
        dfmean_rslt <- compute_edf(family=family, gmat_lst=gmat_lst[ps_shps_tosim], dmat_lst=dmat_lst[ps_shps_tosim],
                                   mu0_lst=NULL, np_ix_lst=np_ix_lst[ps_shps_tosim],
                                   weights=NULL, nsim=nsim, wt.iter=wt.iter,
                                   parallel=parallel, pnt=pnt, shps_to_sims=shps_to_sims)
        
        dfmean_caic <- dfmean_rslt$dfmean_caic
        dfmean_cbic <- dfmean_rslt$dfmean_cbic
        
        dfmean_cbic_vals[ps_shps_tosim, ix] <- dfmean_cbic
        
        if (infocrit == "CBIC") {
          cic_vals_max[ps_shps_tosim] <- llh_vec[ps_shps_tosim] + a * log(n) * (dfmean_cbic) / n^(degree)
          cic_vals_edf0[ps_shps_tosim] <- llh_vec[ps_shps_tosim] + a * log(n) * (dfmean_cbic) / n^(degree)
        }
        
        if (infocrit == "CAIC") {
          #new
          cic_vals_max[ps_shps_tosim] <- llh_vec[ps_shps_tosim] + 2 * (dfmean_cbic) / n
          cic_vals_edf0[ps_shps_tosim] <- llh_vec[ps_shps_tosim] +  2 * (dfmean_cbic) / n
        }
      }
      
      #new
      if(any(cic_vals_max < cic_st)){
        if (wt.iter & infocrit == "CBIC") {
          if(length(cic_vals_edf0) == 5){
            cic_vals[1:length(cic_vals_edf0), ix] <- cic_vals_edf0
          } else {
            #tree
            cic_vals[1:2, ix] <- cic_vals_edf0[1:2]
            cic_vals[9, ix] <- cic_vals_edf0[3]
          }
        } else {
          if(length(cic_vals_edf0) == 9){
            cic_vals[, ix] <- cic_vals_edf0
          } else if(length(cic_vals_edf0) == 5){ #new: set4
            cic_vals[1:length(cic_vals_edf0), ix] <- cic_vals_edf0
          } else {
            #tree:
            cic_vals[1:2, ix] <- cic_vals_edf0[1:2]
            cic_vals[9, ix] <- cic_vals_edf0[3]
          }
        }
        
        ps_shpi <- which(cic_vals_edf0 == min(cic_vals_edf0))
        
        if(length(ps_shpi) > 1) {
          ps_shpi <- rev(ps_shpi)[1]
        }
        
        shpi <- shps[ps_shpi]
        etahatkeep <- etahat_lst[[ps_shpi]]
        coefskeep <- bhat_lst[[ps_shpi]]
        mskeep <- ms_lst[[ps_shpi]]
        kts <- kts_lst[[ps_shpi]] #will be NULL for ordinal
        dd <- dd_lst[[ps_shpi]]
        dm <- dm_lst[[ps_shpi]]
        ps <- ps_vec[ps_shpi]
        
        if(shpi != 18) {
          #new:
          if(shpi == 2 & attr(x, 'type') == 'tree'){
            ddkeep_x <- NULL
            ddkeep_z <- dd
            np_add <- ncol(dd)
          } else {
            ddkeep_x <- dd
            ddkeep_z <- NULL
            np_add <- 0
            if (shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5) {
              ddkeep_z <- cbind(x - mean(x))
              np_add <- 1
            }
          }
          
          if(pnt){
            dmat_keep <- dm #new dmat, use dmat for all dmat's combined as block diag
            ps_keep <- ps
          }
          ktskeep <- kts
        }
        
        if(shpi == 18) {
          ddkeep_x <- NULL
          if(pnt){
            dmat_keep <- NULL
            ps_keep <- 0
          }
          ktskeep <- NULL
          #cat('include linear', '\n')
          ddkeep_z <- cbind(x - mean(x))
          np_add <- 1
        }
        cic_st_next_step_0 <- cic_vals_edf0[ps_shpi]
        cic_st_next_step <- cic_st_next_step_0
        cic_st <- cic_st_next_step_0
      }
    }
  }
  rslt <- list(cic_vals = cic_vals, ddkeep_x = ddkeep_x, ktskeep = ktskeep,
               ddkeep_z = ddkeep_z, cic_st = cic_st, cic_st_next_step = cic_st_next_step,
               etahatkeep = etahatkeep, coefskeep = coefskeep, mskeep = mskeep,
               rm_id = rm_id, dmat_keep = dmat_keep, ps_keep = ps_keep,
               dfmean_cbic_vals = dfmean_cbic_vals)
  return (rslt)
}

#----------------------------------------------------------------------------------------------------------
#backward step used in classo or stepwise
#----------------------------------------------------------------------------------------------------------
backward_cic <- function(xmat, y, varlst, varlst_z, ddkeep_x_all, 
                         ddkeep_x_lst, knots_lst=NULL, boundary.knots = NULL,
                         ddkeep_z_all, ddkeep_z_lst, mskeep_lst, mod_tried, nrep, totx,
                         oldcic = NULL, infocrit = "CAIC", dfmean_old, dfmeanmax_old, penalty_old, oldcoefs,
                         oldetahat, active, type = "stepwise", nsim = 100, wt.iter = FALSE,
                         family = gaussian(), parallel = FALSE, space = "E", pnt = FALSE, 
                         dmat_keep_lst = list(NULL), 
                         ps_keep_vec = NULL, a=0.5, 
                         add_edfu=0, add_kts=0, shps_all)
{
  cicfamily <- CicFamily(family)
  llh.fun <- cicfamily$llh.fun
  #new: use log link in gamma
  linkfun <- cicfamily$linkfun
  etahat.fun <- cicfamily$etahat.fun
  gr.fun <- cicfamily$gr.fun
  wt.fun <- cicfamily$wt.fun
  zvec.fun <- cicfamily$zvec.fun
  muhat.fun <- cicfamily$muhat.fun
  ysim.fun <- cicfamily$ysim.fun
  deriv.fun <- cicfamily$deriv.fun
  dev.fun <- cicfamily$dev.fun
  chs <- which(mod_tried[nrep-1, 1:totx] != 0) #already removed z (shp = 0)
  changed <- FALSE
  rm_chi <- mod_tried[nrep-1, 1:totx]
  shps <- 9:17
  #new: if we use cbic and wt.iter then we only have 4 options
  if(wt.iter & infocrit == "CBIC"){
    shps <- c(9:12, 17)
  }
  
  ddkeep_x_all_0 <- ddkeep_x_all
  #cat ('dim of x0', dim(ddkeep_x_all), '\n')
  ddkeep_z_all_0 <- ddkeep_z_all
  dmat_keep_lst_0 <- dmat_keep_lst
  dmat_keep_all <- dmat_keep_all_0 <- dmat <- dm <- NULL
  
  #new: avoid computing coefskeep again
  varlst_0 <- varlst
  varlst_z_0 <- varlst_z
  
  #new! order dmat_keep_lst to get dmat_keep_all
  #follow the same order of ddkeep_x_all
  #cannot use this: may add some x twice if the x was out of the model in the previous step: 2 5 3 1 6 2
  seq_added <- mod_tried[1:nrep, totx + 3]
  #test!
  seq_added <- seq_added[!duplicated(seq_added, fromLast = TRUE)]
  if(pnt){
    dmat_keep_all <- Matrix::bdiag(Filter(Negate(is.null), dmat_keep_lst[seq_added]))
    dmat_keep_all <- as.matrix(dmat_keep_all)
    if(prod(dim(dmat_keep_all)) == 0){
      dmat_keep_all <- NULL
    }
    dmat_keep_all_0 <- dmat_keep_all
  }
  
  cic_st_next_step <- oldcic
  cic_st <- oldcic
  n <- NROW(xmat)
  one <- cbind(rep(1, n))
  #dfmeankeep <- dfmean_old
  etahatkeep <- oldetahat
  coefskeep <- oldcoefs
  back_inactive <- NULL
  back_inactive_check <- rep(FALSE, length(chs))
  
  totx <- NCOL(xmat)
  coefskeep_x_lst <- vector("list", length = totx)
  coefskeep_z_lst <- vector("list", length = totx)
  
  #xtxkeep <- oldxtx
  #facekeep <- oldface
  for(chi in chs){
    #print (varlst_0)
    #print (chi)
    #-----------------------------------------------------------------------------------------
    #change the following to include shapes loop
    #go through all un tried shapes for the previously chosen x's to see if cic is smaller
    #-----------------------------------------------------------------------------------------
    x <- xmat[, chi, drop=TRUE]
    #----------------------------------------------------
    #won't check again for categorical predictor ?
    #rm_chi[chi] records the chosen shape
    if (rm_chi[chi] == 1 & inherits(x, "factor")) {
      #should update varlst_0 and varlst_z_0 here
      varlst_0 <- varlst
      varlst_z_0 <- varlst_z
      shps_not_tried <- 0
      # ddkeep_z_all track changing z columns
      # ddkeep_z_all_0 track final z columns whenever cic_st is decreased
      ddkeep_z_all <- ddkeep_z_all_0[, -which(varlst_z == chi), drop=FALSE]
      gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all_0)
      np_ix <- NCOL(ddkeep_z_all) + 1
      
      #new
      if(pnt){
        dmat <- dmat_keep_all_0
      }
      
      dfmean_rslt <- compute_edf(family=family, gmat_lst=list(gmat), 
                                 dmat_lst = list(dmat), 
                                 mu0_lst=NULL, np_ix_lst=list(np_ix),
                                 weights=NULL, nsim=nsim, wt.iter=wt.iter,
                                 parallel=parallel, pnt=pnt)
      dfmean_caic <- dfmean_rslt$dfmean_caic #val without z; val with z is oldcic
      dfmean_cbic <- dfmean_rslt$dfmean_cbic
      
      nc <- ncol(gmat)
      amat <- diag(nc)
      ps_zero <- 1:np_ix
      amat[ps_zero, ps_zero] <- 0
      bvec <- cbind(rep(0, nc))
      xtx <- crossprod(gmat)
      
      if(pnt & !is.null(dmat)){
        dm_zero <- matrix(0, nrow = nrow(dmat), ncol = np_ix)
        dmat <- cbind(dm_zero, dmat)
        dv <- crossprod(dmat)
        xtx <- xtx + dv
      }
      
      if(np_ix > 0 & NCOL(gmat) > np_ix) {# z and x
        zvec <- y
        if(!wt.iter){
          #test
          if(!pnt | pnt & all(dmat == 0)){
            dsend <- gmat[, -c(1:np_ix), drop = FALSE]
            zsend <- gmat[, 1:np_ix, drop = FALSE]
            ans <- coneB(zvec, dsend, zsend)
            bhat <- ans$coefs
            bhat = round(bhat, 10)
            etahat <- gmat %*% bhat
          } else {
            xty <- crossprod(gmat, zvec)
            # calculate eigenvalues
            ev <- eigen(xtx, only.values = TRUE)$values
            # check if the smallest eigenvalue is too close to zero or negative
            if (min(ev) < 1e-8) {
              diag(xtx) <- diag(xtx) + 1e-8
            }
            ans <- solve.QP(xtx, xty, amat, bvec)
            #ans = solve.QP(xtx, xty, t(amat), bvec)
            bhat <- ans$solution
            bhat = round(bhat, 10)
            etahat <- gmat %*% bhat
          }
        } else {
          gmat0 <- gmat
          weights <- rep(1, n)
          etahat <- etahat.fun(n, y, fml = family$family)
          gr <- gr.fun(y, etahat, weights, fml = family$family)
          wt <- wt.fun(y, etahat, n, weights, fml = family$family)
          cvec <- wt * etahat - gr
          zvec <- zvec.fun(cvec, wt, y, fml = family$family)
          muhat <- family$linkinv(etahat)
          #sm <- 1e-5
          sm <- 1e-7
          if (family$family == "binomial") {
            mdiff <- abs(max(muhat) - 1) > sm
          } else {mdiff <- TRUE}
          
          diff <- 1
          niter <- 0
          
          while (diff > sm  & mdiff & niter < 20)
          {
            oldmu <- muhat
            niter <- niter + 1
            gr <- gr.fun(y, etahat, weights=NULL, fml = family$family)
            wt <- wt.fun(y, etahat, n, weights=NULL, fml = family$family)
            
            cvec <- wt * etahat - gr
            zvec <- zvec.fun(cvec, wt, y, fml = family$family)
            gmat <- sweep(gmat0, MARGIN = 1, sqrt(wt), FUN = '*')
            
            #test
            if(!pnt | pnt & all(dmat == 0)){
              dsend <- gmat[, -c(1:np_ix), drop = FALSE]
              zsend <- gmat[, 1:np_ix, drop = FALSE]
              ansi <- coneB(zvec, dsend, zsend)
              bhat <- ansi$coefs
            }
            
            if(pnt & !is.null(dmat) & !all(dmat == 0)){
              xtx <- crossprod(gmat)
              xtx <- xtx + dv
              xty <- crossprod(gmat, zvec)
              # calculate eigenvalues
              ev <- eigen(xtx, only.values = TRUE)$values
              # check if the smallest eigenvalue is too close to zero or negative
              if (min(ev) < 1e-8) {
                diag(xtx) <- diag(xtx) + 1e-8
              }
              ansi <- solve.QP(xtx, xty, amat, bvec)
              #ansi = solve.QP(xtx, xty, t(amat), bvec)
              bhat <- ansi$solution
            }
            
            bhat = round(bhat, 10)
            etahat <- gmat0 %*% bhat
            muhat <- family$linkinv(etahat)
            mdiff <- abs(max(muhat) - 1)
            if (family$family == "binomial") {
              mdiff <- abs(max(muhat) - 1) > sm
            } else {mdiff <- TRUE}
            
            diff <- mean((muhat - oldmu)^2)
          }
        }
      } else if (np_ix > 0 & NCOL(gmat) == np_ix) { #all z #shpi = 18 and no constrained x has been chosen
        zvec <- y
        if(!wt.iter){
          #bhat <- solve(crossprod(gmat), t(gmat)) %*% zvec
          bhat <- ginv(crossprod(gmat)) %*% t(gmat) %*% zvec
          bhat = round(bhat, 10)
          etahat <- gmat %*% bhat
        } else {
          m0 <- stats::glm(y ~ gmat - 1, family = family)
          etahat <- m0$linear.predictors
          bhat <- coef(m0)
          muhat <- family$linkinv(etahat)
        }
      }
      
      # #mu0 is used for CIC simulation later
      rss <- sum((y - etahat)^2)
      if(wt.iter){
        llh <- llh.fun(y, muhat, etahat, phihat=NULL, n, weights=NULL, fml = family$family)
      } else {llh <- log(rss)}
      
      if (infocrit == "CBIC") {
        #cic_val <- llh + log(n) * (dfmean_cbic) / n
        #new
        #cic_val <- llh + a * log(n) * (dfmean_cbic) / n^(6/7)
        
        all_shps <- shps_all
        all_shps[chi] <- 0
        if(any(all_shps %in% c(9, 10))){ #smooth incr/decr
          degree <- 6/7
        } else if (any(all_shps > 10 & all_shps < 17)){ #smooth conv/conc
          degree <- 6/7#8/9 #6/7#
        } else if (any(all_shps > 0 & all_shps < 9)) { #ord or parametric
          degree <- 6/7#1 #6/7#
        } else if (all(all_shps %in% c(0, 17))) { #won't happen
          stop('check!')
        }
        cic_val <- llh + a * log(n) * (dfmean_cbic) / n^(degree)
      }
      
      if (infocrit == "CAIC") {
        cic_val <- llh + 2 * (dfmean_caic) / n
      }
      
      # #new: compare without z vs with z; cic_val[1] is with z
      bool <- (cic_val - cic_st > 0) & (cic_val - cic_st) < 1e-8
      if(cic_val < cic_st | bool){
        #cat('z changed in backward!', 'nrep:', nrep, '\n')
        changed <- TRUE
        new_shp <- shps_not_tried #will always = 0
        mod_tried[nrep, chi] <- new_shp
        #new: check more
        mod_tried[nrep, totx + 2] <- cic_val
        etahatkeep <- etahat
        coefskeep <- bhat
        
        ddkeep_z <- NULL
        ddkeep_z_lst[chi] <- list(NULL)
        #new:
        varlst_z_0 <- varlst_z_0[-which(varlst_z_0 == chi)]
        #new:
        #coefskeep_z_lst[chi] <- list(NULL)
        
        mskeep_lst[chi] <- list(NULL)
        np_add <- 0
        knots_lst[chi] <- list(NULL)
        ddkeep_x <- NULL
        ddkeep_x_lst[chi] <- list(NULL) #is already null bc it is factor?
        
        #new: need to update ddkeep_z_all_0 and ddkeep_x_all_0 here in case >= 2 predictors' shapes are changed
        seq_added <- mod_tried[1:nrep, totx + 3]
        ddkeep_z_lst_ord <- ddkeep_z_lst[seq_added]
        ddkeep_z_all_0 <- do.call(base::cbind, ddkeep_z_lst_ord)
        
        obs <- 1:totx
        obs_ord <- obs[seq_added]
        
        #new
        coefs_nointercept <- coefskeep[-1]
        #print (varlst_z_0)
        if(length(varlst_z_0) > 0){
          zid <- 1:length(varlst_z_0)
          zcoefs <- coefs_nointercept[zid]
          xcoefs <- coefs_nointercept[-zid]
        }
        if(length(varlst_z_0) == 0){
          zcoefs <- NULL
          xcoefs <- coefs_nointercept
        }
        
        #new: we should update ddkeep_x_all and ddkeep_z_all here
        ddkeep_x_all_temp <- do.call(base::cbind, ddkeep_x_lst)
        ddkeep_z_all_temp <- do.call(base::cbind, ddkeep_z_lst)
        accx <- NCOL(ddkeep_x_all_temp)
        accz <- NCOL(ddkeep_z_all_temp)
        
        coefskeep_x_lst <- vector("list", length = totx)
        coefskeep_z_lst <- vector("list", length = totx)
        
        if(accx > 0){
          #new
          #print (varlst_0)
          ux_added <- unique(varlst_0) 
          for(x_added in ux_added){
            coefskeep_x_lst[[x_added]] <- xcoefs[which(varlst_0 %in% x_added)]
          }
          
          #order to be put back in coefskeep
          ord_x <- order(varlst_0)
          xcoefs <- xcoefs[ord_x]
        }
        
        if(accz > 0){
          varlst_z_lst <- lapply(1:length(obs_ord), function(i, zlst) {if(NCOL(zlst[[i]]) >= 1) {rep(obs_ord[i], NCOL(zlst[[i]]))}}, zlst = ddkeep_z_lst_ord)
          varlst_z <- do.call(base::c, varlst_z_lst)
          
          #new
          uz_added <- unique(varlst_z_0)
          for(z_added in uz_added){
            coefskeep_z_lst[[z_added]] <- zcoefs[which(varlst_z_0 %in% z_added)]
          }
          
          ord_z <- order(varlst_z_0)
          zcoefs <- zcoefs[ord_z]
        } else {varlst_z <- NULL}
        
        #new: put b0 back
        coefskeep <- c(coefskeep[1], zcoefs, xcoefs)
        
        cic_st <- cic_val
        if (infocrit == "CBIC"){
          #cic_st_next_step_0 <- llh + log(n) * (dfmean_cbic) / n
          
          #new
          #cic_st_next_step_0 <- llh + a * log(n) * (dfmean_cbic) / n^(6/7)
          cic_st_next_step_0 <- llh + a * log(n) * (dfmean_cbic) / n^(degree)
        }
        if (infocrit == "CAIC"){
          cic_st_next_step_0 <- llh + 2 * (dfmean_caic) / n
        }
        
        if(cic_st_next_step_0 < oldcic){
          cic_st_next_step <- cic_st_next_step_0
        }
        #new: put it back to active set
        back_inactive_check[which(chs %in% chi)] <- TRUE
      }
      #}
    } else {
      #shps <- 9:18
      #should update varlst_0 and varlst_z_0 here
      varlst_0 <- varlst
      varlst_z_0 <- varlst_z
      
      if(attr(x, 'type') %in% c('num', 'numeric')){
        shps <- 9:17
      }
      if(attr(x, 'type') %in% c('ord', 'ordinal')){
        shps <- c(1:8, 0)
      }
      if(attr(x, 'type') %in% c('tree')){
        shps <- c(1:2, 0)
      }
      #new: if we use cbic and wt.iter then we only have 4 options
      if(wt.iter & infocrit == "CBIC"){
        if(attr(x, 'type') %in% c('num', 'numeric')){
          shps <- c(9:12, 17)
        }
        if(attr(x, 'type') %in% c('ord', 'ordinal')){
          shps <- c(1:4, 0)
        }
        if(attr(x, 'type') %in% c('tree')){
          shps <- c(1:2, 0)
        }
      }
      shps_not_tried <- shps[!shps %in% rm_chi[chi]]
      #new: speed up the code
      if(attr(x, 'type') %in% c('num', 'numeric')){
        shps_not_tried <- 17
      } else {
        shps_not_tried <- 0
      }
      #cic_st = mod_tried[nrep, totx + 2]
      #test more!
      nobs <- length(shps_not_tried) + 1 #include the chosen shape too
      obs <- 1:nobs
      gmat_lst <- vector("list", length = nobs) #accumulative
      np_ix_lst <- vector("list", length = nobs)
      #mu0_lst <- vector("list", length = nobs)
      ms_lst <- vector("list", length = nobs)
      dd_lst <- vector("list", length = nobs) #for new
      kts_lst <- vector("list", length = nobs)
      
      dmat_lst <- vector("list", length = nobs) #accumulative
      dm_lst <- vector("list", length = nobs) #for new
      ps_vec <- rep(0, nobs)
      
      #new
      #st_new_vec <- rep(0, nobs)
      #ed_new_vec <- rep(0, nobs)
      
      #add the chosen shp first
      if(rm_chi[chi] %in% c(17, 0)){
        dd_lst[1] <- list(NULL)
        #dmat_lst[1] <- list(NULL)
      } else if (rm_chi[chi] %in% c(18) | rm_chi[chi] == 2 & attr(x, 'type') == 'tree'){
        if(rm_chi[chi] %in% c(18)){
          dd_lst[[1]] <- cbind(x - mean(x))
        } else {
          dd_lst[[1]] <- ddkeep_z_all_0[, which(varlst_z == chi), drop=FALSE]
        }
        #dmat_lst[1] <- list(NULL)
      } else {
        dd_lst[[1]] <- ddkeep_x_all_0[, which(varlst == chi), drop=FALSE]
        if(pnt){
          dm_lst[[1]] <- dmat_keep_lst_0[[chi]]
        }
        #dmat_lst[1] <- dmat_keep_all_0[, which(varlst == chi), drop=FALSE]
      }
      gmat_lst[[1]] <- cbind(one, ddkeep_z_all_0, ddkeep_x_all_0)
      np_ix_lst[[1]] <- NCOL(ddkeep_z_all_0) + 1
      
      #new:
      if(pnt & !is.null(dmat_keep_all_0)){
        dmat_lst[[1]] <- dmat_keep_all_0
      }
      
      #record gmat, np_ix, etc
      for(i in 2:nobs)
      {
        shpi <- shps_not_tried[i-1]
        if(rm_chi[chi] %in% c(17, 18, 0) | rm_chi[chi] == 2 & attr(x, 'type') == 'tree'){
          ddkeep_x_all <- ddkeep_x_all_0
          if(pnt){
            dmat_keep_all <- dmat_keep_all_0
          }
        } else {
          ddkeep_x_all <- ddkeep_x_all_0[, -which(varlst == chi), drop=FALSE]
          if(pnt){
            dmat_keep_all <- Matrix::bdiag(Filter(Negate(is.null), dmat_keep_lst_0[seq_added[-which(seq_added %in% chi)]]))
            dmat_keep_all <- as.matrix(dmat_keep_all)
            if(prod(dim(dmat_keep_all)) == 0){
              dmat_keep_all <- NULL
            }
          }
        }
        #if the chosen shape is conv or conc, then need to remove from ddkeep_z_all
        if(rm_chi[chi] > 10 & rm_chi[chi] < 13 | rm_chi[chi] == 18 | rm_chi[chi] > 2 & rm_chi[chi] < 5 | rm_chi[chi] == 2 & attr(x, 'type') == 'tree'){
          ddkeep_z_all <- ddkeep_z_all_0[, -which(varlst_z == chi), drop=FALSE]
        } else {
          ddkeep_z_all <- ddkeep_z_all_0
        }
        
        if (shpi == 17 | shpi == 0) {
          ms <- NULL
          dd <- NULL
          kts <- NULL
          gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all)
          np_ix <- NCOL(ddkeep_z_all) + 1
          if(pnt){
            dmat <- dmat_keep_all
            #wrong! shp not tried is 17, then cic_val is not oldcic
            #cic_val <- oldcic
            #not_skip_cic <- FALSE
            dm <- NULL
            ps <- 0
          }
        } else if (shpi == 18) {
          #dd <- cbind(x)
          ms <- mean(x)
          dd <- cbind(x - mean(x))
          kts <- NULL
          gmat <- cbind(one, ddkeep_z_all, dd, ddkeep_x_all)
          np_ix <- NCOL(ddkeep_z_all) + 2
          if(pnt){
            dmat <- dmat_keep_all
            dm <- NULL
            ps <- 0
          }
          #not_skip_cic <- TRUE
        } else {
          
          if(shpi <= 8){
            
            if(attr(x, 'type') %in% c('ord')){
              dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = 0, space = space)
              ms <- dd_ans$ms
              dd <- t(dd_ans$amat)
              if(shpi > 2 & shpi < 5){
                gmat <- cbind(one, ddkeep_z_all, x-mean(x), ddkeep_x_all, dd)
                np_ix <- NCOL(ddkeep_z_all) + 2
              } else {
                gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all, dd)
                np_ix <- NCOL(ddkeep_z_all) + 1
              }
              dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
            }
            
            if(attr(x, 'type') %in% c('tree')){
              dd <- t(tree.fun(x))
              ms <- apply(dd, 2, mean)
              dd <- apply(dd, 2, function(e) e - mean(e))
              if(shpi == 1){
                gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all, dd)
                np_ix <- NCOL(ddkeep_z_all) + 1
                dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
              }
              if(shpi == 2){
                gmat <- cbind(one, ddkeep_z_all, dd, ddkeep_x_all)
                np_ix <- NCOL(ddkeep_z_all) + 1 + NCOL(dd)
                dm <- NULL
              }
            }
            
            ps <- 0
            dmat <- dmat_keep_all
            if(!is.null(dmat)){
              dm_lst_tmp <- vector('list', length = 2)
              dm_lst_tmp[[1]] <- dmat
              dm_lst_tmp[[2]] <- dm
              dmat <- as.matrix(bdiag(dm_lst_tmp))
            } else {
              dmat <- dm
            }
          }
          
          if(shpi > 8){
            mult <- ifelse(pnt, 2, 1)
            spl_degree_i <- ifelse(shpi %in% c(9, 10), 2, 3) #use 2 for i-spl, 3 for c-spl
            n1 <- length(unique(x))
            #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/5)) + 5, trunc(n1^(1/7)) + 5)
            #new
            #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/5)) + 5, trunc(n1^(1/7)) + 5)
            nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/7)) + 6 + add_kts, trunc(n1^(1/9)) + 6 + add_kts)
            edfu <- ifelse(spl_degree_i == 2, nkts/mult, nkts/mult + 1) + add_edfu
            kts <- seq.int(min(x), max(x), length.out = nkts)
            
            #new: boundary knots
            kts0 <- boundary.knots[[chi]]
            if(length(kts0) >= 2){
              x_rng <- range(x, na.rm = TRUE)
              kts0_rng <- range(kts0, na.rm = TRUE)
              
              tol <- 1e-10
              reset_lwr <- kts0_rng[1] > x_rng[1] + tol
              reset_upp <- kts0_rng[2] < x_rng[2] - tol
              
              if(reset_lwr){
                kts_lwr <- x_rng[1]
              } else {
                kts_lwr <- kts0_rng[1]
              }
              if(reset_upp){
                kts_upp <- x_rng[2]
              } else {
                kts_upp <- kts0_rng[2]
              }
              
              kts[1] <- kts_lwr
              kts[nkts] <- kts_upp
            }
            
            dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = kts, space = space)
            #dd is orth to 1
            ms <- dd_ans$ms
            dd <- t(dd_ans$amat)
            
            if(shpi > 10 & shpi < 13){
              gmat <- cbind(one, ddkeep_z_all, x-mean(x), ddkeep_x_all, dd)
              np_ix <- NCOL(ddkeep_z_all) + 2
            } else {
              gmat <- cbind(one, ddkeep_z_all, ddkeep_x_all, dd)
              np_ix <- NCOL(ddkeep_z_all) + 1
            }
            
            if(pnt){
              nc <- NCOL(dd)
              dm <- makedmat_1D_ic(m = nc, knots = kts, shpi = shpi) #different from makedmat_1D for testpar
              qvi <- crossprod(dd)
              dvi <- crossprod(dm)
              #ps <- 1e-6
              #9, 10, 11, 12, 13, 14, 15, 16
              #if(shpi %in% c(9, 11, 13)){
              #ps <- uniroot(f = .search_ps, qv0 = qvi, dv = dvi, dd = dd,
              #              edfu = edfu, interval = c(1e-10, 1e+10), tol = .Machine$double.eps^0.25)$root
              #max_log_pen <- log(NROW(dd))
              log_root_res <- uniroot(f = function(log_pen,...) .search_ps_shapeselect(exp(log_pen),...),
                                      interval = c(-8, 8), 
                                      qv0 = qvi, dv = dvi, dd = dd, edfu = edfu,
                                      extendInt = "yes",
                                      tol = .Machine$double.eps^0.25)
              # convert back to the original scale
              #temp
              if(abs(log_root_res$root) > 8) {
                ps <- floor(exp(8))
              } else {
                ps <- exp(log_root_res$root)
              }
              #ps <- exp(log_root_res$root)
              
              ps_shpi <- which(shps_not_tried %in% shpi) + 1 #1st is for shp chosen
              #print (ps_shpi)
              ps_vec[ps_shpi] <- ps
              #} else if (shpi == 10){
              #   ps <- ps_vec[2]
              #   ps_shpi <- which(shps_not_tried %in% shpi) + 1 #1st is for shp chosen
              #   ps_vec[ps_shpi] <- ps
              #} else if (shpi == 12) {
              #   ps <- ps_vec[4]
              #   ps_shpi <- which(shps_not_tried %in% shpi) + 1 #1st is for shp chosen
              #   ps_vec[ps_shpi] <- ps
              #} else {
              #   ps <- ps_vec[6]
              #   ps_shpi <- which(shps_not_tried %in% shpi) + 1 #1st is for shp chosen
              #   ps_vec[ps_shpi] <- ps
              #}
              
              dm <- sqrt(ps) * dm
              dmat <- dmat_keep_all
              if(!is.null(dmat)){
                dm_lst_tmp <- vector('list', length = 2)
                dm_lst_tmp[[1]] <- dmat
                dm_lst_tmp[[2]] <- dm
                dmat <- as.matrix(bdiag(dm_lst_tmp))
              } else {
                dmat <- dm
              }
            }
            #not_skip_cic <- TRUE
          }
        }
        
        if(shpi == 17 | shpi == 0) {
          dd_lst[i] <- list(NULL)
          ms_lst[i] <- list(NULL)
          kts_lst[i] <- list(NULL)
        } else {
          if(shpi == 18) {
            #dd_lst[i] <- list(NULL)
            dd_lst[[i]] <- dd
            ms_lst[[i]] <- ms
            kts_lst[i] <- list(NULL)
          } else {
            dd_lst[[i]] <- dd
            ms_lst[[i]] <- ms
            #NULL for ordinal
            if(shpi > 8){
              kts_lst[[i]] <- kts
            }
          }
        }
        
        if(pnt & !is.null(dm)){
          dm_lst[[i]] <- dm
        }
        
        gmat_lst[[i]] <- gmat
        np_ix_lst[[i]] <- np_ix
        #mu0_lst[[i]] <- mu0
        
        if(pnt & !is.null(dmat)){
          dmat_lst[[i]] <- dmat
        }
      }
      
      shp17_id <- NULL
      if(any(shps_not_tried %in% c(17, 0))){
        #shp17_id = which(shps_not_tried %in% 17)
        shp17_id <- which(shps_not_tried %in% c(17, 0)) + 1 #1st one is for chosen shp
      }
      
      dfmean_rslt <- compute_edf(family=family, gmat_lst=gmat_lst, dmat_lst=dmat_lst,
                                 mu0_lst=NULL, np_ix_lst=np_ix_lst,
                                 weights=NULL, shp17_id=shp17_id, nsim=nsim, wt.iter=wt.iter, 
                                 parallel=parallel, pnt=pnt)
      dfmean_caic <- dfmean_rslt$dfmean_caic #if 17 is in shps not tried, the caic has two vals, the 2nd one is for 17
      dfmean_cbic <- dfmean_rslt$dfmean_cbic
      #print (dfmean_cbic)
      #go through all shapes to compute rss and cic/bic
      for(shpi in shps_not_tried){
        #print (shpi)
        ps_shpi <- which(shps_not_tried %in% shpi) + 1 #1st is for shp chosen
        #ps_shpi = which(9:18 %in% shpi)
        #new
        #st_new <- st_new_vec[ps_shpi]
        #ed_new <- ed_new_vec[ps_shpi]
        
        np_ix <- np_ix_lst[[ps_shpi]]
        gmat <- gmat_lst[[ps_shpi]]
        ms <- ms_lst[[ps_shpi]]
        if(shpi %in% c(17, 0)) {
          dd <- NULL
          dm <- NULL
        } else {
          dd <- dd_lst[[ps_shpi]]
          if(pnt){
            dm <- dm_lst[[ps_shpi]]
          }
        }
        kts <- kts_lst[[ps_shpi]]
        
        if(pnt){
          dmat <- dmat_lst[[ps_shpi]]
          ps <- ps_vec[ps_shpi]
        }
        
        nc <- ncol(gmat)
        amat <- diag(nc)
        #amat[1:np_ix, 1:np_ix] = 0
        ps_zero <- 1:np_ix
        amat[ps_zero, ps_zero] <- 0
        bvec <- cbind(rep(0, nc))
        xtx <- crossprod(gmat)
        if(pnt & !is.null(dmat)){
          dm_zero <- matrix(0, nrow = nrow(dmat), ncol = np_ix)
          dmat <- cbind(dm_zero, dmat)
          dv <- crossprod(dmat)
          xtx <- xtx + dv
        }
        
        if(np_ix > 0 & NCOL(gmat) > np_ix) {# z and x
          zvec <- y
          if(!wt.iter){
            #test
            if(!pnt | pnt & all(dmat == 0)){
              dsend <- gmat[, -c(1:np_ix), drop = FALSE]
              zsend <- gmat[, 1:np_ix, drop = FALSE]
              ans <- coneB(zvec, dsend, zsend)
              bhat <- ans$coefs
              bhat = round(bhat, 10)
              etahat <- gmat %*% bhat
              face <- ans$face
            } else {
              xty <- crossprod(gmat, zvec)
              #test
              # if(qr(xtx)$rank < NCOL(xtx)){
              #   xtx <- xtx + diag(1e-8, NCOL(xtx))
              # }
              # calculate eigenvalues
              ev <- eigen(xtx, only.values = TRUE)$values
              # check if the smallest eigenvalue is too close to zero or negative
              if (min(ev) < 1e-8) {
                diag(xtx) <- diag(xtx) + 1e-8
              }
              ans <- solve.QP(xtx, xty, amat, bvec)
              #ans = solve.QP(xtx, xty, t(amat), bvec)
              bhat <- ans$solution
              bhat = round(bhat, 10)
              etahat <- gmat %*% bhat
              face <- ans$iact
            }
          } else {
            gmat0 <- gmat
            weights <- rep(1, n)
            etahat <- etahat.fun(n, y, fml = family$family)
            gr <- gr.fun(y, etahat, weights, fml = family$family)
            wt <- wt.fun(y, etahat, n, weights, fml = family$family)
            cvec <- wt * etahat - gr
            zvec <- zvec.fun(cvec, wt, y, fml = family$family)
            
            muhat <- family$linkinv(etahat)
            #sm <- 1e-5
            sm <- 1e-7
            if (family$family == "binomial") {
              mdiff <- abs(max(muhat) - 1) > sm
            } else {mdiff <- TRUE}
            
            diff <- 1
            niter <- 0
            
            while (diff > sm  & mdiff & niter < 20)
            {
              oldmu <- muhat
              niter <- niter + 1
              gr <- gr.fun(y, etahat, weights=NULL, fml = family$family)
              wt <- wt.fun(y, etahat, n, weights=NULL, fml = family$family)
              
              cvec <- wt * etahat - gr
              zvec <- zvec.fun(cvec, wt, y, fml = family$family)
              gmat <- sweep(gmat0, MARGIN = 1, sqrt(wt), FUN = '*')
              
              #test
              if(!pnt | pnt & all(dmat == 0)){
                dsend <- gmat[, -c(1:np_ix), drop = FALSE]
                zsend <- gmat[, 1:np_ix, drop = FALSE]
                ansi <- coneB(zvec, dsend, zsend)
                bhat <- ansi$coefs
              }
              
              if(pnt & !is.null(dmat) & !all(dmat == 0)){
                xtx <- crossprod(gmat)
                xtx <- xtx + dv
                xty <- crossprod(gmat, zvec)
                #test
                # if(qr(xtx)$rank < NCOL(xtx)){
                #   xtx <- xtx + diag(1e-8, NCOL(xtx))
                # }
                # calculate eigenvalues
                ev <- eigen(xtx, only.values = TRUE)$values
                # check if the smallest eigenvalue is too close to zero or negative
                if (min(ev) < 1e-8) {
                  diag(xtx) <- diag(xtx) + 1e-8
                }
                ansi <- solve.QP(xtx, xty, amat, bvec)
                #ansi = solve.QP(xtx, xty, t(amat), bvec)
                bhat <- ansi$solution
                face <- ansi$iact
              }
              
              bhat = round(bhat, 10)
              etahat <- gmat0 %*% bhat
              muhat <- family$linkinv(etahat)
              mdiff <- abs(max(muhat) - 1)
              if (family$family == "binomial") {
                mdiff <- abs(max(muhat) - 1) > sm
              } else {mdiff <- TRUE}
              
              diff <- mean((muhat - oldmu)^2)
            }
          }
          #new
          # coef_xc_new <- NULL
          # coef_x_new <- NULL
          # if(shpi > 10 & shpi < 13){
          #   coef_xc_new <- bhat[np_ix]
          #   coef_x_new <- bhat[st_new:ed_new]
          # } else {
          #   #coef_xc_new <- NULL
          #   coef_x_new <- bhat[st_new:ed_new]
          # }
          
        } else if (np_ix > 0 & NCOL(gmat) == np_ix) { #all z #shpi = 18 and no constrained x has been chosen
          zvec <- y
          if(!wt.iter){
            #bhat <- solve(crossprod(gmat), t(gmat)) %*% zvec
            bhat <- ginv(crossprod(gmat)) %*% t(gmat) %*% zvec
            bhat = round(bhat, 10)
            etahat <- gmat %*% bhat
          } else {
            m0 <- stats::glm(y ~ gmat - 1, family = family)
            etahat <- m0$linear.predictors
            bhat <- coef(m0)
            muhat <- family$linkinv(etahat)
          }
        }
        
        rss <- sum((y - etahat)^2)
        if(wt.iter){
          llh <- llh.fun(y, muhat, etahat, phihat=NULL, n, weights=NULL, fml = family$family)
          #test!
        } else {llh <- log(rss)}
        
        ps_shpi_bic <- which(shps_not_tried %in% shpi) + 1 #1st is for chosen shp
        
        if (infocrit == "CBIC" & shpi %in% c(17, 0)) {
          #cic_val <- llh + log(n) * (dfmean_caic[2]) / n
          #new
          #cic_val <- llh + a * log(n) * (dfmean_caic[2]) / n^(6/7)
          
          all_shps <- shps_all
          all_shps[chi] <- shpi
          if(any(all_shps %in% c(9, 10))){ #smooth incr/decr
            degree <- 6/7
          } else if (any(all_shps > 10 & all_shps < 17)){ #smooth conv/conc
            degree <- 6/7 #8/9#6/7#
          } else if (any(all_shps > 0 & all_shps < 9)) { #ord or parametric; e.g., all_shps = c(8,17)
            degree <- 6/7 #1# 6/7#
          } else if (all(all_shps %in% c(0, 17))) {
            stop('check!')
          }
          cic_val <- llh + a * log(n) * (dfmean_caic[2]) / n^(degree)
          
        }
        if (infocrit == "CBIC" & shpi != 17 & shpi != 0) {
          #cic_val <- llh + log(n) * (dfmean_cbic[ps_shpi_bic]) / n
          #new:
          #cic_val <- llh + a * log(n) * (dfmean_cbic[ps_shpi_bic]) / n^(6/7)
          
          all_shps <- shps_all
          all_shps[chi] <- shpi
          if(any(all_shps %in% c(9, 10))){ #smooth incr/decr
            degree <- 6/7
          } else if (any(all_shps > 10 & all_shps < 17)){ #smooth conv/conc
            degree <- 6/7# 8/9#6/7#
          } else if (any(all_shps > 0 & all_shps < 9)) { #ord or parametric
            degree <- 6/7# 1#6/7#
          } else {
            stop('check!')
          }
          cic_val <- llh + a * log(n) * (dfmean_cbic[ps_shpi_bic]) / n^(degree)
        }
        if (infocrit == "CAIC" & shpi %in% c(17, 0)) {
          cic_val <- llh + 2 * (dfmean_caic[2]) / n
        }
        if (infocrit == "CAIC" & shpi != 17 & shpi != 0) {
          #test
          cic_val <- llh + 2 * (dfmean_cbic[ps_shpi_bic]) / n
        }
        
        #test more!
        bool <- (cic_val - cic_st > 0) & (cic_val - cic_st) < 1e-8
        #new: if etahat = oldetahat, then we should remove x even though cic_val > cic_st?
        bool2 <- isTRUE(all.equal(etahat, etahatkeep, tolerance = 1e-8))
        bool <- bool | bool2
        if(cic_val < cic_st | bool){
          #cat(paste0('x', chi), 'changed in backward!', 'nrep:', nrep, '\n')
          changed <- TRUE
          new_shp <- shpi
          #new:
          if(new_shp == 17){
            new_shp <- 0
          }
          #new: put chi back to active set
          if(shpi == 17 | shpi == 0) {
            back_inactive_check[which(chs %in% chi)] <- TRUE
            #cat('x', chi, 'put back in inactive set!', 'nrep:', nrep, '\n')
          }
          mod_tried[nrep, chi] <- new_shp
          #new: check more
          mod_tried[nrep, totx+2] <- cic_val
          etahatkeep <- etahat
          #print (head(etahatkeep))
          coefskeep <- bhat
          
          if(shpi != 18) {
            #new:
            if(shpi != 2 | attr(x, 'type') != 'tree'){
              ddkeep_x <- dd
            } else {
              ddkeep_z <- dd
              ddkeep_x <- NULL
            }
            
            if(pnt){
              dmat_keep <- dm
              ps_keep <- ps
            }
            ktskeep <- kts
            if(shpi == 17 | shpi == 0) {
              ddkeep_x_lst[chi] <- list(NULL)
              #new
              if(any(varlst_0 == chi)){
                varlst_0 <- varlst_0[-which(varlst_0 == chi)]
              }
              
              #test more
              ddkeep_z_lst[chi] <- list(NULL)
              #new
              if(any(varlst_z_0 == chi)){
                varlst_z_0 <- varlst_z_0[-which(varlst_z_0 == chi)]
              }
              
              if(pnt){
                dmat_keep_lst[chi] <- list(NULL)
                ps_keep_vec[chi] <- 0
              }
              
              mskeep_lst[chi] <- list(NULL)
              kts_lst[chi] <- list(NULL)
            } else {
              if(shpi == 2 & attr(x, 'type') == 'tree'){
                ddkeep_x_lst[chi] <- list(NULL)
                if(any(varlst_0 == chi)){
                  varlst_0 <- varlst_0[-which(varlst_0 == chi)]
                }
                if(pnt){
                  dmat_keep_lst[chi] <- list(NULL)
                  ps_keep_vec[chi] <- 0
                }
                mskeep_lst[[chi]] <- ms
              } else {
                ddkeep_x_lst[[chi]] <- dd
                #new
                if(any(varlst_0 == chi)){
                  varlst_0 <- varlst_0[-which(varlst_0 == chi)]
                }
                varlst_0 <- c(varlst_0, rep(chi, ncol(dd)))
                #new
                #coefskeep_x_lst[[chi]] <- coef_x_new
                
                if(pnt){
                  dmat_keep_lst[[chi]] <- dm
                  ps_keep_vec[chi] <- ps
                }
                mskeep_lst[[chi]] <- ms
                #NULL for ordinal
                if(shpi > 8){
                  knots_lst[[chi]] <- kts
                }
              }
            }
            
            if(shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5 | shpi == 2 & attr(x, 'type') == 'tree'){
              if(shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5){
                ddkeep_z <- cbind(x - mean(x))
                ddkeep_z_lst[[chi]] <- cbind(x - mean(x))
                #new
                if(any(varlst_z_0 == chi)){
                  varlst_z_0 <- varlst_z_0[-which(varlst_z_0 == chi)]
                }
                varlst_z_0 <- c(varlst_z_0, rep(chi, 1))
                np_add <- 1
              } else {
                ddkeep_z <- dd
                ddkeep_z_lst[[chi]] <- dd
                #new
                if(any(varlst_z_0 == chi)){
                  varlst_z_0 <- varlst_z_0[-which(varlst_z_0 == chi)]
                }
                varlst_z_0 <- c(varlst_z_0, rep(chi, NCOL(dd)))
                np_add <- NCOL(dd)
              }
            } else {
              ddkeep_z <- NULL
              ddkeep_z_lst[chi] <- list(NULL)
              #new
              if(any(varlst_z_0 == chi)){
                varlst_z_0 <- varlst_z_0[-which(varlst_z_0 == chi)]
              }
              #new
              #coefskeep_z_lst[chi] <- list(NULL)
              np_add <- 0
            }
            #varlst_lst[[chi]] <- rep(chi, NCOL(dd))
          }
          
          if (shpi == 18) {
            ddkeep_z <- cbind(x - mean(x))
            ddkeep_z_lst[[chi]] <- cbind(x - mean(x))
            #new
            if(any(varlst_z_0 == chi)){
              varlst_z_0 <- varlst_z_0[-which(varlst_z_0 == chi)]
            }
            varlst_z_0 <- c(varlst_z_0, rep(chi, 1))
            
            mskeep_lst[[chi]] <- mean(x)
            np_add <- 1
            knots_lst[chi] <- list(NULL)
            ddkeep_x <- NULL
            ddkeep_x_lst[chi] <- list(NULL)
            #new
            if(any(varlst_0 == chi)){
              varlst_0 <- varlst_0[-which(varlst_0 == chi)]
            }
            if(pnt){
              dmat_keep <- NULL
              ps_keep <- 0
              dmat_keep_lst[chi] <- list(NULL)
              ps_keep_vec[chi] <- 0
            }
          }
          
          #new: need to update ddkeep_z_all_0 and ddkeep_x_all_0 here in case >= 2 predictors' shapes are changed
          seq_added <- mod_tried[1:nrep, totx + 3]
          # ddkeep_x_lst_ord will only store dd's for chosen x
          ddkeep_x_lst_ord <- ddkeep_x_lst[seq_added]
          ddkeep_z_lst_ord <- ddkeep_z_lst[seq_added]
          
          if(pnt){
            dmat_keep_lst_0 <- dmat_keep_lst
            #dmat_keep_all_0 <- as.matrix(bdiag(Filter(Negate(is.null), dmat_keep_lst_0)))
            #new!
            #dmat_keep_all_0 <- as.matrix(bdiag(dmat_keep_lst_0[seq_added]))
            dmat_keep_all_0 <- as.matrix(bdiag(Filter(Negate(is.null), dmat_keep_lst_0[seq_added])))
          }
          
          ddkeep_x_all_0 <- do.call(base::cbind, ddkeep_x_lst_ord)
          ddkeep_z_all_0 <- do.call(base::cbind, ddkeep_z_lst_ord)
          
          obs <- 1:totx
          obs_ord <- obs[seq_added]
          
          #new: not in the order of x1,x2...; the chi one is at the end.
          coefs_nointercept <- coefskeep[-1]
          #print (varlst_z_0)
          if(length(varlst_z_0) > 0){
            zid <- 1:length(varlst_z_0)
            zcoefs <- coefs_nointercept[zid]
            xcoefs <- coefs_nointercept[-zid]
          }
          
          if(length(varlst_z_0) == 0){
            zcoefs <- NULL
            xcoefs <- coefs_nointercept
          }
          
          #cat('chi:', chi, '\n')
          #cat('shpi:', shpi, '\n')
          #print (varlst_0)
          #print (zcoefs)
          #new: we should update ddkeep_x_all and ddkeep_z_all here
          ddkeep_x_all_temp <- do.call(base::cbind, ddkeep_x_lst)
          ddkeep_z_all_temp <- do.call(base::cbind, ddkeep_z_lst)
          accx <- accz <- 0
          accx <- NCOL(ddkeep_x_all_temp)
          accz <- NCOL(ddkeep_z_all_temp)
          
          totx <- NCOL(xmat)
          #new
          coefskeep_x_lst <- vector("list", length = totx)
          coefskeep_z_lst <- vector("list", length = totx)
          
          if(accx > 0){
            varlst_lst <- lapply(1:length(obs_ord), function(i, xlst) {if(NCOL(xlst[[i]]) >= 1) {rep(obs_ord[i], NCOL(xlst[[i]]))}}, xlst = ddkeep_x_lst_ord)
            #print (varlst_lst)
            varlst <- do.call(base::c, varlst_lst)
            
            #new
            #print (varlst_0)
            ux_added <- unique(varlst_0) 
            for(x_added in ux_added){
              coefskeep_x_lst[[x_added]] <- xcoefs[which(varlst_0 %in% x_added)]
            }
            
            #test!?
            ord_x <- order(varlst_0)
            xcoefs <- xcoefs[ord_x]
          } else {varlst <- NULL}
          
          if(accz > 0){
            varlst_z_lst <- lapply(1:length(obs_ord), function(i, zlst) {if(NCOL(zlst[[i]]) >= 1) {rep(obs_ord[i], NCOL(zlst[[i]]))}}, zlst = ddkeep_z_lst_ord)
            varlst_z <- do.call(base::c, varlst_z_lst)
            
            #new
            uz_added <- unique(varlst_z_0)
            for(z_added in uz_added){
              coefskeep_z_lst[[z_added]] <- zcoefs[which(varlst_z_0 %in% z_added)]
            }
            #test!? correct; match order of columns in xmat
            ord_z <- order(varlst_z_0)
            zcoefs <- zcoefs[ord_z]
          } else {varlst_z <- NULL}
          
          #new
          # totx <- NCOL(xmat)
          # coefskeep_x_lst <- vector("list", length = totx)
          # coefskeep_z_lst <- vector("list", length = totx)
          # if(length(varlst_0) > 0){
          #   ux_added <- unique(varlst_0) 
          #   for(x_added in ux_added){
          #     coefskeep_x_lst[[x_added]] <- xcoefs[which(varlst_0 %in% x_added)]
          #   }
          # }
          
          # if(length(varlst_z_0) > 0){
          #   uz_added <- unique(varlst_z_0)
          #   for(z_added in uz_added){
          #     coefskeep_z_lst[[z_added]] <- zcoefs[which(varlst_z_0 %in% z_added)]
          #   }
          # }
          
          #new: put b0 back
          coefskeep <- c(coefskeep[1], zcoefs, xcoefs)
          
          cic_st <- cic_val
          
          if (infocrit == "CBIC") {
            #degree is determined earlier
            cic_st_next_step_0 <- llh + a * log(n) * (dfmean_cbic[ps_shpi_bic]) / n^(degree)
          }
          
          if (infocrit == "CAIC") {
            #cic_val <- llh + 2 * (dfmean_caic) / n
            cic_st_next_step_0 <- llh +  2 * (dfmean_cbic[ps_shpi_bic]) / n
          }
          
          if(cic_st_next_step_0 < oldcic){
            cic_st_next_step <- cic_st_next_step_0
          }
        }
      }
    }
  }
  #new:
  if(any(back_inactive_check)){
    back_inactive <- chs[back_inactive_check]
  }
  rslt = list(ddkeep_x_lst = ddkeep_x_lst, ddkeep_z_lst = ddkeep_z_lst, 
              dmat_keep_lst = dmat_keep_lst, ps_keep_vec = ps_keep_vec,
              knots_lst = knots_lst,
              mskeep_lst = mskeep_lst,
              mod_tried = mod_tried, changed = changed,
              cic_st = cic_st, cic_st_next_step = cic_st_next_step,
              back_inactive = back_inactive,
              etahatkeep = etahatkeep, coefskeep = coefskeep, 
              coefskeep_x_lst = coefskeep_x_lst, coefskeep_z_lst = coefskeep_z_lst)
  return (rslt)
}

#--------------------------------------------------------------------
.compute_edf0 = function(isim, gmat_lst, mu0_lst=NULL, np_ix_lst,
                         dmat_lst=NULL, weights=NULL,
                         wt.iter=FALSE, family=gaussian(), pnt=FALSE,
                         xtx_lst=NULL, uinv_lst=NULL,
                         xtxinv_lst=NULL, qv0_lst=NULL, dv_lst=NULL,
                         amat_lst=NULL, bvec_lst=NULL,
                         atil_lst=NULL, edfi_base_vec=NULL, shps_to_sims=NULL)
{
  cicfamily <- CicFamily(family)
  llh.fun <- cicfamily$llh.fun
  #new: use log link in gamma
  linkfun <- cicfamily$linkfun
  etahat.fun <- cicfamily$etahat.fun
  gr.fun <- cicfamily$gr.fun
  wt.fun <- cicfamily$wt.fun
  zvec.fun <- cicfamily$zvec.fun
  muhat.fun <- cicfamily$muhat.fun
  ysim.fun <- cicfamily$ysim.fun
  deriv.fun <- cicfamily$deriv.fun
  dev.fun <- cicfamily$dev.fun
  
  nobs <- length(gmat_lst)
  #9:16, 18 #17 is in the list in the backward step
  obs <- 1:nobs #9  or 10 shapes, including linear, will be fit in each isim?
  eps <- 1e-8
  edfs <- rep(0, nobs)
  
  # ysim.fun <- function(n, mu0 = NULL, fml = object$family, shp0 = NULL) {
  #   if (fml == "binomial") {
  #     ysim <- 1:n*0
  #     ysim[runif(n) < .5] <- 1
  #   }
  #   if (fml == "poisson") {
  #     if (!is.null(mu0)) {
  #       ysim <- rpois(n, mu0)
  #     }
  #   }
  #   if (fml == "gaussian") {
  #     ysim <- rnorm(n)
  #   }
  #   if (fml == "Gamma") {
  #     ysim <- rgamma(n, shape=1)
  #   }
  #   ysim
  # }
  #new:
  if(nobs > 1) {
    if(is.null(shps_to_sims)){
      if(nobs == 9) {
        obs <- c(1, 3, 5, 9)
      }
      if(nobs == 8) {
        obs <- c(1, 3, 5)
      }
      #new for bin and cbic
      if(nobs == 4) {
        obs <- c(1, 3)
      }
      if(nobs == 5) { #backward step for bin and cbic
        obs <- c(1, 3, 5)
      }
      #if(nobs == 10) {
      #  obs <- c(1, 3, 5, 9, 10)
      #}
    } else {
      if(nobs == 2 | nobs == 3){
        obs <- 1:nobs
      } else {
        mono_ps <- mono_ps_other <- mono_ps_1 <- NULL
        conv_ps <- conv_ps_other <- conv_ps_1 <- NULL
        if(any(shps_to_sims %in% c(1, 2, 9, 10))){
          mono_ps <- sort(which(shps_to_sims %in% c(1, 2, 9, 10)))
          mono_ps_1 <- mono_ps[1]
          if(length(mono_ps) > 1) {
            mono_ps_other <- mono_ps[-1]
          }
        }
        
        if(any(shps_to_sims > 2 & shps_to_sims < 9) | any(shps_to_sims > 10 & shps_to_sims < 17)){
          conv_ps <- sort(which(shps_to_sims > 2 & shps_to_sims < 9 | shps_to_sims > 10 & shps_to_sims < 17))
          conv_ps_1 <- conv_ps[1]
          if(length(conv_ps) > 1){
            conv_ps_other <- conv_ps[-1]
          }
        }
        obs <- c(mono_ps_1, conv_ps_1)
      }
    }
  }
  
  for(i in obs){
    #print(i)
    gmat <- gmat_lst[[i]]
    #mu0 <- mu0_lst[[i]]
    #mu0 <- rep(1, NROW(gmat))
    np_ix <- np_ix_lst[[i]]
    
    #new
    amat <- amat_lst[[i]]
    bvec <- bvec_lst[[i]]
    xtx <- xtx_lst[[i]]
    xtxinv <- xtxinv_lst[[i]]
    uinv <- uinv_lst[[i]]
    qv0 <- qv0_lst[[i]]
    atil <- atil_lst[[i]]
    edfi_base <- edfi_base_vec[i]
    
    #new
    if(pnt){
      dmat <- dmat_lst[[i]] #ignored np_ix
      #kts <- knots_lst[[i]] #used for bin
      dv <- dv_lst[[i]] #considered np_ix
    } else {
      dmat <- NULL
    }
    if(NCOL(gmat) == np_ix) {
      #print ('no sim')
      #if(!wt.iter){
      edfi <- np_ix
      #} else {
      #  edfi <- np_ix + 1 #include intercept and remove later
      #}
    }
    
    if (NCOL(gmat) > np_ix) {
      n <- NROW(gmat)
      #test
      #mu0 <- rep(0, n)
      #if(!wt.iter){
      #  ysim <- mu0 + rnorm(n)
      #ysim <- rnorm(n)
      #} else {
      ysim <- ysim.fun(n, mu0 = NULL, fml = family$family, shp0 = NULL)
      #}
      
      #-------------------------------------------
      #use qprog
      xtx <- qv0
      if(pnt & !is.null(dmat)){
        #test
        dmat_zero <- matrix(0, nrow = nrow(dmat), ncol = np_ix)
        dmat2 <- cbind(dmat_zero, dmat)
        dv <- crossprod(dmat2)
        xtx <- qv0 + dv
      }
      #-------------------------------------------
      #if(np_ix > 0){
      if(!wt.iter){
        #ans <- suppressMessages(try(coneB(ysim, dsend, zsend), silent = TRUE))
        #-------------------------------------------
        #use qprog
        #gmat0 = gmat
        if (!pnt | pnt & all(dmat == 0)) {
          zvec <- ysim
          dsend <- gmat[, -c(1:np_ix), drop = FALSE]
          zsend <- gmat[, 1:np_ix, drop = FALSE]
          #test!
          ans <- suppressMessages(try(coneB(zvec, dsend, zsend), silent = TRUE))
          if (inherits(ans, "try-error")) {
            break
          } else {
            face <- ans$face
            bhat <- ans$coefs
            bhat = round(bhat, 10)
          }
        } else {
          xty <- crossprod(gmat, ysim)
          #test
          #amat is symmetric, no need to transpose it 
          ans <- try(solve.QP(xtx, xty, amat, bvec), silent = TRUE)
          if (!inherits(ans, "try-error")) {
            bhat <- ans$solution
            bhat = round(bhat, 10)
            face <- ans$iact
          }
        }
      } else {
        gmat0 <- gmat
        #weights <- rep(1, n)
        etahat <- etahat.fun(n, ysim, fml = family$family)
        gr <- gr.fun(ysim, etahat, weights=NULL, fml = family$family)
        wt <- wt.fun(ysim, etahat, n, weights=NULL, fml = family$family)
        cvec <- wt * etahat - gr
        zvec <- zvec.fun(cvec, wt, ysim, fml = family$family)
        #muhat <- binomial(link = "logit")$linkinv(etahat)
        muhat <- family$linkinv(etahat)
        sm <- 1e-5
        #sm <- 1e-7
        
        if (family$family == "binomial") {
          mdiff <- abs(max(muhat) - 1) > sm
        } else {mdiff <- TRUE}
        
        diff <- 1
        niter <- 0
        
        while (diff > sm & mdiff & niter < 10)
        {
          oldmu <- muhat
          niter <- niter + 1
          gr <- gr.fun(ysim, etahat, weights=NULL, fml = family$family)
          wt <- wt.fun(ysim, etahat, n, weights=NULL, fml = family$family)
          #ch = muhat > 1e-5 & muhat < 1-(1e-5)
          #wt[ch] = muhat[ch]*(1-muhat[ch])
          #wt[!ch] = 1e-5
          cvec <- wt * etahat - gr
          zvec <- zvec.fun(cvec, wt, ysim, fml = family$family)
          gmat <- sweep(gmat0, MARGIN = 1, sqrt(wt), FUN = '*')
          #dsend <- gmat[, -c(1:np_ix), drop = FALSE]
          #zsend <- gmat[, 1:np_ix, drop = FALSE]
          #test!
          #ans <- suppressMessages(try(coneB(zvec, dsend, cbind(onewt, zsend)), silent = TRUE))
          
          if(!pnt | pnt & all(dmat == 0)){
            dsend <- gmat[, -c(1:np_ix), drop = FALSE]
            zsend <- gmat[, 1:np_ix, drop = FALSE]
            #test!
            ans <- suppressMessages(try(coneB(zvec, dsend, zsend), silent = TRUE))
            if (inherits(ans, "try-error")) {
              break
            } else {
              face <- ans$face
              bhat <- ans$coefs
            }
          }
          
          if(pnt & !is.null(dmat) & !all(dmat == 0)){
            xtx <- crossprod(gmat)
            xtx <- xtx + dv
            #test!
            xty <- crossprod(gmat, zvec)
            ans <- try(solve.QP(xtx, xty, amat, bvec), silent = TRUE)
            if (inherits(ans, "try-error")) {
              break
            } else {
              bhat <- ans$solution
              face <- ans$iact
            }
          }
          
          bhat = round(bhat, 10)
          etahat <- gmat0 %*% bhat
          muhat <- family$linkinv(etahat)
          
          mdiff <- abs(max(muhat) - 1)
          if (family$family == "binomial") {
            mdiff <- abs(max(muhat) - 1) > sm
          } else {mdiff <- TRUE}
          
          diff <- mean((muhat - oldmu)^2)
        }
      }
      #ans <- tryCatch({coneB(ysim, dsend, zsend)}, error = function(e) {e})
      #}
      
      if(inherits(ans, "try-error")) {
        edfi <- NULL
      } else {
        #bhat <- ans$solution
        #face <- ans$iact
        if(!pnt | pnt & all(dmat == 0)){
          edfi <- sum(abs(na.omit(bhat)) > eps)
          #edfi <- sum(abs(na.omit(ans$solution)) > eps)
        } else {
          #new: if !wt.iter, use pre-computed edfi_base; otherwise, use pre-computed xtxinv
          if(wt.iter | !wt.iter){
            umat <- chol(xtx)
            uinv <- backsolve(umat, diag(ncol(xtx)))
            atil <- amat %*% uinv
            xtxinv <- tcrossprod(uinv)
            bigp_base <- xtxinv %*% crossprod(gmat)
            edfi_base <- sum(diag(bigp_base))
          }
          
          if (all(face == 0) | length(face) == 0) {
            edfi <- edfi_base
          } else {
            smat <- -t(atil[face, ,drop = FALSE])
            #new: avoid solve error
            qr_smat <- qr(smat)
            smat <- qr.Q(qr_smat)
            pmat_polar <- tcrossprod(smat)
            if(!wt.iter){
              bigp_polar <- uinv %*% tcrossprod(pmat_polar, uinv) %*% qv0
            }else{
              bigp_polar <- uinv %*% tcrossprod(pmat_polar, uinv) %*% crossprod(gmat)
            }
            edfi <- edfi_base - sum(diag(bigp_polar))
          }
        }
      }
    }
    if(!is.null(edfi)){
      edfs[i] <- edfi
    }
  }
  #new:
  if(nobs > 1){
    
    if(nobs >= 4){
      if(is.null(shps_to_sims)){
        edfs[2] <- edfs[1]
        edfs[4] <- edfs[3]
        if(nobs > 5){
          edfs[6:8] <- edfs[5]
        }
      } else {
        if(!is.null(mono_ps_other)){
          edfs[mono_ps_other] <- edfs[mono_ps_1]
        }
        if(!is.null(conv_ps_other)){
          edfs[conv_ps_other] <- edfs[conv_ps_1]
        }
      }
    }
  }
  return (edfs)
}

#--------------------------------------------------------------------
#simulate edf by lapply or mclapply
#--------------------------------------------------------------------
compute_edf <- function(family = gaussian(link="identity"), gmat_lst, dmat_lst=NULL,
                        mu0_lst=NULL, np_ix_lst, weights=NULL,
                        shp17_id=NULL, nsim=50, wt.iter=FALSE,
                        parallel=FALSE, pnt=FALSE, shps_to_sims=NULL){
  #not faster
  #totcores = detectCores()
  #new:
  #library(RhpcBLASctl)
  # This is the "Magic Fix" for matrix math in parallel?
  #blas_set_num_threads(1) 
  #omp_set_num_threads(1)
  
  ng <- length(gmat_lst)
  xtx_lst <- vector('list', length = ng)
  xtxinv_lst <- vector('list', length = ng)
  uinv_lst <- vector('list', length = ng)
  amat_lst <- vector('list', length = ng)
  bvec_lst <- vector('list', length = ng)
  qv0_lst <- vector('list', length = ng)
  dv_lst <- vector('list', length = ng)
  atil_lst <- vector('list', length = ng)
  edfi_base_vec <- numeric(ng)
  
  for(i in 1:ng){
    gmat <- gmat_lst[[i]]
    np_ix <- np_ix_lst[[i]]
    dmat <- dmat_lst[[i]]
    
    nc <- ncol(gmat)
    amat <- diag(nc)
    ps_zero <- 1:np_ix
    amat[ps_zero, ps_zero] <- 0
    bvec <- cbind(rep(0, nc))
    qv0 <- crossprod(gmat)
    xtx <- qv0
    if(pnt & !is.null(dmat)){
      #test
      dmat_zero <- matrix(0, nrow = nrow(dmat), ncol = np_ix)
      dmat2 <- cbind(dmat_zero, dmat)
      dv <- crossprod(dmat2)
      dv_lst[[i]] <- dv #need this for binomial
      
      #test, compute this term later 
      #xtx <- qv0 + dv
      #xtx <- qv0
    }
    
    #test
    #qv0 <- Matrix(qv0, sparse = TRUE)
    qv0_lst[[i]] <- qv0
    
    #test
    #xtx <- Matrix(xtx, sparse = TRUE)
    # calculate eigenvalues
    ev <- eigen(xtx, only.values = TRUE)$values
    # check if the smallest eigenvalue is too close to zero or negative
    if (min(ev) < 1e-8) {
      diag(xtx) <- diag(xtx) + 1e-8
    }
    xtx_lst[[i]] <- xtx
    
    amat_lst[[i]] <- amat
    bvec_lst[[i]] <- bvec
    
    umat <- chol(xtx)
    uinv <- backsolve(umat, diag(ncol(xtx)))
    atil <- amat %*% uinv
    
    uinv_lst[[i]] <- uinv
    atil_lst[[i]] <- atil
    
    xtxinv <- tcrossprod(uinv)
    xtxinv_lst[[i]] <- xtxinv
    if(!wt.iter){
      bigp_base <- xtxinv %*% qv0
      edfi_base <- sum(diag(bigp_base))
      edfi_base_vec[i] <- edfi_base
    }
  }
  if(parallel){
    cores <- min(8, parallel::detectCores(logical = FALSE) - 1)
    #test!
    #cl <- makeCluster(cores)
    #edfs <- parLapply(cl, 1:nsim, function(isim){
    #  .compute_edf0(isim, gmat_lst, mu0_lst=NULL, np_ix_lst, dmat_lst, weights=NULL, wt.iter, family, pnt)
    #})
    #stopCluster(cl)
    edfs <- suppressWarnings(parallel::mclapply(1:nsim,
                                                .compute_edf0, gmat_lst, mu0_lst=NULL,
                                                np_ix_lst, dmat_lst,
                                                weights=NULL, wt.iter,
                                                family, pnt, xtx_lst, uinv_lst,
                                                xtxinv_lst, qv0_lst, dv_lst,
                                                amat_lst, bvec_lst, atil_lst,
                                                edfi_base_vec, shps_to_sims, mc.cores = cores))
    
    if(all(sapply(edfs, is.null))){
      edfs <- parallel::mclapply(1:nsim, .compute_edf0, gmat_lst, mu0_lst=NULL,
                                 np_ix_lst, dmat_lst,
                                 weights=NULL, wt.iter,
                                 family, pnt, xtx_lst, uinv_lst,
                                 xtxinv_lst, qv0_lst, dv_lst,
                                 amat_lst, bvec_lst, atil_lst,
                                 edfi_base_vec, shps_to_sims, mc.cores = 1)
    }
    
  }else{
    edfs <- lapply(1:nsim, .compute_edf0, gmat_lst, mu0_lst=NULL,
                   np_ix_lst, dmat_lst,
                   weights=NULL, wt.iter, family, pnt, xtx_lst, uinv_lst,
                   xtxinv_lst, qv0_lst, dv_lst, amat_lst, bvec_lst,
                   atil_lst, edfi_base_vec, shps_to_sims)
  }
  edfs <- Filter(Negate(is.null), edfs)
  edfs <- simplify2array(edfs)
  #print (edfs)
  #edfs is a # of shapes by nsim matrix
  #row averages will be used for CBIC
  #col maxes' average will be used for CAIC
  
  nobs <- length(gmat_lst)
  if(nobs > 1){
    if(is.null(shp17_id)){
      dfmean_cbic <- rowMeans(edfs)
      dfmean_caic <- max(dfmean_cbic)
    } else {
      dfmean_cbic <- rowMeans(edfs)
      dfmean_caic_no17 <- max(dfmean_cbic[-shp17_id])
      dfmean_caic_17 <- dfmean_cbic[shp17_id]
      dfmean_caic <- c(dfmean_caic_no17, dfmean_caic_17)
    }
    #save(edfs, file='edfs.Rda')
    rslt <- list(dfmean_caic = dfmean_caic, dfmean_cbic = dfmean_cbic)
  } else {
    dfmean <- mean(edfs)
    rslt <- list(dfmean_caic = dfmean, dfmean_cbic = dfmean)
  }
  return (rslt)
}


#--------------------------------
#visualization
#--------------------------------
#plot.shapeselect <- function(object, pages = 1, ci = FALSE, mn = TRUE, mu = FALSE,...)
plot.shapeselect <- function(x,...)
{

  object <- x
  extras <- list(...)
  
  if(exists("pages", where = extras)){
    pages <- extras$pages
  } else {
    pages <- 1
  }
  
  if(exists("ci", where = extras)){
    ci <- extras$ci
  } else {
    ci <- FALSE
  }
  
  if(exists("mn", where = extras)){
    mn <- extras$mn 
  } else {
    mn <- TRUE
  }
  
  if(exists("mu", where = extras)){
    mu <- extras$mu
  } else {
    mu <- FALSE
  }
  
  layout_for_cgam <- function(n_plots, pages = 1) {
    ppp <- ceiling(n_plots / max(pages, 1L))
    nr <- floor(sqrt(ppp))
    nc <- ceiling(ppp / nr)
    return (c(nr, nc))
  }
  
  xmat <- object$xmat
  p <- ncol(xmat)
  types <- object$types
  
  for(ix in 1:p){
    if(inherits(xmat[,ix], "factor") & types[ix] %in% "tree"){
      xmat[,ix] <- as.numeric(as.character(xmat[,ix]))
    }
    if(inherits(xmat[,ix], "factor") & types[ix] %in% "ord"){
      xmat[,ix] <- as.numeric(as.character(xmat[,ix]))
    }
  }
  
  family <- object$family
  etacomps <- object$etacomps
  #theme_use <- ggthemes::theme_solarized()
  #theme_use <- ggthemes::theme_economist()
  shps <- object$shps
  xnms <- object$xnms
  not_selected <- xnms[which(shps %in% c(0, 17))]
  
  zmat <- object$zmat
  np <- object$d0
  capk <- object$capk
  varlist <- object$varlst #used varlst in classo
  varlist <- sort(varlist) #follow xmat order, not the x's added order
  varlist_z <- object$varlst_z
  varlist_z <- sort(varlist_z)
  #dd <- object$dd #include 1
  dd <- t(object$bigmat) #only contains basis functions of selected predictors
  #if(onlySelect){
  #kept <- !shps %in% c(0, 1, 17) #remove z
  kept <- !shps %in% c(0, 17) 
  xmat <- xmat[, kept, drop = FALSE]
  types <- types[kept]
  etacomps <- etacomps[kept, ,drop = FALSE]
  shps <- shps[kept]
  xnms <- xnms[kept]
  #}
  p <- ncol(xmat)
  if(pages == 1){
    oldpar <- par(mfrow = layout_for_cgam(n_plots = p, pages = pages))
  }
  if(ci){
    cov_beta <- vcov.cgam(object)
  } else {
    cov_beta <- NULL
  }
  extras <- list(...)
  
  uvar_x <- unique(varlist)
  uvar_z <- unique(varlist_z)
  
  ncon <- 1 #for conv/conc
  ziter <- 1 #for z
  #types <- object$types
  shapes_char <- object$shapes[kept]
  #p <- ncol(xmat) #need to define it again! only use chosen xs
  for (k in 1:p) {
    x_k <- xmat[, k]
    xnm_k <- xnms[k]
    
    #shp_k <- ShapeToChar_shapeselect(shps[k], types[k])
    shp_k <- shapes_char[k]
    eta_k <- etacomps[k, ]
    
    if(mu){
      eta_k <- family$linkinv(eta_k)
    }
    
    #ps_k <- which(varlist == k)
    if(is.factor(x_k)){
      #x_k <- as.numeric(levels(x_k))[x_k] 
      ps_k <- which(varlist_z == uvar_z[ziter]) + 1 #consider constant
      ziter <- ziter + 1
    } else {
      ps_k <- which(varlist == uvar_x[k])
      ps_k <- ps_k + np
    }
    
    #range_eta <- range(eta_k, na.rm = TRUE)
    #if(!is.null(cov_beta) & !xnm_k %in% not_selected){
    if(!is.null(cov_beta)){
      if(shps[k] > 2 & shps[k] < 5 | shps[k] > 10 & shps[k] < 13){
        ps_k <- c(1 + capk + ncon, ps_k) #constant, nominal, conv/conc
        ncon <- ncon + 1
      }
      cov_k <- cov_beta[ps_k, ps_k]
      dd_k <- dd[, ps_k, drop = FALSE]
      sd_eta_k <- diag(dd_k %*% cov_k %*% t(dd_k))^.5
      lwr_eta_k <- eta_k - 2 * sd_eta_k 
      upp_eta_k <- eta_k + 2 * sd_eta_k 
      
      if(mu){
        lwr_eta_k <- family$linkinv(lwr_eta_k)
        upp_eta_k <- family$linkinv(upp_eta_k)
      }
      
      range_eta <- range(lwr_eta_k, upp_eta_k, na.rm = TRUE)
    } else {
      range_eta <- range(eta_k, na.rm = TRUE)
    }
    padding <- 0.05 * diff(range_eta)
    ylim_use <- c(range_eta[1] - padding, range_eta[2] + padding)
    #--------------------------
    #vis each component 
    
    ord_k <- order(x_k)
    tp_k <- ifelse(is.factor(x_k), 'p', 'l')
    xaxt_k <- ifelse(is.factor(x_k), 'n', 's')
    x_k_tmp <- x_k
    if(is.factor(x_k_tmp)){
      #x_k_tmp <- as.numeric(levels(x_k))[x_k] 
      x_k_tmp <- as.integer(x_k) #- 1
    }
    if(mn){
      plot(x_k_tmp[ord_k], eta_k[ord_k], xlab = bquote(~ .(xnm_k)), ylab = bquote(hat(eta)[.(k)]),
           ylim = ylim_use, 
           main = bquote(~ hat(eta)[.(k)] ~ "vs" ~ .(xnm_k) ~ "(" * .(shp_k) * ")"),
           #main = bquote(~ .(shp_k)),
           xaxt = xaxt_k,
           type = tp_k, lwd = 1)
    }else{
      plot(x_k_tmp[ord_k], eta_k[ord_k], xlab = bquote(~ .(xnm_k)), ylab = "",
           ylim = ylim_use, 
           #main = bquote(~ hat(eta)[.(k)] ~ "vs" ~ .(xnm_k) ~ "(" * .(shp_k) * ")"),
           main = NULL,
           xaxt = xaxt_k,
           type = tp_k, lwd = 1)
    }    
    
    if(!exists("rug", where = extras)){
      rug <- TRUE
    } else {
      rug <- extras$rug
    }
    
    if(rug){
      rug(x_k, ticksize = 0.03, side = 1)
    }
    
    #if(!is.null(cov_beta) & !xnm_k %in% not_selected){
    if(!is.null(cov_beta) & !is.factor(x_k)){
      lines(x_k[ord_k], lwr_eta_k[ord_k], col = 1, lty = 2)
      lines(x_k[ord_k], upp_eta_k[ord_k], col = 1, lty = 2)
    }
    #print (attr(x_k, 'type') == 'tree')
    if(is.factor(x_k) | types[k] == 'tree'){
      points(x_k_tmp, eta_k, pch=20)
      if(is.factor(x_k)){
        lines(x_k_tmp[ord_k], eta_k[ord_k], lty = 2)
        axis(1, at = unique(x_k_tmp), labels = unique(levels(x_k)))  # custom
      }
    }
    
    if(pages > 1){
      readline("Press <Enter> to continue...")
    }
  }
  if(pages == 1){
    on.exit(par(oldpar))  # restores original settings
  }
}

#---------------------------------------------------------------------------------------------------------------------------
predict.shapeselect <- function(object, newdata, interval = c("none", "confidence", "prediction"), 
                                type = c("response", "link"), level = 0.95, n.mix = 100,...)
{
  interval <- match.arg(interval)
  type <- match.arg(type)
  
  family <- object$family
  space <- object$space
  pnt <- object$pnt
  center <- object$center
  standardize <- object$standardize
  boundary.knots <- object$boundary.knots
  #cicfamily <- CicFamily(family)
  #muhat.fun <- cicfamily$muhat.fun
  if (!inherits(object, "classo")) {
    warning("calling predict.classo(<fake-classo-object>) ...")
  }
  if (missing(newdata) | is.null(newdata)) {
    etahat <- object$etahat
    muhat <- object$muhat
    #muhat <- muhat.fun(etahat, fml = family$family)
    #muhat <- family$linkinv(etahat)
    ans <- list(fit = muhat, etahat = etahat)
    return (ans)
  }
  if (!is.data.frame(newdata)) {
    stop ("newdata must be a data frame!")
  }
  if(interval != "none"){
    cmat <- vcov.cgam(object)
  } else {
    cmat <- NULL
  }
  #xmatpr = object$xmat
  #match newdata by the xnms
  knots_lst <- object$knots_lst
  mskeep_lst <- object$mskeep_lst
  cvec <- object$cvec
  
  sdx <- object$sdx
  mx <- object$mx
  coefskeep_withone <- object$coefskeep_withone #will lead to wrong prediction when we have conc/conv; bigmat is not ordered the same as cbind(one, ddkeepz_all, ddkeepx_all)
  coefs_x_lst <- object$coefskeep_x_lst
  coefs_z_lst <- object$coefskeep_z_lst
  
  #temp
  xnms <- object$xnms
  shps <- object$shps
  seq_added <- object$seq_added
  
  xnms_kept <- xnms[which(!shps %in% c(0))]
  #if(!all(colnames(newdata) %in% xnms_kept) | !all(xnms_kept %in% colnames(newdata))){
  if(!all(names(newdata) %in% xnms)){
    stop ("newdata must include at least all selected predictors!")
  } else {
    if(NCOL(newdata) == length(xnms_kept)){
      newdata <- newdata[, xnms_kept, drop = F]
    } else if(NCOL(newdata) == length(xnms)){
      newdata <- newdata[, xnms, drop = F]
    }
  }
  
  #######################
  #local helper function#
  #######################
  my_line <- function(xp = NULL, y, x, end, start) {
    slope <- NULL
    intercept <- NULL
    yp <- NULL
    slope <- (y[end] - y[start]) / (x[end] - x[start])
    intercept <- y[end] - slope * x[end]
    yp <- intercept + slope * xp
    ans <- new.env()
    ans$slope <- slope
    ans$intercept <- intercept
    ans$yp <- yp
    ans
  }
  
  if(any(shps > 0)){
    xmat <- object$xmat #not scaled
    types <- object$types
    #xmat0 <- xmat
    n <- nrow(xmat)
    p <- length(xnms)
    for(i in 1:p){
      if(is.numeric(xmat[,i]) & types[i] %in% c('num', 'numeric')){
        #if(center){
        xmat[,i] <- (xmat[,i] - mx[i]) / sdx[i]
        newdata[,i] <- (newdata[,i] - mx[i]) / sdx[i]
        #} else {
        #  xmat[,i] <- xmat[,i] / sdx[i]
        #  newdata[,i] <- newdata[,i] / sdx[i]
        #}
        
        #new: scale user-defined boundary knots
        # if(is.list(boundary.knots)){
        #   kts0 <- boundary.knots[[i]]
        # }
        # 
        # if(length(kts0) >= 2){
        #   #if(center){
        #     kts0 <- (kts0 - mx[i]) / sdx[i]
        #   #} else {
        #   #  kts0 <- kts0 / sdx[i]
        #   #}
        # }
        #xmat[,i] <- xmat[,i] / sdx[i]
        #xmat[,ix] <- scale(xmat[,ix], center=center)
      }
      #test 
      if(inherits(xmat[,i], "factor") & types[i] %in% "tree"){
        xmat[,i] <- as.numeric(as.character(xmat[,i]))
        newdata[,i] <- as.numeric(as.character(newdata[,i]))
      }
      if(inherits(xmat[,i], "factor") & types[i] %in% "ord"){
        xmat[,i] <- as.numeric(as.character(xmat[,i]))
        newdata[,i] <- as.numeric(as.character(newdata[,i]))
      }
    }
    #print (newdata)
    #create xmatpr for new data
    zmat <- NULL
    xcmat <- NULL #x's when shp is conv or conc
    delta <- NULL
    ddkeep_x_lst <- vector("list", length = p)
    ddkeep_z_lst <- vector("list", length = p)
    one <- rep(1, nrow(newdata))
    varlst <- varlst_z <- NULL
    #shps_added <- shps[seq_added]
    #xnms_added <- xnms[seq_added]
    for(ix in 1:p){
      #print (ix)
      type_ix <- types[ix]
      shpi <- shps[ix]
      if(shpi %in% c(0, 17)){
        next 
      }
      x <- newdata[, ix]
      ms <- mskeep_lst[[ix]]
      if(inherits(x, "factor")){
        #dd <- model.matrix(~ x)[, -1, drop=FALSE]
        #dd <- apply(dd, 2, function(e) e - mean(e))
        lvli <- levels(x)
        xui <- unique(x)
        ztbi <- levels(xmat[, ix])
        dd <- NULL
        zcoefid <- NULL
        if (!any(lvli %in% ztbi)) {
          stop ("new factor level must be among factor levels in the fit!")
        } else {
          #baseline: ztbi[1]
          #if(any(lvli != ztbi[1])){
          klvls <- length(ztbi)
          rn <- nrow(newdata)
          #remove baseline?
          dd <- matrix(0, nrow = rn, ncol = (klvls - 1))
          #new: test!
          #if(nrow(ms) == 1) {
          for (i1 in 1:rn) {
            if (x[i1] != ztbi[1]) {
              id_col <- which(ztbi %in% x[i1]) - 1
              #print (id_col)
              dd[i1, id_col] <- 1
              dd[i1, id_col] <- dd[i1, id_col] - ms[id_col] #ms[i1, id_col]
            }
          }
          #temp
          for(ic in 1:ncol(dd)){
            ps_base <- which(dd[, ic] == 0)
            dd[ps_base, ic] <- dd[ps_base, ic] - ms[ic] #ms[ps_base, ic]
          }
          #}
          
          # if(nrow(ms) > 1) {
          #   oldx <- xmat[, ix]
          #   ddold <- model.matrix(~ oldx)[, -1, drop=FALSE]
          #   #xs <- sort(oldx)
          #   #ord <- order(oldx)
          #   nr <- nrow(dd)
          #   #nc <- length(x)
          #   nc <- ncol(dd)
          #   ms2 <- matrix(0, nrow = nr, ncol = nc)
          #   for (i1 in 1:nc) {
          #     xs <- sort(ddold[, i1])
          #     ord <- order(ddold[, i1])
          #     for (i2 in 1:nr) {
          #       ms2[i2, i1] <- my_line(xp = x[i1], y = ms[i2, ][ord], x = xs, end = n, start = 1)$yp
          #     }
          #   }
          #   
          # }
          
          zmat <- cbind(zmat, dd) 
          
          #test
          ddkeep_z_lst[[ix]] <- dd
          varlst_z <- c(varlst_z, rep(ix, NCOL(dd)))
        }
        #}
      } else {
        oldx <- xmat[, ix]
        #if(any(x > max(oldx)) | any(x < min(oldx))){
        #  stop ("No extrapolation is allowed in classo prediction!")
        #}
        shpi <- shps[ix]
        ms <- mskeep_lst[[ix]]
        
        if(type_ix %in% c("num", "numeric")){
          kts <- knots_lst[[ix]]
          # if(type_ix %in% c("ord", "ordinal")){
          #   kts <- 0
          # }
          #print (kts)
          dd_ans <- suppressWarnings(
            makedelta(x, sh = shpi, numknots = 0, knots = kts, space = space, suppre = TRUE, interp = TRUE))
          dd <- dd_ans$amat
          #cat('finished delta')
          
          if(shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5){
            ms <- t(ms)
          }
          #need to make basis orthogonal to one or [1, x]
          if (shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5) {
            xs <- sort(oldx)
            ord <- order(oldx)
            nr <- nrow(dd)
            nc <- length(x)
            ms2 <- matrix(0, nrow = nr, ncol = nc)
            for (i1 in 1:nc) {
              for (i2 in 1:nr) {
                ms2[i2, i1] <- my_line(xp = x[i1], y = ms[i2, ][ord], x = xs, end = n, start = 1)$yp
              }
            }
            dd <- dd - ms2
          } else {
            dd <- dd - ms
          }
          delta <- cbind(delta, t(dd))
        }
        
        if(type_ix %in% c("ord", "ordinal")){
          dd <- pred_del(oldx, shpi, x, ms)  
        }
        
        if(type_ix %in% c("tree")){
          #easier than predict.cgam
          uoldx <- unique(oldx)
          #pl <- uoldx[1]
          if (is.numeric(oldx)) {
            if (0 %in% oldx) {
              pl <- 0
            } else {pl <- min(oldx)}
          } else {
            pl <- uoldx[1]
          }
          
          nx <- uoldx[uoldx != pl]
          is_pl <- x %in% pl
          #use the order of unique oldx to determine the treatment level of new: x
          dd <- matrix(0, nrow = length(nx), ncol = length(x))
          #if new x only has placebo, then just center each row of dd by ms
          for(i in 1:nrow(dd)){
            xi <- nx[i]
            dd[i, ] <- dd[i, ] - ms[which(nx %in% xi)] 
          }
          
          if(any(!is_pl)){
            for(i in 1:nrow(dd)){
              xi <- nx[i]
              if(any(x %in% xi)){
                dd[i, which(x %in% xi)] <- 1 - ms[which(nx %in% xi)] 
              }
            }
          }
        }
        
        if(shpi == 2 & type_ix %in% c("tree")) {
          ddkeep_z_lst[[ix]] <- t(dd)
          varlst_z <- c(varlst_z, rep(ix, NROW(dd)))
        } else {
          ddkeep_x_lst[[ix]] <- t(dd)
          varlst <- c(varlst, rep(ix, NROW(dd)))
        }
        
        if(shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5){
          #zmat <- cbind(zmat, x - mean(oldx))
          xcmat <- cbind(xcmat, x - mean(oldx))
          ddkeep_z_lst[[ix]] <- x - mean(oldx)
          varlst_z <- c(varlst_z, rep(ix, 1))
        }
      }
    }
    
    #temp
    # if(shrinketa){
    #   shps_added <- shps[seq_added]
    #   xnms_added <- xnms[seq_added]
    #   for(i in 1:length(seq_added)){
    #     ps <- seq_added[i]
    #     xnmi <- xnms_added[i]
    #     shpi <- shps_added[i]
    #     if(shpi == 0){
    #       next
    #     }
    #     if(shpi == 1 | shpi == 2 & type_ix %in% c("tree")) {
    #       coefs_z_lst[[ps]] <- coefs_z_lst[[ps]] * cvec[ps]
    #     } else {
    #       coefs_x_lst[[ps]] <- coefs_x_lst[[ps]] * cvec[ps]
    #       if(shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5){
    #         coefs_z_lst[[ps]] <- coefs_z_lst[[ps]] * cvec[ps]
    #       }
    #     }
    #   }
    # }
    #xmatpr <- cbind(one, zmat, xcmat, delta)
    #print (dim(xmatpr))
    #get new fit
    #fit <- xmatpr %*% coefskeep_withone[1:ncol(xmatpr)]
    #coefsuse <- c(onecoef, zcoefs_nox, xccoefs, xcoefs)
    
    #ddkeep_x_lst_ord <- ddkeep_x_lst[seq_added]
    #ddkeep_z_lst_ord <- ddkeep_z_lst[seq_added]
    delta <- do.call(base::cbind, ddkeep_x_lst)
    zmat <- do.call(base::cbind, ddkeep_z_lst) #include conv/conc
    xmatpr <- cbind(one, zmat, delta)
    
    #coefs_x_lst_ord <- coefs_x_lst[seq_added]
    #coefs_z_lst_ord <- coefs_z_lst[seq_added]
    coefsuse_z <- do.call(base::c, coefs_z_lst)
    coefsuse_x <- do.call(base::c, coefs_x_lst)
    
    coefsuse <- c(coefskeep_withone[1], coefsuse_z, coefsuse_x)
    fit <- xmatpr %*% coefsuse
    #fit <- switch(type, link = fit, response = object$family$linkinv(fit))
    
    etacomps <- matrix(0, nrow = p, ncol = length(fit))
    etacomps_lwr <- etacomps_upp <- etacomps
    #shps_added <- shps[seq_added]
    #xnms_added <- xnms[seq_added]
    if(!is.null(cmat)){
      cmat_noone <- cmat[2:nrow(cmat),2:nrow(cmat)]
      if(NCOL(zmat) > 0){
        ps_z <- 1:NCOL(zmat)
        cmat_z <- cmat_noone[ps_z, ps_z] #include conv/conc
        cmat_x <- cmat_noone[-ps_z, -ps_z]
      } else {
        cmat_x <- cmat_noone
      }
    } else {
      cmat_z <- cmat_x <- NULL
    }
    
    for(k in 1:p){
      #ps <- seq_added[i]
      # xnmi <- xnms_added[i]
      # shpi <- shps_added[i]
      shpi <- shps[k]
      if(shpi == 0){#flat
        next
      } else {
        if (shpi == 1 | shpi == 2) {#need to use 1 for incr later
          #etacomps[ps, ] <- ddkeep_z_lst[[ps]] %*% coefs_z_lst[[ps]]
          if(types[k] %in% c('nom', 'nominal') | shpi == 2 & types[k] %in% 'tree'){
            etacomps[k, ] <- ddkeep_z_lst[[k]] %*% coefs_z_lst[[k]]
            if(!is.null(cmat_z)){
              ps_k <- which(varlst_z == k)
              dd_k <- zmat[, ps_k, drop = FALSE]
              cov_k <- cmat_z[ps_k, ps_k]
            }
          }
          if(types[k] %in% c('ord', 'ordinal') | shpi == 1 & types[k] %in% 'tree'){
            etacomps[k, ] <- ddkeep_x_lst[[k]] %*% coefs_x_lst[[k]]
            if(!is.null(cmat_x)){
              ps_k <- which(varlst == k)
              dd_k <- delta[, ps_k, drop = FALSE]
              cov_k <- cmat_x[ps_k, ps_k]
            }
          }
        } else if (shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5) {
          etacomps[k, ] <- ddkeep_x_lst[[k]] %*% coefs_x_lst[[k]]
          etacomps[k, ] <- etacomps[k, ] + ddkeep_z_lst[[k]] * coefs_z_lst[[k]]
          if(!is.null(cmat)){
            ps_k <- c(which(varlst_z == k), which(varlst == k))
            dd_k <- cbind(zmat[, which(varlst_z == k), drop = FALSE], delta[, which(varlst == k), drop = FALSE])
            cov_k <- cmat_noone[ps_k, ps_k]
          }
        } else {
          etacomps[k, ] <- ddkeep_x_lst[[k]] %*% coefs_x_lst[[k]]
          if(!is.null(cmat_x)){
            ps_k <- which(varlst == k)
            dd_k <- delta[, ps_k, drop = FALSE]
            cov_k <- cmat_x[ps_k, ps_k]
          }
        }
        
        if(!is.null(cmat)){
          eta_k <- etacomps[k, ]
          #ps_k <- which(varlst == k)
          #dd_k <- xmat[obs, ps_k, drop = FALSE]
          #cov_k <- cov_beta[ps_k, ps_k]
          sd_eta_k <- diag(dd_k %*% cov_k %*% t(dd_k))^.5
          mult <- qnorm((1 - level)/2, lower.tail=FALSE)
          lwr_eta_k <- eta_k - mult * sd_eta_k
          upp_eta_k <- eta_k + mult * sd_eta_k
          etacomps_lwr[k, ] <- lwr_eta_k
          etacomps_upp[k, ] <- upp_eta_k
        }
      }
    }
    
    if(interval != "none"){
      coveta <- diag(xmatpr %*% cmat %*% t(xmatpr))
      se.fit <- sqrt(coveta)
      # new: C.I. level
      mult <- qnorm((1 - level)/2, lower.tail=FALSE)
      upp <- fit + mult*se.fit 
      lwr <- fit - mult*se.fit 
      fit <- switch(type, link = fit, response = object$family$linkinv(fit))
      lwr <- switch(type, link = lwr, response = object$family$linkinv(lwr))
      upp <- switch(type, link = upp, response = object$family$linkinv(upp))
      #se.fit <- switch(type, link = se.fit, response = object$family$linkinv(se.fit))
      ans <- list(fit = fit, lwr = lwr, upp = upp, se.fit = se.fit, xmatpr = xmatpr, newddx = delta, newddz = zmat,
                  etacomps = etacomps, etacomps_lwr = etacomps_lwr, etacomps_upp = etacomps_upp)
    } else {
      #print (type)
      fit <- switch(type, link = fit, response = object$family$linkinv(fit))
      ans <- list(fit = fit, xmatpr = xmatpr, newddx = delta, newddz = zmat, etacomps = etacomps)
    }
  }
  
  
  if(all(shps == 0)){
    y <- object$y
    m0 <- glm(y ~ 1, family = family$family)
    fit <- predict(m0, newdata)
    xmatpr <- predict(m0, newdata, type = 'terms')
    ans <- list(fit = fit, xmatpr = xmatpr)
  } 
  
  return (ans)
}

#############################
#predict delta for ordinal x#
#############################
# pred_del = function(x, sh, xp, ms) {
#   n = length(xp)
#   #x = (x - min(x)) / (max(x) - min(x))
#   xu = sort(unique(x))
#   n1 = length(xu)
#   sigma = NULL
#   #######################
#   #local helper function#
#   #######################
#   my_line = function(xp = NULL, y, x, end, start) {
#     slope = NULL
#     intercept = NULL
#     yp = NULL
#     slope = (y[end] - y[start]) / (x[end] - x[start])
#     intercept = y[end] - slope * x[end]
#     yp = intercept + slope * xp
#     ans = new.env()
#     ans$slope = slope
#     ans$intercept = intercept
#     ans$yp = yp
#     ans
#   }
#   #  increasing or decreasing
#   if (sh < 3) {
#     sigma = matrix(0, nrow = n1 - 1, ncol = n)
#     for (i in 1: (n1 - 1)) {
#       sigma[i, xp > xu[i]] = 1
#     }
#     if (sh == 2) {sigma = -sigma}
#     for (i in 1:(n1 - 1)) {sigma[i, ] = sigma[i, ] - ms[i]}
#   }
#   if (sh == 3 | sh == 4) {
#     #  convex or concave
#     sigma = matrix(0, nrow = n1 - 2, ncol = n)
#     #for (i in 1: (n1 - 2)) {
#     #	sigma[i, x > xu[i]] = x[x > xu[i]] - xu[i]
#     #}
#     for (i in 1: (n1 - 2)) {
#       sigma[i, xp > xu[i+1]] = xp[xp > xu[i+1]] - xu[i+1]
#     }
#     if (sh == 4) {sigma = -sigma}
#     #xm = cbind(1:n*0+1, xp)
#     #xpx = solve(t(xm) %*% xm)
#     #pm = xm %*% xpx %*% t(xm)
#     #sigma = sigma - sigma %*% t(pm)
#     xs = sort(x)
#     ord = order(x)
#     nx = length(x)
#     obs = 1:nx
#     m = nrow(ms)
#     ms0 = matrix(0, nrow = m, ncol = n)
#     for (i1 in 1:n) {
#       for (i2 in 1:m) {
#         ms0[i2, i1] = my_line(xp = xp[i1], y = ms[i2, ][ord], x = xs, end = nx, start = 1)$yp
#       }
#     }
#     sigma = sigma - ms0
#   }
#   if (sh > 4 & sh < 9) {
#     sigma = matrix(0, nrow = n1 - 1, ncol = n)
#     if (sh == 5) { ### increasing convex
#       for (i in 1:(n1 - 1)) {
#         sigma[i, xp > xu[i]] = (xp[xp > xu[i]] - xu[i]) / (max(x) - xu[i])
#       }
#       for (i in 1:(n1 - 1)) {sigma[i,] = sigma[i,] - ms[i]}
#     } else if (sh == 6) {  ## decreasing convex
#       for (i in 1:(n1 - 1)) {
#         sigma[i, xp < xu[i + 1]] = (xp[xp < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
#       }
#       for (i in 1:(n1 - 1)) {sigma[i,] = sigma[i,] - ms[i]}
#       #print (ms)
#     } else if (sh == 7) { ## increasing concave
#       for (i in 1:(n1 - 1)) {
#         sigma[i, xp < xu[i + 1]] = (xp[xp < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
#       }
#       for (i in 1:(n1 - 1)) {sigma[i,] = -sigma[i,] + ms[i]}
#     } else if (sh == 8) {## decreasing concave
#       for (i in 1:(n1 - 1)) {
#         sigma[i, xp > xu[i]] = (xp[xp > xu[i]] - xu[i]) / (max(x) - xu[i])
#       }
#       for (i in 1:(n1 - 1)) {sigma[i,] = -sigma[i,] + ms[i]}
#     }
#   }
#   return (sigma)
# }

#--------------------------------------------------------------------------------
# zvec.fun <- function(cvec = NULL, wt = NULL, y, sm = 1e-7, fml = object$family) {
#   n <- length(y)
#   if (fml == "gaussian") {
#     #zvec <- y
#     zvec <- wt^(1/2) * y
#   }
#   if (fml == "poisson") {
#     #zvec <- cvec / wt
#     zvec <- cvec / sqrt(wt)
#   }
#   if (fml == "binomial") {
#     zvec = 1:n*0
#     zvec[wt == 0] <- 1 / sm
#     zvec[wt > 0] <- cvec[wt > 0] / sqrt(wt[wt > 0])
#   }
#   if (fml == "Gamma") {
#     zvec <- cvec / sqrt(wt)
#   }
#   zvec
# }

# wt.fun <- function(y, etahat = NULL, n = NULL, weights = NULL, fml = object$family){
#   if (is.null(weights)) {
#     weights <- 1:n*0 + 1
#   }
#   w <- weights
#   if (fml == "poisson") {
#     wt <-  w * exp(etahat)
#   }
#   if (fml == "binomial") {
#     if (all(etahat == 0)){
#       #wt <- 1:n*0 + 1/4
#       wt <- w * (1:n*0 + 1/4)
#     } else {
#       wt <- 1:n*0
#       wt <- (1 - 1 / (1 + exp(etahat))) * (1 / (1 + exp(etahat))) #pi_ij(1 - pi_ij)
#       wt <- c(wt)
#       # for (i in 1:n) {
#       #   if (etahat[i] > 100) {
#       #     wt[i] <- 0
#       #   } else {
#       #     wt[i] <- w[i] * exp(etahat[i]) / ((1 + exp(etahat[i]))^2)
#       #   }
#       # }
#     }
#   }
#   if (fml == "gaussian") {
#     wt <- w # (1:n*0 + 1) / w
#   }
#   if (fml == "Gamma") {
#     wt <- w * y * exp(-etahat)
#   }
#   wt <- as.vector(wt)
#   wt
# }
# 
# etahat.fun <- function(n, y, fml = object$family){
#   if (fml == "poisson") {
#     etahat <- 1:n*0 + log(mean(y))
#   }
#   if (fml == "binomial") {
#     etahat <- 1:n*0
#   }
#   if (fml == "Gamma") {
#     etahat <- 1:n*0 + log(mean(y))
#   }
#   etahat
# }

# gr.fun <- function(y, etahat = NULL, weights = NULL, fml = object$family){
#   n <- length(y)
#   if (is.null(weights)) {
#     weights <- 1:n*0 + 1
#   }
#   w <- weights
#   if (fml == "poisson") {
#     gr <- w * (exp(etahat) -  y)
#   }
#   if (fml == "binomial") {
#     if (all(etahat == 0)) {
#       gr <- w * (1/2 - y)
#     } else {
#       gr <- 1:n*0
#       #ignore weights
#       gr <- 1 - 1 / (1 + exp(etahat)) - y
#       
#       # for (i in 1:n) {
#       #   if (etahat[i] > 100) {
#       #     gr[i] <- w[i] * (1 - y[i])
#       #   } else {gr[i] <- w[i] * (exp(etahat[i]) / (1 + exp(etahat[i])) - y[i])}
#       # }
#     }
#   }
#   if (fml == "Gamma") {
#     gr <- w * (1 - y * exp(-etahat))
#   }
#   gr
# }

# ysim.fun <- function(n, mu0 = NULL, fml = object$family, shp0 = NULL) {
#   if (fml == "binomial") {
#     ysim <- 1:n*0
#     ysim[runif(n) < .5] <- 1
#   }
#   if (fml == "poisson") {
#     if (!is.null(mu0)) {
#       ysim <- rpois(n, mu0)
#     }
#   }
#   if (fml == "gaussian") {
#     ysim <- rnorm(n)
#   }
#   if (fml == "Gamma") {
#     ysim <- rgamma(n, shape=1)
#   }
#   ysim
# }

# llh.fun <- function(y, muhat = NULL, etahat = NULL, phihat = NULL, n = NULL, weights = NULL, fml = object$family){
#   #sm <- 1e-7
#   sm <- 1e-5
#   if (is.null(weights)) {
#     weights <- 1:n*0 + 1
#   }
#   w <- weights
#   #new: avoid Inf
#   if (fml == "poisson") {
#     llh <- 2 * sum(w[w!=0] * (muhat[w!=0] - y[w!=0] * etahat[w!=0])) / n
#   }
#   if (fml == "binomial") {
#     llh <- 0
#     if (all(0 <= y) & all(y <= 1)) {
#       for (i in 1:n) {
#         if (muhat[i] > 0 & muhat[i] < 1) {
#           llh <- llh + w[i] * (y[i] * log(muhat[i]) + (1 - y[i]) * log(1 - muhat[i]))
#         }
#       }
#       #ch <- muhat > 0 & muhat < 1
#       #llh <- sum(w[ch] * (log(1 - muhat[ch]) - y[ch] * etahat[ch]))
#       llh <- (-2/n) * llh
#     } else {
#       stop ("y values must be 0 <= y <= 1!")
#     }
#   }
#   if (fml == "gaussian") {
#     if (all(w == 1)) {
#       llh <- log(sum((y - etahat)^2))
#     } else {
#       llh <- log(sum(w[w!=0] * (y[w!=0] - etahat[w!=0])^2)) - sum(log(w[w!=0])) / n
#     }
#   }
#   if (fml == "Gamma") {
#     vuhat <- 1 / phihat
#     #print (vuhat)
#     #llh <- 2 * sum(w[w!=0] * (etahat[w!=0] + y[w!=0] * exp(-etahat[w!=0]))) / n
#     llh <- 2 / n * (vuhat * sum(w[w!=0] * (etahat[w!=0] + y[w!=0] * exp(-etahat[w!=0]))) + n * (log(gamma(vuhat)) - vuhat * log(vuhat)) - (vuhat-1) * sum(log(y)))
#     #print (llh)
#   }
#   llh
# }

.search_ps_shapeselect = function(pen, qv0, dv, dd, edfu, covmat_inv=NULL) {
  # 1. Calculate qv
  qv = qv0 + pen * dv
  # 2. Stronger Regularization
  # If 1e-10 failed, your matrix is extremely sensitive. 
  # Let's use a slightly larger nudge or a relative nudge.
  diag(qv) <- diag(qv) + (max(diag(qv)) * 1e-12) + 1e-8
  #diag(qv) <- diag(qv) + 1e-10
  # 3. Use tryCatch for the whole block to catch "Singular" errors
  val <- tryCatch({
    umat = chol(qv)
    # Efficiently calculate the trace: 
    # instead of solve(qv, t(dd)), use backsolve/forsolve with the Cholesky factor
    # This is numerically much more stable
    # Let A = dd %*% qv_inv %*% t(dd). We want trace(A).
    tmp <- backsolve(umat, t(dd), transpose = TRUE)
    trace_val <- sum(tmp^2) 
    
    return(trace_val - .8 * edfu)
    
  }, error = function(e) {
    # If the math breaks, return a value that tells uniroot to keep looking
    # Usually a large positive value if we are looking for a root from above
    return(1e10) 
  })
  
  return(val)
}

#--------------------------------------------------------------------
#cross validation to find the c vector
#remember to acknowledge
#--------------------------------------------------------------------
create_folds <- function(n_data, nfolds = 5){
  folds <- rep(list(numeric(0)), nfolds)
  available_samples <- 1:n_data
  fold_size <- floor(n_data/nfolds)
  if(fold_size == n_data/nfolds){
    for(fold_id in 1:nfolds){
      folds[[fold_id]] <- sample(available_samples, fold_size)
      available_samples <- available_samples[!(available_samples %in% folds[[fold_id]])]
    }
  } else{
    for(fold_id in 1:(nfolds-1)){
      folds[[fold_id]] <- sample(available_samples, min(fold_size + rbinom(1, 1, 0.5), length(available_samples)))
      available_samples <- available_samples[!(available_samples %in% folds[[fold_id]])]
    }
    folds[[nfolds]] <- available_samples
  }
  return(folds)
}

#------------------------------------------------------------------------------------------
cv.classo = function(fmat, y, nfolds = 5, type = c("glmnet", "solve.qp"), family=gaussian()) {
  #sdxs = apply(fmat, 2, sd) #used later to go back to the original scale
  #mxs = apply(fmat, 2, mean)
  #fmat = scale(fmat, center = FALSE, scale = TRUE) #scale new x
  type = match.arg(type)
  
  #temp:
  if(NCOL(fmat) == 1) {
    #type = "solve.qp"
    cvec = 1
    rslt = list(cvec = cvec)
    #return (rslt)
  } else {
    if(type == "solve.qp") {
      n = length(y)
      #vec = 1:n
      #shuffle the observations
      #vec_sample = sample(vec, size = n, replace = FALSE) 
      #split this data into k folds with roughly equal sample sizes
      #folds = split(vec_sample, 1:nfolds) 
      folds = create_folds(n, nfolds)
      
      p = NCOL(fmat)
      amat = diag(p)
      #new: each c value <= 1
      amat = rbind(amat, -diag(p))
      amat = rbind(amat, c(rep(-1, p)))
      
      #amat = cbind(0, amat)
      #what's the steps for s?
      s.max = 1*p
      #s.max = 1.5*p
      s.min = 0.01
      s.step = (s.max - s.min)/49
      s.grid = seq(from=s.min, to=s.max, by=s.step)
      ns = length(s.grid)
      
      #record all 5-fold cv error for s on the s.grid
      cv_errors_all_s = vector("numeric", length = ns) 
      cv_errors_mat = matrix(NA, nrow = ns, ncol = nfolds)
      ck_mat = matrix(NA, nrow = ns, ncol = (p))
      
      #y is already centered
      df = cbind(y, fmat)
      for(is in 1:ns){
        #print (is)
        si = s.grid[is]
        test_errors_s = vector("numeric", nfolds) #initiate an error vector
        for(i in 1:nfolds){
          test_set_id = folds[[i]] #take out the ith fold as the test set ID
          test_set = df[test_set_id, ,drop = FALSE] #specify the test set
          train_set = df[-test_set_id, ,drop = FALSE]  #specify the training set
          xs_train = train_set[, -1, drop = FALSE] #put y at the 1st column
          xs_test = test_set[, -1, drop = FALSE] #put y at the 1st column
          y_train = train_set[, 1, drop = FALSE]
          y_test = test_set[, 1, drop = FALSE]
          #apply the trained model to the test set and get yhat_test
          #trained model: for si, get a set of ck's, fk's are already estimated in step 1
          qv = crossprod(xs_train)
          cv = crossprod(xs_train, y_train)
          #new: each c value >= 0 and <=1
          bvec = c(rep(0, p), rep(-1, p), -si)
          ansi = solve.QP(qv, cv, t(amat), bvec)
          ck_vec = ansi$solution
          
          ck_mat[is, ] = ck_vec #no intercept
          yhat_test = xs_test %*% ck_vec
          test_errors_s[i] = mean((y_test - yhat_test)^2) #compute the test error for the ith iteration
        }
        cv_errors_all_s[is] = mean(test_errors_s) # sum(k-fold cv)/k
        cv_errors_mat[is, ] = test_errors_s
      }
      
      #print (cv_errors_all_s)
      #print (which.min(cv_errors_all_s))
      #if there is a tie, then choose the smallest s
      best.s_ps = which(cv_errors_all_s == min(cv_errors_all_s))[1]
      best.s = s.grid[best.s_ps]
      cvec_kfold = ck_mat[best.s_ps, ,drop = TRUE] 
      
      #cvec_kfold = (cvec_kfold / sdxs) |> round(5)
      
      #final model: use all data points
      qv = crossprod(fmat)
      cv = crossprod(fmat, y)
      #best.s is not cv.glmnet$lambda.min
      #new: each c value >= 0 and <=1
      bvec = c(rep(0, p), rep(-1, p), -best.s)
      ans.qp = solve.QP(qv, cv, t(amat), bvec)
      cvec = ans.qp$solution
      cvec = cvec |> round(7)
      #back to the original scale
      #cvec = (cvec / sdxs) |> round(5)
      #for(ic in seq(NCOL(ck_mat))) {ck_mat[,ic] = ck_mat[,ic]/sdxs[ic]}
      rslt = list(ck_mat = ck_mat, cvec = cvec, best.s_ps = best.s_ps, best.s = best.s)
      #return (rslt)
    }
    
    #y is not centered.
    if(type == "glmnet") {
      cv_lasso = glmnet::cv.glmnet(fmat, y, alpha=1, intercept=TRUE, lower.limits=0, 
                                   standardize=FALSE, nfolds=nfolds, family=family, type.measure = "deviance")
      
      #cv_lasso = glmnet::cv.glmnet(fmat, y, alpha=1, intercept=TRUE, 
      #                             standardize=FALSE, nfolds=nfolds, family=family, type.measure = "deviance")
      
      cvec = coef(cv_lasso, s = "lambda.min")[-1]
      rslt = list(cvec = cvec)#, nng.mod = nng.mod, cv_lasso = cv_lasso)
      #return (rslt)
    }
  }
  return (rslt)
}

##############################
#tranform shapes back to char
##############################
ShapeToChar_shapeselect = function(shp, tag = "x") {
  #if (max(shp) > 16 | min(shp) < 0) {
  #	stop ('No such a shape! A shape value must be between 0 and 16.')
  #}
  if (tag == "x") {
    if (shp == 0) {
      #shp = 18
      shp = 19 #this should include z too...need to add tag later
    }
    switch(shp,
           cs1 = {ch = 'incr'},
           #need to change later:
           #cs1 = {ch = 'in the model'},
           cs2 = {ch = 'decreasing'},
           cs3 = {ch = 'convex'},
           cs4 = {ch = 'concave'},
           cs5 = {ch = 'increasing convex'},
           cs6 = {ch = 'decreasing convex'},
           cs7 = {ch = 'increasing concave'},
           cs8 = {ch = 'decreasing concave'},
           cs9 = {ch = 'smooth increasing'},
           cs10 = {ch = 'smooth deceasing'},
           cs11 = {ch = 'smooth convex'},
           cs12 = {ch = 'smooth concave'},
           cs13 = {ch = 'smooth increasing convex'},
           cs14 = {ch = 'smooth increasing concave'},
           cs15 = {ch = 'smooth decreasing convex'},
           cs16 = {ch = 'smooth decreasing concave'},
           #cs17 = {ch = 's'},
           cs17 = {ch = 'flat'},
           cs18 = {ch = 'linear'},
           cs19 = {ch = 'out of the model'} 
           #{print ('No such a shape')}
    )
  } else if (tag == "z") {
    if (shp == 0) {
      shp = 2
    }
    switch(shp,
           cs1 = {ch = 'in the model'},
           cs2 = {ch = 'out of the model'}
    )
  } else if (tag == "tree") {
    if (shp == 0) {
      shp = 3
    }
    switch(shp,
           cs1 = {ch = 'tree'},
           cs2 = {ch = 'unordered'},
           cs3 = {ch = 'out of the model'}
    )
  }
  return (ch)
}

#--------------------------------------------------------------------
#print.classo: show the final choices
#--------------------------------------------------------------------
print.shapeselect = function(x,...)
{
  print(x$family)
  cat("Call:\n")
  #print(object$mod_final)
  if(inherits(x, "shapeselectAll")){
    cat("Top models:\n")
    print(head(x$mod_tried))
  } else {
    print(x$mod_tried)
  }
  cat("\n")
  invisible(x)
}

#--------------------------------
#all subsets
#--------------------------------
shapeselectAll <- function(xmat, y, family = gaussian, nsim = 100,
                           eps = 1e-8, space = 'E', infocrit = c("CBIC", "CAIC"), pnt = TRUE, parallel = TRUE,
                           message = FALSE, a = .5, add_edfu = 0, decrease_edfu = 0, types = NULL, center = FALSE)
{
  cicfamily <- CicFamily(family)
  llh.fun <- cicfamily$llh.fun
  #new: use log link in gamma
  #linkfun <- cicfamily$linkfun
  etahat.fun <- cicfamily$etahat.fun
  gr.fun <- cicfamily$gr.fun
  wt.fun <- cicfamily$wt.fun
  zvec.fun <- cicfamily$zvec.fun
  #muhat.fun <- cicfamily$muhat.fun
  #ysim.fun <- cicfamily$ysim.fun
  #deriv.fun <- cicfamily$deriv.fun
  #dev.fun <- cicfamily$dev.fun
  
  infocrit <- match.arg(infocrit)
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family))
    stop("'family' not recognized!")
  fml <- family$family
  if(fml == "gaussian"){
    wt.iter <- FALSE
  } else {wt.iter <- TRUE}
  xnms <- colnames(xmat)
  if(is.null(xnms)){xnms <- paste('x', 1:NCOL(xmat), sep='')}
  #new: use types to check if x is for spline or ordinal basis or just nominal
  p <- NCOL(xmat)
  for(ix in 1:p){
    attr(xmat[, ix], "type") <- types[ix]
  }
  n <- length(y)
  #--------------------------------
  #step 0: project y onto 1
  #--------------------------------
  one <- rep(1, n)
  etahat <- muhat <- mean(y)
  rss <- sum((y - etahat)^2)
  
  #scale numeric x first:
  p <- NCOL(xmat)
  sdx <- mx <- rep(0, p)
  for(ix in 1:p){
    if(is.numeric(xmat[,ix]) & types[ix] %in% c("num", "numeric")){
      #xmat[,ix] <- scale(xmat[,ix], center=F)
      #sdx[ix] <- attr(xmat[,ix], "scaled:scale")
      
      xmat[,ix] <- scale(xmat[,ix], center=center)
      #sdx <- c(sdx, attr(xmat[,ix], "scaled:scale"))
      sdx[ix] <- attr(xmat[,ix], "scaled:scale")
      if(center){
        mx[ix] <- attr(xmat[, ix], 'scaled:center')
      }
    }
    if(inherits(xmat[,ix], "factor") & types[ix] %in% "tree"){
      xmat[,ix] <- as.numeric(as.character(xmat[,ix]))
    }
    attr(xmat[, ix], "type") <- types[ix]
  }
  
  #need a larger value to start cic/bic?
  if(infocrit %in% c("CAIC", "CIC")) {
    if(!wt.iter){
      cic <- log(rss) +  2*1/n
    } else {
      llh <- llh.fun(y=y, muhat=rep(1/2, n), etahat=rep(0, n), n=n, fml=fml)
      #test!
      #llh <- 2 * (llh * n / 2) #2 times neg-loglike
      cic <- llh + 2*1/n
    }
  }
  
  if(infocrit == "CBIC") {
    if(!wt.iter){
      #cic <- n*log(rss) + log(n) * 1
      cic <- log(rss) + a * log(n) * 1 / n^(6/7)
      #cic <- log(rss) + log(n) * 1 / n^(6/7)
    } else {
      #m0 <- glm(y ~ 1, family=family)
      #cic <- n*log(m0$deviance) + log(n) * 1 
      #test
      #cic <- cic / n
      llh <- llh.fun(y=y, muhat=rep(1/2, n), etahat=rep(0, n), n=n, fml=fml)
      #cic <- llh + log(n) * 1 / n
      #cic <- llh + a * log(n) * 1 / n^(6/7)
      cic <- llh + a * log(n) * (1+add_edfu) / n^(6/7)
      #print (cic)
    }
  }
  #----------------------------------------------------------------------
  #step 1: list all shape combs
  #----------------------------------------------------------------------
  totx <- NCOL(xmat)
  tot_combs <- 1
  shps_lst <- vector("list", length = totx)
  for(ic in 1:totx){
    if(inherits(xmat[[ic]], 'factor')){
      acc <- 2 #include flat
      #shps_lst[[ic]] <- c(0, 1)
      #new
      #if(attr(xmat[[ic]], 'type') %in% c('nom', 'nominal')){
      shps_lst[[ic]] <- c(0, 1)
      #}
    } else {
      #acc <- 9 #include flat: 17
      if(attr(xmat[[ic]], 'type') %in% c('ord', 'ordinal')){
        shps_lst[[ic]] <- c(1:8, 0)
      }
      if(attr(xmat[[ic]], 'type') %in% c('tree')){
        shps_lst[[ic]] <- c(1:2, 0)
      }
      if(attr(xmat[[ic]], 'type') %in% c('num', 'numeric')){
        shps_lst[[ic]] <- 9:17
      }
      acc <- length(shps_lst[[ic]])
      #new: if we use cbic and wt.iter then we only have 4 options
      if(wt.iter & infocrit == "CBIC"){
        if(attr(xmat[[ic]], 'type') %in% c('ord', 'ordinal')){
          shps_lst[[ic]] <- c(1:4, 0)
        }
        if(attr(xmat[[ic]], 'type') %in% c('num', 'numeric')){
          shps_lst[[ic]] <- c(9:12, 17)
        }
        if(attr(xmat[[ic]], 'type') %in% c('tree')){
          shps_lst[[ic]] <- c(1:2, 0)
        }
        acc <- length(shps_lst[[ic]])
      }
    }
    tot_combs <- tot_combs * acc
  }
  
  msg <- paste("We will have a total of", tot_combs, "models.", "Continue?")
  if(message){
    if(tot_combs > 160){
      #library(svDialogs)
      res <- dlgMessage(msg, "yesno")$res
      if (res == "yes") {
        cat("Computing CBIC for each model....", "\n")
      } else {
        stop("You can call shapeselect instead.", "\n")
      }
    }
  }
  
  mod_tried <- matrix(0, nrow = (tot_combs+1), ncol = (totx + 5))
  colnames(mod_tried) <- c(xnms, infocrit, "edf0", "-2/n*loglike", "penalty", "penalty for spline")
  #the 1st row is for etahat = ybar
  mod_tried[1, totx+1] <- cic
  mod_tried[1, totx+2] <- 1
  
  #----------------------------------------------------------------------
  #step 2: make delta, compute edf and cic
  #----------------------------------------------------------------------
  shp_combs <- expand.grid(shps_lst)
  colnames(shp_combs) <- xnms
  zmat_all <- varlst_all <- varlst_z_all <- delta_all <- coefs_lst <- vector("list", length = tot_combs)
  for(icomb in 1:tot_combs){
    #print (icomb)
    shp_comb_i <- as.numeric(shp_combs[icomb, ])
    if(all(shp_comb_i %in% c(0, 17))){
      mod_tried[icomb+1, ] <- mod_tried[1, ,drop = TRUE] 
    } else {
      non_flat_ps <- which(!shp_comb_i %in% c(0, 17))
      #new: for ordinal
      #if(all(shp_comb_i[non_flat_ps] < 9)){
      #  pnt = FALSE
      #}
      gmat <- NULL
      np_ix <- 1
      zmat <- NULL
      delta <- NULL
      dmat <- NULL
      varlst <- NULL
      varlst_z <- NULL
      dm_lst <- vector("list", length = length(non_flat_ps))
      for(ix in non_flat_ps){
        x <- xmat[, ix]
        if(inherits(x, "factor") & attr(x, 'type') != 'tree'){
          #new:
          #if(attr(x, 'type') %in% c('nom', 'nominal')){
          dd <- model.matrix(~ x)[, -1, drop=FALSE]
          ms <- apply(dd, 2, mean)
          dd <- apply(dd, 2, function(e) e - mean(e))
          zmat <- cbind(zmat, dd) #no intercept
          varlst_z <- c(varlst_z, rep(ix, ncol(dd)))
          #}
          
        } else {
          shpi <- shp_comb_i[ix]
          #new: include ordinal
          if(shpi <= 8){
            if(attr(x, 'type') %in% c('ord', 'ordinal')){
              dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = 0, space = space)
              ms <- dd_ans$ms
              dd <- t(dd_ans$amat)
              ps <- 0
              dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
            }
            
            #new
            if(attr(x, 'type') %in% c('tree')){ 
              #ms <- NULL
              dd <- t(tree.fun(x))
              ms <- apply(dd, 2, mean)
              dd <- apply(dd, 2, function(e) e - mean(e))
              #gmat <- cbind(one, dd)
              if(shpi == 1){
                dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
                ps <- 0
              }
              
              if(shpi == 2){
                dm <- NULL
                ps <- 0
              }
            }
            
            if(shpi != 2 | attr(x, 'type') != 'tree'){
              delta <- cbind(delta, dd)
            }
            
            if(shpi != 2 | attr(x, 'type') != 'tree'){
              varlst <- c(varlst, rep(ix, ncol(dd)))
            }
            
            if(shpi > 2 & shpi < 5 | shpi == 2 & attr(x, 'type') == 'tree'){
              if(shpi > 2 & shpi < 5){
                zmat <- cbind(zmat, x - mean(x)) #remove mean(x)?
                varlst_z <- c(varlst_z, ix)
              }
              if(shpi == 2 & attr(x, 'type') == 'tree'){
                zmat <- cbind(zmat, dd)
                varlst_z <- c(varlst_z, rep(ix, ncol(dd)))
              }
            }
            
            ps <- 0
            if(shpi != 2 | attr(x, 'type') != 'tree'){
              dm <- matrix(0, nrow=ncol(dd), ncol=ncol(dd))
              dm_lst[[ix]] <- dm #no penalty
            } 
          }
          
          if(shpi > 8){
            mult <- ifelse(pnt, 2, 1)
            spl_degree_i <- ifelse(shpi %in% c(9, 10), 2, 3) #use 2 for i-spl, 3 for c-spl
            n1 <- length(unique(x))
            #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/5)) + 5, trunc(n1^(1/7)) + 5)
            nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/7)) + 6, trunc(n1^(1/9)) + 6)
            #new
            #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/5)) + 5 - decrease_edfu, 
            #                      trunc(n1^(1/7)) + 5 + add_edfu)
            edfu <- ifelse(spl_degree_i == 2, nkts/mult, nkts/mult + 1)
            #print (c(nkts, edfu))
            kts <- seq.int(min(x), max(x), length.out = nkts)
            dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = kts, space = space)
            #dd is orth to 1
            ms <- dd_ans$ms
            dd <- t(dd_ans$amat)
            #kts <- dd_ans$knots
            delta <- cbind(delta, dd) #splines excluding conv/conc identity term
            varlst <- c(varlst, rep(ix, ncol(dd)))
            
            if(shpi > 10 & shpi < 13){
              zmat <- cbind(zmat, x - mean(x))
              varlst_z <- c(varlst_z, ix)
            } 
            
            if(pnt){
              nc <- NCOL(dd)
              dm <- makedmat_1D_ic(m = nc, knots = kts, shpi = shpi) 
              qvi <- crossprod(dd)
              dvi <- crossprod(dm)
              ps <- uniroot(f = .search_ps_shapeselect, qv0 = qvi, dv = dvi, dd = dd,
                            edfu = edfu, interval = c(1e-10, 1e+10), tol = .Machine$double.eps^0.25)$root
              
              #print (ps)
              #ps <- .5
              #test!
              #if(shpi > 10) {ps = 0}
              #if(shpi < 11) {ps = 1}
              dm <- sqrt(ps) * dm
              dm_lst[[ix]] <- dm
            }
          }
        }
      }
      #---------------
      #compute etahat 
      #---------------
      gmat <- cbind(one, zmat, delta)
      
      if(!is.null(zmat)){
        zmat_all[[icomb]] <- zmat 
        varlst_z_all[[icomb]] <- varlst_z
      } 
      
      if(!is.null(delta)){
        delta_all[[icomb]] <- delta
        varlst_all[[icomb]] <- varlst
      }
      
      gmat0 <- gmat
      np_ix <- 1 + NCOL(zmat)
      dmat <- as.matrix(bdiag(Filter(Negate(is.null), dm_lst)))
      #new: change dmat to be 0 if it is 1 by 1 0
      # if(dim(dmat)[1] == 1 & dim(dmat)[2] == 1){
      #   dmat <- dmat[1, 1]
      # }
      
      qv0 <- crossprod(gmat)
      xtx <- qv0
      nc <- ncol(gmat)
      amat <- diag(nc)
      ps_zero <- 1:np_ix
      amat[ps_zero, ps_zero] = 0
      bvec <- cbind(rep(0, nc))
      
      if(pnt & !is.null(dmat)){
        dmat_zero = matrix(0, nrow = nrow(dmat), ncol = np_ix)
        dmat2 = cbind(dmat_zero, dmat)
        dv = crossprod(dmat2)
        print (dim(qv0))
        print (dim(dv))
        xtx = qv0 + dv
      }
      
      if(!wt.iter){
        xty = crossprod(gmat, y)
        # if(!pnt){
        #   dsend <- gmat[, -c(1:np_ix), drop = FALSE]
        #   zsend <- gmat[, 1:np_ix, drop = FALSE]
        #   ans <- suppressMessages(try(coneB(y, dsend, zsend), silent = TRUE))
        # }
        #if(pnt & !is.null(dmat)){
        ans = suppressMessages(try(solve.QP(xtx, xty, t(amat), bvec), silent = TRUE))
        #}
        if (!inherits(ans, "try-error")) {
          
          #if(pnt){
          bhat = ans$solution
          face = ans$iact
          #}
          
          # if(!pnt){
          #   face = ans$face
          #   bhat = ans$coefs
          # }
          
          etahat = gmat %*% bhat
          muhat = etahat
        }
      } else {
        etahat <- etahat.fun(n, y, fml = family$family)
        gr <- gr.fun(y, etahat, weights=NULL, fml = family$family)
        wt <- wt.fun(y, etahat, n, weights=NULL, fml = family$family)
        cvec <- wt * etahat - gr
        zvec <- zvec.fun(cvec, wt, y, fml = family$family)
        muhat <- family$linkinv(etahat)
        #sm <- 1e-8
        sm <- 1e-5
        
        if (family$family == "binomial") {
          mdiff <- abs(max(muhat) - 1) > sm
        } else {mdiff <- TRUE}
        
        diff <- 1
        niter <- 0
        
        while (diff > sm & mdiff & niter < 10) 
        {
          oldmu <- muhat
          niter <- niter + 1
          gr <- gr.fun(y, etahat, weights=NULL, fml = family$family)
          wt <- wt.fun(y, etahat, n, weights=NULL, fml = family$family)
          
          cvec <- wt * etahat - gr
          zvec <- zvec.fun(cvec, wt, y, fml = family$family)
          gmat <- sweep(gmat0, MARGIN = 1, sqrt(wt), FUN = '*')
          
          if(!pnt | pnt & all(dmat == 0)){
            dsend <- gmat[, -c(1:np_ix), drop = FALSE]
            zsend <- gmat[, 1:np_ix, drop = FALSE]
            ans <- suppressMessages(try(coneB(zvec, dsend, zsend), silent = TRUE))
            if (inherits(ans, "try-error")) {
              break
            } else {
              face <- ans$face
              bhat <- ans$coefs
            }
          }
          
          if(pnt & !is.null(dmat)){
            xtx = crossprod(gmat)
            xtx = xtx + dv
            xty = crossprod(gmat, zvec)
            ans = suppressMessages(try(solve.QP(xtx, xty, t(amat), bvec), silent = TRUE))
            if (inherits(ans, "try-error")) {
              break
            } else {
              bhat = ans$solution
              face = ans$iact
            }
          }
          
          etahat <- gmat0 %*% bhat
          muhat <- family$linkinv(etahat)
          
          mdiff <- abs(max(muhat) - 1)
          if (family$family == "binomial") {
            mdiff <- abs(max(muhat) - 1) > sm
          } else {mdiff <- TRUE}
          
          diff <- mean((muhat - oldmu)^2)
        }
      }
      
      coefs_lst[[icomb]] <- bhat
      #---------------
      #compute cic
      #---------------
      cores = min(8, parallel::detectCores(logical = FALSE) - 1)
      if(parallel){
        edfs = parallel::mclapply(1:nsim, 
                                  .compute_edf0_fix_shps, gmat=gmat0, np_ix, dmat, 
                                  pnt, qv0, amat, bvec, dv,
                                  family, ps, mc.cores = (cores))
      }else{
        edfs = lapply(1:nsim,
                      .compute_edf0_fix_shps, gmat=gmat0, np_ix, dmat, 
                      pnt, qv0, amat, bvec, dv, family, ps)
      }
      edfs = Filter(Negate(is.null), edfs) 
      edfs = simplify2array(edfs)
      #cat('shpi comb', icomb, '\n')
      #print (edfs)
      if(length(edfs) > 0){
        dfmean <- mean(edfs, na.rm = TRUE)
      }
      #temp
      if(length(edfs) == 0) {
        dfmean <- np_ix
      }
      
      if(wt.iter){
        llh <- llh.fun(y, muhat, etahat, phihat=NULL, n, weights=NULL, fml=family$family)
        #test!
        #llh <- llh * n / 2
        #llh <- 2 * (llh * n / 2) #2 times neg-loglike
      } else {
        rss <- sum((y - etahat)^2)
        llh <- log(rss)
      }
      
      if (infocrit == "CBIC") {
        #cic_val <- n * llh + log(n) * (dfmean)
        #test
        #cic_val <- llh + log(n) * (dfmean) / n
        #new
        #cic_val <- llh + a * log(n) * (dfmean) / n^(6/7)
        if(any(non_flat_ps %in% c(9, 10))){ #smooth incr/decr
          degree <- 6/7
        } else if (any(non_flat_ps > 10 & non_flat_ps < 17)){ #smooth conv/conc
          degree <- 6/7#8/9
        } else if (all(non_flat_ps > 0 & non_flat_ps < 9)) { #ord or parametric
          degree <- 6/7#1
        } else {
          stop('check!')
        }
        cic_val <- llh + a * log(n) * (dfmean) / n^degree
        
        #new
        # if(shpi %in% c(9, 10)){
        #   bic_p <- 3
        # }
        # if(shpi %in% c(11, 12, 13, 14, 15, 16)){
        #   bic_p <- 3 
        # }
        # if(shpi %in% 18){
        #   bic_p <- 2
        # }
        # cic_val <- llh + a * log(n) * (dfmean) * n^(-1) #* n^(-2bic_p/(2*bic_p+1))
      } 
      
      if (infocrit == "CAIC") {
        cic_val <- llh +  2 * (dfmean) / n
      }
      
      mod_tried[icomb + 1, 1:totx] <- shp_comb_i
      mod_tried[icomb + 1, totx + 1] <- cic_val
      mod_tried[icomb + 1, totx + 2] <- dfmean
      mod_tried[icomb + 1, totx + 3] <- llh
      mod_tried[icomb + 1, totx + 4] <- cic_val - llh
      mod_tried[icomb + 1, totx + 5] <- ps
    }
  }
  #----------------------------------------------------------------------
  #step 3: choose the best one
  #----------------------------------------------------------------------
  #mod_tried <- unique(mod_tried)
  cics <- as.vector(mod_tried[, totx+1])
  optim_id <- which.min(cics)[1]
  optim_shp_comb <- mod_tried[optim_id, 1:totx]
  optim_shp_comb <- as.vector(optim_shp_comb)
  ordcic <- order(cics)
  mod_tried_unordered <- mod_tried
  mod_tried <- mod_tried[ordcic, ]
  #----------------------------------------------------------------------
  #step 4: compute etacomps for the best 
  #----------------------------------------------------------------------
  mod_best <- mod_tried[1, ,drop = FALSE]
  n <- length(y)
  etacomps <- matrix(0, nrow = totx, ncol = n)
  etahatkeep <- ifelse(family$family == 'gaussian', rep(mean(y), n), rep(0, n))
  muhatkeep <- family$linkinv(etahatkeep)
  if(optim_id > 1) {
    optim_id_noflat <- optim_id-1
    coefskeep <- coefs_lst[[optim_id_noflat]]
    ddkeep_z_all <- zmat_all[[optim_id_noflat]]
    ddkeep_x_all <- delta_all[[optim_id_noflat]]
    varlst <- varlst_all[[optim_id_noflat]]
    varlst_z <- varlst_z_all[[optim_id_noflat]]
    coefskeep0 <- coefskeep
    coefskeep <- coefskeep[-1] #remove intercept
    capk <- length(varlst_z)
    if(capk > 0){
      zcoefs <- coefskeep[1:capk]
      xcoefs <- coefskeep[-c(1:capk)]
    } else {
      zcoefs <- NULL
      xcoefs <- coefskeep
    }
    #new
    etahatkeep <- cbind(one, ddkeep_z_all, ddkeep_x_all) %*% coefskeep0
    muhatkeep <- family$linkinv(etahatkeep)
    for(i in 1:totx){
      ps <- i
      xi <- xmat[, i]
      shpi <- optim_shp_comb[i]
      #if(shpi == 17 | ps %in% x_rm){#flat
      if(shpi == 17 | shpi == 0){#flat
        next
      } else {
        #one <- cbind(rep(1, n))
        #xi = as.numeric(xi)
        if(shpi == 18) {
          etacomps[ps, ] <- ddkeep_z_all[, varlst_z %in% ps, drop = FALSE] %*% zcoefs[which(varlst_z %in% ps)]
        } else if (shpi == 1 | shpi == 2) {#need to use 1 for incr later
          if(attr(xi, 'type') %in% c('nom', 'nominal') | attr(xi, 'type') %in% 'tree' & shpi == 2){
            etacomps[ps, ] <- (ddkeep_z_all[, which(varlst_z %in% ps), drop = FALSE]) %*% zcoefs[which(varlst_z %in% ps)]
          }
          if(attr(xi, 'type') %in% c('ord', 'ordinal') | attr(xi, 'type') %in% 'tree' & shpi == 1){
            etacomps[ps, ] <- (ddkeep_x_all[, which(varlst %in% ps), drop = FALSE]) %*% xcoefs[which(varlst %in% ps)] 
          }
        } else if (shpi > 10 & shpi < 13 | shpi > 2 & shpi < 5) {
          etacomps[ps, ] <- (ddkeep_x_all[, which(varlst %in% ps), drop = FALSE]) %*% xcoefs[which(varlst %in% ps)] 
          etacomps[ps, ] <- etacomps[ps, ] + (xi-mean(xi)) * zcoefs[which(varlst_z %in% ps)]
        } else {
          etacomps[ps, ] <- (ddkeep_x_all[, which(varlst %in% ps), drop = FALSE]) %*% xcoefs[which(varlst %in% ps)]
        }
      }
    }
  }
  #----------------------------------------------------------------------
  #step 5: Lasso shrink? no need 
  #----------------------------------------------------------------------
  #----------------------------------------------------------------------
  #return result
  #----------------------------------------------------------------------
  if(any(optim_shp_comb == 17)) {
    ps_17 = which(optim_shp_comb == 17)
    optim_shp_comb[ps_17] = 0
  }
  tags <- rep('x', p)
  tags[which(types %in% c('nom', 'nominal'))] <- 'z'
  tags[which(types %in% 'tree')] <- 'tree'
  #shapes in character
  shapes <- mapply(ShapeToChar_Classo, optim_shp_comb, tags)
  rslt <- list(mod_tried = mod_tried, etacomps = etacomps, mod_best = mod_best, shps = optim_shp_comb, shapes=shapes,
               family = family, etahat = etahatkeep, muhat = muhatkeep)
  class(rslt) <- c("shapeselectAll", "shapeselect")
  return (rslt)
}

#####################################
#simulate edf for a fixed shape comb#
#####################################
.compute_edf0_fix_shps <- function(isim, gmat, np_ix, dmat = NULL, pnt = FALSE, 
                                   qv0 = NULL, amat = NULL, bvec = NULL, dv = NULL,
                                   family = gaussian(link = 'identity'), ps=0)
{
  cicfamily <- CicFamily(family)
  llh.fun <- cicfamily$llh.fun
  #new: use log link in gamma
  linkfun <- cicfamily$linkfun
  etahat.fun <- cicfamily$etahat.fun
  gr.fun <- cicfamily$gr.fun
  wt.fun <- cicfamily$wt.fun
  zvec.fun <- cicfamily$zvec.fun
  muhat.fun <- cicfamily$muhat.fun
  ysim.fun <- cicfamily$ysim.fun
  deriv.fun <- cicfamily$deriv.fun
  dev.fun <- cicfamily$dev.fun
  
  eps = 1e-8
  wt.iter = ifelse(family$family == 'gaussian', FALSE, TRUE)
  if(NCOL(gmat) == np_ix) {
    edfi <- np_ix 
  } 
  if (NCOL(gmat) > np_ix) {
    n <- NROW(gmat)
    #ysim <- ysim.fun(n, mu0 = NULL, fml = family$family, shp0 = NULL)
    if(family$family == 'binomial'){
      ysim <- rbinom(n, size = 1, prob = .5)
    }else{
      ysim <- ysim.fun(n, mu0 = NULL, fml = family$family, shp0 = NULL)
    }
    #-------------------------------------------
    #use qprog
    #nc = ncol(gmat)
    #amat = diag(nc)
    #ps_zero = 1:np_ix
    #amat[ps_zero, ps_zero] = 0
    #bvec = cbind(rep(0, nc))
    gmat0 = gmat
    #qv0 = crossprod(gmat0)
    xtx = qv0
    
    #use this to test edfi
    #if(round(ps, 8) == 0){
    #   pnt = FALSE
    #}
    
    if(pnt & !is.null(dmat)){
      #dmat_zero = matrix(0, nrow = nrow(dmat), ncol = np_ix)
      #dmat2 = cbind(dmat_zero, dmat)
      #dv = crossprod(dmat2)
      xtx = qv0 + dv
    }
    #-------------------------------------------
    if(!wt.iter){
      xty = crossprod(gmat, ysim)
      #if(pnt){
      ans = suppressMessages(try(solve.QP(xtx, xty, t(amat), bvec), silent = TRUE))
      #}
      # if(!pnt){
      #   dsend <- gmat[, -c(1:np_ix), drop = FALSE]
      #   zsend <- gmat[, 1:np_ix, drop = FALSE]
      #   ans <- suppressMessages(try(coneB(ysim, dsend, zsend), silent = TRUE))
      # }
      if (!inherits(ans, "try-error")) {
        #if(pnt){
        bhat = ans$solution
        face = ans$iact
        #}
        # if(!pnt){
        #   face <- ans$face
        #   bhat <- ans$coefs
        # }
      }
    } else {
      etahat <- etahat.fun(n, ysim, fml = family$family)
      gr <- gr.fun(ysim, etahat, weights=NULL, fml = family$family)
      wt <- wt.fun(ysim, etahat, n, weights=NULL, fml = family$family)
      cvec <- wt * etahat - gr
      zvec <- zvec.fun(cvec, wt, ysim, fml = family$family)
      muhat <- family$linkinv(etahat)
      #sm <- 1e-8
      sm <- 1e-5
      
      if (family$family == "binomial") {
        mdiff <- abs(max(muhat) - 1) > sm
      } else {mdiff <- TRUE}
      
      diff <- 1
      niter <- 0
      
      while (diff > sm & mdiff & niter < 10) 
      {
        oldmu <- muhat
        niter <- niter + 1
        gr <- gr.fun(ysim, etahat, weights=NULL, fml = family$family)
        wt <- wt.fun(ysim, etahat, n, weights=NULL, fml = family$family)
        cvec <- wt * etahat - gr
        zvec <- zvec.fun(cvec, wt, ysim, fml = family$family)
        gmat <- sweep(gmat0, MARGIN = 1, sqrt(wt), FUN = '*')
        
        if(!pnt | pnt & all(dmat == 0)){
          dsend <- gmat[, -c(1:np_ix), drop = FALSE]
          zsend <- gmat[, 1:np_ix, drop = FALSE]
          #test!
          ans <- suppressMessages(try(coneB(zvec, dsend, zsend), silent = TRUE))
          if (inherits(ans, "try-error")) {
            break
          } else {
            face <- ans$face
            bhat <- ans$coefs
          }
        }
        
        if(pnt & !is.null(dmat)){
          xtx = crossprod(gmat)
          xtx = xtx + dv
          #test!
          xty = crossprod(gmat, zvec)
          ans = suppressMessages(try(solve.QP(xtx, xty, amat, bvec), silent = TRUE))
          #ans = suppressMessages(try(solve.QP(xtx, xty, t(amat), bvec), silent = TRUE))
          if (inherits(ans, "try-error")) {
            break
          } else {
            bhat = ans$solution
            face = ans$iact
          }
        }
        
        etahat <- gmat0 %*% bhat
        muhat <- family$linkinv(etahat)
        
        mdiff <- abs(max(muhat) - 1)
        if (family$family == "binomial") {
          mdiff <- abs(max(muhat) - 1) > sm
        } else {mdiff <- TRUE}
        
        diff <- mean((muhat - oldmu)^2)
      }
    }
    
    if(inherits(ans, "try-error")) {
      edfi <- NULL
    } else {
      if(!pnt | pnt & all(dmat == 0)){
        #test
        bhat <- round(bhat, 10)
        edfi <- sum(abs(na.omit(bhat)) > eps)
      } else {
        #edfi <- sum(abs(na.omit(bhat)) > eps)
        umat = chol(xtx)
        uinv = backsolve(umat, diag(ncol(xtx)))
        atil = amat %*% uinv
        xtxinv = tcrossprod(uinv)
        if(!wt.iter){
          bigp_base = xtxinv %*% qv0
        }else{
          bigp_base = xtxinv %*% crossprod(gmat)
        }
        edfi_base = sum(diag(bigp_base))
        
        if (all(face == 0) | length(face) == 0) {
          edfi = edfi_base
        } else {
          smat = -t(atil[face, ,drop = FALSE])
          #new: avoid solve error
          qr_smat = qr(smat)
          smat = qr.Q(qr_smat)
          #smat = qr.Q(qr_smat)[, 1:(qr_smat$rank), drop = FALSE]
          pmat_polar = tcrossprod(smat)
          if(!wt.iter){
            bigp_polar = uinv %*% tcrossprod(pmat_polar, uinv) %*% qv0
          }else{
            bigp_polar = uinv %*% tcrossprod(pmat_polar, uinv) %*% crossprod(gmat)
          }
          edfi = edfi_base - sum(diag(bigp_polar))
        }
      }
    }
  }
  return (edfi)
}


























































