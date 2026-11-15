#' Shape selection for constrained additive models
#'
#' Fits an additive model by selecting variables with shape or order constraints
#' for each predictor using Cone Information Bayesian Criterion (CBIC).
#' The function can perform forward selection or stepwise (forward/backward) updates
#' over a set of candidate shapes appropriate to each predictor type.
#'
#' @param xmat A data frame or matrix of predictors. Numeric predictors should
#'   be used with `types = "num"`, and factor predictors with
#'   `types = "nom"`. Column names are used in the output when available.
#'   Numeric predictors are fitted using shape-constrained smooth functions,
#'   whereas nominal predictors are included in or excluded from the model.
#'
#' @param y A response vector with length equal to `nrow(xmat)`. For a Gaussian
#'   model, `y` should be numeric; for a binomial model, it should be binary.
#'
#' @param types A character vector of length `ncol(xmat)` specifying the type
#'   of each predictor. Currently supported values are `"num"` for numeric
#'   predictors and `"nom"` for nominal categorical predictors.
#'
#' @param direction The direction of the model search. `"forward"` performs
#'   forward selection of predictors and their shapes. `"both"` performs
#'   forward selection followed by backward revisions, including the removal
#'   of predictors. The default is `"both"`.
#'
#' @param family A family object or a character string naming a family
#'   function. It specifies the error distribution and link function, as in
#'   `glm()`. Currently supported families are Gaussian and binomial.
#'
#' @param boundary.knots An optional list of boundary knots for the numeric
#'   predictors. The default is `NULL`. For a numeric predictor, the
#'   corresponding element should be a numeric vector containing at least two
#'   finite values, such as the values returned by
#'   `range(x, na.rm = TRUE)`. For a nominal predictor, the corresponding
#'   element should be `NULL`.
#'
#' @param nsim The number of simulations used to estimate the mean effective
#'   degrees of freedom in the CBIC calculation at each selection step.
#'
#' @param pnt Logical; if `TRUE`, include penalty-related calculations during
#'   model fitting. The default is `TRUE`.
#'  
#' @param parallel Logical. Whether to enable parallel computation. The default is FALSE. 
#'   If parallel is TRUE, then the default number of cores = 2. 
#'   To use a different number of cores for example, to use 4 cores, run: options(cgam.cores = 4L)
#'
#' @param user_edfus An optional numeric vector of three positive finite
#'   starting effective degrees of freedom. If `NULL`, starting values are
#'   determined automatically and tuned before variable selection. This
#'   argument is primarily intended for package maintenance and development.
#'
#' @param scale.method The scaling method for numeric predictors. `"minmax"`
#'   scales each predictor to the unit interval, `"zscore"` standardizes each
#'   predictor using its mean and standard deviation, and `"none"` performs no
#'   scaling. The default is `"minmax"`.
#'
#' @param set4 Logical; if `TRUE`, consider only the four basic smooth shapes
#'   for numeric predictors: increasing, decreasing, convex, and concave.
#'   If `FALSE`, also consider combined monotonicity-and-curvature shapes.
#'
#' @param tune The method used to tune the effective degrees of freedom of
#'   selected numeric predictors. Choices are `"cv5"`, `"cv10"`, `"gcv"`,
#'   and `"none"`. The default is `"cv5"`.
#'
#' @param cv_seed An integer random seed used for cross-validation when
#'   `tune` is not `"none"`. The default is `123`.
#'   
#' @return An object of class `"shapeselect"` containing the fitted model and
#'   the selection results. Important components include:
#'   \describe{
#'     \item{shps}{Numeric shape code for each predictor; zero denotes an
#'       unselected predictor.}
#'     \item{shapes}{Character labels describing the selected shape for each
#'       predictor.}
#'     \item{active}{Indices of the selected predictors.}
#'     \item{active_shps}{Shape codes for the selected predictors.}
#'     \item{mod_tried}{Data frame recording the selection path and criterion
#'       values.}
#'     \item{etahat}{Fitted linear predictors.}
#'     \item{muhat}{Fitted mean responses on the response scale.}
#'     \item{etacomps}{Data frame of fitted component contributions, with one row per original predictor.}
#'     \item{coefskeep}{Estimated coefficients for the final model.}
#'   }
#'
#' For numeric predictors, the shape codes are 9 (smooth increasing), 10
#' (smooth decreasing), 11 (smooth convex), 12 (smooth concave), 13 (smooth
#' increasing convex), 14 (smooth increasing concave), 15 (smooth decreasing
#' convex), 16 (smooth decreasing concave), and 19 (smooth unconstrained).
#' When `set4 = TRUE`, only codes 9--12 and 19 are considered. Code 0 denotes
#' a numeric predictor that is not selected; the corresponding label in
#' `shapes` is `"out of the model"`. Code 17 represents a flat shape when it
#' occurs in the fitted object.
#'
#' For nominal predictors, code 1 denotes that the predictor is included in
#' the model, and code 0 denotes that it is not selected. Their corresponding
#' labels in `shapes` are `"in the model"` and `"out of the model"`,
#' respectively.
#' 
#' @details
#' The algorithm begins with an intercept-only fit, then evaluates each inactive predictor under
#' a set of allowable shapes (depending on `types`). At each step it selects the predictor/shape
#' combination that improves the chosen information criterion. 
#'
#' @examples
#' # Example 1: gaussian response, two numeric predictors, unpenalized
#' library(MASS)
#' data("Rubber", package = "MASS")
#' y <- Rubber[, 1]
#' xmat <- Rubber[, -1]
#' # the sample size = 30. For a small sample, it is better to use pnt = FALSE
#' ans <- shapeselect(xmat, y, types = c("num", "num"), pnt = FALSE)
#' # print out the chosen predictors and their shapes
#' print (ans)
#' # plot chosen components
#' plot(ans)
#' # make a prediction
#' xtest <- data.frame(hard = mean(Rubber$hard), tens = mean(Rubber$tens))
#' pfit <- predict(ans, xtest, type = "response")$fit
#' 
#' # Example 2: binomial response, one nominal predictor and six numeric predictors, penalized 
#' # each numeric predictor has 4 shape options
#' data(Pima.tr, package = "MASS")
#'
#' df <- Pima.tr |> dplyr::filter(age <= 70) |>
#'   dplyr::mutate(npreg_grp = cut(npreg, breaks = c(-Inf, 0, 1, 5, 9, Inf),
#'       labels = c("0", "1", "2", "3", "4"), ordered_result = FALSE)) |>
#'   dplyr::select(-npreg) |>
#'   dplyr::select(npreg_grp, dplyr::everything())
#'
#' y <- as.integer(df$type == "Yes")
#' xmat <- df[, 1:7, drop = FALSE]
#' 
#' # Default nsim = 500
#' \dontrun{
#' ans <- shapeselect(xmat, y, types = c("nom", rep("num", 6)), family = "binomial",
#'   nsim = 200, pnt = TRUE, set4 = TRUE)
#'
#' print(ans)
#' plot(ans)
#' }
#' @export
#-------------------------------------------------------------------------------------------------------
#main function
#-------------------------------------------------------------------------------------------------------
shapeselect <- function(xmat, y, types = NULL, direction = c("both", "forward"), family = gaussian(link = "identity"), 
                        boundary.knots = NULL, nsim = 500, pnt = TRUE, parallel = FALSE,
                        user_edfus = NULL, scale.method = c("minmax", "zscore", "none"), set4 = FALSE, 
                        tune = c("cv5", "cv10", "gcv", "none"), cv_seed = 123){ 
  
  direction <- match.arg(direction)
  #new
  scale.method <- match.arg(scale.method)
  tune <- match.arg(tune)
  
  #new: force xmat to be a data frame
  if (is.matrix(xmat)) {
    warning("xmat is a matrix; factor information may be lost.", "please provide xmat as a data.frame.")
    xmat <- as.data.frame(xmat)
  }
  
  n <- length(y)
  p <- ncol(xmat)
  if (n != nrow(xmat)) {
    stop("'length(y)' must equal 'nrow(xmat)'.")
  }
  if (is.null(types) || length(types) != p) {
    stop("'types' must be provided as one of the 2: 'nom' or 'num' and have length equal to ncol(xmat).")
  }
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family))
    stop("'family' not recognized!")

  wt.iter <- (family$family == "binomial")
  xnms <- colnames(xmat) %||% paste0("x", 1:p)
  xmat_0 <- xmat #used xmat_0 later in plot
  #pre-process xmat
  #ans_process <- process_xmat(xmat, types, boundary.knots = boundary.knots, standardize = standardize, center = center, 
  #                            add_kts = add_kts, pnt = pnt)
  ans_process <- process_xmat(xmat, types, boundary.knots = boundary.knots, scale.method = scale.method,
                              add_kts = 0, pnt = pnt)
  xmat <- ans_process$xmat
  boundary.knots <- ans_process$boundary.knots
  m_total <- ans_process$m_total
  mx <- ans_process$mx
  sdx <- ans_process$sdx
  names(mx) <- names(sdx) <- xnms
  add_kts <- ans_process$add_kts
  pnt <- ans_process$pnt
  
  #initialize global model state
  active_0 <- integer(0)
  active_shps_0 <- integer(0)
  inactive_0 <- seq_len(p)
  #compute step 0 cbic
  #allow cbic only
  infocrit <- "CBIC"
  cic_new <- calculate_initial_cic(y, family, a = 1, infocrit)
  etahatkeep <- if(wt.iter) rep(0, n) else rep(mean(y), n)
  
  #set up the history matrix to keep track of each predictor's adding/deleting, adding shapes
  mod_tried <- data.frame(matrix(0, nrow = 0, ncol = (p + 5)))
  xnms_with_id <- paste0(xnms, '|', 1:p)
  colnames(mod_tried) <- c(xnms_with_id, "# of x's", infocrit, "penalty", "new add", "action")
  #current model state: 0 means not currently in the model
  current_shps <- rep(0, p)
  
  #pre-compute objects that will be used a lot!
  #params <- list(xmat = xmat, types = types, boundary.knots = boundary.knots, pnt = pnt, set4 = set4, 
  #               wt.iter = wt.iter, infocrit = infocrit, add_kts = add_kts, add_edfu = add_edfu, 
  #               space = space, edfu_plan = edfu_plan, step = 1) #step=1 can be removed.
  #new: define starting edfus to depend on n
  
  # validate user-supplied starting edfus
  if (!is.null(user_edfus) && (!is.numeric(user_edfus) || length(user_edfus) != 3L || anyNA(user_edfus) ||
                               any(!is.finite(user_edfus)) || any(user_edfus <= 0))) {
    warning(
      "'user_edfus' must be a numeric vector containing three positive, ",
      "finite values. Using automatically selected starting EDFs instead."
    )
    user_edfus <- NULL
  }
  
  # defaults when initial EDF tuning is skipped
  initial_r2 <- NA_real_
  edfu_scaler <- rep(1, 3)
  tuned_edfus <- NULL
  weights <- NULL
  space <- "E"
  # null or vector with 3 values; if user provides edfus, then avoid the initial tuning
  start_edfus <- user_edfus 
  # not tune edfus if unpenalized
  if(pnt && is.null(user_edfus)){
    edfu_m <- 2 * n^(1/7) - 1 #match the one used in 
    user_edfus <- c(edfu_m, edfu_m + 1, edfu_m + 7)
    start_edfus <- user_edfus 
    
    #at this step: user_edfus will affect the penalty matrix, but tune_start_edfu will 
    #use the one without absorbing the penalty term
    precomps_initial <- precompute_lists(xmat = xmat, types = types, boundary.knots = boundary.knots, pnt = pnt, set4 = set4,
                                         wt.iter = wt.iter, infocrit = infocrit, add_kts = 0, add_edfu = 0, 
                                         user_edfus = user_edfus, space = space, edfu_plan = 1)
    
    #new:
    user_edfus <- tune_start_edfu(precomps_initial, y, user_edfus = user_edfus, family = family, 
                                  r2_gam = NULL)
    initial_r2 <- attr(user_edfus, "initial_r2")
    edfu_scaler <- attr(user_edfus, "edfu_scaler") 
    tuned_edfus <- user_edfus 
    #new: binomial weights in x'wx
    weights <- attr(user_edfus, "weights")
  }
  
  #update the penalty terms and penalty matrix in precomps using the tuned user_edfus
  precomps <- precompute_lists(xmat = xmat, types = types, boundary.knots = boundary.knots, pnt = pnt, set4 = set4,
                               wt.iter = wt.iter, infocrit = infocrit, add_kts = 0, add_edfu = 0, 
                               user_edfus = user_edfus, space = space, edfu_plan = 1, weights = weights)
  
  #initialize objects that will be kept
  #_0 kept copy
  varlst_0 <- varlst_s_0 <- varlst_z_0 <- NULL #track 19, other non-flat shapes for x, nom
  ddkeep_x <- ddkeep_z <- NULL
  ddkeep_x_all_0 <- dv_keep_all_0 <- ddkeep_z_all_0 <-  NULL 
  #dv_keep_lst is used in the backward step
  dv_keep_lst_0 <- vector("list", length = p)
  #for plot and predict
  ps_keep <- mskeep <- ktskeep <- gmatkeep <- coefskeep <- NULL 
  check <- length(inactive_0) >= 1
  #new
  edf0keep <- 1
  #main stepwise loop
  nreps <- 0
  #new: num_active_nom used in backward step to compute minedf0 (temp)
  num_active_nom <- 0
  while(check) {
    success_fwd <- FALSE
    success_bwd <- FALSE
    #print (mod_tried)
    cic_old <- cic_new
    #oldetahat <- etahatkeep
    if(length(inactive_0) >= 1) {
      #forward step 
      ans_fwd <- forward_cic(xmat, y, cic_old, active_0 = active_0, active_shps_0 = active_shps_0, 
                             nsim = nsim, family = family, oldetahat = etahatkeep,
                             parallel = parallel, pnt = pnt, a = 1, set4 = set4, precomps = precomps, 
                             ddkeep_x_all_0 = ddkeep_x_all_0, 
                             ddkeep_z_all_0 = ddkeep_z_all_0, 
                             dv_keep_all_0 = dv_keep_all_0, 
                             dv_keep_lst_0 = dv_keep_lst_0,
                             varlst_s_0 = varlst_s_0, varlst_0 = varlst_0, 
                             varlst_z_0 = varlst_z_0, oldcoefs = coefskeep, 
                             oldgmat = gmatkeep, oldedf0 = edf0keep, 
                             oldnum_active_nom = num_active_nom)
      #cat('finished fwd', '\n')
      #keep these
      cic_new <- ans_fwd$cic_new
      picked_x <- ans_fwd$picked_x
      picked_shp <- ans_fwd$picked_shp
      #new: used to update num_active_nom
      num_active_nom <- ans_fwd$num_active_nom
      ps_keep <- ans_fwd$ps_keep
      success_fwd <- !is.null(picked_x) && !is.null(picked_shp)
      
      if(success_fwd){
        current_shps[picked_x] <- picked_shp
        mod_tried <- append_mod_tried(mod_tried, current_shps, cic = cic_new, ps_keep = ps_keep,
                                      xid = picked_x, action = "add", xnms_with_id = xnms_with_id)
        active_0 <- c(active_0, picked_x)
        active_shps_0 <- c(active_shps_0, picked_shp) #used for backward
        inactive_0 <- inactive_0[inactive_0 != picked_x]
        #put the following in a function:
        ddkeep_x <- ans_fwd$ddkeep_x
        ddkeep_z <- ans_fwd$ddkeep_z
        if (!is.null(ddkeep_x)) {
          ddkeep_x_all_0 <- cbind(ddkeep_x_all_0, ddkeep_x)
          varlst_0 <- c(varlst_0, rep(picked_x, ncol(ddkeep_x)))
          if (picked_shp == 19) {
            varlst_s_0 <- c(varlst_s_0, rep(1, ncol(ddkeep_x)))
          } else {
            varlst_s_0 <- c(varlst_s_0, rep(0, ncol(ddkeep_x)))
          }
          if (pnt) {
            dmat_keep <- ans_fwd$dmat_keep
            dv_keep <- ans_fwd$dv_keep
            #dv_keep_lst_0[[picked_x]] <- dv_keep #ignore np
            dv_keep_lst_0[picked_x] <- list(dv_keep) #dv_keep might be NULL; pnt = T but we choose a z
            dv_keep_all_0 <- combine_penalties(fixed_dv = dv_keep_all_0, new_dmat = dmat_keep) #ignore np
          }
        }
        if (!is.null(ddkeep_z)) {
          ddkeep_z_all_0 <- cbind(ddkeep_z_all_0, ddkeep_z)
          varlst_z_0 <- c(varlst_z_0, rep(picked_x, ncol(ddkeep_z)))
        }
        
        etahatkeep <- ans_fwd$etahatkeep
        gmatkeep <- ans_fwd$gmatkeep
        coefskeep <- ans_fwd$coefskeep
        edf0keep <- ans_fwd$edf0keep
      }
    } else {
      picked_x <- picked_shp <- NULL
    }
    
    #backward step
    if (direction == "both" && length(active_0) >= 2){
      #test
      if(!is.null(picked_shp)){ #if the previous forward step succeeded, then skip checking the most recent predictor
        check_set <- active_0[-length(active_0)]
        check_shps <- active_shps_0[-length(active_0)]
      } else {
        check_set <- active_0
        check_shps <- active_shps_0
      }
      
      #check_set <- active_0[-length(active_0)]
      #check_shps <- active_shps_0[-length(active_0)]
      #the cic_old in the back ward step is the most recent one returned by th forward step
      #not assign cic_new here; if it is not improved, then this assignment will stop the checking
      ans_bwd <- backward_cic(xmat, y, cic_old = cic_new, active_0 = active_0, 
                              active_shps_0 = active_shps_0, check_set = check_set, check_shps = check_shps,
                              nsim = nsim, family = family, oldetahat = etahatkeep, parallel = parallel, pnt = pnt,
                              a = 1, precomps = precomps, ddkeep_x_all_0 = ddkeep_x_all_0, ddkeep_z_all_0 = ddkeep_z_all_0,
                              dv_keep_lst_0 = dv_keep_lst_0, varlst_s_0 = varlst_s_0, varlst_0 = varlst_0,
                              varlst_z_0 = varlst_z_0, dv_keep_all_0 = dv_keep_all_0, 
                              xnms_with_id = xnms_with_id, mod_tried = mod_tried,
                              current_shps = current_shps, oldcoefs = coefskeep, 
                              oldgmat = gmatkeep, oldedf0 = edf0keep, oldnum_active_nom = num_active_nom)
      #cat('finished bwd', '\n')
      success_bwd <- ans_bwd$has_changed
      if(success_bwd){
        cic_new <- ans_bwd$cic_new
        active_0 <- ans_bwd$active_0
        #new
        num_active_nom <- ans_bwd$num_active_nom
        active_shps_0 <- ans_bwd$active_shps_0
        inactive_0 <- ans_bwd$inactive_0
        etahatkeep <- ans_bwd$etahatkeep
        gmatkeep <- ans_bwd$gmatkeep
        coefskeep <- ans_bwd$coefskeep
        ddkeep_x <- ans_bwd$ddkeep_x
        ddkeep_z <- ans_bwd$ddkeep_z
        ddkeep_x_all_0 <- ans_bwd$ddkeep_x_all_0
        ddkeep_z_all_0 <- ans_bwd$ddkeep_z_all_0 
        varlst_0 <- ans_bwd$varlst_0
        varlst_s_0 <- ans_bwd$varlst_s_0
        varlst_z_0 <- ans_bwd$varlst_z_0
        dv_keep_all_0 <- ans_bwd$dv_keep_all_0
        #print (dim(dv_keep_all_0))
        #print (dim(gmatkeep))
        dv_keep_lst_0 <- ans_bwd$dv_keep_lst_0
        mod_tried <- ans_bwd$mod_tried
        current_shps <- ans_bwd$current_shps
        edf0keep <- ans_bwd$edf0keep
      }
    }
    #check stopping criterion
    #check <- cic_new < cic_old && length(inactive_0) > 0 
    nreps <- nreps + 1
    #tol <- 1e-8
    #check <- cic_old - cic_new > tol && nreps < 100 
    check <- (success_fwd || success_bwd) && nreps < 100 
  }
  #cat('finish choosing', '\n')
  # final assembly 
  shps <- if (nrow(mod_tried) > 0) {as.numeric(mod_tried[nrow(mod_tried), 1:p])} else {rep(0, p)}
  capk <- .nc_safe(ddkeep_z_all_0)
  shapes <- shape_tags_to_char(shps, types)
  #new: tune edfus for the kept predictors by gcv -> get new coefskeep
  #need to update mod_tried' penalty too
  edfu_cv <- edfu_start <- ps_cv <- cv_min <- cv_start <- NULL
  a_cv <- rep(1, sum(types[active_0] == "num"))
  
  if (!pnt && tune != "none") {
    #message("'tune' is ignored when pnt = FALSE; using tune = 'none'.")
    tune <- "none"
  }
  
  if(tune != "none" && any(types[active_0] == "num") && pnt){
    cv_rslt <- tune_edfu(active=active_0, active_shps=active_shps_0, gmatkeep=gmatkeep, 
                         varlst_s=varlst_s_0, np_ix=(1+capk), precomps=precomps, y=y, 
                         types_active=types[active_0], wt.iter=wt.iter, 
                         xnms_active=xnms[active_0], tune=tune, cv_seed=cv_seed, class_cut=0.5, muhatkeep=family$linkinv(etahatkeep))
    #cat('finished tuning', '\n')
    #all terms below follow active_0 order not xmat order
    coefskeep <- cv_rslt$coefskeep
    etahatkeep <- cv_rslt$etahatkeep #muhat is updated in the 'fit' object
    ps_cv <- cv_rslt$ps
    #mod_tried$penalty <- ps_gcv #not work if there is any "change"
    edfu_cv <- cv_rslt$edfu_cv
    edfu_start <- cv_rslt$edfu_start
    cv_min <- cv_rslt$cv_min
    cv_start <- cv_rslt$cv_start
    a_cv <- cv_rslt$a_cv
    active_tuned <- cv_rslt$active_tuned
    mod_tried <- update_cv_penalties(mod_tried, active_tuned=active_tuned, ps_cv=ps_cv,
                                     edfu_start=edfu_start, edfu_cv=edfu_cv)
    #shps <- if (nrow(mod_tried) > 0) {as.numeric(mod_tried[nrow(mod_tried), 1:p])} else {rep(0, p)}
    #shapes <- shape_tags_to_char(shps, types)
  }
  etacomps <- assemble_shapechoice(active_0, active_shps_0, xmat, types, coefskeep, 
                                   varlst_z_0, varlst_0, ddkeep_z_all_0, ddkeep_x_all_0)
  etacomps <- as.data.frame(etacomps)
  rownames(etacomps) <- xnms
  rownames(mod_tried) <- NULL
  #useful to users: etahat, muhat, etacomps for plot, shps, shapes, active
  fit <- list(xmat_0 = xmat_0, xmat = xmat, xnms = xnms, mod_tried = mod_tried, 
              etahat = as.numeric(etahatkeep), etacomps = etacomps, bigmat = gmatkeep, coefskeep = coefskeep, 
              muhat = family$linkinv(etahatkeep), types = types, family = family, pnt = pnt,
              ddkeep_x_all = ddkeep_x_all_0, ddkeep_z_all = ddkeep_z_all_0, np = (1+capk), capk = capk, 
              varlst_z = varlst_z_0, varlst = varlst_0, varlst_s = varlst_s_0,
              shapes = shapes, shps = shps, active = active_0, boundary.knots = boundary.knots,
              active_shps = active_shps_0, precomps = precomps, sdx = sdx, mx = mx, m_total = m_total,
              edf0keep = edf0keep, edfu_cv = edfu_cv,  edfu_start = edfu_start, ps_cv = ps_cv, tuned_edfus = tuned_edfus,
              cv_min = cv_min, cv_start = cv_start, a_cv = a_cv, start_edfus = start_edfus,
              initial_r2 = initial_r2, edfu_scaler = edfu_scaler)
  class(fit) <- "shapeselect"
  return(fit)
}

#-----------------------------------------------------------------
#new: add a column showing the tuned penalty for numeric active xs
#-----------------------------------------------------------------
update_cv_penalties <- function(mod_tried, active_tuned, ps_cv, edfu_start, edfu_cv) {
  mod_tried$penalty_cv <- NA_real_
  mod_tried$edfu_start <- NA_real_
  mod_tried$edfu_cv <- NA_real_
  for (i in seq_along(active_tuned)) {
    xid <- active_tuned[i]
    rows_x <- which(mod_tried$`new add` == xid)
    if (length(rows_x) > 0L) {
      last_row <- tail(rows_x, 1)
      mod_tried$penalty_cv[last_row] <- ps_cv[i]
      mod_tried$edfu_start[last_row] <- edfu_start[i]
      mod_tried$edfu_cv[last_row] <- edfu_cv[i]
    }
  }
  return (mod_tried)
}

#------------------------------------------------------------
#tune starting edfu according to r2 and explained deviance
#only used for pnt = T
#------------------------------------------------------------
tune_start_edfu <- function(precomps, y, user_edfus = c(3, 3, 10), 
                            family = gaussian(link = "identity"),
                            r2_gam = NULL) 
{
  if (!is.numeric(user_edfus) || length(user_edfus) == 0L) {
    stop("user_edfus must be a nonempty numeric vector.")
  }
  
  if(is.null(r2_gam)){
    #new: if family is binomial, let the initial weight to be mu(1-mu) where mu is simple average
    if(family$family == "binomial"){
      mu0 <- rep(mean(y), length(y))
      wt0 <- mu0*(1-mu0)
    } else {
      wt0 <- NULL
    }
    
    # use the precomps to get an intial gmat and dmat
    p <- length(precomps)
    dd_lst <- dv_lst <- vector("list", length = p)
    zmat <- xmat <- edfu_vec <- lam_vec <- NULL
    
    # not follow active set order, active = NULL, follow x1,x2,...,xp in xmat  
    for(ix in seq_len(p)){
      pc_ix <- precomps[[ix]]
      type_ix <- pc_ix$type_ix
      
      if(type_ix == "num"){
        dv_ix <- pc_ix$dv_lst_raw$`9` #convex spl with identity
        dv_lst[[ix]] <- dv_ix
        edfu_vec <- c(edfu_vec, 1*user_edfus[1]) 
      }
      
      if(type_ix == "num"){
        dd_ix <- pc_ix$dd_lst$`9` #ignore shapes now, just use the basis for increasing
        dd_lst[[ix]] <- dd_ix
        xmat <- cbind(xmat, dd_ix)
      } else if (type_ix == "nom" || type_ix == "tree" || type_ix == "ord"){
        dd_z <- pc_ix$dd_lst$`1`
        zmat <- cbind(zmat, dd_z)
      }
    }
    
    current_gmat <- cbind(1, zmat, xmat)
    
    if (!is.null(edfu_vec)){ #at least 1 numeric x
      dv_lst <- Filter(Negate(is.null), dv_lst)
      dd_lst <- Filter(Negate(is.null), dd_lst)
      p_num <- length(dd_lst)
      for(ix in seq_len(p_num)){
        dd <- dd_lst[[ix]]
        dv <- dv_lst[[ix]]
        qv <- crossprod(dd)
        edfu <- edfu_vec[ix]
        
        # cap:
        # nkts <- ncol(dd) - 1 #when we used 13
        nkts <- ncol(dd)
        #edfu <- min(max(edfu, 1.5), nkts - 1)
        edfu <- min(max(edfu, 2), nkts - 1)
        
        log_lam <- tryCatch({
          log_lo <- -10
          log_hi <- 18
          f_lo <- try(.search_ps_sel(exp(log_lo), qv0 = qv, dv = dv, dd = dd, edfu = edfu, weights = wt0), silent = TRUE)
          
          #if the lower endpoint fails, move it inward
          while (inherits(f_lo, "try-error") && log_lo < log_hi) {
            log_lo <- log_lo + 1
            f_lo <- try(.search_ps_sel(exp(log_lo), qv0 = qv, dv = dv, dd = dd,
                                   edfu = edfu, weights = wt0), silent = TRUE)
          }
          
          f_hi <- try(.search_ps_sel(exp(log_hi), qv0 = qv, dv = dv, dd = dd, edfu = edfu, weights = wt0), silent = TRUE)
          
          # if the upper endpoint fails, move it inward
          while (inherits(f_hi, "try-error") && log_hi > log_lo) {
            log_hi <- log_hi - 1
            f_hi <- try(.search_ps_sel(exp(log_hi), qv0 = qv, dv = dv, dd = dd, edfu = edfu, weights = wt0), silent = TRUE)
          }
          
          # check that two valid endpoints were found
          if (inherits(f_lo, "try-error") || inherits(f_hi, "try-error")) {
            stop("could not find two valid penalty endpoints.")
          }
          
          if (f_lo <= 0) {
            #root is at or below the smallest usable log penalty
            log_lo
          } else if (f_hi >= 0) {
            # root is at or above the largest usable log penalty
            log_hi
          } else {
            # f_lo > 0 and f_hi < 0, so the root is bracketed
            log_root_res <- try(
              uniroot(f = function(log_pen, ...) {.search_ps_sel(exp(log_pen), ...)},
                      interval = c(log_lo, log_hi), qv0 = qv, dv = dv, dd = dd, edfu = edfu, weights = wt0, 
                      extendInt = "no", tol = .Machine$double.eps^0.25), silent = TRUE)
            
            if (!inherits(log_root_res, "try-error")) {
              log_root_res$root
            } else {
              # both endpoints are valid and bracket the root.
              # use the endpoint whose function value is closest to zero.
              if (abs(f_lo) <= abs(f_hi)) {
                log_lam_fallback <- log_lo
              } else {
                log_lam_fallback <- log_hi
              }
              cat("unexpected try-error in uniroot\n")
              cat("edfu:", edfu, "\n")
              cat("usable log interval:", log_lo, log_hi, "\n")
              #cat("endpoint values:", f_lo, f_hi, "\n")
              
              log_lam_fallback
            }
          }
          
        }, error = function(err) {
          
          #cat("Failed to find log_lam\n")
          #cat("edfu:", edfu, "\n")
          #cat("Error:", conditionMessage(err), "\n")
          
          NA_real_
        })
        
        lam <- exp(log_lam)
        lam_vec <- c(lam_vec, lam)
      }
      
      dv_lst_tmp <- Map(function(D, lam) lam * D, dv_lst, lam_vec)
      current_dv <- Matrix::bdiag(dv_lst_tmp) |> as.matrix()
      
    } else {
      current_dv <- NULL
    }
    
    if(!is.null(zmat)){
      np_ix <- 1L + .nc_safe(zmat) 
    } else {
      np_ix <- 1L
    }
    
    if(!is.null(xmat)){
      varlst_s <- rep(0, .nc_safe(xmat))
    } else {
      varlst_s <- NULL
    }
    
    # call the fit engine
    # new: if initial_tune, then do unconstrained fit
    initial_fit <- fit_model_engine(y, current_gmat, current_dv, family, 
                                    pnt = TRUE, np_ix, varlst_s, is_first_step = FALSE,
                                    compute_gcv = FALSE, initial_tune = TRUE) 
    
    # obtain r2 or deviance explained
    etahat <- initial_fit$etahat
    # new: get w for binomial
    if(family$family == "binomial"){
      muhat <- family$linkinv(etahat)
      wvec <- muhat * (1 - muhat)
      wvec <- round(wvec, 10)
    } else {
      wvec <- rep(1, length(etahat))
    }
    r2 <- compute_dev(y, etahat, family)
    
  } else {
    r2 <- r2_gam
    #should not be used, just a placeholder
    wvec <- rep(1, length(y))
  }
  
  if (!is.finite(r2)) {
    r2 <- 0
  }
  
  # more flexible
  r2_use <- min(1, max(0, r2)) 
  #------------------------------------------------------------------------------------------#
  #new: use linear, not piece wise linear, compute edfu directly, no need to compute sc first#
  #------------------------------------------------------------------------------------------#
  user_edfus_scaled <- rep(0, 3)
  n <- length(y)
  edfu_m <- 2 * n^(1/7)
  if (family$family == 'gaussian') {
    user_edfus_scaled[1] <- (0.5 + r2_use / 2) * (edfu_m - 1)
    user_edfus_scaled[2] <- (0.5 + r2_use / 2) * edfu_m
    user_edfus_scaled[3] <- (0.75 + r2_use / 1.6) * (edfu_m + 6)
  } else if (family$family == 'binomial') {
    edfu_m <- 1.4 * n^(1/7)
    user_edfus_scaled[1] <- (0.5 + r2_use / 4) * (edfu_m - 1)
    user_edfus_scaled[2] <- (0.5 + r2_use / 4) * (edfu_m)
    user_edfus_scaled[3] <- (0.75 + r2_use / 1.6) * (edfu_m + 6)
  }
  
  sc <- user_edfus_scaled / user_edfus
  # attach useful information without changing the numeric result
  attr(user_edfus_scaled, "initial_r2") <- r2_use
  attr(user_edfus_scaled, "edfu_scaler") <- sc
  attr(user_edfus_scaled, "weights") <- wvec
  #attr(user_edfus_scaled, "pen") <- lam_vec
  return (user_edfus_scaled)
}

#--------------------------------------------------------------------------------
# computer either explain deviance or r2 used to scale the starting edfus
#--------------------------------------------------------------------------------
compute_dev <- function(y, etahat, family) {
  #no weight, if y is binomial then it is 0/1
  wt <- rep(1, length(y))
  muhat <- family$linkinv(etahat)
  
  if (family$family == "binomial") {
    dev <- sum(family$dev.resids(y, muhat, wt = wt))
    mu0 <- rep(mean(y), length(y))
    dev0 <- sum(family$dev.resids(y, mu0, wt = wt))
    
    if (!is.finite(dev0) || dev0 <= 0) {
      return(0)
    }
    
    return(1 - dev / dev0)
  } else {
    ssto <- sum((y - mean(y))^2)
    sse <- sum((y - muhat)^2)
    
    if (!is.finite(ssto) || ssto <= 0) {
      return(0)
    }
    
    return(1 - sse / ssto)
  }
}

#--------------------------------------------------------------------------------
#tune edfu after shape and variable selection
#new: add muhatkeep to compute the weight for IRLS in x'wx
#--------------------------------------------------------------------------------
tune_edfu <- function(active, active_shps, gmatkeep, varlst_s, np_ix, precomps, y,
                      types_active, wt.iter = FALSE, xnms_active,
                      a_vec = seq(.1, .4, by = .1),
                      edfu_grid = NULL,
                      tune = c("cv5", "cv10", "gcv"),
                      cv_seed = 123, class_cut = 0.5, muhatkeep,...) {
  
  tune <- match.arg(tune)
  
  if(wt.iter){
    family <- binomial(link = 'logit')
  } else {
    family <- gaussian(link = "identity")
  }
  
  #num_id <- which(types_active == "num")
  
  ## ---------------------------------------------------------
  ## keep numeric predictors
  ## ---------------------------------------------------------
  
  #kept <- types_active == "num"
  kept <- types_active %in% c("num", "ord") | (types_active == "tree" & active_shps == 1)
  
  if (!any(types_active == "num")) {
    stop("No numeric predictor is chosen. no edfus to tune")
  }
  
  active_tuned <- active[kept]
  active_shps <- active_shps[kept]
  # follow the order of active: when active is c(2,1),
  # precomps should follow this order instead of 1,2,3...
  precomps <- precomps[active_tuned]
  xnms_tuned <- xnms_active[kept]
  #new:
  types_active_kept <- types_active[kept]
  
  ## ------------------------------------------------------------------------------------------
  ## start with the edfu used by the kept num predictors
  ## follow the order in active, not follow the column order in xmat
  ## new: no need to record dd for tree, shp=1 or ord bc there is no penalty
  ## need to record their dmat/dv, which is zero, otherwise there will be a xtx+dv dimension error
  ## ------------------------------------------------------------------------------------------
  p <- length(precomps)
  edfu_start <- numeric(p)
  nkts <- numeric(p)
  dd_lst <- dv_lst <- vector("list", length = p)
  
  for(ix in seq_len(p)){
    shpi_ix <- active_shps[ix]
    shp_name <- as.character(active_shps[ix])
    pc_ix <- precomps[[ix]]
    
    dv_ix <- pc_ix$dv_lst_raw[[shp_name]]
    dv_lst[[ix]] <- dv_ix
    
    dd_ix <- pc_ix$dd_lst[[shp_name]]
    
    # Needed only when EDF is tuned
    if (types_active_kept[ix] == "num") {
      dd_lst[[ix]] <- dd_ix
    }
    
    edfu_vec <- pc_ix$edfu_vec
    edfu_ix <- edfu_vec[shp_name]
    
    # use this one so when the edfu strategy changes in
    # precompute_lists(), it still works
    nkts[ix] <- length(pc_ix$kts_lst[[shp_name]])
    edfu_start[ix] <- edfu_ix
  }
  
  ## ----------------------------------------------------------------------
  ## start with the flexible edfu and then tune it
  ## ----------------------------------------------------------------------
  # new: let user-defined edfu_grid be an option; assume it is provided as a vector
  if(!is.null(edfu_grid)){ #not in use
    if (!is.numeric(edfu_grid) || length(edfu_grid) < 2L || any(!is.finite(edfu_grid)) ||any(edfu_grid <= 0)) {
      stop("edfu_grid must be a numeric vector of at least two positive finite values.")
    }
    edfu_grid_mat <- matrix(0, nrow = length(edfu_grid), ncol = p)
    for(ix in seq_len(p)){
      edfu_grid_mat[, ix] <- edfu_grid
    }
    edfu_grid <- edfu_grid_mat
  } else {
    if (!is.numeric(a_vec) || length(a_vec) < 2L || any(!is.finite(a_vec)) ||any(a_vec <= 0)) {
      stop("a_vec must be a numeric vector of at least two positive finite values.")
    }
    # start_frac <- edfu_start / nkts
    # make sure the current EDF is included in the grid
    # a_vec <- sort(unique(c(a_vec, start_frac)))
    # edfu_grid <- outer(a_vec, nkts)
    
    edfu_grid <- outer(a_vec, nkts)
    # always include each predictor's starting edfu
    edfu_grid <- rbind(edfu_grid, edfu_start)
    rownames(edfu_grid) <- NULL
  }
  
  #new: cap every potential edfu at a minimum of 1.5
  edfu_grid[, types_active_kept == "num"] <- pmax(edfu_grid[, types_active_kept == "num"], 1.5)
  ## -------------------------------------------------------------------
  ## find lambda corresponding to each EDF value
  ## -------------------------------------------------------------------
  # follow the order in active, not xmat
  lams_mat <- matrix(0, nrow = nrow(edfu_grid), ncol = ncol(edfu_grid))
  #new: binomial, dd must be weighted by w^.5
  if(wt.iter){
    weights <- muhatkeep * (1 - muhatkeep)
  } else {
    weights <- NULL
  }
  
  for(ix in seq_len(p)){
    
    #new:
    if(types_active_kept[ix] != "num"){
      next
    }
    
    dd <- dd_lst[[ix]]
    dv <- dv_lst[[ix]]
    qv <- crossprod(dd)
    edfu_vec <- edfu_grid[, ix]
    
    log_lams <- sapply(edfu_vec, function(e) {
      tryCatch({
        log_lo <- -10
        log_hi <- 18
        
        f_lo <- try(.search_ps_sel(exp(log_lo), qv0 = qv, dv = dv, dd = dd, edfu = e, weights = weights), silent = TRUE)
        
        # if the lower endpoint fails, move it inward
        while (inherits(f_lo, "try-error") && log_lo < log_hi) {
          log_lo <- log_lo + 1
          f_lo <- try(.search_ps_sel(exp(log_lo), qv0 = qv, dv = dv, dd = dd, edfu = e, weights = weights), silent = TRUE)
        }
        
        f_hi <- try(.search_ps_sel(exp(log_hi), qv0 = qv, dv = dv, dd = dd, edfu = e, weights = weights), silent = TRUE)
        
        # if the upper endpoint fails, move it inward
        while (inherits(f_hi, "try-error") && log_hi > log_lo) {
          log_hi <- log_hi - 1
          f_hi <- try(.search_ps_sel(exp(log_hi), qv0 = qv, dv = dv, dd = dd, edfu = e, weights = weights), silent = TRUE)
        }
        
        # check that two valid endpoints were found
        if (inherits(f_lo, "try-error") || inherits(f_hi, "try-error")) {
          stop("could not find two valid penalty endpoints.")
        }
        
        if (f_lo <= 0) {
          # root is at or below the smallest usable log penalty
          log_lo
        } else if (f_hi >= 0) {
          # root is at or above the largest usable log penalty
          log_hi
        } else {
          
          # f_lo > 0 and f_hi < 0, so the root is bracketed
          log_root_res <- try(uniroot(f = function(log_pen,...) {.search_ps_sel(exp(log_pen), ...)},
                                      interval = c(log_lo, log_hi), qv0 = qv, dv = dv, dd = dd,
                                      edfu = e, weights = weights, extendInt = "no",
                                      tol = .Machine$double.eps^0.25),
                              silent = TRUE
          )
          
          if (!inherits(log_root_res, "try-error")) {
            log_root_res$root
          } else {
            
            # both endpoints are valid and bracket the root.
            # use the endpoint whose function value is closest to zero.
            if (abs(f_lo) <= abs(f_hi)) {
              log_lam_fallback <- log_lo
            } else {
              log_lam_fallback <- log_hi
            }
            
            #cat("unexpected try-error in uniroot\n")
            #cat("edfu:", e, "\n")
            #cat("usable log interval:", log_lo, log_hi, "\n")
            #cat("endpoint values:", f_lo, f_hi, "\n")
            
            log_lam_fallback
          }
        }
        
      }, error = function(err) {
        
        #cat("Failed to find log_lam\n")
        #cat("edfu:", e, "\n")
        #cat("Error:", conditionMessage(err), "\n")
        
        NA_real_
      })
    })
    
    lams <- exp(log_lams)
    lams_mat[, ix] <- lams
  }
  
  ## -----------------------------------------------------------------------
  ## starting values
  ## -----------------------------------------------------------------------
  start_id <- sapply(seq_len(p), function(ix) {which.min(abs(edfu_grid[, ix] - edfu_start[ix]))})
  current_id <- start_id
  # for each x_j, take out its current lambda
  current_lams <- sapply(seq_len(p), function(j) lams_mat[current_id[j], j])
  
  ## -----------------------------------------------------------------------
  ## helper: create penalty matrix from lambda vector
  ## -----------------------------------------------------------------------
  make_dv <- function(lam_vec) {
    dv_lst_tmp <- Map(function(D, lam) lam * D, dv_lst, lam_vec)
    Matrix::bdiag(dv_lst_tmp) |> as.matrix()
  }
  
  ## -----------------------------------------------------------------------
  ## helper: full-data model fit
  ## -----------------------------------------------------------------------
  fit_for_lams <- function(lam_vec) {
    dv <- make_dv(lam_vec)
    mod <- fit_model_engine(y, current_gmat = gmatkeep, current_dv = dv, family = family,
                            pnt = TRUE, np_ix = np_ix, varlst_s = varlst_s,
                            is_first_step = TRUE, compute_gcv = TRUE)
    return(mod)
  }
  
  ## -----------------------------------------------------------------------
  ## new: create folds if 5- or 10- fold CV is requested
  ## -----------------------------------------------------------------------
  if(tune == "cv5" || tune == "cv10"){
    n <- length(y)
    K <- if(tune == "cv5") 5L else 10L
    #set.seed(cv_seed)
    fold_id <- withr::with_seed(cv_seed, {
      if (wt.iter) {
        fid <- integer(n)
        id0 <- which(y == 0)
        id1 <- which(y == 1)
        fid[id0] <- sample(rep(1:K, length.out = length(id0)))
        fid[id1] <- sample(rep(1:K, length.out = length(id1)))
        fid
      } else {
        sample(rep(1:K, length.out = n))
      }
    })
  }
  
  ## -----------------------------------------------------------------------
  ## new: compute k-fold CV criterion for one lambda vector
  ## -----------------------------------------------------------------------
  cv_for_lams <- function(lam_vec) {
    
    dv <- make_dv(lam_vec)
    
    n <- length(y)
    mse_cv <- numeric(K)
    
    for(ff in seq_len(K)){
      test_id <- which(fold_id == ff)
      train_id <- which(fold_id != ff)
      y_train <- y[train_id]
      y_test <- y[test_id]
      gmat_train <- gmatkeep[train_id, , drop = FALSE]
      gmat_test <- gmatkeep[test_id, , drop = FALSE]
      
      ## fit model using training observations
      #new: dv should have zero's if type = 'ord' or type = 'tree' & shp = 1
      mod_train <- fit_model_engine(y_train, current_gmat = gmat_train, current_dv = dv,
                                    family = family, pnt = TRUE, np_ix = np_ix, varlst_s = varlst_s,
                                    is_first_step = TRUE, compute_gcv = FALSE)
      
      ## -------------------------------------------------------------
      ## prediction on test set
      ## -------------------------------------------------------------
      eta_test <- as.numeric(gmat_test %*% mod_train$bhat)
      
      if(!wt.iter){
        ## gaussian
        ## test-set MSE
        mse_cv[ff] <-  sum((y_test - eta_test)^2) / length(y_test) #mean((y_test - eta_test)^2)
      }else{
        ## binomial
        ## predicted probability
        prob_test <- family$linkinv(eta_test)
        ## predicted class
        pred_test <- as.integer(prob_test >= class_cut)
        ## misclassification rate
        # mse_cv[ff] <- mean(pred_test != y_test)
        ## new: try brier score
        mse_cv[ff] <- mean((y_test - prob_test)^2)
      }
    }
    ## mean test loss across the 5 folds
    rslt <- mean(mse_cv)
    return (rslt)
  }
  
  ## -----------------------------------------------------------------------
  ## new: general criterion
  ## GCV:
  ## CV5 gaussian: mean test MSE
  ## CV5 binomial: mean test brier score/misclassification rate
  ## -----------------------------------------------------------------------
  score_for_lams <- function(lam_vec) {
    if(tune == "gcv"){
      mod <- fit_for_lams(lam_vec)
      return(mod$gcv)
    }else{
      mse_cv <- cv_for_lams(lam_vec)
      return (mse_cv)
    }
  }
  
  ## -----------------------------------------------------------------------
  ## criterion at the original EDF values
  ## -----------------------------------------------------------------------
  score_start <- score_for_lams(lam_vec = current_lams)
  score_current <- score_start
  #cat('score before tuning: ', score_current, '\n')
  #cat('lams before tuning: ', current_lams, '\n')
  
  ## optional history
  max_pass <- 5L
  score_path <- vector("list", max_pass)
  id_path <- vector("list", max_pass)
  
  ## -----------------------------------------------------------------------
  ## one-by-one tuning
  ## -----------------------------------------------------------------------
  for(pass in seq_len(max_pass)){
    changed <- FALSE
    score_pass <- vector("list", p)
    
    for(j in rev(seq_len(p))){
      #new
      if(types_active_kept[j] != "num"){
        next
      }
      score_j <- rep(NA_real_, nrow(lams_mat))
      for(k in seq_len(nrow(lams_mat))){
        # current solution for all coordinates,
        # except coordinate j is replaced by candidate k
        lam_try <- current_lams
        lam_try[j] <- lams_mat[k, j]
        
        #if(k == current_id[j]){
        #  score_j[k] <- score_current
        #}else{
        score_j[k] <- score_for_lams(lam_try)
        #}
      }
      
      k_best <- which.min(score_j)
      score_best <- score_j[k_best]
      
      # if(score_current != 0){
      #   improved <- (score_current - score_best) / abs(score_current)
      # }else{
      #   improved <- score_current - score_best
      # }
      
      improved <- (score_current - score_best)
      tol <- 1e-8 * max(1, abs(score_current))
      
      #if(improved > 1e-2){
      if(improved > tol){
        changed <- TRUE
        current_id[j] <- k_best
        current_lams[j] <- lams_mat[k_best, j]
        score_current <- score_best
        #cat('score current: ', score_current, '\n')
        #cat('current lams: ', current_lams, '\n')
      }
      score_pass[[j]] <- score_j
    }
    
    score_path[[pass]] <- score_pass
    
    id_path[[pass]] <- current_id
    
    ## stop if a complete pass made no changes
    if(!changed){
      break
    }
  }
  
  ## -----------------------------------------------------------------------
  ## final fit
  ## refit the final model using ALL observations
  ## -----------------------------------------------------------------------
  mod_final <- fit_for_lams(lam_vec = current_lams)
  ps <- current_lams
  #a_selected <- a_vec[current_id]
  edfu_selected <- sapply(seq_len(p), function(j) {edfu_grid[current_id[j], j]})
  a_selected <- edfu_selected / nkts
  coefskeep <- mod_final$bhat
  etahatkeep <- mod_final$etahat
  names(a_selected) <- xnms_tuned
  names(edfu_selected) <- xnms_tuned
  names(edfu_start) <- xnms_tuned
  
  ## -----------------------------------------------------------------------
  ## results
  ## -----------------------------------------------------------------------
  rslt <- list(ps = ps, edfu_cv = edfu_selected, edfu_start = edfu_start, criterion = tune,
               cv_min = score_current, cv_start = score_start,  a_cv = a_selected,
               coefskeep = coefskeep, etahatkeep = as.numeric(etahatkeep),
               lams_mat = lams_mat, score_path = score_path,
               id_path = id_path, active_tuned = active_tuned)
  
  return(rslt)
}

#----------------------------------------------------------------------------------------------------------
#if any x is added to the active set, or if any shape changes for predictors
#already in the active set: update gmat, np_ix, varlst_s according to ddkeep_x, ddkeep_z, etc.
#----------------------------------------------------------------------------------------------------------
forward_cic <- function(xmat, y, cic_old, active_0 = integer(0), active_shps_0 = integer(0), nsim = 400, 
                        family = gaussian(link = 'identity'), oldetahat = NULL, parallel = TRUE, pnt = TRUE, 
                        a = 1, set4 = FALSE, precomps = NULL, ddkeep_x_all_0 = NULL, ddkeep_z_all_0 = NULL, 
                        dv_keep_all_0 = NULL, dv_keep_lst_0 = NULL, varlst_s_0 = NULL, varlst_0 = NULL, 
                        varlst_z_0 = NULL, oldcoefs = NULL, oldgmat = NULL, oldedf0 = NULL, oldnum_active_nom = 0)
{
  n <- nrow(xmat)
  p <- ncol(xmat)
  is_first_step <- length(active_0) == 0L
  inactive_0 <- setdiff(seq_len(p), active_0)
  xmat_inact <- xmat[, inactive_0, drop = FALSE]
  #only use inactive predictors
  capl <- ncol(xmat_inact)  
  etahatkeep <- oldetahat
  coefskeep <- oldcoefs
  gmatkeep <- oldgmat
  #new
  edf0keep <- oldedf0
  ktskeep <- mskeep <- NULL
  ddkeep_x <- ddkeep_z <- dmat_keep <- dv_keep <- NULL
  precomps_inact <- precomps[inactive_0]
  picked_x <- picked_shp <- NULL
  ps_keep <- 0
  #new
  num_active_nom <- oldnum_active_nom
  
  cic_new <- cic_old #the cic from the previous step; if no improvement, then just return the old cic
  cic_st <- cic_old #cic_st will be updated when going through predictors
  nxc <- sum(active_shps_0 %in% c(3, 4, 11, 12))
  #if(is_first_step) runall <- FALSE else runall <- TRUE
  runall <- TRUE
  for(ix in seq_len(capl)){
    #print (ix)
    x <- xmat_inact[, ix]
    pc <- precomps_inact[[ix]] 
    type_ix <- pc$type_ix
    shps <- pc$shps
    dmat_lst <- pc$dmat_lst
    ms_lst <- pc$ms_lst 
    kts_lst <- pc$kts_lst 
    dd_lst <- pc$dd_lst
    ps_vec <- pc$ps_vec
    #c(shps, dmat_lst, ms_lst, kts_lst, dd_lst, ps_vec) %<-% pc[c("shps", "dmat_lst", "ms_lst", "kts_lst", "dd_lst", "ps_vec")]
    if(is_first_step){
      gmat_lst <- pc$gmat_lst 
      edfu_un <- pc$edfu_un
      dv_lst <- pc$dv_lst
      np_ix_lst <- pc$np_ix_lst
      varlst_s_lst <- pc$varlst_s_lst
      #c(gmat_lst, edfu_un, dv_lst, np_ix_lst, varlst_s_lst) %<-% pc[c("gmat_lst", "edfu_un", "dv_lst", "np_ix_lst", "varlst_s_lst")]
    } else {
      comps_update <- update_precomp_lists(x = x, type_ix = type_ix, ddkeep_x_all = ddkeep_x_all_0, 
                                           ddkeep_z_all = ddkeep_z_all_0, 
                                           dv_keep_all = dv_keep_all_0, 
                                           varlst_s_keep = varlst_s_0, 
                                           pc = pc, pnt = pnt)
      gmat_lst <- comps_update$gmat_lst
      dv_lst <- comps_update$dv_lst 
      np_ix_lst <- comps_update$np_ix_lst
      varlst_s_lst <- comps_update$varlst_s_lst
      edfu_un <- NULL #check more
    }
    nobs <- length(shps)
    bhat_lst <- etahat_lst <- vector("list", length = nobs)
    llh_vec <- cic_lwrs <- numeric(nobs)
    #new: modify lower bound by the number of active xs so far
    num_active <- length(active_0)
    #minedf0 <- edf0keep #+ 1
    #new: from Mary's code
    minedf0 <- edf0keep # + .4
    #print (edf0keep)
    for(ps_shpi in seq_along(shps)){
      shpi <- shps[ps_shpi]
      #print (shpi)
      gmat <- gmat_lst[[ps_shpi]]
      dv <- dv_lst[[ps_shpi]]
      #dmat <- NULL
      np_ix <- np_ix_lst[[ps_shpi]]
      varlst_s <- varlst_s_lst[[ps_shpi]]
      ans <- fit_model_engine(y, current_gmat = gmat, current_dv = dv, 
                              family = family, pnt = pnt, np_ix = np_ix, 
                              varlst_s = varlst_s, 
                              is_first_step = is_first_step)
      bhat_lst[[ps_shpi]] <- ans$bhat
      etahat_lst[[ps_shpi]] <- ans$etahat
      llh <- calculate_llh(y, ans$etahat, family)
      llh_vec[ps_shpi] <- llh
      #cic_lwrs[ps_shpi] <- llh + 1 * log(n) * np_ix / n^(6/7)
      #cic_lwrs[ps_shpi] <- llh + 1 * log(n) * num_active / n^(6/7)
      cic_lwrs[ps_shpi] <- llh + log(n) * minedf0 / n^(6/7)
    }
    
    ps_to_sim <- shps_to_sim <- NULL
    #use relative tol here
    tol <- 1e-5
    improved_lwr <- (cic_st - cic_lwrs) / abs(cic_st)
    if(any(improved_lwr > tol)){
      ps_to_sim <- which(improved_lwr > tol)
      shps_to_sim <- shps[ps_to_sim]
      #new
      if((!19 %in% shps_to_sim) && is_first_step){
        edfu_un <- NULL
      }
      # if(family$family == "binomial"){
      #   edfu_un <- NULL
      # }
      edf_args <- list(gmat_lst = gmat_lst[ps_to_sim], dv_lst = dv_lst[ps_to_sim],
                       varlst_s_lst = varlst_s_lst[ps_to_sim], np_ix_lst = np_ix_lst[ps_to_sim], 
                       nsim = nsim, parallel = parallel, pnt = pnt, shps_to_sim = shps_to_sim, type_ix = type_ix,
                       edfu_un = edfu_un, runall = runall, family = family, is_first_step = is_first_step)
      df_rslt <- do.call(compute_edf, edf_args)
      dfmean_cbic <- df_rslt$dfmean_cbic
      
      capz <- 0 
      if(!is.null(ddkeep_z_all_0) && (NCOL(ddkeep_z_all_0) - nxc) > 0){ #not adjust for conv/conc
        #not include 1
        capz <- capz + NCOL(ddkeep_z_all_0) - nxc
      }
      #check more
      if(type_ix == "nom" && 1 %in% shps_to_sim){
        capz <- capz + ncol(dd_lst[[which(shps == 1)]])
      }
      if(type_ix == "tree" && 2 %in% shps_to_sim){
        capz <- capz + ncol(dd_lst[[which(shps == 2)]])
      }
      
      #test
      if(capz > 0){
        dfmean_cbic <- dfmean_cbic - capz + capz * n^(-1/7)
      }
      
      cic_vals <- llh_vec[ps_to_sim] + log(n) * dfmean_cbic / n^(6/7)
      #use relative tol here
      tol <- 1e-5
      improved <- (cic_st - cic_vals) / abs(cic_st)
      if(any(improved > tol)){
        best_local <- which.min(cic_vals)
        picked_shp <- shps_to_sim[best_local]
        pskeep <- match(picked_shp, shps)
        etahatkeep <- etahat_lst[[pskeep]] 
        coefskeep <- bhat_lst[[pskeep]]
        gmatkeep <- gmat_lst[[pskeep]]
        mskeep <- ms_lst[[pskeep]] 
        ktskeep <- kts_lst[[pskeep]] 
        dd <- dd_lst[[pskeep]] 
        #new: keep dmat; not accumulative
        dmat <- dmat_lst[[pskeep]]
        ps <- ps_vec[pskeep]
        
        if(picked_shp == 1 && type_ix == "nom"){
          ddkeep_x <- NULL
          ddkeep_z <- dd
          #new
          num_active_nom <- num_active_nom + 1
        } else if(picked_shp == 2 && type_ix == "tree"){#nom case in tree shape
          ddkeep_x <- NULL
          ddkeep_z <- dd
        } else {
          ddkeep_x <- dd
          ddkeep_z <- NULL
          if ((picked_shp > 10 && picked_shp < 13) || (picked_shp > 2 && picked_shp < 5)) {
            ddkeep_z <- cbind(x - mean(x))
          } 
        }
        
        if (pnt) {
          dmat_keep <- dmat
          #if we are checking nom, both dmat and dv_keep should be NULL
          dv_keep   <- if (!is.null(dmat)) crossprod(dmat) else NULL
          ps_keep   <- ps
        }
        
        cic_st <- cic_vals[best_local]
        #new:
        edf0keep <- dfmean_cbic[best_local]
        cic_new <- cic_st
        picked_x <- inactive_0[ix]
        #print (edf0keep)
        # state_commit <- commit_forward_state(ddkeep_x=ddkeep_x, ddkeep_z=ddkeep_z, 
        #                                      ddkeep_x_all_0 = ddkeep_x_all_0, ddkeep_z_all_0=ddkeep_z_all_0,
        #                                  varlst_0=varlst_0, varlst_s_0=varlst_s_0, 
        #                                  varlst_z_0=varlst_z_0, picked_x=picked_x, picked_shp=picked_shp, 
        #                                  dv_keep_all_0=dv_keep_all_0, dv_keep_lst_0=dv_keep_lst_0, 
        #                                  dmat_keep=dmat_keep, dv_keep=dv_keep, pnt=pnt)
        # 
        # ddkeep_x_all_0 <- state_commit$ddkeep_x_all_0
        # ddkeep_z_all_0 <- state_commit$ddkeep_z_all_0
        # varlst_0 <- state_commit$varlst_0
        # varlst_s_0 <- state_commit$varlst_s_0
        # varlst_z_0 <- state_commit$varlst_z_0
        # dv_keep_all_0 <- state_commit$dv_keep_all_0 
        # dv_keep_lst_0 <- state_commit$dv_keep_lst_0
      }
    }
  }
  rslt <- list(ddkeep_x = ddkeep_x, ddkeep_z = ddkeep_z, dmat_keep = dmat_keep, dv_keep = dv_keep, 
               ps_keep = ps_keep, mskeep = mskeep, ktskeep = ktskeep, edf0keep = edf0keep,
               etahatkeep = etahatkeep, coefskeep = coefskeep, gmatkeep = gmatkeep,
               cic_new = cic_new, picked_x = picked_x, picked_shp = picked_shp, num_active_nom = num_active_nom)
  return (rslt)
}

#-----------------------------------------------------------------------------------------------------
#idea is very similar to forward cic
#_0 is for current model; they will change if any new cic < current cic (will be updated in this function)
#the _all, _lst, varlst_ should all be re-set to the _0's whenever start checking a new predictor
#-----------------------------------------------------------------------------------------------------
backward_cic <- function(xmat, y, cic_old, active_0 = integer(0), active_shps_0 = integer(0),
                         check_set = integer(0), check_shps = integer(0), nsim = 400, 
                         family = gaussian(link = "identity"), 
                         oldetahat = NULL, parallel = FALSE, pnt = FALSE, a = 1, precomps = NULL, 
                         ddkeep_x_all_0 = NULL, ddkeep_z_all_0 = NULL, dv_keep_lst_0 = NULL, 
                         varlst_s_0 = NULL, varlst_0 = NULL, varlst_z_0 = NULL, dv_keep_all_0 = NULL, 
                         xnms_with_id = NULL, mod_tried = NULL, current_shps = NULL, 
                         oldcoefs = NULL, oldgmat = NULL, oldedf0 = NULL, oldnum_active_nom = 0)
{
  has_changed <- FALSE
  #norem <- TRUE #new
  n <- nrow(xmat)
  p <- ncol(xmat)
  inactive_0 <- setdiff(seq_len(p), active_0)
  xmat_check <- xmat[, check_set, drop = FALSE]
  #only use predictors in active set except for the most recent one
  capl <- NCOL(xmat_check)  
  etahatkeep <- oldetahat
  ktskeep <- mskeep <- NULL
  coefskeep <- oldcoefs
  gmatkeep <- oldgmat
  ddkeep_x <- ddkeep_z <- dmat_keep <- NULL
  precomps_check <- precomps[check_set]
  picked_x <- picked_shp <- removed_x <- NULL
  ps_keep <- 0
  #new:
  edf0keep <- oldedf0
  #new
  num_active_nom <- oldnum_active_nom
  
  cic_new <- cic_old #the cic from the previous step; if no improvement, then just return the old cic
  cic_st <- cic_old #cic_st will be updated when going through predictors
  runall <- TRUE
  for(ix in seq_len(capl)){
    #make a copy of what the currect model has: _0; if any(cic_vals < cic_st): current model will be updated, and those _0's will be modified
    ddkeep_x_all <- ddkeep_x_all_0
    ddkeep_z_all <- ddkeep_z_all_0
    dv_keep_lst <- dv_keep_lst_0 
    varlst_s <- varlst_s_0 #track 19 #will update as ix goes to capl
    varlst  <- varlst_0 #track num's other shapes
    varlst_z  <- varlst_z_0 #track z, nom, tree 2, ord 3,4 and num 11, 12
    active <- active_0
    inactive <- inactive_0
    active_shps <- active_shps_0
    
    x <- xmat_check[, ix]
    #the predictor we are now checking
    xid <- check_set[ix]
    pc <- precomps_check[[ix]] 
    type_ix <- pc$type_ix
    #new: add a flat shape and everything it needs
    pc <- update_pc_with_flat(pc)
    shps <- pc$shps
    dmat_lst <- pc$dmat_lst
    ms_lst <- pc$ms_lst 
    kts_lst <- pc$kts_lst 
    dd_lst <- pc$dd_lst
    ps_vec <- pc$ps_vec
    #c(shps, dmat_lst, ms_lst, kts_lst, dd_lst, ps_vec) %<-% pc[c("shps", "dmat_lst", "ms_lst", "kts_lst", "dd_lst", "ps_vec")]
    #-----------------------------------
    #update shps to be those not tried
    #-----------------------------------
    shp_chosen <- check_shps[ix]
    #current_shp_ps <- which(shps == shp_chosen)
    kept_shps <- shps != shp_chosen
    shps <- shps[kept_shps]
    pos <- match(xid, active)
    if (is.na(pos)) stop("Backward step error: xid not found in active set.")
    #--------------------------------------------------------------------------------------------------------
    #update active, inactive, active_shps, in this current model: remove xid and its shape
    #add it at the end of the tracking vectors if the new picked_shp != 17 or 0
    #--------------------------------------------------------------------------------------------------------
    active <- active[-pos]
    inactive <- c(inactive, xid)
    active_shps <- active_shps[-pos]
    nxc <- sum(active_shps %in% c(3, 4, 11, 12))
    #-------------------------------------------------------------------------------------------------------------
    #update ddkeep_x_all, ddkeep_z_all
    #varlst_s_0, dv_keep_all
    #they will remove the chosen predictor's dd_x, dd_z, vals in varlst_s_0 and reduce dv_keep_all
    #-------------------------------------------------------------------------------------------------------------
    cand_state <- remove_predictor_from_state(xid = xid, type_ix = type_ix, shp_chosen = shp_chosen,
                                              ddkeep_x_all = ddkeep_x_all, ddkeep_z_all = ddkeep_z_all, dv_keep_lst = dv_keep_lst,
                                              varlst_s = varlst_s, varlst = varlst, varlst_z = varlst_z, pnt = pnt)
    ddkeep_x_all <- cand_state$ddkeep_x_all
    ddkeep_z_all <- cand_state$ddkeep_z_all
    dv_keep_all  <- cand_state$dv_keep_all
    dv_keep_lst <- cand_state$dv_keep_lst
    varlst_s     <- cand_state$varlst_s
    varlst       <- cand_state$varlst
    varlst_z     <- cand_state$varlst_z
    #---------------------------------------
    #update lists used for computing
    #---------------------------------------
    comps_update <- update_precomp_lists(x = x, type_ix = type_ix, ddkeep_x_all = ddkeep_x_all,  
                                         ddkeep_z_all = ddkeep_z_all, 
                                         dv_keep_all = dv_keep_all, 
                                         varlst_s_keep = varlst_s, pc = pc, 
                                         pnt = pnt, skip_shape = shp_chosen)
    gmat_lst <- comps_update$gmat_lst
    dv_lst <- comps_update$dv_lst 
    np_ix_lst <- comps_update$np_ix_lst
    varlst_s_lst <- comps_update$varlst_s_lst
    keep_id <- comps_update$keep_id
    
    #remove the chosen one, to be used when new cic is smaller
    dmat_lst <- dmat_lst[keep_id]
    ms_lst  <- ms_lst[keep_id]
    kts_lst <- kts_lst[keep_id]
    dd_lst <- dd_lst[keep_id]
    ps_vec <- ps_vec[keep_id]
    #edfu_un <- NULL #check more
    
    nobs <- length(shps)
    bhat_lst <- etahat_lst <- vector("list", length = nobs)
    llh_vec <- cic_lwrs <- numeric(nobs)
    names(bhat_lst) <- names(etahat_lst) <- as.character(shps)
    names(llh_vec) <- names(cic_lwrs) <- as.character(shps)
    #new: update cic lower bounds by the number of active xs
    num_active <- length(active)
    #new: from Mary's code
    #minedf0 = dim(zmod)[2]*(1-n^(-1/7))+sum(curmod[1:nc]>0)+1	
    capz <- 0 
    if(!is.null(ddkeep_z_all) && (NCOL(ddkeep_z_all) - nxc) > 0){
      #not include 1
      capz <- capz + NCOL(ddkeep_z_all) - nxc
    }
    #test: 1 is for the intercept
    minedf0 = capz*(1-n^(-1/7)) + (num_active-num_active_nom) + 1
    for(ps_shpi in seq_along(shps)){
      gmat <- gmat_lst[[ps_shpi]]
      dv <- dv_lst[[ps_shpi]]
      np_ix <- np_ix_lst[[ps_shpi]]
      #new name to avoid overwriting
      varlst_s_use <- varlst_s_lst[[ps_shpi]]
      ans <- fit_model_engine(y, current_gmat = gmat, current_dv = dv, family = family, pnt = pnt, np_ix = np_ix, 
                              varlst_s = varlst_s_use, is_first_step = FALSE)
      bhat_lst[[ps_shpi]] <- ans$bhat
      etahat_lst[[ps_shpi]] <- ans$etahat
      llh <- calculate_llh(y, ans$etahat, family)
      llh_vec[ps_shpi] <- llh
      #cic_lwrs[ps_shpi] <- llh + 1 * log(n) * np_ix / n^(6/7)
      #cic_lwrs[ps_shpi] <- llh + 1 * log(n) * num_active / n^(6/7)
      #new: from Mary's code
      cic_lwrs[ps_shpi] <- llh + log(n) * minedf0 / n^(6/7)
    }
    
    ps_to_sim <- shps_to_sim <- NULL
    #use relative tolerance here
    tol <- 1e-5
    improved_lwr <- (cic_st - cic_lwrs) / abs(cic_st)
    if(any(improved_lwr > tol)){
      ps_to_sim <- which(improved_lwr > tol)
      #ps_to_sim <- which(cic_lwrs < cic_st) #ok to use this here because it will only be run after the threshold using tol is passed
      shps_to_sim <- shps[ps_to_sim]
      edf_args <- list(gmat_lst = gmat_lst[ps_to_sim], dv_lst = dv_lst[ps_to_sim],
                       varlst_s_lst = varlst_s_lst[ps_to_sim], np_ix_lst = np_ix_lst[ps_to_sim], 
                       nsim = nsim, parallel = parallel, pnt = pnt, shps_to_sim = shps_to_sim, type_ix = type_ix,
                       edfu_un = NULL, runall = runall, family = family, is_first_step = FALSE)
      df_rslt <- do.call(compute_edf, edf_args)
      dfmean_cbic <- df_rslt$dfmean_cbic
      
      capz <- 0 
      if(!is.null(ddkeep_z_all) && (NCOL(ddkeep_z_all) - nxc) > 0){
        #not include 1
        capz <- capz + NCOL(ddkeep_z_all) - nxc
      }
      
      #check more: won't happen for z in the backward step; the only shape to check is 0
      # if(type_ix == "nom" && all(shps_to_sim == 1)){
      #   capz <- NCOL(ddkeep_z_all) + NCOL(dd_lst[[1]]) 
      # }
      
      if(type_ix == "tree" && 2 %in% shps_to_sim){
        capz <- capz + ncol(dd_lst[[which(shps == 2)]])
      }
      
      #test
      if(capz > 0){
        dfmean_cbic <- dfmean_cbic - capz + capz * n^(-1/7)
      }
      
      cic_vals <- llh_vec[ps_to_sim] + log(n) * dfmean_cbic / n^(6/7)
      
      #print (dim(ddkeep_x_all_0))
      #use relative tolerance here
      tol <- 1e-5
      improved <- (cic_st - cic_vals) / abs(cic_st)
      if(any(improved > tol)){
        has_changed <- TRUE
        #cat('better model: ', ix, '\n')
        best_local <- which.min(cic_vals)
        picked_shp <- shps_to_sim[best_local]
        #picked_shp_to_cat <- if(picked_shp == 17) 0 else picked_shp
        #cat(xnms_with_id[xid], 'changed!', 'new shp = ', if(picked_shp == 17) 0 else picked_shp , '\n')
        #pskeep <- which(shps == picked_shp)
        pskeep <- match(picked_shp, shps)
        etahatkeep <- etahat_lst[[pskeep]] 
        coefskeep <- bhat_lst[[pskeep]]
        gmatkeep <- gmat_lst[[pskeep]]
        mskeep <- ms_lst[[pskeep]] 
        ktskeep <- kts_lst[[pskeep]] 
        dd <- dd_lst[[pskeep]] 
        #new: keep dmat; not accumulative
        dmat <- dmat_lst[[pskeep]]
        ps <- ps_vec[pskeep]
        
        state_commit <- commit_backward_state(picked_shp = picked_shp, xid = xid,
                                              x = x, type_ix = type_ix, dd = dd, dmat = dmat, ddkeep_x_all = ddkeep_x_all, 
                                              ddkeep_z_all = ddkeep_z_all, dv_keep_lst = dv_keep_lst,
                                              varlst_s = varlst_s, varlst = varlst, varlst_z = varlst_z, pnt = pnt)
        
        ddkeep_x <- state_commit$ddkeep_x
        ddkeep_z <- state_commit$ddkeep_z
        ddkeep_x_all_0 <- state_commit$ddkeep_x_all_0
        ddkeep_z_all_0 <- state_commit$ddkeep_z_all_0
        varlst_s_0 <- state_commit$varlst_s_0
        varlst_0 <- state_commit$varlst_0
        varlst_z_0 <- state_commit$varlst_z_0
        dv_keep_lst_0 <- state_commit$dv_keep_lst_0
        dv_keep_all_0 <-  state_commit$dv_keep_all_0
        
        if(pnt){
          dmat_keep <- dmat
          ps_keep <- ps
        }
        
        cic_st <- cic_vals[best_local]
        #new
        edf0keep <- dfmean_cbic[best_local]
        
        cic_new <- cic_st
        #picked_x <- inactive[ix]
        #pos <- match(xid, active_0)
        #if (is.na(pos)) stop("Internal error: xid not found in active set.")
        if (picked_shp != 0 && picked_shp != 17) {
          picked_x <- xid
          removed_x <- NULL
          action <- "change"
          #cat(action, '\n')
          #active_shps[pos] <- picked_shp
          active_0 <- c(active, xid)
          inactive_0 <- inactive[inactive != xid]
          active_shps_0 <- c(active_shps, picked_shp)
          #norem <- TRUE
        } else {
          removed_x <- xid
          picked_x <- NULL
          action <- "remove"
          #cat(action, '\n')
          inactive_0 <- inactive
          #keep <- active != xid
          active_0 <- active
          active_shps_0 <- active_shps
          #new
          if(type_ix == "nom"){
            num_active_nom <- num_active_nom - 1
          }
          #norem <- FALSE
        }
        
        current_shps[xid] <- ifelse(picked_shp == 17, 0, picked_shp)
        mod_tried <- append_mod_tried(mod_tried, current_shps, cic = cic_new, ps_keep = ps_keep,
                                      xid = xid, action = action, xnms_with_id = xnms_with_id)
      }
    }
  }
  
  rslt <- list(ddkeep_x = ddkeep_x, ddkeep_z = ddkeep_z, dmat_keep = dmat_keep, ps_keep = ps_keep, mskeep = mskeep, ktskeep = ktskeep,
               etahatkeep = etahatkeep, coefskeep = coefskeep, gmatkeep = gmatkeep,
               cic_new = cic_new, removed_x = removed_x,
               picked_x = picked_x, picked_shp = picked_shp, mod_tried = mod_tried, active_0 = active_0, 
               active_shps_0 = active_shps_0,
               inactive_0 = inactive_0, current_shps = current_shps, 
               ddkeep_x_all_0 = ddkeep_x_all_0,  ddkeep_z_all_0 = ddkeep_z_all_0,
               dv_keep_all_0 = dv_keep_all_0, dv_keep_lst_0 = dv_keep_lst_0, 
               varlst_z_0 = varlst_z_0, varlst_s_0 = varlst_s_0, varlst_0 = varlst_0, 
               has_changed = has_changed, edf0keep = edf0keep, num_active_nom = num_active_nom)
  return (rslt)
}

#--------------------------------------------------------------------
remove_predictor_from_state <- function(xid, type_ix, shp_chosen, ddkeep_x_all, ddkeep_z_all, dv_keep_lst,
                                        varlst_s, varlst, varlst_z, pnt = FALSE){
  
  
  has_z <- (type_ix == "nom" || (type_ix == "tree" && shp_chosen == 2) ||(type_ix == "ord" && shp_chosen %in% c(3, 4)) 
            ||(type_ix == "num" && shp_chosen %in% c(11, 12)))
  has_x <- (type_ix %in% c("num", "ord") || (type_ix == "tree" && shp_chosen == 1))
  
  if (has_z) {
    keep_z <- varlst_z != xid
    if (!is.null(ddkeep_z_all)) {
      ddkeep_z_all <- ddkeep_z_all[, keep_z, drop = FALSE]
    }
    varlst_z <- varlst_z[keep_z]
  } 
  
  if(has_x) {
    keep_x <- varlst != xid
    if (!is.null(ddkeep_x_all)) {
      ddkeep_x_all <- ddkeep_x_all[, keep_x, drop = FALSE]
    }
    varlst_s <- varlst_s[keep_x]
    varlst <- varlst[keep_x]
    
    if(pnt){
      dv_keep_lst[xid] <- list(NULL)
    }
  }
  
  dv_keep_all <- NULL
  if(pnt){
    u_keep_x <- unique(varlst)
    blocks <- dv_keep_lst[u_keep_x]
    # check for NULL blocks
    if (any(vapply(blocks, is.null, logical(1)))) {
      warning(
        "NULL dv detected in dv_keep_lst during backward step",
        "check backward update."
      )
    }
    dv_keep_all <- as.matrix(Matrix::bdiag(blocks))
  }
  
  list(
    ddkeep_x_all = ddkeep_x_all,
    ddkeep_z_all = ddkeep_z_all,
    dv_keep_lst  = dv_keep_lst,
    dv_keep_all  = dv_keep_all,
    varlst_s     = varlst_s,
    varlst       = varlst,
    varlst_z     = varlst_z
  )
}

#--------------------------------------------------------------------
commit_backward_state <- function(picked_shp, xid, x, type_ix, dd, dmat,
                                  ddkeep_x_all, ddkeep_z_all, dv_keep_lst,
                                  varlst_s, varlst, varlst_z, pnt = FALSE) {
  ddkeep_x <- NULL
  ddkeep_z <- NULL
  dv_keep_all_0 <- NULL
  
  ddkeep_x_all_0 <- ddkeep_x_all
  ddkeep_z_all_0 <- ddkeep_z_all
  varlst_s_0 <- varlst_s
  varlst_0 <- varlst
  varlst_z_0 <- varlst_z
  dv_keep_lst_0 <- dv_keep_lst
  
  removed <- (picked_shp == 0 || picked_shp == 17) 
  #the first case will never happen: kept for logic completeness
  #one could have both z and x: conv/conc
  has_z <- ((picked_shp == 1 && type_ix == "nom") || (picked_shp == 2 && type_ix == "tree") ||(picked_shp %in% c(3, 4) && type_ix == "ord") ||
              (picked_shp %in% c(11, 12) && type_ix == "num") )
  
  has_x <- ((picked_shp != 17 && type_ix == "num") ||(picked_shp != 0  && type_ix == "ord") ||
              (picked_shp == 1  && type_ix == "tree"))
  
  if (removed) {
    # predictor removed entirely
    ddkeep_x <- NULL
    ddkeep_z <- NULL
    
    if(pnt){
      dv_keep_lst_0[xid] <- list(NULL)
      u_keep_x <- unique(varlst_0)
      blocks <- dv_keep_lst_0[u_keep_x]
      # check for NULL blocks
      if (any(vapply(blocks, is.null, logical(1)))) {
        warning("NULL dv detected in dv_keep_lst during backward step. check backward update.")
      }
      dv_keep_all_0 <- as.matrix(Matrix::bdiag(blocks))
    }
    
  } 
  
  if (has_z) {
    # nominal predictor kept as z-part
    if(type_ix == "nom" || type_ix == "tree"){
      ddkeep_z <- dd
      ddkeep_z_all_0 <- cbind(ddkeep_z_all, ddkeep_z)
      varlst_z_0 <- c(varlst_z, rep(xid, ncol(ddkeep_z)))
      #new: if type_ix is 'tree' and the new shape is 2, we need to remove the zero dmat and create dv_keep_all_0 again
      if (type_ix == "tree") {
        
        if(pnt){
          dv_keep_lst_0[xid] <- list(NULL)
          u_keep_x <- unique(varlst_0)
          blocks <- dv_keep_lst_0[u_keep_x]
          # check for NULL blocks
          if (any(vapply(blocks, is.null, logical(1)))) {
            warning("NULL dv detected in dv_keep_lst during backward step. check backward update.")
          }
          dv_keep_all_0 <- as.matrix(Matrix::bdiag(blocks))
        }
        
      }
    }
    
    if(type_ix == "ord" || type_ix == "num"){
      # some numeric shapes also contribute a z-part
      ddkeep_z <- matrix(x - mean(x), ncol = 1)
      ddkeep_z_all_0 <- cbind(ddkeep_z_all, ddkeep_z)
      varlst_z_0 <- c(varlst_z, xid)
    }
  }
  
  if (has_x) {
    # predictor kept as x-part
    ddkeep_x <- dd
    ddkeep_x_all_0 <- cbind(ddkeep_x_all, ddkeep_x)
    varlst_0 <- c(varlst, rep(xid, ncol(ddkeep_x)))
    if (picked_shp == 19) {
      varlst_s_0 <- c(varlst_s, rep(1, ncol(ddkeep_x)))
    } else {
      varlst_s_0 <- c(varlst_s, rep(0, ncol(ddkeep_x)))
    }
    
    #check
    if(pnt){
      #if we have nom, then dmat is NULL, crossprod will fail
      #dv_keep_lst_0[[xid]] <- if (!is.null(dmat)) crossprod(dmat) else NULL
      dv_keep_lst_0[xid] <- list(if (!is.null(dmat)) crossprod(dmat) else NULL)
      
      u_keep_x <- unique(varlst_0)
      blocks <- dv_keep_lst_0[u_keep_x]
      # check for NULL blocks
      if (any(vapply(blocks, is.null, logical(1)))) {
        warning("NULL dv detected in dv_keep_lst during backward step. check backward update.")
      }
      dv_keep_all_0 <- as.matrix(Matrix::bdiag(blocks))
    }
  }
  
  list(
    ddkeep_x = ddkeep_x,
    ddkeep_z = ddkeep_z,
    ddkeep_x_all_0 = ddkeep_x_all_0,
    ddkeep_z_all_0 = ddkeep_z_all_0,
    dv_keep_lst_0 = dv_keep_lst_0,
    dv_keep_all_0 = dv_keep_all_0, 
    varlst_s_0 = varlst_s_0,
    varlst_0 = varlst_0,
    varlst_z_0 = varlst_z_0
  )
}

#-----------------------------------------------------------------------------------------------
#prediction are not in use now
#add confidence later; just use vcov.cgam
#check the predict function in classo_use.R for how to add c.i. bands for each comp
#-----------------------------------------------------------------------------------------------
predict.shapeselect <- function(object, newdata, interval = c("none"), 
                                type = c("response", "link"),...)
{
  interval <- match.arg(interval)
  type <- match.arg(type)
  
  if (!inherits(object, "shapeselect")) {
    warning("calling predict.shapeselect(<fake-shapeselect-object>) ...")
  }
  if (missing(newdata) || is.null(newdata)) {
    etahat <- object$etahat
    muhat <- object$muhat
    ans <- list(fit = muhat, etahat = etahat)
    return (ans)
  }
  if (!is.data.frame(newdata)) {
    stop ("newdata must be a data frame!")
  }
  if(is.null(object$active_shps)){
    stop("no predictor is chosen. don't use the predict function!")
  }
  
  active <- object$active
  active_shps <- object$active_shps
  xnms <- object$xnms
  xnms_kept <- xnms[active]
  
  if (!all(xnms_kept %in% names(newdata))) {
    stop("newdata must include all selected predictors!")
  }
  
  newdata <- newdata[, xnms_kept, drop = FALSE]
  
  if(any(active_shps > 0)){
    family <- object$family
    coefskeep <- object$coefskeep #with one
    pnt <- object$pnt
    #---------------------------
    #only focus on the chosen xs
    #---------------------------
    sdx <- object$sdx[active]
    mx <- object$mx[active]
    pcs <- object$precomps[active]
    xmat <- object$xmat[, active, drop = FALSE] #chosen ones; processed
    types <- object$types[active]
    
    varlst <- object$varlst
    varlst_z <- object$varlst_z
    p <- ncol(xmat)
    n <- nrow(xmat)
    
    #idx <- types == "num" & vapply(newdata, is.numeric, logical(1))
    #newdata[idx] <- lapply(which(idx), function(i) {(newdata[[i]] - mx[i]) / sdx[i]})
    for(i in 1:p){
      if(types[i] == "num"){
        newdata[,i] <- (newdata[,i] - mx[i]) / sdx[i]
      }
      #test more
      if(inherits(xmat[,i], "factor") && types[i] == "tree"){
        xmat[,i] <- as.numeric(as.character(xmat[,i]))
        #newdata[,i] <- as.numeric(as.character(newdata[,i]))
      }
      if(inherits(newdata[,i], "factor") && types[i] == "tree"){
        #xmat[,i] <- as.numeric(as.character(xmat[,i]))
        newdata[,i] <- as.numeric(as.character(newdata[,i]))
      }
      if(inherits(xmat[,i], "factor") && types[i] == "ord"){
        xmat[,i] <- as.numeric(as.character(xmat[,i]))
        #newdata[,i] <- as.numeric(as.character(newdata[,i]))
      }
      if(inherits(newdata[,i], "factor") && types[i] == "ord"){
        #xmat[,i] <- as.numeric(as.character(xmat[,i]))
        newdata[,i] <- as.numeric(as.character(newdata[,i]))
      }
    }
    
    #get etacomps for newdata first -> accumulative sum is etahat
    #easier to check
    dd_lst <- vector("list", length = p)
    names(dd_lst) <- xnms_kept
    new_ddkeep_z_all <- new_ddkeep_x_all <- NULL
    #--------------------------------------
    #create basis functions for newdata
    #only focus on the chosen ones
    #--------------------------------------
    if(!is.null(varlst_z)){
      uvar_z <- unique(varlst_z)
      for(zid in uvar_z){
        pos <- match(zid, active)
        pc <- pcs[[pos]]
        type_ix <- pc$type_ix
        x <- newdata[, pos]
        shpi <- active_shps[pos]
        #nom shape == 1
        #tree shape == 2
        #ord/num shape = 3,4,11,12
        oldx <- xmat[, pos]
        if(type_ix %in% c("num", "ord") && shpi %in% c(3, 4, 11, 12)){
          dd <- cbind(x - mean(oldx))
          #dd <- cbind(x) #already center and scaled -> should be changed bc min-max is a new option
        } else if (type_ix == "nom" || (type_ix == "tree" && shpi == 2)){
          dd <- create_newbasis(x, oldx, pc, shpi)
        }
        new_ddkeep_z_all <- cbind(new_ddkeep_z_all, dd)
        dd_lst[[pos]] <- dd
      }
    }
    
    if(!is.null(varlst)){
      uvar_x <- unique(varlst)
      for(xid in uvar_x){
        pos <- match(xid, active)
        oldx <- xmat[, pos]
        pc <- pcs[[pos]]
        x <- newdata[, pos]
        shpi <- active_shps[pos]
        dd <- create_newbasis(x, oldx, pc, shpi)
        new_ddkeep_x_all <- cbind(new_ddkeep_x_all, dd)
        dd_lst[[pos]] <- dd
      }
    }
    
    #etacomps in this part follow the order of active set; not x1 to xp
    etacomps <- assemble_shapechoice(active, active_shps, xmat = newdata, types, coefskeep, varlst_z,
                                     varlst, ddkeep_z_all = new_ddkeep_z_all, 
                                     ddkeep_x_all = new_ddkeep_x_all, use_predict = TRUE)
    #order it back to be the same as what the object has: x1 to xp
    #etacomps <- etacomps[order(active),]
    fit <- colSums(etacomps)
    fit <- fit + coefskeep[1] #add the intercept back
    fit <- switch(type, link = fit, response = object$family$linkinv(fit))
    
    #easier for visualization
    etacomps <- t(etacomps)
    etacomps <- as.data.frame(etacomps)
    colnames(etacomps) <- xnms_kept
    
    ans <- list(fit = fit, etacomps = etacomps, dd_lst = dd_lst,
                new_ddkeep_z_all = new_ddkeep_z_all, new_ddkeep_x_all = new_ddkeep_x_all)
    return (ans)
  }
}

#-----------------------------------------------------------------------------------------------
create_newbasis <- function(x, oldx, pc, shpi){
  n <- length(oldx)
  type_ix <- pc$type_ix
  shps <- pc$shps
  shp_pos <- match(shpi, shps)
  ms <- pc$ms_lst[[shp_pos]]
  dd <- NULL
  
  #######################
  #local helper function#
  #######################
  my_line = function(xp = NULL, y, x, end, start) {
    slope = NULL
    intercept = NULL
    yp = NULL
    slope = (y[end] - y[start]) / (x[end] - x[start])
    intercept = y[end] - slope * x[end]
    yp = intercept + slope * xp
    ans = new.env()
    ans$slope = slope
    ans$intercept = intercept
    ans$yp = yp
    ans
  }

  if(type_ix == "nom") {
    if (!inherits(oldx, "factor")) {
      stop("oldx must be a factor for nominal predictors.")
    }
    
    x <- factor(as.character(x), levels = levels(oldx))
    
    if (any(is.na(x) & !is.na(as.character(x)))) {
      stop("All factor levels in newdata must be among factor levels in the fit!")
    }
    
    # x_chr <- as.character(x)
    # allowed <- levels(oldx)
    # if (!all(stats::na.omit(unique(x_chr)) %in% allowed)) {
    #   stop("All factor levels in newdata must be among factor levels in the fit!")
    # }
    # x <- factor(x_chr, levels = allowed)
    
    # if(!inherits(x, "factor")){
    #   x <- as.factor(x)
    # }
    
    lvli <- levels(x) #new x
    #xui <- unique(x)
    ztbi <- levels(oldx) #old x
    
    #dd <- NULL
    if (!all(lvli %in% ztbi)) {
      stop("All factor levels in newdata must be among factor levels in the fit!")
    } else {
      #baseline: ztbi[1]
      klvls <- length(ztbi)
      rn <- length(x)
      #remove baseline?
      dd <- matrix(0, nrow = rn, ncol = (klvls - 1))
      for (i1 in 1:rn) {
        if (x[i1] != ztbi[1]) {
          id_col <- which(ztbi %in% x[i1]) - 1
          #print (id_col)
          dd[i1, id_col] <- 1
          dd[i1, id_col] <- dd[i1, id_col] - ms[id_col] #ms[i1, id_col]
        }
      }
      for(ic in 1:ncol(dd)){
        ps_base <- which(dd[, ic] == 0)
        dd[ps_base, ic] <- dd[ps_base, ic] - ms[ic] #ms[ps_base, ic]
      }
    }
  } else if(type_ix == "ord") {
    #shps <- pc$shps
    #shp_pos <- match(shpi, shps)
    #ms <- pc$ms_lst[[shp_pos]]
    dd <- pred_del(oldx, shpi, x, ms)  
    dd <- t(dd)
  } else if(type_ix == "tree") {
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
    #seems not good?
    #dd <- matrix(0, nrow = length(nx), ncol = length(x))
    ##if new x only has placebo, then just center each row of dd by ms
    # for(i in 1:nrow(dd)){
    #   xi <- nx[i]
    #   dd[i, ] <- dd[i, ] - ms[which(nx %in% xi)]
    # }
    # if(any(!is_pl)){
    #   for(i in 1:nrow(dd)){
    #     xi <- nx[i]
    #     if(any(x %in% xi)){
    #       dd[i, which(x %in% xi)] <- 1 - ms[which(nx %in% xi)]
    #     }
    #   }
    # }
    
    #test
    dd <- tree.fun(x, pl = pl)
    #if(any(!is_pl)){
    for(i in 1:nrow(dd)){
      xi <- nx[i]
      is_xi <- x %in% xi
      not_xi <- !is_xi
      ms_xi <- ms[which(nx %in% xi)]
      if(any(is_xi)){
        dd[i, is_xi] <- 1 - ms_xi
      }
      if(any(not_xi)){
        dd[i, not_xi] <- 0 - ms_xi
      }
    }
    #}
    
    dd <- t(dd)
  } else if(type_ix == "num") {
    kts <- pc$kts_lst[[shp_pos]]
    #suppress the boundary knots warning
    dd_ans <- suppressWarnings(makedelta(x, sh = ifelse(shpi != 19, shpi, shpi - 2), 
                                         knots = kts, space = "E", suppre = TRUE, interp = TRUE))
    dd <- dd_ans$amat
    if(shpi > 10 && shpi < 13 || shpi > 2 && shpi < 5){
      ms <- t(ms)
    }
    #need to make basis orthogonal to one or [1, x]
    if (shpi > 10 && shpi < 13 || shpi > 2 && shpi < 5) {
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
    dd <- t(dd)
  } else {
    stop('type of x is not allowed!')
  }
  return (dd)
}

# ------------------------------------------------------------------------------------------------------------------------
#get etacomps: each smooth component or nominal component
#when use_predict = F, xconv/xconc are already centered
#when use_predict = T, xconv/xconc are not yet centered and should be centered in predict.shapeselect
# ------------------------------------------------------------------------------------------------------------------------
assemble_shapechoice <- function(active, active_shps, xmat, types, 
                                 coefskeep, varlst_z, varlst, ddkeep_z_all, 
                                 ddkeep_x_all, use_predict = FALSE){
  n <- nrow(xmat)
  p <- ncol(xmat)
  nz <- length(varlst_z)
  coefskeep_no_one <- coefskeep[-1]
  if (nz > 0) {
    # first element is intercept
    coefskeep_z <- coefskeep_no_one[1:nz]
    coefskeep_x <- coefskeep_no_one[-c(1:nz)]
  } else {
    coefskeep_z <- NULL
    coefskeep_x <- coefskeep_no_one
  }
  #when use_predict = FALSE, each component will be computed
  #in predict.shapeselect, only kept predictors will be computed
  etacomps <- matrix(0, nrow = p, ncol = n)
  for(i in seq_along(active)){
    ps <- active[i]
    out_idx <- if (use_predict) i else ps
    type_ix <- if (use_predict) types[i] else types[ps]
    
    xi <- xmat[, out_idx]
    shpi <- active_shps[i]
    is_z <- ((shpi == 1 && type_ix == "nom") || (shpi == 2 && type_ix == "tree"))
    is_x_conv <- (shpi > 10 && shpi < 13 && type_ix == "num" || shpi > 2 && shpi < 5  && type_ix == "ord")
    is_x_no_conv <- (type_ix == "num" && !shpi %in% c(11, 12) || type_ix == "ord" && !shpi %in% c(3, 4) || type_ix == "tree" && shpi == 1)
    if(is_z){
      idz <- which(varlst_z == ps)
      etacomps[out_idx, ] <- ddkeep_z_all[, idz, drop = FALSE] %*% coefskeep_z[idz]
    } else if(is_x_conv){
      idx <- which(varlst == ps)
      idz <- which(varlst_z == ps)
      etacomps[out_idx, ] <- ddkeep_x_all[, idx, drop = FALSE] %*% coefskeep_x[idx]
      etacomps[out_idx, ] <- etacomps[out_idx, ] + ddkeep_z_all[, idz, drop = FALSE] %*% coefskeep_z[idz]
    } else if(is_x_no_conv){
      idx <- which(varlst == ps)
      etacomps[out_idx, ] <- ddkeep_x_all[, idx, drop = FALSE] %*% coefskeep_x[idx]
    }
  }
  etacomps
}

# ------------------------------------------------------------
#translate numeric shapes into characters
# ------------------------------------------------------------
shape_tags_to_char <- function(shps, types) {
  p <- length(types)
  
  tags <- rep("x", p)
  tags[types == "nom"]  <- "z"
  tags[types == "tree"] <- "tree"
  
  shapes <- mapply(ShapeToChar_shapeselect, shps, tags)
  
  return (shapes)
}

# ------------------------------------------------------------
#safe NCOL
# ------------------------------------------------------------
.nc_safe <- function(z) {if (is.null(z)) 0L else NCOL(z)}

#--------------------------------------------------------------------
#print.shapeselect: show the final choices
#--------------------------------------------------------------------
print.shapeselect <- function(x,...)
{
  mod_tried_print <- x$mod_tried
  mod_tried_print[c('penalty','edfu_start')] <- list(NULL, NULL)
  if ("penalty_cv" %in% names(mod_tried_print)) {
    mod_tried_print$penalty_cv <- NULL
  }
  print(mod_tried_print)

  mod_final <- paste(x$xnms[x$active], x$shapes[x$active], sep = " = ")
  cat("Final Model: \n")
  name_width <- max(nchar(x$xnms[x$active]))
  cat(sprintf("  %-*s\n", name_width + 3, mod_final), sep = "")
  #cat(paste(strwrap(mod_final, width = getOption("width") - 2, initial = "  ", prefix = "  "),
  #    collapse = "\n"), "\n", sep = "")
  #cat("Final model: ", paste(x$xnms[x$active], x$shapes[x$active], sep = " = ", collapse = "; "), "\n", sep = "")
  cat("\n")
  invisible(x)
}


# ------------------------------------------------------------
#update mod_tried
# ------------------------------------------------------------
append_mod_tried <- function(mod_tried, current_shps, cic, ps_keep, xid, action, xnms_with_id) {
  infocrit <- "CBIC"
  
  row <- c(
    as.list(current_shps),
    list(
      `# of x's` = sum(!current_shps %in% c(0, 17)),
      CBIC = cic,
      penalty = ps_keep,
      `new add` = xid,
      action = action
    )
  )
  
  names(row)[seq_along(xnms_with_id)] <- xnms_with_id
  rbind(mod_tried, as.data.frame(row, check.names = FALSE))
}

# ----------------------------------------------------------------------------------
#standardize xmat's numeric predictors, record their center and sd for prediction
#compute total number of possible edges roughly to adjust pnt 
# ----------------------------------------------------------------------------------
process_xmat <- function(xmat, types, boundary.knots = NULL,
                         scale.method = "minmax", add_kts = NULL, pnt = TRUE) {
  n <- nrow(xmat)
  p <- ncol(xmat)
  mult <- 1
  m_total <- 0
  
  mx  <- rep(0, p) #easier for the predict function
  sdx <- rep(1, p) #easier for the predict function
  
  for (ix in seq_len(p)) {
    #print (ix)
    if (is.numeric(xmat[, ix]) && types[ix] == "num") {
      kts0 <- NULL
      
      if (is.list(boundary.knots)) {
        kts0 <- boundary.knots[[ix]]
        
        if (!is.null(kts0) &&
            (!is.numeric(kts0) || any(!is.finite(kts0)) || length(kts0) < 2)) {
          warning(
            "Invalid 'boundary.knots': must be a numeric vector of at least two finite values. Using range(x) instead.",
            call. = FALSE
          )
          boundary.knots[ix] <- list(NULL)
          kts0 <- NULL
        }
      }
      
      #if (standardize) {
      if (scale.method == "zscore") {
        x_scaled <- scale(xmat[, ix], center = TRUE, scale = TRUE)
        #xmat[[ix]] <- x_scaled #works only for data.frame
        xmat[, ix] <- x_scaled #works for both matrix and data.frame
        sdx[ix] <- attr(x_scaled, "scaled:scale")
        
        #if (center) {
        mx[ix] <- attr(x_scaled, "scaled:center")
        #}
        
        if (!is.null(kts0) && length(kts0) >= 2) {
          #if (center) {
          kts0 <- kts0 - mx[ix]
          #}
          kts0 <- kts0 / sdx[ix]
          boundary.knots[[ix]] <- kts0
        }
      }
      
      # new: minmax
      if(scale.method == "minmax") {
        
        mx[ix] <- min(xmat[, ix], na.rm = TRUE)
        sdx[ix] <- max(xmat[, ix], na.rm = TRUE) -  min(xmat[, ix], na.rm = TRUE)
        x_scaled <- (xmat[, ix] - mx[ix]) / sdx[ix]
        xmat[, ix] <- x_scaled
        
        #skip for now
        if (!is.null(kts0) && length(kts0) >= 2) {
          kts0 <- kts0 - mx[ix]
          kts0 <- kts0 / sdx[ix]
          boundary.knots[[ix]] <- kts0
        }
      }
      
      #new: none
      if(scale.method == "none"){
        mx[ix] <- 0
        sdx[ix] <- 1
      }
      
      n1 <- length(unique(xmat[, ix]))
      #m_total <- m_total + mult * (trunc(n1^(1/7)) + 6)
      #new: use the edfu_plan 1
      m_total <- m_total + round(4 * n1^(1/9) + 10)
      # print (m_total)
      # not sure
      # if(p > n){
      #   cat('number of features > sample size!', 'use unpenalized estimator!', '\n')
      #   pnt <- FALSE
      # }
    }
    
    # if (inherits(xmat[, ix], "factor") && types[ix] %in% c("tree", "ord")) {
    #   xmat[, ix] <- as.numeric(as.character(xmat[, ix]))
    # }
    
    if (inherits(xmat[, ix], "factor") && types[ix] %in% c("tree", "ord")) {
      #xmat[, ix] <- as.numeric(as.character(xmat[, ix]))
      col <- xmat[, ix]
      xmat[, ix] <- as.numeric(levels(col))[col]
    }
    
    attr(xmat[, ix], "type") <- types[ix]
  }
  
  #new rule for add_kts
  # if(is.null(add_kts)){
  #   add_kts <- if (pnt) 0 else -2
  # }
  
  list(
    xmat = xmat,
    boundary.knots = boundary.knots,
    m_total = m_total,
    mx = mx,
    sdx = sdx,
    add_kts = add_kts,
    pnt = pnt
  )
}

# ------------------------------------------------------------
calculate_initial_cic <- function(y, family = gaussian(link = "identity"), 
                                  a = 1, infocrit = "CBIC") {
  wt.iter <- (family$family != "gaussian")
  n <- length(y)
  if (!wt.iter) {
    rss <- sum((y - mean(y))^2)
    # Degree = 1 at step 0
    degree <- 6/7
    return(log(rss/n) + log(n) / n^degree) #test
    #return(log(rss) + a * log(n) * 1 / n^degree)
  } else {
    # For GLMs, use the likelihood function
    llh <- llh.fun(y = y, muhat = rep(.5, n), etahat = rep(0, n), n = n, fml = family$family)
    #return(llh + a * log(n) * 1 / n)
    #new 10/7/2026
    degree <- 6/7
    return(llh + log(n) / n^degree)
  }
}

# ------------------------------------------------------------
#no_constraint: call solve
#with constraint: call solve.QP, pnt=T: xtx absorts dv or pnt=F
#no iteration
# ------------------------------------------------------------
solve_engine <- function(xtx, xty, amat = NULL, bvec = NULL, no_constraint = FALSE){
  if(no_constraint){
    ans <- try(solve(xtx, xty), silent = TRUE)
    if(inherits(ans, 'try-error')){
      #cat('see error in solve in the engine', '\n')
      ev <- eigen(xtx, only.values = TRUE)$values
      # check if the smallest eigenvalue is too close to zero or negative
      #new:
      tol <- 1e-8 * max(1, max(abs(ev)))
      if (min(ev) < tol) {
        diag(xtx) <- diag(xtx) + (tol - min(ev))
      }
      ans <- solve(xtx, xty)
    } 
    bhat <- ans
    face <- 1:length(bhat)
  } else {
    ans <- try(solve.QP(xtx, xty, amat, bvec), silent = TRUE)
    if(inherits(ans, 'try-error')){
      #cat('see error in solve.QP in the engine', '\n')
      ev <- eigen(xtx, only.values = TRUE)$values
      # check if the smallest eigenvalue is too close to zero or negative
      # if (min(ev) < 1e-8) {
      #   diag(xtx) <- diag(xtx) + 1e-8
      # }
      #new:
      tol <- 1e-8 * max(1, max(abs(ev)))
      if (min(ev) < tol) {
        diag(xtx) <- diag(xtx) + (tol - min(ev))
      }
      
      ans <- solve.QP(xtx, xty, amat, bvec)
    }
    bhat <- ans$solution
    face <- ans$iact
  }
  return (list(bhat = bhat, face = face))
}

# ------------------------------------------------------------
# LS or IRLS, pnt or not pnt
# ------------------------------------------------------------
fit_model_engine <- function(y, current_gmat, current_dv, family, 
                             pnt = TRUE, np_ix, varlst_s, is_first_step = FALSE,
                             compute_gcv = FALSE, initial_tune = FALSE) {
  wt.iter <- (family$family == "binomial")
  # This essentially performs the optimization:
  # min ||y - Xb||^2 s.t. Ab >= 0 (for constrained shapes)
  # or uses IRLS for GLM families
  gmat <- current_gmat
  dv <- current_dv
  nc <- ncol(gmat)
  n <- length(y)
  ps_zero <- 1:np_ix
  
  #new: for linear space and unconstrained x comps
  #move varlst_s by np_ix if np_ix > 0 s.t. it can track all columns in gmat
  if(any(varlst_s == 1)){
    if(is_first_step) { #at the first step, the varlst_s used in forward is only for smooth x
      ps_zero <- c(ps_zero, np_ix + which(varlst_s == 1))
    } else {
      ps_zero <- c(ps_zero, which(varlst_s == 1))
    }
  }
  
  #if(length(ps_zero) == nc) { #every column in gmat is either linear or unconstrained
  #new:
  if(length(ps_zero) == nc || initial_tune){
    no_constraint <- TRUE
    amat <- NULL
    bvec <- NULL
  } else {
    no_constraint <- FALSE
    amat <- diag(nc)
    amat[ps_zero, ps_zero] <- 0
    bvec <- rep(0, nc)
  }
  
  qv0 <- crossprod(gmat)
  xtx <- qv0
  #test! avoid linear, shp = 18
  if(pnt){
    is_empty_matrix <- function(x) {
      is.matrix(x) && identical(dim(x), c(0L, 0L))
    }
    if(!is.null(dv) && !is_empty_matrix(dv)){ #dv = NULL or <0 x 0 matrix> means nom/tree(2) only
      dv <- combine_penalties(fixed_dv = dv, np_ix = np_ix)
      xtx <- xtx + dv
    } else {
      dv <- matrix(0, nc, nc)
    }
  } else {
    dv <- matrix(0, nc, nc)
  }
  
  if (!wt.iter) {
    # gaussian case: Constrained Least Squares
    # If dmat is provided, it's a penalized fit; if not, it's a pure projection
    xty <- crossprod(gmat, y)
    ans <- solve_engine(xtx, xty, amat=amat, bvec=bvec, no_constraint=no_constraint)
    bhat <- ans$bhat
    bhat <- round(bhat, 10)
    etahat <- gmat %*% bhat
    muhat <- etahat
    face <- ans$face
  } else {
    # GLM Case: IRLS loop
    # Implementation follows standard GLM logic but calls 
    # constrained solver at each weight update
    gmat0 <- gmat
    weights <- rep(1, n)
    etahat <- etahat.fun(n, y, fml = family$family)
    gr <- gr.fun(y, etahat, weights, fml = family$family)
    wt <- wt.fun(y, etahat, n, weights, fml = family$family)
    cvec <- wt * etahat - gr
    zvec <- zvec.fun(cvec, wt, y, fml = family$family)
    muhat <- family$linkinv(etahat)
    sm <- 1e-6
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
      xtx <- crossprod(gmat)
      xtx <- xtx + dv #dv is not changed; pnt=F, then dv=0
      xty <- crossprod(gmat, zvec)
      ans <- solve_engine(xtx, xty, amat=amat, bvec=bvec, no_constraint=no_constraint)
      bhat <- ans$bhat
      #bhat <- round(bhat, 10)
      face <- ans$face
      etahat <- gmat0 %*% bhat
      muhat <- family$linkinv(etahat)
      mdiff <- abs(max(muhat) - 1)
      if (family$family == "binomial") {
        mdiff <- abs(max(muhat) - 1) > sm
      } else {mdiff <- TRUE}
      diff <- mean((muhat - oldmu)^2)
    }
  }
  
  #new: for tuning edfu
  gcv <- NULL
  if(compute_gcv){
    llh <- llh.fun(y, muhat, etahat, phihat = NULL, n, weights = NULL, fml = family$family) * n / 2
    sse <- sum((y - muhat)^2)
    dev <- ifelse(wt.iter, llh, sse)
    umat <- try(chol(xtx), silent = TRUE)
    if(inherits(umat, 'try-error')){
      ev <- eigen(xtx, only.values = TRUE)$values
      # check if the smallest eigenvalue is too close to zero or negative
      # if (min(ev) < 1e-8) {
      #   diag(xtx) <- diag(xtx) + 1e-8
      # }
      #new:
      tol <- 1e-8 * max(1, max(abs(ev)))
      if (min(ev) < tol) {
        diag(xtx) <- diag(xtx) + (tol - min(ev))
      }
      umat <- chol(xtx)
    }
    uinv <- backsolve(umat, diag(ncol(xtx)))
    if(!no_constraint){
      atil <- amat %*% uinv
    }
    xtxinv <- tcrossprod(uinv)
    if(!wt.iter){
      bigp_base <- xtxinv %*% qv0
    } else {
      bigp_base <- xtxinv %*% crossprod(gmat)
    }
    edf_base <- sum(diag(bigp_base))
    
    if (all(face == 0) || length(face) == 0 || no_constraint) {
      edf <- edf_base
    } else {
      smat <- -t(atil[face, , drop = FALSE])
      qr_smat <- qr(smat)
      smat <- qr.Q(qr_smat)
      pmat_polar <- tcrossprod(smat)
      
      if (!wt.iter) {
        bigp_polar <- uinv %*% tcrossprod(pmat_polar, uinv) %*% qv0
      } else {
        bigp_polar <- uinv %*% tcrossprod(pmat_polar, uinv) %*% crossprod(gmat)
      }
      edf <- edf_base - sum(diag(bigp_polar))
    }
    gcv <- dev / (1 - edf / n)^2
  }
  #add: xtx, amat, face, qv0 for compute_edf0
  #gmat is not the same after the while loop for binomial; needs to be returned for edf0 compute
  fit <- list(bhat = bhat, etahat = etahat, muhat = muhat, xtx = xtx, amat = amat,
              face = face, qv0 = qv0, no_constraint = no_constraint, gmat = gmat, gcv = gcv)
  return(fit)
}

# ------------------------------------------------------------
#used when finished each step and update the dmat_lst for all?
# ------------------------------------------------------------
combine_penalties <- function(fixed_dv, new_dmat = NULL, np_ix = 0) {
  no_fixed_dv <- is.null(fixed_dv) ||
    (is.matrix(fixed_dv) && identical(dim(fixed_dv), c(0L, 0L)))
  if (no_fixed_dv && !is.null(new_dmat)) {
    dv <- crossprod(new_dmat)
  } else if (!no_fixed_dv && !is.null(new_dmat) && np_ix == 0) {
    p1 <- ncol(fixed_dv)
    p2 <- ncol(new_dmat)
    dv <- matrix(0, p1 + p2, p1 + p2)
    dv[1:p1, 1:p1] <- fixed_dv
    dv[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- crossprod(new_dmat)
  } else if (!no_fixed_dv && is.null(new_dmat) && np_ix > 0) {
    p <- ncol(fixed_dv)
    dv <- matrix(0, np_ix + p, np_ix + p)
    dv[(np_ix + 1):(np_ix + p), (np_ix + 1):(np_ix + p)] <- fixed_dv
  } else if (!no_fixed_dv && !is.null(new_dmat) && np_ix > 0) {
    p1 <- ncol(fixed_dv)
    p2 <- ncol(new_dmat)
    dv <- matrix(0, np_ix + p1 + p2, np_ix + p1 + p2)
    dv[(np_ix + 1):(np_ix + p1), (np_ix + 1):(np_ix + p1)] <- fixed_dv
    dv[(np_ix + p1 + 1):(np_ix + p1 + p2),
       (np_ix + p1 + 1):(np_ix + p1 + p2)] <- crossprod(new_dmat)
  } else if (!no_fixed_dv && is.null(new_dmat) && np_ix == 0) {
    dv <- fixed_dv
  } else {
    stop("no fixed_dv or np_ix or new_dmat")
  }
  dv
}

# ------------------------------------------------------------
calculate_llh <- function(y, etahat, family) {
  wt.iter <- (family$family == "binomial")
  n <- length(y)
  if (!wt.iter) {
    rss <- sum((y - etahat)^2)
    #llh <- log(rss)
    llh <- log(rss/n)
  } else {
    llh <- llh.fun(y, muhat = family$linkinv(etahat), etahat, phihat = NULL, n, weights = NULL, fml = family$family)
  }
  llh 
}

# ------------------------------------------------------------
# candidate shapes for one predictor
# ------------------------------------------------------------
get_candidate_shapes <- function(x, type_ix, set4 = FALSE, wt.iter = FALSE, step = 1) {
  #step = 1: forward 
  #step = 2: backward
  # if (wt.iter) {
  #   set4 <- TRUE
  # }
  #if (inherits(x, "factor") & identical(type_ix, "nom")) {
  if(type_ix == "nom") {
    #if(step == 1) 
    shps <- 1
    #if(step == 2) shps <- c(1, 0)
    #return("nom")
  } else if (type_ix == "num") {
    #if(step == 1){
    shps <- if (set4) c(9:12, 19) else c(9:16, 19)
    #} 
    #if(step == 2){
    #  shps <- if (set4) c(9:12, 19, 17) else c(9:16, 19, 17)
    #}  
  } else if (type_ix == "ord") {
    #if(step == 1){
    shps <- if (set4) 1:4 else 1:8
    #}
    #if(step == 2){
    #  shps <- if (set4) c(1:4, 0) else c(1:8, 0)
    #}  
  } else if (type_ix == "tree") {
    #if(step == 1){
    shps <- 1:2
    #} 
    #if(step == 2){
    #  shps <- c(1, 2, 0)
    #} 
  } else {
    stop("Unknown predictor type: ", type_ix)
  }
  shps
}

# ------------------------------------------------------------
# main precomputation routine 
# add user_edfus
# ------------------------------------------------------------
precompute_lists <- function(xmat, types, boundary.knots, pnt = TRUE, set4 = FALSE, 
                             wt.iter = FALSE, infocrit = "CBIC", add_kts = 0, add_edfu = 0, 
                             user_edfus = NULL, space = "E", edfu_plan = 1, weights = NULL){
  #list2env(params, envir = environment())
  
  # expected names in params:
  # xmat, types, boundary.knots, pnt, set4, wt.iter, infocrit,
  # add_kts, add_edfu, space, edfu_plan
  
  #default: forward step
  # if(is.null(step)){
  #   step <- 1
  # }
  # 
  capl <- NCOL(xmat)
  n <- NROW(xmat)
  one <- matrix(1, nrow = n, ncol = 1)
  
  if (is.null(boundary.knots)) {
    boundary.knots <- vector("list", capl)
  }
  
  if (length(boundary.knots) < capl) {
    stop("boundary.knots must be NULL or a list of length at least NCOL(xmat).")
  }
  
  # helper: safe access to boundary knots
  .get_boundary_knots_i <- function(boundary.knots, ix) {
    if (length(boundary.knots) < ix) return(numeric(0))
    kts0 <- boundary.knots[[ix]]
    if (is.null(kts0)) numeric(0) else kts0
  }
  
  # helper: compute edfu_un for shape 19
  .compute_edfu_un <- function(gmat, dm, np_ix) {
    qv0 <- crossprod(gmat)
    dv <- combine_penalties(fixed_dv = crossprod(dm), np_ix = np_ix)
    xtx <- qv0 + dv
    umat <- chol(xtx)
    qvinv <- chol2inv(umat)
    sum(diag(qvinv %*% qv0))
  }
  
  #.nc_safe <- function(z) {if (is.null(z)) 0L else NCOL(z)}
  
  precomp <- vector("list", capl)
  #new: keep track of the 1st continuous predictor and use its ps_vec for future continuous predictors
  first_num_idx <- NULL
  first_num_ps_vec <- NULL
  for (ix in seq_len(capl)) {
    #print (ix)
    x <- xmat[, ix]
    type_ix <- types[ix]
    #new:
    #no need to include shp=17 to compute the elements
    shps <- get_candidate_shapes(x = x, type_ix = type_ix, set4 = set4, wt.iter = wt.iter)
    nobs <- length(shps)
    
    gmat_lst <- vector("list", length = nobs)
    np_ix_lst <- vector("list", length = nobs)
    ms_lst <- vector("list", length = nobs)
    dd_lst <- vector("list", length = nobs)
    kts_lst <- vector("list", length = nobs)
    dmat_lst <- vector("list", length = nobs)
    dv_lst <- vector("list", length = nobs)
    #new: used in tuning
    dv_lst_raw <- vector("list", length = nobs)
    ps_vec <- rep(0, length = nobs)
    edfu_vec <- rep(0, length = nobs)
    #edfu_vec_unpnt <- rep(0, length = nobs)
    varlst_s_lst <- vector("list", length = nobs) #varlst_s has the same length of all x comps' edges, not ncol(gmat)
    edfu_un <- NULL
    
    names(gmat_lst) <- names(np_ix_lst) <- names(ms_lst) <- names(dd_lst) <-
      names(kts_lst) <- names(dv_lst) <- names(dv_lst_raw) <- names(dmat_lst) <- names(varlst_s_lst) <- as.character(shps)
    names(ps_vec) <- names(edfu_vec) <- as.character(shps)
    
    for (i in seq_len(nobs)) {
      shpi <- shps[i]
      #print (shpi)
      #print (type_ix)
      dd <- ms <- kts <- gmat <- np_ix <- dmat <- dv <- dv_raw <- varlst_s <- NULL
      ps <- 0
      
      if(type_ix == "nom" && shpi == 1){#for nom z
        #new: some user forgets to code x as factor
        if (!is.factor(x)) {
          #message("This predictor is specified as type 'nom' but is not a factor; converting it to a factor.")
          x <- factor(x)
        }
        dd <- model.matrix(~ x)[, -1, drop = FALSE]
        ms <- colMeans(dd)
        dd <- scale(dd, center = TRUE, scale = FALSE)
        newz <- dd
        newx <- NULL
        
        gmat <- cbind(one, newz)
        #not use: just for conceptual correctness
        #dmat <- matrix(0, nrow = (ncol(dd)-1), ncol = ncol(dd))
        #penalty matrix 
        #dv <- matrix(0, nrow = ncol(dd), ncol = ncol(dd))
        ps <- 0
        np_ix <- 1 + ncol(newz)
      } else if (type_ix %in% c("ord", "tree") && shpi <= 8) {#for x: ord, num, tree
        if (type_ix == "ord") {
          dd_ans <- makedelta(x, sh = shpi, numknots = 0, knots = 0, space = space)
          #centered terms used for prediction
          ms <- dd_ans$ms
          #edges
          dd <- t(dd_ans$amat)
          #penalty matrix 
          dmat <- matrix(0, nrow = (ncol(dd) - 1), ncol = ncol(dd))
          dv <- matrix(0, nrow = ncol(dd), ncol = ncol(dd))
          #new: dv_raw is for gcv tuning, not used for nom or ord but we need it as a retuend obj
          dv_raw <- dv
          ps <- 0
          
          if (shpi > 2 && shpi < 5) {
            newz <- x - mean(x)
            #newx <- dd
            #gmat <- cbind(one, x - mean(x), dd)
            #np_ix <- 2
          } else {
            newz <- NULL
            #newx <- dd
            #gmat <- cbind(one, dd)
            #np_ix <- 1
          }
          newx <- dd
          gmat <- cbind(one, newz, newx)
          np_ix <- 1 + .nc_safe(newz)
        }
        
        if (type_ix == "tree") {
          dd <- t(tree.fun(x))
          ms <- colMeans(dd)
          dd <- scale(dd, center = TRUE, scale = FALSE)
          #gmat <- cbind(one, dd)
          
          if (shpi == 1) {#tree
            newx <- dd
            newz <- NULL
            dmat <- matrix(0, nrow = (ncol(dd) - 1), ncol = ncol(dd))
            dv <- matrix(0, nrow = ncol(dd), ncol = ncol(dd))
            dv_raw <- dv
            ps <- 0
            #np_ix <- 1
          }
          
          if (shpi == 2) {#nominal
            newx <- NULL
            newz <- dd
            #penalty matrix 
            #dmat <- matrix(0, nrow = ncol(dd), ncol = ncol(dd))
            ps <- 0
            #np_ix <- 1 + ncol(dd)
          }
          
          #not use: dmat, just for conceptual correctness
          #dmat <- matrix(0, nrow = (ncol(dd)-1), ncol = ncol(dd))
          #dv <- matrix(0, nrow = ncol(dd), ncol = ncol(dd))
          gmat <- cbind(one, newz, newx)
          np_ix <- 1 + .nc_safe(newz)
        }
      } else if (type_ix == "num" && shpi > 8) {      #for type_ix = 'num'
        mult <- ifelse(pnt, 2, 1)
        spl_degree_i <- ifelse(shpi %in% c(9, 10), 2, 3)
        n1 <- length(unique(x))
        #test
        #n1 <- length(x)
        #nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/7)) + 6 + add_kts, trunc(n1^(1/9)) + 6 + add_kts)
        #test nkts
        if(edfu_plan == 1){
          #nkts <- mult * ifelse(spl_degree_i == 2, .5 * round(n1^(1/7)*6), .5 * round(n1^(1/9)*6 + 6))
          if(shpi %in% c(9, 10)) {
            nkts <- round(6 * n1^(1/7))
          } else if (shpi %in% c(11, 12)){
            nkts <- round(4 * n1^(1/9) + 10)
          } else {
            nkts <- round(4 * n1^(1/9) + 10)
          }
          #temp
          if(!pnt) nkts <- round(nkts/2)
        } else {
          nkts <- mult * ifelse(spl_degree_i == 2, trunc(n1^(1/7)) + 6, trunc(n1^(1/9)) + 6)
        }
        
        # if (spl_degree_i == 2) {
        #   edfu <- ifelse(edfu_plan == 1, nkts / mult, n1^(1/7) + 3)
        # } else if (spl_degree_i == 3 && shpi %in% c(11, 12)) {
        #   edfu <- ifelse(edfu_plan == 1, nkts / mult + 1, n1^(1/9) + 4)
        # } else if (spl_degree_i == 3 && shpi > 12) {
        #   if (shpi != 19) {#13, etc
        #     #test edfu:
        #     #test_edfu <- mult * (trunc(n1^(1/9)) + 6) - add_edfu
        #     #test_edfu <- nkts - add_edfu
        #     #edfu <- ifelse(edfu_plan == 1, nkts / mult + 1 + add_edfu, test_edfu + add_edfu)
        #     edfu <- ifelse(edfu_plan == 1, nkts / mult + 1 + add_edfu, n1^(1/9) + 5 + add_edfu)
        #   } else {
        #     edfu <- ifelse(edfu_plan == 1, nkts / mult + 1, n1^(1/9) + 5)
        #   }
        # }
        
        #test
        #the else part is not used
        #call tune_start_edfu and then store the result in user_edfus
        if (spl_degree_i == 2) {
          edfu <- if(!is.null(user_edfus)) user_edfus[1] else nkts * .5 - 3 #  0.4 * nkts - 1  #  nkts * .4  # 0.4 * nkts - 1  # nkts * .5 - 1.5 # 0.4 * nkts - 1 # nkts * .5 - 1.5 # nkts * .5 - 1.5 # 0.4 * nkts - 1# nkts * .5 - 1.5 # 0.4 * nkts # round(nkts / 2) #0.4 * nkts  #nkts / 2 - 2
        } else if (spl_degree_i == 3 && shpi %in% c(11, 12)) {
          edfu <- if(!is.null(user_edfus)) user_edfus[2] else nkts * .5 - 3 # 0.3 * nkts - 1 # nkts * .4   # 0.3 * nkts - 1 # nkts * .5 - 1 # 0.3 * nkts - 1 # nkts * .5 - 1 # nkts * .5 - 1 # 0.3 * nkts - 1# nkts * .5 - 1 # 0.3 * nkts # round(nkts / 2) #0.3 * nkts #nkts / 2 
        } else if (spl_degree_i == 3 && shpi > 12) {
          if (shpi != 19) {#13, etc
            edfu <- if(!is.null(user_edfus)) user_edfus[3] else nkts * .5 + 3 # nkts * .5 + 3 #nkts * .5 + 3   # 0.3 * nkts + 3 # nkts * .5 + 3 # 0.3 * nkts + 3 # nkts * .5 + 3 #nkts * .5 + 3 #0.3 * nkts + 3 # nkts * .5 + 3 # 0.3 * nkts + 3 #round(nkts / 2) + 1 #0.3 * nkts + 3 #nkts / 2 + 1 + add_edfu #- 3
          } else {
            #user_edfus[3] is too large when the real shape is 19 -> 19 not chosen
            edfu <- if(!is.null(user_edfus)) user_edfus[2] else nkts * .5 + 3 # nkts * .4  # 0.3 * nkts - 1 # nkts * .5 - 3 # 0.3 * nkts - 1# nkts * .5 - 1 # 0.3 * nkts #round(nkts / 2) # 0.3 * nkts #nkts / 2 + 1 #- 3
          }
        }
        
        if(pnt){
          edfu_vec[i] <- edfu
        }
       
        kts <- seq.int(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = nkts)
        if(anyNA(x)){
          message('some predictor has NA!')
        }
        kts0 <- .get_boundary_knots_i(boundary.knots, ix)
        if (length(kts0) >= 2) {
          x_rng <- range(x, na.rm = TRUE)
          kts0_rng <- range(kts0, na.rm = TRUE)
          tol <- 1e-10
          reset_lwr <- kts0_rng[1] > x_rng[1] + tol
          reset_upp <- kts0_rng[2] < x_rng[2] - tol
          
          kts[1] <- if (reset_lwr) x_rng[1] else kts0_rng[1]
          kts[nkts] <- if (reset_upp) x_rng[2] else kts0_rng[2]
        }
        
        dd_ans <- makedelta(x, sh = ifelse(shpi != 19, shpi, (shpi-2)), numknots = 0,
                            knots = kts, space = space)
        
        ms <- dd_ans$ms
        dd <- t(dd_ans$amat)
        
        if (shpi > 10 && shpi < 13) {
          newz <- x - mean(x)
          #gmat <- cbind(one, x - mean(x), dd)
          #np_ix <- 2
        } else {
          newz <- NULL
          #gmat <- cbind(one, dd)
          #np_ix <- 1
        }
        newx <- dd
        gmat <- cbind(one, newz, newx)
        np_ix <- 1 + .nc_safe(newz)
        #varlst_s <- rep(0, ncol(newx))
        
        #new: record the unpenalized edfu if pnt=F
        if(!pnt){
          edfu_vec[i] <- ncol(dd)
        }
        
        #new: as an approx, if we do minmax scaling then use the same ps for each continuous predictor
        if (pnt) {
          nc <- NCOL(dd)
          #check!
          dm0 <- makedmat_1D_ic(m = nc, knots = kts, shpi = shpi)
          qvi <- crossprod(dd)
          dvi <- crossprod(dm0)
          
          #new: cap edfu at 1.5 (smallest) and nkts-1 (largest)
          #edfu <- min(max(edfu, 1.5), nkts - 1)
          edfu <- min(max(edfu, 2), nkts - 1)
          
          #new: 1st continuous x, compute ps
          #if(tune_start){
          #  ps <- 1
          #  ps_vec[i] <- ps
          #} else if(!is.null(user_ps)){
          #if(!is.null(user_ps)){
          #  ps <- user_ps[i]
          #  ps_vec[i] <- ps
          #} else if(is.null(first_num_idx)){
          
          if(is.null(first_num_idx)){
            if (shpi %in% c(9, 11, 13, 19)) { 
              #new: compute ps for 19 because when set4 = T, ps for 13 is not computed but user_edfu for 19 is assigned as that for 13
              #-18 might be too small to cause chol() error in .search_ps_sel
              #replace it with -10
              log_lo <- -10
              log_hi <- 18
              
              #new
              f_lo <- try(.search_ps_sel(exp(log_lo), qv0 = qvi, dv = dvi,
                                     dd = dd, edfu = edfu, weights = weights), silent = TRUE)
              
              #new: if the lower endpoint fails, move it inward
              while (inherits(f_lo, "try-error") && log_lo < log_hi) {
                log_lo <- log_lo + 1
                f_lo <- try(.search_ps_sel(exp(log_lo), qv0 = qvi, dv = dvi,
                                       dd = dd, edfu = edfu, weights = weights), silent = TRUE)
              }
              
              f_hi <- try(.search_ps_sel(exp(log_hi), qv0 = qvi, dv = dvi,
                                     dd = dd, edfu = edfu, weights = weights), silent = TRUE)
              
              #new: if the upper endpoint fails, move it inward
              while (inherits(f_hi, "try-error") && log_hi > log_lo) {
                log_hi <- log_hi - 1
                f_hi <- try(.search_ps_sel(exp(log_hi), qv0 = qvi, dv = dvi,
                                       dd = dd, edfu = edfu, weights = weights),silent = TRUE)
              }
              
              if (inherits(f_lo, "try-error") || inherits(f_hi, "try-error")) {
                stop("could not find two valid penalty endpoints.")
              } else if (f_lo <= 0) {
                # the target edfu requires a penalty below the usable range.
                # use the smallest penalty that was successfully evaluated.
                ps <- exp(log_lo)
              } else if (f_hi >= 0) {
                # the target edfu requires a penalty above the searched range.
                # use the largest penalty that was successfully evaluated.
                ps <- exp(log_hi)
              } else {
                # f_lo > 0 and f_hi < 0, the root can be found within the interval
                log_root_res <- try(uniroot(f = function(log_pen, ...) {
                  .search_ps_sel(exp(log_pen), ...)
                }, interval = c(log_lo, log_hi), qv0 = qvi,
                dv = dvi, dd = dd, edfu = edfu, weights = weights, extendInt = "no",
                tol = .Machine$double.eps^0.25), silent = TRUE)
                
                if (!inherits(log_root_res, "try-error")) {
                  ps <- exp(log_root_res$root)
                } else {
                  # both endpoints are valid and contains the root, so this is an
                  # unexpected numerical failure. use the endpoint closest to the root.
                  if (abs(f_lo) <= abs(f_hi)) {
                    ps <- exp(log_lo)
                  } else {
                    ps <- exp(log_hi)
                  }
                  #cat("unexpected try-error in uniroot\n")
                  #cat("edfu:", edfu, "\n")
                  #cat("nkts:", nkts, "\n")
                  #cat("usable log interval:", log_lo, log_hi, "\n")
                  #cat("endpoint values:", f_lo, f_hi, "\n")
                  #print(log_root_res)
                }
              }
              
              ps_vec[i] <- ps
            
            } else if (shpi == 10) {
              # borrows ps from shape 9
              ps <- ps_vec[which(shps %in% 9)]
              ps_vec[i] <- ps
            } else if (shpi == 12) {
              # borrows ps from shape 11
              ps <- ps_vec[which(shps %in% 11)]
              ps_vec[i] <- ps
            } else {
              # borrows ps from shape 13
              ps <- ps_vec[which(shps %in% 13)]
              ps_vec[i] <- ps
            }
          } else {
            ps <- as.numeric(first_num_ps_vec[as.character(shpi)])
            ps_vec[i] <- ps
          }
          
          #raw
          dv_raw <- crossprod(dm0)
          dv <- ps * crossprod(dm0)
          dmat <- sqrt(ps) * dm0
        }
      }
      
      # for future QP logic: 1 marks unconstrained smooth basis columns
      newx_smooth_track <- if (shpi == 19) rep(1, .nc_safe(newx)) else rep(0, .nc_safe(newx))
      varlst_s <- newx_smooth_track
      
      dd_lst[[i]] <- dd
      ms_lst[[i]] <- ms
      if(!is.null(kts)){
        kts_lst[[i]] <- kts
      }
      gmat_lst[[i]] <- gmat
      np_ix_lst[[i]] <- np_ix
      if(!is.null(dv)){
        dv_lst[[i]] <- dv #absored ps; 19 is different from 13
        #new: not absorb ps
        dv_lst_raw[[i]] <- dv_raw
      }
      
      if(!is.null(dmat)){
        dmat_lst[[i]] <- dmat
      }
      
      #dm_lst[[i]] <- dm #no need; same as dmat
      ps_vec[i] <- ps
      varlst_s_lst[[i]] <- varlst_s #same length as ncol(newx), only shp=19 columns = 1
      
      if (shpi == 19 && !wt.iter) { #skip computing edfu_un for binomial
        if (pnt && !is.null(dmat)) {
          edfu_un <- .compute_edfu_un(gmat, dmat, np_ix)
        } else {
          edfu_un <- NCOL(gmat)
        }
      }
    }
    
    #not use: we switch to data driven approach and each x should have diff ps
    #if (type_ix == "num" && is.null(first_num_idx)) {
    #  first_num_idx <- ix
    #  first_num_ps_vec <- ps_vec
    #}
    
    precomp[[ix]] <- list(shps = shps, gmat_lst = gmat_lst, np_ix_lst = np_ix_lst, 
                          ms_lst = ms_lst, dd_lst = dd_lst, kts_lst = kts_lst,
                          dmat_lst = dmat_lst, dv_lst = dv_lst, dv_lst_raw = dv_lst_raw,
                          ps_vec = ps_vec, varlst_s_lst = varlst_s_lst, 
                          edfu_un = edfu_un, edfu_vec = edfu_vec, 
                          type_ix = type_ix)
  }
  
  precomp
}


# ------------------------------------------------------------
#update pc with shp = 0 or 17
#for the backward step only
# ------------------------------------------------------------
update_pc_with_flat <- function(pc){
  shps <- pc$shps
  gmat_lst <- pc$gmat_lst
  np_ix_lst <- pc$np_ix_lst
  dd_lst <- pc$dd_lst
  ms_lst <- pc$ms_lst
  kts_lst <- pc$kts_lst
  dmat_lst <- pc$dmat_lst
  dv_lst <- pc$dv_lst
  ps_vec <- pc$ps_vec
  varlst_s_lst <- pc$varlst_s_lst
  edfu_un <- pc$edfu_un
  type_ix <- pc$type_ix
  
  if(type_ix == "nom") {#add 0 
    shps <- c(shps, 0)
  } else if (type_ix == "tree" || type_ix == "ord") {#add 17
    shps <- c(shps, 0)
  } else if (type_ix == "num") {#add 17
    shps <- c(shps, 17)
  }
  
  gmat_lst <- append(gmat_lst, list(NULL))
  np_ix_lst <- append(np_ix_lst, list(NULL))
  ms_lst <- append(ms_lst, list(NULL))
  dd_lst <- append(dd_lst, list(NULL))
  kts_lst <- append(kts_lst, list(NULL))
  dv_lst <- append(dv_lst, list(NULL))
  dmat_lst <- append(dmat_lst, list(NULL))
  varlst_s_lst <- append(varlst_s_lst, list(NULL))
  ps_vec <- c(ps_vec, 0)
  
  names(gmat_lst) <- names(np_ix_lst) <- names(ms_lst) <- names(dd_lst) <-
    names(kts_lst) <- names(dv_lst) <- names(dmat_lst) <- names(varlst_s_lst) <- as.character(shps)
  names(ps_vec) <- as.character(shps)
  
  pc$shps <- shps
  pc$gmat_lst <- gmat_lst
  pc$np_ix_lst <- np_ix_lst
  pc$dd_lst <- dd_lst
  pc$ms_lst <- ms_lst
  pc$kts_lst <- kts_lst
  pc$dmat_lst <- dmat_lst
  pc$dv_lst <- dv_lst
  pc$ps_vec <- ps_vec
  pc$varlst_s_lst <- varlst_s_lst
  #pc$edfu_un <- edfu_un
  #pc$type_ix <- type_ix
  return (pc)
}

# ------------------------------------------------------------
# update precomputed candidate lists after current model changes
# ------------------------------------------------------------
update_precomp_lists <- function(x, type_ix, ddkeep_x_all, ddkeep_z_all, dv_keep_all=NULL, varlst_s_keep,
                                 pc, pnt = TRUE, skip_shape = NULL) {
  
  #for the backward:
  shps_all <- pc$shps
  keep_id <- seq_along(shps_all)
  if (!is.null(skip_shape)) {
    keep_id <- keep_id[shps_all != skip_shape]
  }
  
  shps <- shps_all[keep_id]
  #shps <- pc$shps
  #type_ix <- pc$type
  dd_lst <- pc$dd_lst[keep_id]
  dm_lst <- pc$dmat_lst[keep_id] # individual dmat absorbing sqrt(ps)
  #ms_lst <- pc$ms_lst
  #kts_lst <- pc$kts_lst
  #ps_vec <- pc$ps_vec
  
  n <- length(x)
  one <- matrix(1, nrow = n, ncol = 1)
  nobs <- length(shps)
  gmat_lst <- vector("list", nobs)
  #new!
  dv_lst <- vector("list", nobs)
  np_ix_lst <- vector("list", nobs)
  varlst_s_lst <- vector("list", nobs)
  
  #.nc <- function(z) if (is.null(z)) 0L else NCOL(z)
  
  for (i in seq_along(shps)) {
    shpi <- shps[i]
    #print (shpi)
    #print (type_ix)
    dd <- dd_lst[[i]]
    #ms <- ms_lst[[i]]
    #kts <- kts_lst[[i]]
    #ps <- ps_vec[i]
    dm <- dm_lst[[i]]
    gmat <- dmat <- np_ix <- varlst_s <- NULL #dv = NULL means pnt=F
    dv <- NULL
    #newz <- NULL
    #newx <- NULL
    #newx_smooth_track <- NULL
    
    if(type_ix == "nom"){#nom
      if(shpi == 1){
        newz <- dd
        newx <- NULL
        # if(!is.null(dv_keep_all)){
        #   dv <- combine_penalties(fixed_dv = dv_keep_all, np_ix = .nc(newz))
        # }
        dv <- dv_keep_all
        #dmat <- NULL
        #ps <- 0
      } else if (shpi == 0) {
        #do nothing when shpi = 0
        newz <- NULL
        newx <- NULL
        dv <- dv_keep_all
      }
      newx_smooth_track <- NULL
    } else if (type_ix %in% c("ord", "tree") && shpi <= 8) {#ord and tree
      
      if(shpi == 0) {
        #do nothing 
        newz <- NULL
        newx <- NULL
        dv <- dv_keep_all
        newx_smooth_track <- NULL
      } else {
        if (type_ix == "ord") {
          if (shpi > 2 && shpi < 5) {
            newz <-  x - mean(x)
          } else {
            newz <- NULL
          }
          newx <- dd
        }
        if (type_ix == "tree") {
          if (shpi == 1) {
            newz <- NULL
            newx <- dd
          }
          if (shpi == 2) {
            newx <- NULL
            newz <- dd
          }
        }
        # if (!is.null(dmat_keep_all) && !is.null(dm)) {
        #   dmat <- combine_penalties(fixed_dmat = dmat_keep_all, new_dmat = dm) 
        # } else if (!is.null(dmat_keep_all)) {
        #   dmat <- dmat_keep_all
        # } else {
        #   dmat <- dm
        # }
        #check more: should be the same as use np_ix?
        #if(!is.null(dv_keep_all)){
        if (pnt) {
          dv <- combine_penalties(fixed_dv = dv_keep_all, new_dmat = dm)
        } else {
          dv <- dv_keep_all
        }
        #}
        #dv <- dv_keep_all
        newx_smooth_track <- rep(0, .nc_safe(newx))
      }
      
    } else if (type_ix == "num" && shpi > 8) {#num
      
      if(shpi == 17) {
        #do nothing; equivalent to shpi = 0 when type_ix = 'nom' or 'ord' or 'tree
        newz <- NULL
        newx <- NULL
        dv <- dv_keep_all
        newx_smooth_track <- NULL
      } else {
        if (shpi > 10 && shpi < 13) {
          newz <- x - mean(x)
        } else {
          newz <- NULL
        }
        newx <- dd
        if (pnt) {
          dv <- combine_penalties(fixed_dv = dv_keep_all, new_dmat = dm)
          # if (!is.null(dmat_keep_all)) {
          #   #dmat <- as.matrix(Matrix::bdiag(dmat_keep_all, dm))
          #   dmat <- combine_penalties(fixed_dmat = dmat_keep_all, new_dmat = dm) 
          # } else {
          #   dmat <- dm
          # }
        } else {
          dv <- dv_keep_all
        }
        newx_smooth_track <- if (shpi == 19) rep(1, .nc_safe(newx)) else rep(0, .nc_safe(newx))
      }
    }
    
    gmat <- cbind(one, ddkeep_z_all, newz, ddkeep_x_all, newx)
    np_ix <- 1L + .nc_safe(ddkeep_z_all) + .nc_safe(newz)
    #didn't update dv for nom or non-continuous x, do the update here, just add a zero block, will also add a zero block for continuous x
    #dv <- combine_penalties(fixed_dv = dv, np_ix = np_ix)
    
    #varlst_s <- c(rep(0, np_ix), varlst_s_0, newx_smooth_track)
    varlst_s <- c(rep(0, np_ix), varlst_s_keep, newx_smooth_track)
    
    gmat_lst[[i]] <- gmat
    
    if(!is.null(dv)){
      dv_lst[[i]] <- dv
    }
    np_ix_lst[[i]] <- np_ix
    varlst_s_lst[[i]] <- varlst_s
  }
  
  list(gmat_lst = gmat_lst, dv_lst = dv_lst, np_ix_lst = np_ix_lst, 
       varlst_s_lst = varlst_s_lst, keep_id = keep_id)
}

# ---------------------------------------------------------------------------------------------------------------------
.compute_edf0 <- function(isim, gmat_lst = NULL, dv_lst = NULL, varlst_s_lst = NULL, np_ix_lst = NULL, 
                          pnt = FALSE, shps_to_sim = NULL, type_ix = NULL, runall = TRUE,  
                          family = gaussian(link = 'identity'), is_first_step = FALSE, uinv_lst = NULL) {
  wt.iter <- (family$family == "binomial")
  nobs <- length(gmat_lst)
  edfs <- rep(0, nobs)
  eps <- 1e-8
  obs <- get_obs_indices(nobs, shps_to_sim, type_ix, runall)
  
  for (i in obs) {
    gmat <- gmat_lst[[i]]
    np_ix <- np_ix_lst[[i]]
    varlst_s <- varlst_s_lst[[i]]
    
    #new: speed up the sims; used in !wt.iter case only
    uinv <- uinv_lst[[i]]
    
    if (pnt) {
      dv <- dv_lst[[i]]  # ignored np_ix?
    } else {
      dv <- NULL
    }
    
    #check more 
    #if (NCOL(gmat) == (np_ix + sum(varlst_s))) { #without constr
    #  edfi <- np_ix + sum(varlst_s)
    #}
    
    #if (NCOL(gmat) > (np_ix + sum(varlst_s))) { #with constr
    n <- NROW(gmat)
    ysim <- ysim.fun(n, mu0 = NULL, fml = family$family, shp0 = NULL)
    ans <- fit_model_engine(ysim, current_gmat = gmat, current_dv = dv, family = family, pnt = pnt, np_ix = np_ix, 
                            varlst_s = varlst_s, is_first_step = is_first_step)
    bhat <- ans$bhat
    xtx <- ans$xtx
    amat <- ans$amat
    face <- ans$face
    qv0 <- ans$qv0
    no_constraint <- ans$no_constraint
    if (inherits(ans, "try-error")) {
      edfi <- NULL
    } else {
      if (!pnt || (pnt && all(dv == 0))) {
        edfi <- np_ix + sum(abs(bhat[-c(1:np_ix)]) > eps)
      } else {
        #pre-computed for gaussian; only !wt.iter can use the pre-computed uinv; binomial needs update
        
        if(wt.iter){
          umat <- try(chol(xtx), silent = TRUE)
          if(inherits(umat, 'try-error')){
            ev <- eigen(xtx, only.values = TRUE)$values
            # check if the smallest eigenvalue is too close to zero or negative
            #if (min(ev) < 1e-8) {
            #  diag(xtx) <- diag(xtx) + 1e-8
            #}
            
            # new
            tol <- 1e-8 * max(1, max(abs(ev)))
            if (min(ev) < tol) {
              diag(xtx) <- diag(xtx) + (tol - min(ev))
            }
            
            umat <- chol(xtx)
          }
          uinv <- backsolve(umat, diag(ncol(xtx)))
        }
        
        if(!no_constraint){
          atil <- amat %*% uinv
        }
        xtxinv <- tcrossprod(uinv)
        if(!wt.iter){
          bigp_base <- xtxinv %*% qv0
        } else {
          gmat <- ans$gmat #for binomial; gmat is updated with wt
          bigp_base <- xtxinv %*% crossprod(gmat)
        }
        edfi_base <- sum(diag(bigp_base))
        
        if (all(face == 0) || length(face) == 0 || no_constraint) {
          edfi <- edfi_base
        } else {
          smat <- -t(atil[face, , drop = FALSE])
          qr_smat <- qr(smat)
          smat <- qr.Q(qr_smat)
          pmat_polar <- tcrossprod(smat)
          
          if (!wt.iter) {
            bigp_polar <- uinv %*% tcrossprod(pmat_polar, uinv) %*% qv0
          } else {
            bigp_polar <- uinv %*% tcrossprod(pmat_polar, uinv) %*% crossprod(gmat)
          }
          edfi <- edfi_base - sum(diag(bigp_polar))
        }
      }
    }
    #}
    
    if (!is.null(edfi)) {
      edfs[i] <- edfi
    }
  }
  
  if (!runall) {
    if (!is.null(shps_to_sim) && type_ix != "tree") {
      #check!
      copy_group_edf <- function(edfs, shps_to_sim, ids) {
        pos <- which(shps_to_sim %in% ids)
        if (length(pos) <= 1) return(edfs)
        rep_pos <- pos[1]
        other_pos <- pos[-1]
        edfs[other_pos] <- edfs[rep_pos]
        edfs
      }
      edfs <- copy_group_edf(edfs, shps_to_sim, c(1, 2, 9, 10))   # Mono
      edfs <- copy_group_edf(edfs, shps_to_sim, c(3, 4, 11, 12))  # Conv
      edfs <- copy_group_edf(edfs, shps_to_sim, c(5, 8, 13, 16))  # Combo 1
      edfs <- copy_group_edf(edfs, shps_to_sim, c(6, 7, 14, 15))  # Combo 2
    } else if (is.null(shps_to_sim)) {
      if (nobs >= 4) {
        edfs[2] <- edfs[1]
        edfs[4] <- edfs[3]
        if (nobs > 5) {
          edfs[6:8] <- edfs[5]
        }
      }
    }
  }
  
  edfs
}

# ----------------------------------------------------------------------------------------------------------
compute_edf <- function(gmat_lst = NULL, dv_lst = NULL, varlst_s_lst = NULL,
                        np_ix_lst = NULL, nsim = 200, parallel = FALSE, 
                        pnt = FALSE, shps_to_sim = NULL, type_ix = NULL, edfu_un = NULL, 
                        runall = TRUE, family = gaussian(link = "identity"), 
                        is_first_step = FALSE) {
  if(length(shps_to_sim) == 1 && shps_to_sim == 19){#check more
    dfmean_caic <- edfu_un
    dfmean_cbic <- edfu_un
  } else {
    ng <- length(gmat_lst)
    use_id <- if (!is.null(edfu_un)) which(shps_to_sim != 19) else 1:ng
    #use_id <- 1:ng
    # if (skipcurrent) {
    #   use_id <- use_id[use_id != 1]
    # }
    
    #new: could add more lsts to speed up the sims
    uinv_lst <- vector("list", length = ng)
    
    is_empty_matrix <- function(x) {
      is.matrix(x) && identical(dim(x), c(0L, 0L))
    }
    
    safe_chol <- function(x, eps = 1e-8) {
      out <- try(chol(x), silent = TRUE)
      if (!inherits(out, "try-error")) {
        return(out)
      }
      ev <- eigen(x, only.values = TRUE)$values
      #if (min(ev) < eps) {
      #  diag(x) <- diag(x) + eps
      #}
      tol <- eps * max(1, max(abs(ev)))
      if (min(ev) < tol) {
        diag(x) <- diag(x) + (tol - min(ev))
      }
      chol(x)
    }
    
    uinv_lst <- lapply(seq_len(ng), function(i) {
      gmat <- gmat_lst[[i]]
      dv   <- dv_lst[[i]]
      np_ix <- np_ix_lst[[i]]
      xtx <- crossprod(gmat)
      if(!is.null(dv) && !is_empty_matrix(dv)){
        dv <- combine_penalties(fixed_dv = dv, np_ix = np_ix)
        xtx <- xtx + dv
      }
      umat <- safe_chol(xtx)
      backsolve(umat, diag(ncol(xtx)))
    })
    
    # for(i in 1:ng){
    #   gmat <- gmat_lst[[i]]
    #   dv <- dv_lst[[i]]
    #   xtx <- crossprod(gmat)
    #   nc <- ncol(gmat)
    #   if(pnt){
    #     if(!is.null(dv) && !is_empty_matrix(dv)){ 
    #       xtx <- xtx + dv
    #     }
    #   }
    #   umat <- try(chol(xtx), silent = TRUE)
    #   if(inherits(umat, 'try-error')){
    #     ev <- eigen(xtx, only.values = TRUE)$values
    #     if (min(ev) < 1e-8) {
    #       diag(xtx) <- diag(xtx) + 1e-8
    #     }
    #     umat <- chol(xtx)
    #   }
    #   uinv <- backsolve(umat, diag(ncol(xtx)))
    #   uinv_lst[[i]] <- uinv
    # }
    
    worker <- function(isim) {
      .compute_edf0(isim = isim, gmat_lst = gmat_lst[use_id], dv_lst = dv_lst[use_id],
                    varlst_s_lst = varlst_s_lst[use_id],
                    np_ix_lst = np_ix_lst[use_id], pnt = pnt, 
                    shps_to_sim = shps_to_sim[use_id],
                    type_ix = type_ix, runall = runall, family = family,
                    is_first_step = is_first_step, uinv_lst = uinv_lst)
    }
    
    if (parallel) {
      # let the user choose ncores? 
      ncores <- parallel::detectCores(logical = FALSE)
      cores_avail <- if (is.na(ncores)) 1L else ncores
      #cores <- max(1L, min(4L, cores_avail - 1L))
      #new: let the user maker changes about how many cores
      user_cores <- .get_cgam_option("cores")
      
      if (!is.numeric(user_cores) || user_cores < 1L) {
        stop("'cgam.cores' must be a positive integer.")
      }
      
      user_cores <- as.integer(user_cores)
      cores <- min(user_cores, cores_avail)
      
      # new: allow windows and Mac
      is_windows = .Platform$OS.type == "windows"
      if(is_windows){
        # windows: use parLapply with a PSOCK cluster
        cl <- parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl))  # ensure cleanup
        edfs <- parallel::parLapply(cl, seq_len(nsim), worker)
      } else {
        edfs <- suppressWarnings(parallel::mclapply(seq_len(nsim), worker, mc.cores = cores))
      }
    
      #if the parallel code not work, then try lapply
      if (all(vapply(edfs, is.null, logical(1)))) {
        edfs <- lapply(seq_len(nsim), worker)
      }
    } else {
      edfs <- lapply(seq_len(nsim), worker)
    }
    
    edfs <- Filter(Negate(is.null), edfs)
    edfs <- simplify2array(edfs)
    #nobs <- length(gmat_lst[use_id])
    
    dfmean_cbic <- if (inherits(edfs, "matrix")) rowMeans(edfs) else mean(edfs)
    #check more
    #if (!is.null(edfu_un) && length(use_id) < ng) {
    if (!is.null(edfu_un)) {
      dfmean_cbic <- c(dfmean_cbic, edfu_un)
    } #else if (!is.null(edfu_un) && length(use_id) == ng) {
    #stop('check edfu_un!')
    #}
    
    #if (is.null(shp17_id)) {
    dfmean_caic <- max(dfmean_cbic)
    # } else {
    #   # if (skipcurrent) {
    #   #   shp17_id <- shp17_id - 1
    #   # }
    #   dfmean_caic_no17 <- max(dfmean_cbic[-shp17_id])
    #   dfmean_caic_17 <- dfmean_cbic[shp17_id]
    #   dfmean_caic <- c(dfmean_caic_no17, dfmean_caic_17)
    # }
  }
  rslt <- list(dfmean_caic = dfmean_caic, dfmean_cbic = dfmean_cbic)
  rslt
}

# ------------------------------------------------------------
get_obs_indices <- function(nobs, shps_to_sim, type_ix, runall) {
  if (runall) return(1:nobs)
  if (!is.null(type_ix) && type_ix == 'tree') return(1:nobs)
  
  # grouping logic for shape simulations
  get_first <- function(ids) {
    pos <- which(shps_to_sim %in% ids)
    if (length(pos) > 0) pos[1] else NULL
  }
  
  obs <- c(
    get_first(c(1, 2, 9, 10)),   # Mono
    get_first(c(3, 4, 11, 12)),  # Conv
    get_first(c(5, 8, 13, 16)),  # Combo 1
    get_first(c(6, 7, 14, 15)),  # Combo 2
    which(shps_to_sim == 19),   # Unconstr
    which(shps_to_sim %in% c(0, 17)) # Flat
  )
  return(as.vector(na.omit(obs)))
}

#--------------------------------------------------------------------------------
zvec.fun <- function(cvec = NULL, wt = NULL, y, sm = 1e-7, fml = object$family) {
  n <- length(y)
  if (fml == "gaussian") {
    #zvec <- y
    zvec <- wt^(1/2) * y
  }
  if (fml == "poisson") {
    #zvec <- cvec / wt
    zvec <- cvec / sqrt(wt)
  }
  if (fml == "binomial") {
    zvec = 1:n*0
    zvec[wt == 0] <- 1 / sm
    zvec[wt > 0] <- cvec[wt > 0] / sqrt(wt[wt > 0])
  }
  if (fml == "Gamma") {
    zvec <- cvec / sqrt(wt)
  }
  zvec
}

wt.fun <- function(y, etahat = NULL, n = NULL, weights = NULL, fml = object$family){
  if (is.null(weights)) {
    weights <- 1:n*0 + 1
  }
  w <- weights
  if (fml == "poisson") {
    wt <-  w * exp(etahat)
  }
  if (fml == "binomial") {
    if (all(etahat == 0)){
      #wt <- 1:n*0 + 1/4
      wt <- w * (1:n*0 + 1/4)
    } else {
      wt <- 1:n*0
      wt <- (1 - 1 / (1 + exp(etahat))) * (1 / (1 + exp(etahat))) #pi_ij(1 - pi_ij)
      wt <- c(wt)
      # for (i in 1:n) {
      #   if (etahat[i] > 100) {
      #     wt[i] <- 0
      #   } else {
      #     wt[i] <- w[i] * exp(etahat[i]) / ((1 + exp(etahat[i]))^2)
      #   }
      # }
    }
  }
  if (fml == "gaussian") {
    wt <- w # (1:n*0 + 1) / w
  }
  if (fml == "Gamma") {
    wt <- w * y * exp(-etahat)
  }
  wt <- as.vector(wt)
  wt
}

etahat.fun <- function(n, y, fml = object$family){
  if (fml == "poisson") {
    etahat <- 1:n*0 + log(mean(y))
  }
  if (fml == "binomial") {
    etahat <- 1:n*0
  }
  if (fml == "Gamma") {
    etahat <- 1:n*0 + log(mean(y))
  }
  etahat
}

gr.fun <- function(y, etahat = NULL, weights = NULL, fml = object$family){
  n <- length(y)
  if (is.null(weights)) {
    weights <- 1:n*0 + 1
  }
  w <- weights
  if (fml == "poisson") {
    gr <- w * (exp(etahat) -  y)
  }
  if (fml == "binomial") {
    if (all(etahat == 0)) {
      gr <- w * (1/2 - y)
    } else {
      gr <- 1:n*0
      #ignore weights
      gr <- 1 - 1 / (1 + exp(etahat)) - y
      
      # for (i in 1:n) {
      #   if (etahat[i] > 100) {
      #     gr[i] <- w[i] * (1 - y[i])
      #   } else {gr[i] <- w[i] * (exp(etahat[i]) / (1 + exp(etahat[i])) - y[i])}
      # }
    }
  }
  if (fml == "Gamma") {
    gr <- w * (1 - y * exp(-etahat))
  }
  gr
}

ysim.fun <- function(n, mu0 = NULL, fml = object$family, shp0 = NULL) {
  if (fml == "binomial") {
    ysim <- 1:n*0
    ysim[runif(n) < .5] <- 1
  }
  if (fml == "poisson") {
    if (!is.null(mu0)) {
      ysim <- rpois(n, mu0)
    }
  }
  if (fml == "gaussian") {
    ysim <- rnorm(n)
  }
  if (fml == "Gamma") {
    ysim <- rgamma(n, shape=1)
  }
  ysim
}

llh.fun <- function(y, muhat = NULL, etahat = NULL, phihat = NULL, n = NULL, weights = NULL, fml = object$family){
  #sm <- 1e-7
  sm <- 1e-5
  if (is.null(weights)) {
    weights <- 1:n*0 + 1
  }
  w <- weights
  #new: avoid Inf
  if (fml == "poisson") {
    llh <- 2 * sum(w[w!=0] * (muhat[w!=0] - y[w!=0] * etahat[w!=0])) / n
  }
  if (fml == "binomial") {
    llh <- 0
    if (all(0 <= y) & all(y <= 1)) {
      for (i in 1:n) {
        if (muhat[i] > 0 & muhat[i] < 1) {
          llh <- llh + w[i] * (y[i] * log(muhat[i]) + (1 - y[i]) * log(1 - muhat[i]))
        }
      }
      #ch <- muhat > 0 & muhat < 1
      #llh <- sum(w[ch] * (log(1 - muhat[ch]) - y[ch] * etahat[ch]))
      llh <- (-2/n) * llh
    } else {
      stop ("y values must be 0 <= y <= 1!")
    }
  }
  if (fml == "gaussian") {
    if (all(w == 1)) {
      llh <- log(sum((y - etahat)^2))
    } else {
      llh <- log(sum(w[w!=0] * (y[w!=0] - etahat[w!=0])^2)) - sum(log(w[w!=0])) / n
    }
  }
  if (fml == "Gamma") {
    vuhat <- 1 / phihat
    #print (vuhat)
    #llh <- 2 * sum(w[w!=0] * (etahat[w!=0] + y[w!=0] * exp(-etahat[w!=0]))) / n
    llh <- 2 / n * (vuhat * sum(w[w!=0] * (etahat[w!=0] + y[w!=0] * exp(-etahat[w!=0]))) + n * (log(gamma(vuhat)) - vuhat * log(vuhat)) - (vuhat-1) * sum(log(y)))
    #print (llh)
  }
  llh
}

#------------------------------------------------------
#visualize each chosen predictor's etahat or muhat
#------------------------------------------------------
plot.shapeselect <- function(x,...) {
  
  object <- x
  pages <- 1
  mn <- TRUE
  mu <- FALSE
  #use_ggplot <- FALSE
  
  layout_for_cgam <- function(n_plots, pages = 1L) {
    ppp <- ceiling(n_plots / max(pages, 1L))
    nr <- floor(sqrt(ppp))
    nc <- ceiling(ppp / nr)
    c(nr, nc)
  }
  
  extras <- list(...)
  rug_on <- if (is.null(extras$rug)) TRUE else isTRUE(extras$rug)

  xmat <- object$xmat_0
  p <- ncol(xmat)
  types <- object$types
  
  xmat <- as.data.frame(xmat)
  xmat[] <- lapply(seq_len(p), function(ix) {
    #col <- xmat[[ix]]
    col <- xmat[,ix]
    if (inherits(col, "factor") && types[ix] %in% c("tree", "ord")) {#assume tree/ord only use integers
      return(as.numeric(levels(col))[col])
    }
    col
  })
  
  family <- object$family
  etacomps <- object$etacomps
  shps <- object$shps
  xnms <- object$xnms
  
  zmat <- object$ddkeep_z_all
  np <- object$np
  varlist <- sort(object$varlst)
  varlist_z <- sort(object$varlst_z)
  dd <- object$bigmat
  
  kept <- !shps %in% c(0, 17)
  xmat <- xmat[, kept, drop = FALSE]
  types <- types[kept]
  etacomps <- etacomps[kept, , drop = FALSE]
  shps <- shps[kept]
  xnms <- xnms[kept]
  shapes_char <- object$shapes[kept]
  
  p <- ncol(xmat)
  if (p == 0L) {
    warning("No selected predictors to plot.")
    return(invisible(NULL))
  }
  
  uvar_x <- unique(varlist)
  uvar_z <- unique(varlist_z)
  
  ## ------------------------------------------------------------------
  ## base graphics path
  ## ------------------------------------------------------------------
  if (pages == 1L) {
    oldpar <- par(mfrow = layout_for_cgam(n_plots = p, pages = pages))
    on.exit(par(oldpar), add = TRUE)
  }
  
  ziter <- 1L
  for (k in seq_len(p)) {
    #x_k <- xmat[[k]]
    x_k <- xmat[,k]
    xnm_k <- xnms[k]
    shp_k <- shapes_char[k]
    eta_k <- etacomps[k, ]
    type_k <- types[k]
    
    if (mu) {
      eta_k <- family$linkinv(as.numeric(eta_k))
    }
    
    #if (is.factor(x_k)) {
    if(type_k == "nom" || (type_k == "tree" && shp_k == "unordered")){  
      ps_k <- which(varlist_z == uvar_z[ziter]) + 1L
      ziter <- ziter + 1L
    } else {
      ps_k <- which(varlist == uvar_x[k]) + np
    }
    
    range_eta <- range(eta_k, na.rm = TRUE)
    d_eta <- diff(range_eta)
    padding <- if (is.finite(d_eta) && d_eta > 0) 0.05 * d_eta else 0.05
    ylim_use <- c(range_eta[1] - padding, range_eta[2] + padding)
    
    ord_k <- order(x_k)
    tp_k <- if (is.factor(x_k)) "p" else "l"
    xaxt_k <- if (is.factor(x_k)) "n" else "s"
    
    x_k_tmp <- x_k
    if (is.factor(x_k_tmp)) {
      #x_k_tmp <- as.integer(x_k_tmp)
      lvls <- as.numeric(levels(x_k_tmp))
      if(!any(is.na(lvls))){
        x_k_tmp <- as.numeric(levels(x_k_tmp))[x_k_tmp]
      } else {
        x_k_tmp <- as.integer(x_k_tmp)
      }
    }
    
    if (mn) {
      plot(x_k_tmp[ord_k], eta_k[ord_k], xlab = bquote(~ .(xnm_k)),
           ylab = bquote(hat(eta)[.(k)]),  ylim = ylim_use,
           main = bquote(~ hat(eta)[.(k)] ~ "vs" ~ .(xnm_k) ~ "(" * .(shp_k) * ")"),
           xaxt = xaxt_k, type = tp_k, lwd = 1)
    } else {
      plot(x_k_tmp[ord_k], eta_k[ord_k], xlab = bquote(~ .(xnm_k)),
           ylab = "", ylim = ylim_use, main = NULL,
           xaxt = xaxt_k, type = tp_k, lwd = 1)
    }
    
    if (rug_on) {
      rug(x_k_tmp, ticksize = 0.03, side = 1)
    }
    
    if (is.factor(x_k) || types[k] == "tree") {
      points(x_k_tmp, eta_k, pch = 20)
      if (is.factor(x_k)) {
        lines(x_k_tmp[ord_k], eta_k[ord_k], lty = 2)
        axis(1, at = unique(x_k_tmp))
        #axis(1, at = unique(x_k_tmp), labels = lvls) #temp
        #axis(1, at = unique(x_k_tmp), labels = unique(levels(x_k)))
      }
    }
    
    if (pages > 1L) {
      readline("Press <Enter> to continue...")
    }
  }
  
  invisible(NULL)
}

#----------------------------------------------------------------------------------------------------------------------------
#check!
#create penalty matrix
#difference from quadratic or cubic splines
#for i-spl and c-spl
#2nd deri of c-spl is 1st deri of i-spl
#columns are basis at knots
#each row of dmat is: 2nd deriv at the next basis - 2nd deriv at the previous basis, for each  knot, rows are for knots
#only need to compute 9, 11, 13; we only use crossprod(dmat)
#----------------------------------------------------------------------------------------------------------------------------
makedmat_1D_ic <- function(m, knots, shpi = 9) {
  if(shpi < 13) {
    m_use <- m
  } else if (shpi >= 13) {#include x; append a zero column later
    m_use <- m - 1
  }
  # dmat <- matrix(0, nrow = (m_use - 1), ncol = m_use)
  # 
  # nkts <- length(knots)
  # if(shpi %in% c(9, 10)){
  #   vals0 <- iSpline(knots, knots = knots[-c(1, nkts)], Boundary.knots = c(knots[1], knots[nkts]), 
  #                    degree = 1, derivs = 1L, scale = FALSE)
  # }
  # 
  # if(shpi >= 11 && shpi <= 16 || shpi == 19){
  #   vals0 <- cSpline(knots, knots = knots[-c(1, nkts)], Boundary.knots = c(knots[1], knots[nkts]), 
  #                    degree = 1, derivs = 2L, scale = FALSE)
  # }
  # 
  # vals <- diag(vals0)
  # 
  # for (i in 2:m_use) {
  #   dmat[i-1, i-1] <- -1 * vals[i-1]; dmat[i-1, i] <- 1 * vals[i]
  # }
  # 
  # if(shpi >= 13) {
  #   dmat <- cbind(dmat, 0) #last column is for x
  # }
  ##new: use consecutive differences of 2nd derivatives for I-spl; 3rd for C-spl
  nkts <- length(knots)
  midk <- 1:(nkts - 1)
  for(i in 1:(nkts - 1)){midk[i] <- (knots[i] + knots[i+1]) / 2}
  diffmat <- matrix(0, nrow = (nkts - 2), ncol = (nkts - 1))
  for(i in 1:(nkts - 2)){
    diffmat[i, i] <- -1
    diffmat[i, i+1] <- 1
  }
  if(shpi %in% c(9, 10)){
    dmat0 <- iSpline(midk, knots = knots[-c(1, nkts)], Boundary.knots = c(knots[1], knots[nkts]), degree = 1, derivs = 2, scale = FALSE)
    dmat0 <- dmat0 
  }
  if(shpi >= 11 && shpi <= 16 || shpi == 19){
    dmat0 <- cSpline(midk, knots = knots[-c(1, nkts)], Boundary.knots = c(knots[1], knots[nkts]), degree = 1, derivs = 3, scale = FALSE)
    dmat0 <- dmat0 
  }
  ##test
  dmat0 <- dmat0 / max(abs(dmat0))
  dmat <- diffmat %*% dmat0
  if(shpi >= 13) {
    dmat <- cbind(dmat, 0) #last column is for x
  }
  return (dmat)
}

#------------------------------------------
#need a better way to find ps....
#------------------------------------------
.search_ps_sel <- function(pen, qv0, dv, dd, edfu, weights = NULL)
{
  #new:
  if (!is.null(weights)) {
    stopifnot(all(is.finite(weights)), all(weights >= 0))
    
    if (any(weights != 1)) {
      dd <- sweep(dd, MARGIN = 1, sqrt(weights), FUN = '*')
      qv0 <- crossprod(dd)
    }
  }
  
  qv <- qv0 + pen * dv
  umat <- chol(qv)
  qvinv <- chol2inv(umat)
  val <- sum(diag(qvinv %*% qv0)) - 1*edfu
  #val <- sum(diag(qvinv %*% crossprod(dd))) - 1*edfu
  return (val)
}

.compute_edfu_un <- function(gmat, dm, np_ix){
  qv0 <- crossprod(gmat)
  xtx <- qv0
  #check!
  dmat_zero <- matrix(0, nrow = nrow(dm), ncol = np_ix)
  dmat2 <- cbind(dmat_zero, dm)
  dv <- crossprod(dmat2)
  xtx <- qv0 + dv
  umat <- chol(xtx)
  qvinv <- chol2inv(umat)
  edfu_un <- sum(diag(qvinv %*% qv0))
  return (edfu_un)
}

##############################
#tranform shapes back to char
##############################
ShapeToChar_shapeselect <- function(shp, tag = "x") {
  tag <- as.character(tag)
  if (length(shp) != 1 || is.na(shp)) {
    stop("`shp` must be a single non-missing value.")
  }
  if (length(tag) != 1 || is.na(tag)) {
    stop("`tag` must be a single non-missing character value.")
  }
  shp <- as.integer(shp)
  
  shape_map <- list(
    x = c("incr", "decreasing", "convex", "concave",
          "increasing convex", "decreasing convex", "increasing concave", "decreasing concave",
          "smooth increasing", "smooth decreasing", "smooth convex", "smooth concave",
          "smooth increasing convex", "smooth increasing concave", "smooth decreasing convex",
          "smooth decreasing concave", "flat", "linear","smooth unconstrained","out of the model"),
    z = c("in the model", "out of the model"),
    tree = c("tree","unordered", "out of the model")
  )
  
  default_shape <- c(
    x = 20L,
    z = 2L,
    tree = 3L
  )
  
  if (!tag %in% names(shape_map)) {
    stop("`tag` must be one of 'x', 'z', or 'tree'.")
  }
  
  if (shp == 0L) {
    shp <- default_shape[[tag]]
  }
  
  labels <- shape_map[[tag]]
  
  if (shp < 1L || shp > length(labels)) {
    stop(sprintf(
      "Invalid shape code %d for tag '%s'. Valid values are 0 to %d.",
      shp, tag, length(labels)
    ))
  }
  
  labels[[shp]]
}