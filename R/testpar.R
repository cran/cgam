#' Parametric Monotone/Quadratic vs Smooth Constrained Test
#'
#' Performs a hypothesis test comparing a parametric model (e.g., linear, quadratic, warped plane) 
#' to a constrained smooth alternative under shape restrictions.
#'
#' @param formula0 A model formula representing the null hypothesis (e.g., linear or quadratic).
#' @param formula A model formula under the alternative hypothesis, possibly with shape constraints.
#' @param data A data frame or environment containing the variables in the model.
#' @param family A description of the error distribution and link function to be used in the model (e.g., \code{gaussian()}). Gaussian and binomial are now included.
#' @param ps User-defined penalty term. If \code{NULL}, optimal values are estimated.
#' @param edfu User-defined unconstrained effective degrees of freedom. 
#' @param nsim Number of simulations to perform to get the mixture distribution of the test statistic. (default is 1000).
#' @param multicore Logical. Whether to enable parallel computation (default uses global option \code{cgam.parallel}).
#' @param method A character string (default is \code{"testpar.fit"}), currently unused.
#' @param arp Logical. If \code{TRUE}, uses autoregressive structure estimation.
#' @param p Integer vector or value specifying AR(p) order. Ignored if \code{arp = FALSE}.
#' @param space Character. Either \code{"Q"} for quantile-based knots or \code{"E"} for evenly spaced knots. Default is \code{"E"}. 
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing test statistics, fitted values, coefficients, and other relevant objects.
#'
#' @examples
#' 
#' # Example 1: Linear vs. monotone alternative
#' set.seed(123)
#' n <- 100
#' x <- sort(runif(n))
#' y <- 2 * x + rnorm(n)
#' dat <- data.frame(y = y, x = x)
#' ans <- testpar(formula0 = y ~ x, formula = y ~ s.incr(x), multicore = FALSE, nsim = 10)
#' ans$pval
#' summary(ans)
#'\dontrun{
#' # Example 2: Binary response: linear vs. monotone alternative
#' set.seed(123)
#' n <- 100
#' x <- runif(n)
#' eta <- 8 * x - 4
#' mu <- exp(eta) / (1 + exp(eta))
#' y <- rbinom(n, size = 1, prob = mu)
#'
#' ans <- testpar(formula0 = y ~ x, formula = y ~ s.incr(x), 
#' family = binomial(link = "logit"), nsim = 10)
#' summary(ans)
#' 
#' xs = sort(x)
#' ord = order(x)
#'  par(mfrow = c(1,1))
#'  plot(x, y, cex = .2)
#'  lines(xs, binomial(link="logit")$linkinv(ans$etahat0)[ord], col=4, lty=4)
#'  lines(xs, binomial(link="logit")$linkinv(ans$etahat)[ord], col=2, lty=2)
#'  lines(xs, mu[ord])
#'  legend("topleft", legend = c("H0 fit (etahat0)", "H1 fit (etahat)", "True mean"),
#'  col = c(4, 2, 1), lty = c(4, 2, 1), bty = "n", cex = 0.8)
#' }
#' 
#' # Example 3: Bivariate warped-plane surface test
#' set.seed(1)
#' n <- 100
#' x1 <- runif(n)
#' x2 <- runif(n)
#' mu <- 4 * (x1 + x2 - x1 * x2)
#' eps <- rnorm(n)
#' y <- mu + eps
#'
#' # options(cgam.parallel = TRUE) #allow parallel computation
#' ans <- testpar(formula0 = y ~ x1 * x2, formula = y ~ s.incr.incr(x1, x2), nsim = 10)
#' summary(ans)
#' 
#' par(mfrow = c(1, 2)) 
#' persp(ans$k1, ans$k2, ans$etahat0_surf, th=-50, main = 'H0 fit')
#' persp(ans$k1, ans$k2, ans$etahat_surf, th=-50, main = 'H1 fit')
#'
#'\dontrun{
#' # Example 4: AR(1) error, quadratic vs smooth convex
#' set.seed(123)
#' n = 100
#' x = runif(n)
#' mu = 1+x+2*x^2
#' eps = arima.sim(n = n, list(ar = c(.4)), sd = 1)
#' y = mu + eps
#' ans = testpar(y~x^2, y~s.conv(x), arp=TRUE, p=1, nsim=10)
#' ans$pval
#' ans$phi #autoregressive coefficient estimate
#' }
#' \dontrun{
#' # Example 5: Three additive components with shape constraints
#' n <- 100
#' x1 <- runif(n)
#' x2 <- runif(n)
#' x3 <- runif(n)
#' z <- rnorm(n)
#' mu1 <- 3 * x1 + 1
#' mu2 <- x2 + 3 * x2^2
#' mu3 <- -x3^2
#' mu <- mu1 + mu2 + mu3
#' y <- mu + rnorm(n)
#'
#' ans <- testpar(formula0 = y ~ x1 + x2^2 + x3^2 + z, 
#' formula = y ~ s.incr(x1) + s.incr.conv(x2) + s.decr.conc(x3)
#'  + z, family = "gaussian", nsim = 10)
#' ans$pval
#' summary(ans)
#' }
#' @export
testpar <- function(formula0, formula, data, family = gaussian(link = "identity"), 
                    ps = NULL, edfu = NULL, 
                    nsim = 1000, multicore = getOption("cgam.parallel"), method = "testpar.fit",
                    arp = FALSE, p = NULL, space = "E",...) {
  cl <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family))
    stop("'family' not recognized!")
  mf <- match.call(expand.dots = FALSE)
  if(missing(data))
    data <- environment(formula)
  m <- match(c("formula0", "formula", "data"), names(mf), 0L)
  #print (m)
  #----------------------------------------------------------------------------------------------------------------------
  #formula0 is not very useful; we just need it to see if H0 is linear or quadratic, or a linear plane/interaction plane?
  #----------------------------------------------------------------------------------------------------------------------
  form0_expr <- deparse(mf[c(1L, m[1])])
  expr_inside <- sub(".*\\((.*)\\)", "\\1", form0_expr)
  expr_inside <- sub(".*?\\((.*)\\)", "\\1", form0_expr)
  expr_parts <- strsplit(expr_inside, " \\+ ")[[1]]
  expr_parts <- lapply(expr_parts, function(x) {
    x <- trimws(x)
    # Remove outer parentheses ONLY if they wrap the entire string
    if (grepl("^\\(.*\\)$", x)) {
      x <- substr(x, 2, nchar(x) - 1)
    }
    trimws(x)
  })
  
  parametric <- lapply(expr_parts, function(x) {
    if (grepl("\\*", x) && grepl("\\^\\d+", x)) {
      degree <- as.numeric(sub(".*\\^(\\d+).*", "\\1", x))
      if (!is.na(degree) && degree >= 2) {
        warning("Interaction terms (*) with polynomial degrees >= 2 are not allowed in a warped-plane test. Reduced the degree to be 1!", call. = FALSE)
      }
    }
    
    dplyr::case_when(grepl("\\*", x) ~ "warped_plane", 
                     grepl("\\^\\d+", x) &
                       ifelse(grepl("\\^\\d+", x),
                              as.numeric(sub(".*\\^(\\d+).*", "\\1", x)) >= 4, FALSE) ~ "polynomial higher than cubic",
                     grepl("\\^3", x) ~ "cubic", 
                     grepl("\\^2", x) ~ "quadratic", 
                     grepl("\\^1", x) ~ "linear", 
                     grepl("^factor\\(", x) ~ "linear",
                     TRUE ~ "linear")}) |> unlist()
  
  #cat("parametric = ", '\n')
  #print (parametric)
  #new: handle x + I(x^2) + I(x^3); y ~ x + x^2, keep x^2
  #levels <- c("linear", "quadratic", "cubic", "polynomial higher than cubic")
  #if(length(parametric) > 1) {
  #  factor_param <- factor(parametric, levels = levels, ordered = TRUE)
  #  parametric <- as.character(max(factor_param))
  #}
  #cat("parametric = ", '\n')
  #print (parametric)
  form_expr <- deparse(mf[c(1L, m[2])])
  #print (form_expr)
  #remove z here?
  #mf0 <- mf
  #----------------------------------------------------------
  #check attributes in cgam formula; get shapes and x's and y
  #----------------------------------------------------------
  mf <- mf[c(1L, m[-1])] #used formula 
  #print (head(mf))
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  #print (head(mf))
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  #------------------------
  #extract x, y, z, etc...
  #------------------------
  #mf is a data frame and can be treated as a list
  mf_xz <- mf[, -1, drop = FALSE]
  x_or_z <- lapply(mf_xz, function(e) attr(e, "categ"))
  #print (x_or_z)
  #print (names(x_or_z))
  is_z <- sapply(x_or_z, is.null)
  #print (is_z)
  nvars <- nz <- sum(is_z)
  zmat <- x <- shp <- shp_pr <- NULL
  #print (is_z)
  #print (parametric)
  if(any(is_z)){
    #temp!
    #parametric <- parametric[-which(is_z)] #wrong when we have y~x1+x2+x3+x4+z, y~x1+x2+s.incr.incr(x3,x4)+z
    #print (parametric)
    znm_in_form <- names(x_or_z)[which(is_z)]
    z_ps_in_param <- sapply(expr_parts, function(e) e == znm_in_form)
    parametric <- parametric[-which(z_ps_in_param)]
    #print (parametric)
    z <- mf_xz[, is_z, drop = FALSE]
    is_fac <- sapply(z, function(e) is.factor(e))
    #print (which(is_fac))
    #print (class(z))
    #print (head(z[, which(is_fac), drop = FALSE]))
    if(any(is_fac)){
      #test more
      if(sum(!is_fac) > 0){
        zmat <- as.matrix(z[, !is_fac, drop = FALSE])
      } else {zmat <- NULL}
      #print (zmat)
      #add znm later
      #not work:
      #dd <- model.matrix(~ z[, which(is_fac)])[, -1, drop = FALSE]
      #vars_id <- which(is_fac)
      dd <- model.matrix(~ ., data = z[, is_fac, drop = FALSE])[, -1]
      #print (head(dd))
      zmat <- cbind(zmat, dd)
    } else {
      zmat <- as.matrix(z)
    }
    #print (class(zmat))
  }
  # nx <- 0
  # print (parametric)
  if(any(!is_z)){
    #print (parametric)
    x <- mf_xz[, !is_z, drop = FALSE]
    #print (colnames(x))
    add_or_wps <- lapply(x, function(e) attr(e, "categ")) |> simplify2array()
    nx <- length(add_or_wps)
    X <- vector("list", length = nx)
    #print (add_or_wps)
    #print (nx)
    #test! for y~x1+x2+x3*x4, y~s.incr.incr(x1,x2)+s.incr.incr(x3,x4)
    #if(any(add_or_wps == "warp") & any(add_or_wps == "additive")){
    #print (parametric)
    if(any(add_or_wps == "warp")){
      #print (nx)
      parametric2 <- rep("linear", length = nx)
      rp_ps <- which(parametric != "linear")
      #print (rp_ps)
      #rp_ps <- rp_ps[rp_ps <= nx]
      rp_ps0 <- rp_ps
      if(any(rp_ps > nx)){
        rp_ps <- rp_ps - (rp_ps - nx)
      }
      #print (rp_ps)
      parametric2[rp_ps] <- parametric[rp_ps0]
      #parametric2 <- unique(parametric)
      #print (parametric2)
      if(any(parametric == "warped_plane")){
        wp_ps <- which(parametric == "warped_plane")
        #print (wp_ps)
        if(length(wp_ps) == 1){
          if(wp_ps > nx){
            shift <- max(wp_ps) - nx
            parametric2[wp_ps - shift] <- "warped_plane"
          } else {
            #print (wp_ps)
            parametric2[wp_ps] <- "warped_plane"
          }
        }
        if(length(wp_ps) > 1){
          diff_wp_ps <- diff(wp_ps)
          #print (diff_wp_ps)
          if(any(diff_wp_ps > 2)) {
            diff_wp_ps[which(diff_wp_ps > 2)] <- 2
          }
          min_wp_ps <- wp_ps[1]
          #print (min_wp_ps)
          new_wp_ps <- c(min_wp_ps, rep(min_wp_ps, length(wp_ps)-1) + diff_wp_ps)
          #print (new_wp_ps)
          #print (new_wp_ps > nx)
          if(any(new_wp_ps > nx)){
            shift <- max(new_wp_ps) - nx
            new_wp_ps <- new_wp_ps - shift
          }
          parametric2[new_wp_ps] <- "warped_plane"
        }
      }
      parametric <- parametric2
    }
    #print (parametric)
    
    #new: need to match parametric and type first
    type <- rep("monotone", nx)
    for(ix in 1:nx){
      #for(ix in nx:1){
      X[[ix]] <- x[, ix]
      #attr(X[[ix]], "shp") <- attr(x[, ix], "shape")
      shape_ix <- attr(X[[ix]], "shape") 
      if(add_or_wps[ix] == "additive"){
        #if(shape_ix %in% c(11, 12)){
        if(!shape_ix %in% c(9, 10)){
          attr(X[[ix]], "type")  <- "convex"
          #attr(X[[ix]], "parametric") <- "quadratic"
          #test!
          attr(X[[ix]], "parametric") <- parametric[ix]
        } else if(shape_ix %in% c(9, 10)){
          attr(X[[ix]], "type")  <- "monotone"
          #attr(X[[ix]], "parametric") <- "linear"
          attr(X[[ix]], "parametric") <- parametric[ix]
        }
      }
      if(add_or_wps[ix] == "warp"){
        attr(X[[ix]], "type")  <- "monotone_plane"
        #temp
        attr(X[[ix]], "parametric") <- parametric[ix]
        #print (parametric[ix])
      } 
      type[ix] <- attr(X[[ix]], "type")
    }
  }
  # cat("parametric = ", '\n')
  # print (parametric)
  # 
  # cat("type = ", '\n')
  # print (type)
  X <- match_H0_H1(X)
  #----------------------------------------------------------
  #call testpar.fit
  #----------------------------------------------------------
  #print (identical(parent.frame(), environment()))
  #print (ls(parent.frame()))
  #seems to be worse than just call testpar.fit?
  #fit <- eval(call("testpar.fit", x=x, y=y, zmat=zmat, family=family, nsim=nsim, parametric=parametric, type=type, arp=arp, p=p, space=space), envir = environment())
  fit <- testpar.fit(X=X, y=y, zmat=zmat, family=family, ps=ps, edfu=edfu, nsim=nsim, multicore=multicore,
                     GCV=FALSE, lams=NULL, arp=arp, p=p, space=space, parametric=parametric)
  rslt <- structure(c(fit, list(X = X, expr_parts = expr_parts, zmat = zmat, call = cl, formula0 = formula0, formula = formula, terms = mt, data = data,
                                parametric=parametric, type=type, family = family)))
  class(rslt) <- c("testpar", "cgam")
  return(rslt)
}


#----------------------------------------------------------
#check if the user specifies the right H0 and H1 symbols
#----------------------------------------------------------
match_H0_H1 <- function(X){
  nx <- length(X)
  for(ix in 1:nx){
    x <- X[[ix]]
    if(NCOL(x) == 1){
      # if(attr(x, "parametric") == "linear" & attr(x, "type") == "convex") {
      #   attr(x, "parametric") <- "quadratic"
      #   X[[ix]] <- x
      #   warning(paste(attr(x, "nm"), "is defined with dgree = 1 in H0 but the test is for convexity. Not work. Increased the degree to be quadratic!"), call.=FALSE)
      # }
      if(attr(x, "parametric") == "cubic" & attr(x, "type") == "monotone") {
        attr(x, "parametric") <- "quadratic"
        X[[ix]] <- x
        warning(paste(attr(x, "nm"), "is defined with dgree = 3 in H0 but the test is for monotonicity. Not work. Reduced the degree to be quadratic!"), call.=FALSE)
      }
      if(attr(x, "parametric") == "polynomial higher than cubic"){
        attr(x, "parametric") <- "quadratic"
        X[[ix]] <- x
        warning(paste(attr(x, "nm"), "is defined with dgree > 3 in H0. Not work. Reduced the degree to be quadratic!"), call.=FALSE)
      }
    }
    if(NCOL(x) == 2){
      if(attr(x, "parametric") %in% c("quadratic", "cubic", "polynomial higher than cubic")){
        attr(x, "parametric") <- "linear"
        X[[ix]] <- x
        #temp
        warning(paste("H1 is for", attr(x, "type"), "test. Terms in H0 cannot have a degree > 1. All terms are changed to be linear!"), call.=FALSE)
      }
    }
  }
  return(X)
}

#-----------------------------------------------------
#similar to glm.fit
#X is a list with attributes: shp, parametric, type
#-----------------------------------------------------
testpar.fit = function(X, y, zmat = NULL, family = gaussian(link="identity"), 
                       ps = NULL, edfu = NULL, nsim = 1000, multicore = TRUE, 
                       GCV = FALSE, lams = NULL, 
                       arp = FALSE, p = NULL, space = "E",...)
{
  wt.iter = ifelse(family$family == "gaussian", FALSE, TRUE)
  extras = list(...)
  #new:
  #capl = NCOL(x)
  capl = length(X)
  n = length(y)
  k1 = k2 = NULL
  p_optim = 0
  if(arp){
    #if the user didn't provide p, then choose p from 0:3
    if(is.null(p)){
      p = 0:4 
    } else if (p < 0 || p > 4) {
      warning("The order in AR(p) must be an integer >= 0 and <= 4!")
      p = 1
    }
  } else {
    p = 0
  }
  #-------------------
  #get xm, amat, dmat
  #-------------------
  mult = 2
  nz = 0
  if(!is.null(zmat)){
    nz = NCOL(zmat)
  }
  xm0 = NULL #combine dd columnwise
  amat_lst = vector("list", length = capl)
  awmat_lst = vector("list", length = capl)
  dmat_lst = vector("list", length = capl)
  dd_lst = kts_lst = vector("list", length = capl)
  edfu_vec = numeric(capl)
  var_track = var_track_param = NULL
  iadd = ipr = 0
  use_constreg = FALSE
  for(icomp in 1:capl){
    x = X[[icomp]] #|> as.matrix()
    if(NCOL(x) == 1){
      iadd = iadd + 1
      xu = unique(x)
      n1 = length(xu)
      type = attr(x, "type")
      parametric = attr(x, "parametric")
      #new:
      shp = attr(x, "shape")
      if(parametric == "linear"){
        use_constreg = TRUE
      } 
      if(parametric == "quadratic"){
        use_constreg = TRUE
      } 
      if(parametric == "cubic" & type == "convex"){#otherwise give a warning
        use_constreg = TRUE
      }
      #nkts = mult * switch(type, monotone = trunc(n1^(1/5)) + 6, convex = trunc(n1^(1/7)) + 6)
      nkts = mult * switch(type, monotone = trunc(n1^(1/7)) + 6, convex = trunc(n1^(1/9)) + 6)
      if(space == "Q"){
        kts = quantile(xu, probs = seq(0, 1, length = nkts), names = FALSE)
      }
      
      if(space == "E"){
        #this will avoid the inflated test size problem when using s.decr.conv
        kts = 0:(nkts - 1) / (nkts - 1) * (max(x) - min(x)) + min(x)
        #kts = seq.int(min(x), max(x), length.out = nkts) # use this for numerical precision; otherwise may get weird H0 fit; see 1996 in CiscoTL.R
      }
      
      #new: if there is any empty interval; do it again
      #no need to worry about linear or quadratic, which use the clumsy way
      #use clumsy way for cubic as well? faster
      # if(attr(x, "parametric") == "cubic"){
      #   nkts0 = nkts
      #   repeat {
      #     #kts <- seq.int(min(x), max(x), length.out = nkts)
      #     kts = 0:(nkts - 1) / (nkts - 1) * (max(x) - min(x)) + min(x)
      #     bins <- cut(x, breaks = kts, include.lowest = TRUE, right = FALSE)
      #     if (all(table(bins) > 0)) break
      #     nkts <- nkts - 1
      #     if (nkts < nkts0 * 2/3) message("Data is too sparse! Too few knots to offer flexible fit!")
      #   }
      # }
      
      kts_lst[[icomp]] = kts
      #new: make delta here; combine it into xm
      spl_degree = switch(type, monotone = 2L, convex = 3L)
      spl_ord = spl_degree + 1
      #dd_ans = bqspl(x, m=NULL, knots=kts)
      #dd = dd_ans$bmat
      dd = bSpline(x, knots = kts[-c(1, nkts)], degree = spl_degree, intercept = TRUE)
      #dd = bs(x, knots = kts[-c(1, nkts)], degree = spl_degree, intercept = TRUE)
      #xm = cbind(xm, dd)
      #test!
      if(icomp > 1){
        dd = dd[, -1]
      }
      dd_lst[[icomp]] = dd
      var_track = c(var_track, rep(icomp, NCOL(dd)))
      #check more
      #if(is.null(edfu)){
      edfu = switch(type, monotone = (nkts/mult) , convex = (nkts/mult + 1))
      edfu_vec[icomp] = edfu #+ nz
      #}
      shp = attr(x, "shape")
      #parametric = attr(x, "parametric")
      #new: make amat here; add it to amat_lst
      amat = makeamat_1D(shp=shp, spl_ord=spl_ord, kts=kts, x=x, x1=min(x), xn=max(x))
      #check more
      #if(iadd > 1) {
      if(icomp > 1){  
        amat = amat[, -1]
      }
      amat_lst[[icomp]] = amat
      #new: make awmat in case we have plane comps.
      #new: make dmat here; add it to dmat_lst
      nc = NCOL(dd)
      dmat = makedmat_1D(nc, q = spl_ord)
      dmat_lst[[icomp]] = dmat
      #dv_lst[[icomp]] = crossprod(dmat)
      #make xm0 here 
      #xm0 = switch(attr(x, "parametric"), linear = cbind(1, x), quadratic = cbind(1, x, x^2))
      xm0_icomp = switch(parametric, linear = cbind(1, x), quadratic = cbind(1, x, x^2), cubic = cbind(1, x, x^2, x^3))
      #check more
      #if(iadd > 1 |ipr > 1) {
      if(icomp > 1){  
        xm0_icomp = xm0_icomp[, -1]
      }
      xm0 = cbind(xm0, xm0_icomp) 
      var_track_param = c(var_track_param, rep(icomp, NCOL(xm0_icomp)))
      #ddwm_lst[[icomp]] = xm0_icomp 
      #print (shp)
      awmat = makeamat_param(shp = shp, parametric = parametric, x1=min(x), xn=max(x))
      #if(ipr > 1 | iadd > 1){
      if(icomp > 1){  
        awmat = awmat[, -1]
      }
      awmat_lst[[icomp]] = awmat
    }
    
    if(NCOL(x) == 2){
      use_constreg = TRUE
      shp = shp_pr = attr(x, "shape")
      #print (shp_pr)
      parametric = attr(x, "parametric")
      ipr = ipr + 1
      x1 = x[, 1]
      x2 = x[, 2]
      xu1 = unique(x1)
      xu2 = unique(x2)
      m1 = round(5*n^(1/6))
      m2 = round(5*n^(1/6))
      x1sc = (x1 - min(x1)) / (max(x1) - min(x1))
      x2sc = (x2 - min(x2)) / (max(x2) - min(x2))
      if(parametric == "linear"){
        if(space == "Q"){
          k1 = quantile(x1, probs = seq(0, 1, length = m1), names = FALSE)
          k2 = quantile(x2, probs = seq(0, 1, length = m2), names = FALSE)
        }
        if(space == "E"){
          k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) + min(x1)
          k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2) - min(x2)) + min(x2)
        }
      } 
      if(parametric == "warped_plane"){
        if(space == "Q"){
          k1 = quantile(x1sc, probs = seq(0, 1, length = m1), names = FALSE)
          k2 = quantile(x2sc, probs = seq(0, 1, length = m2), names = FALSE)
        }
        if(space == "E"){
          k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1sc) - min(x1sc)) + min(x1sc)
          k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2sc) - min(x2sc)) + min(x2sc)
        }
      }
      #if(is.null(edfu)){
      edfu = 3*(n^(1/3))
      edfu_vec[icomp] = edfu
      #}
      #new: make delta here; combine it into xm
      if(parametric == "linear"){
        sp = space
        dd_ans = makedelta_wps(x1, x2, space = c(sp, sp), k1 = k1, k2 = k2, decreasing = c(FALSE, FALSE))
      }
      if(parametric == "warped_plane"){
        sp = space
        dd_ans = makedelta_wps(x1sc, x2sc, space = c(sp, sp), k1 = k1, k2 = k2, decreasing = c(FALSE, FALSE))
      }
      dd = dd_ans$delta
      if(icomp > 1){
        dd = dd[, -1]
      }
      dd_lst[[icomp]] = dd
      var_track = c(var_track, rep(icomp, NCOL(dd)))
      #new: need to re-define k1 and k2; some cells might be empty
      k1 = dd_ans$k1
      k2 = dd_ans$k2
      kts = list(k1 = k1, k2 = k2)
      kts_lst[[icomp]] = kts 
      m1 = length(k1)
      m2 = length(k2)
      use_constreg = TRUE
      #new: make amat here; add it to amat_lst
      amat = makeamat_nonadd(k1, k2, shp_pr)
      #check more
      if(icomp > 1){
        amat = amat[, -1]
      }
      amat_lst[[icomp]] = amat
      #new: make dmat here; add it to dmat_lst
      dmat = makedmat_2D(k1, k2)
      #check more
      if(icomp > 1){ 
        dmat = dmat[, -1]
      }
      dmat_lst[[icomp]] = dmat
      #dv_lst[[icomp]] = crossprod(dmat)
      #make xm0 here
      #conv.conv or conc.conc not implemented now
      xm0_icomp = switch(parametric, linear = cbind(1, x1, x2), warped_plane = cbind(1, x1sc, x2sc, x1sc*x2sc))
      #check more
      if(icomp > 1){
        xm0_icomp = xm0_icomp[, -1]
      }
      xm0 = cbind(xm0, xm0_icomp)
      var_track_param = c(var_track_param, rep(icomp, NCOL(xm0_icomp)))
      #print (shp)
      #need a new name for makeamat_1D_param; it's not just for 1D
      awmat = makeamat_param(shp, parametric = parametric)
      if(icomp > 1){ 
        awmat = awmat[, -1]
      }
      awmat_lst[[icomp]] = awmat
    } 
    #make xm0 here 
    #xm0 = switch(attr(x, "parametric"), linear = cbind(1, x), quadratic = cbind(1, x, x^2), linear_plane = cbind(1, x1, x2), warped_plane = cbind(1, x1sc, x2sc, x1sc*x2sc))
  }
  amat = Matrix::bdiag(amat_lst) |> as.matrix()
  dmat = Matrix::bdiag(dmat_lst) |> as.matrix()
  #dv = Matrix::bdiag(dv_lst) |> as.matrix()
  
  dd = do.call(cbind, dd_lst)
  #dd = Matrix(dd, sparse = TRUE)
  
  nr_am = NROW(amat)
  #----------------------
  #create bvec for qprog
  #----------------------
  bvec = rep(0, nr_am) #will be different when using constreg
  #----------------------------------------------------------------------------
  #2.get H1 fit: cone \ L satisfying some shape constraint
  #write a new function: fit.alt
  #----------------------------------------------------------------------------
  nc = nc_noz = NCOL(dd)
  if(nz > 0){
    dd = cbind(zmat, dd)
    dmat_zero = matrix(0, nrow = nrow(dmat), ncol = nz)
    dmat = cbind(dmat_zero, dmat)
    amat_zero = matrix(0, nrow = nrow(amat), ncol = nz)
    amat = cbind(amat_zero, amat)
  }
  dv = crossprod(dmat)
  #new:
  qv0 = crossprod(dd)

  use_parallel = isTRUE(getOption("cgam.parallel", FALSE))
  #decide how many cores to use
  cores = 1
  if(use_parallel){
    cores = min(4, getOption("cgam.parallel", parallel::detectCores(logical = FALSE) - 1))
  }
  #detect OS platform
  is_windows = .Platform$OS.type == "windows"
  
  nc_am = NCOL(amat)
  imat = diag(nc_am)
  if(GCV){
    #temp:
    edfu = sum(edfu_vec)
    edfu_grid = c(edfu-1, edfu, edfu + (1:2) * 3)
    if(is.null(lams)){
      lams = sapply(edfu_grid, function(e){uniroot(f = .search_ps, qv0 = qv0, dv = dv, dd = dd, edfu = e, interval = c(1e-10, 2e+2), tol=1e-6)$root})
      lams = c(1e-3, rev(lams))
      if(arp){
        #test!
        lams = 2^(2:8)
        #lams = 2^(1:7)
        lams = lams/2^7/n^(2*spl_ord/(2*spl_ord+1))
        # lams = 2^(1:8)
        # lams = lams/2^8/n^(3/4)
      }
      if(is_windows){
        # Windows: use parLapply with a PSOCK cluster
        cl = parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl))  # ensure cleanup
        #mixture_rslt = parallel::parLapply(cl, 1:n.mix, compute_mixture, ysims=ysims, amat=amat, w=w)
        gcvs_rslt = parallel::parLapply(cl, lams, function(e) fit.hypo(p=p, dd=dd, y=y, amat=amat, bvec=bvec,
                                                                  dv=dv, family=family, arp=arp, ps=e))
      } else {
        # Unix/macOS: use mclapply
        gcvs_rslt = parallel::mclapply(lams, function(e) fit.hypo(p=p, dd=dd, y=y, amat=amat, bvec=bvec,
                                                                  dv=dv, family=family, arp=arp, ps=e), mc.cores = (cores))
      }
      #gcvs_rslt = parallel::mclapply(lams, function(e) fit.hypo(p=p, dd=dd, y=y, amat=amat, bvec=bvec,
      #                                                          dv=dv, family=family, arp=arp, ps=e), mc.cores = (8L))
      #gcvs_rslt = parallel::mclapply(lams, function(e) fit.hypo(p=p, dd=dd, y=y, amat=amat, bvec=bvec,
      #                                                          dv=dv, family=family, arp=arp, ps=e), mc.cores = (8L))
      gcvs = sapply(gcvs_rslt, function(rslti) rslti$gcv, simplify = TRUE)
      edfs = sapply(gcvs_rslt, function(rslti) rslti$edf, simplify = TRUE)
      choice_id = which(gcvs == min(gcvs))[1]
      ps = lams[choice_id] 
      p_optim = p
      c(sse1, etahat, ahatc, face, qv, edf, edfs, gcv, gcvs, sighat, edfu_use, covmat, covmat_inv, phi1, mat1) %<-% gcvs_rslt[[choice_id]][1:15]
    }
    #print (head(etahat))
    etahat = as.matrix(etahat)
    muhat = family$linkinv(etahat)
    #only for H1 fit
  } else {
    #lams = ps
    if(is.null(ps)){
      edfu = edfu_vec
      #if(capl == 1){
        #ps = uniroot(f = .search_ps, qv0 = qv0, dv = dv, dd = dd, edfu = edfu, interval = c(1e-10, 2e+2), tol = .Machine$double.eps^0.32)$root
      #  ps = uniroot(f = .search_ps, qv0 = qv0, dv = dv, dd = dd, edfu = edfu, interval = c(1e-10, 2e+2), tol=1e-6)$root
        #ps = uniroot(f = .search_ps, qv0 = qv0, dv = dv, dd = dd, edfu = edfu, interval = c(1e-10, 2e+2))$root
      #} else if (capl > 1) {
      if (capl >= 1) {
        ps = NULL
        for(ic in 1:capl){
          dd_ic = dd_lst[[ic]]
          qv0_ic = crossprod(dd_ic)
          dm_ic = dmat_lst[[ic]]
          dv_ic = crossprod(dm_ic)
          #psi = uniroot(f = .search_ps, qv0 = qv0_ic, dv = dv_ic, dd = dd_ic, edfu = edfu[ic], interval = c(1e-10, 2e+2), tol = .Machine$double.eps^0.32)$root
          psi = uniroot(f = .search_ps, qv0 = qv0_ic, dv = dv_ic, dd = dd_ic, edfu = edfu[ic], interval = c(1e-10, 2e+2), tol=1e-6)$root
          ps = c(ps, psi)
        }
        #ps = sapply(1:capl, function(ic)uniroot(f = .search_ps, qv0 = crossprod(dd_lst[[ic]]), dv = crossprod(dmat_lst[[ic]]), dd = dd_lst[[ic]], edfu = edfu[ic], interval = c(1e-10, 2e+2), tol = .Machine$double.eps^0.32)$root)
      #}
      }
    } 
    lams = ps
    covmat = covmat_inv = phi = NULL
    #new:
    #dv_lst = vector("list", length = (length(dmat_lst) + sum(nz > 0)))
    if(length(ps) >= 1){
      for(ic in 1:capl){
        #dmat_ic = dmat_lst[[ic]]
        #dv_lst[[ic]] = ps[ic] * crossprod(dmat_ic)
        dmat_lst[[ic]] = sqrt(ps[ic]) * dmat_lst[[ic]]
      }
      #dv = Matrix::bdiag(dv_lst) |> as.matrix()
      #already add zmat in the 1st element of dmat_lst
      dmat = Matrix::bdiag(dmat_lst) |> as.matrix()
      #recreate dmat
      if(nz > 0){
        dmat_zero = matrix(0, nrow = nrow(dmat), ncol = nz)
        dmat = cbind(dmat_zero, dmat)
      }
      dv = crossprod(dmat)
    }
    if(length(p) == 1){ 
      #arp with user-defined order, or !arp and p=0
      #cat('call arp with user-defined order, ps=', ps, '\n')
      ansc = fit.hypo(p=p, dd=dd, y=y, amat=amat, bvec=bvec, dv=dv, wmat=NULL, family=family, ps=ps, arp=arp)
      c(sse1, etahat, ahatc, face, qv, edf, edfs, gcv, gcvs, sighat, edfu_use, covmat, covmat_inv, phi1, mat1) %<-% ansc[1:15]
      p_optim = p
    }else{
      #p = 0:4
      #new: test p=0 first
      ansc = fit.hypo(p=0, dd=dd, y=y, amat=amat, bvec=bvec, dv=dv, family=family, ps=ps, arp=FALSE)
      if(ansc$pval.ts > 0.05){
        c(sse1, etahat, ahatc, face, qv, edf, edfs, gcv, gcvs, sighat, edfu_use, covmat, covmat_inv, phi1, mat1) %<-% ansc[1:15]
        p_optim = 0
      }else{
        if(is_windows){
          # Windows: use parLapply with a PSOCK cluster
          cl = parallel::makeCluster(cores)
          on.exit(parallel::stopCluster(cl))  # ensure cleanup
          aic_rslt = parallel::parLapply(cl, p[-1], fit.hypo, dd=dd, y=y, amat=amat, bvec=bvec, 
                                        dv=dv, family=family, ps=ps, arp=arp)
        } else {
          # Unix/macOS: use mclapply
          aic_rslt = parallel::mclapply(p[-1], fit.hypo, dd=dd, y=y, amat=amat, bvec=bvec, 
                                        dv=dv, family=family, ps=ps, arp=arp, mc.cores=(cores))
        }
        #aic_rslt = parallel::mclapply(p[-1], fit.hypo, dd=dd, y=y, amat=amat, bvec=bvec, 
        #                              dv=dv, family=family, ps=ps, arp=arp, mc.cores=(8L))
        aics = sapply(aic_rslt, function(rslti) rslti$aic, simplify = TRUE)
        #print (aics)
        choice_id = which(aics == min(aics))
        p_optim = (p[-1])[choice_id] 
        c(sse1, etahat, ahatc, face, qv, edf, edfs, gcv, gcvs, sighat, edfu_use, covmat, covmat_inv, phi1, mat1) %<-% (aic_rslt[[choice_id]])[1:15]
      }
    }
    muhat = family$linkinv(etahat)
  }
  #----------------------------------------------------------------------------
  #1.get H0 fit: linear space satisfying some shape constraint
  #use sparse matrix later
  #write a new function: fit.null
  #----------------------------------------------------------------------------
  wmat = x1m0 = x1m = NULL
  #test!
  #use_constreg = FALSE
  used_lm = FALSE
  if(!use_constreg){
    pm0 = xm0 %*% solve(crossprod(xm0), t(xm0))
    #modify dd: keep the columns that give a full column rank; which can pass solve()
    #test!
    qr_obj = qr(dd)
    rk1 = qr_obj$rank
    rk2 = sum(svd(dd)$d > 1e-7)
    rk = rk1 # min(rk2, rk1)
    kept_id = qr_obj$pivot[1:rk]
    #kept_id = 1:rk
    #x1m0 = dd[, kept_id] - pm0 %*% dd[, kept_id]
    x1m0 = dd - pm0 %*% dd
    qr_x1m = qr(x1m0)
    x1m = qr.Q(qr_x1m)[, 1:(qr_x1m$rank), drop = FALSE]
    #test!
    #assume that xm0 is full column rank
    #x1m = qr.Q(qr_x1m, complete = TRUE)[, 1:(rk - ncol(xm0))]
    #xm1tb = crossprod(x1m, dd[, kept_id])
    xm1tb = crossprod(x1m, dd)
    qr_xm1tb = qr(t(xm1tb))
    #new:
    #wmat = qr.Q(qr_xm1tb, complete = TRUE)[, -c(1:rk), drop = FALSE]
    wmat = qr.Q(qr_xm1tb, complete = TRUE)[, -c(1:qr_xm1tb$rank), drop = FALSE]
    ddwm = dd %*% wmat #should match the num of columns in xm0?
    awmat = amat %*% wmat
    
    #new: remove all zero rows
    awmat = round(awmat, 10)
    awmat = awmat[!apply(awmat, 1, function(row) all(row == 0)), , drop = FALSE]
    
    #new: remove redundant rows
    awmat = unique(awmat)
  } else {
    ddwm = xm0
    #avoid solve error for H0
    #ddwm <- scale(xm0, center = FALSE)
    awmat = Matrix::bdiag(awmat_lst) |> as.matrix()
  }
  if(nz > 0){
    #add z before splines
    ddwm = cbind(zmat, ddwm)
    awmat_zero = matrix(0, nrow = nrow(awmat), ncol = nz)
    awmat = cbind(awmat_zero, awmat)
  }
  #ansl = fit.hypo(p=p_optim, dd=dd, y=y, amat=awmat, bvec=rep(0, nrow(awmat)), dv=dv, wmat=wmat, family=family, ps=0, arp=arp)
  #if(!used_lm){
  ansl = fit.hypo(p=p_optim, dd=ddwm, y=y, amat=awmat, bvec=rep(0, nrow(awmat)), 
                  dv=NULL, family=family, ps=0, arp=arp)
  sse0 = ansl$dev
  etahat0 = ansl$etahat
  muhat0 = family$linkinv(etahat0)
  ahatl = ansl$ahat
  qvl = ansl$qv
  #}
  #cat(arp, '\n')
  #----------------------------------------------------------------------------
  #test: H0 vs H1
  #----------------------------------------------------------------------------
  if(wt.iter){
    bval = (sse0 - sse1) / n #sse0 is llh0; sse1 is llh1
  } else {
    bval = (sse0 - sse1) / sse0
  }
  #temp
  #test!
  bval = round(bval, 7)
  sm=1e-7
  #sm = 1e-10
  if (bval > sm & nsim > 0) {
    if (multicore) {
      #tot_cores = parallel::detectCores()
      if(is_windows){
        # Windows: use parLapply with a PSOCK cluster
        cl = parallel::makeCluster(cores)
        on.exit(parallel::stopCluster(cl))  # ensure cleanup
        #mixture_rslt = parallel::parLapply(cl, 1:n.mix, compute_mixture, ysims=ysims, amat=amat, w=w)
        bdist = parallel::parLapply(cl, 1:nsim, .compute_bstat, etahat0=etahat0, n=n, sighat=sighat,
                                   dd=dd, ddwm=ddwm, qv=NULL, qvl=NULL, amat=amat, awmat=awmat,
                                   bvec=bvec, imat=imat, dv=dv, lams=ps, w=NULL, arp=arp, p=p_optim, phi=phi1, 
                                   family=family)
      } else {
        # Unix/macOS: use mclapply
        bdist = parallel::mclapply(1:nsim, .compute_bstat, etahat0=etahat0, n=n, sighat=sighat,
                                   dd=dd, ddwm=ddwm, qv=NULL, qvl=NULL, amat=amat, awmat=awmat,
                                   bvec=bvec, imat=imat, dv=dv, lams=ps, w=NULL, arp=arp, p=p_optim, phi=phi1, 
                                   family=family, mc.cores=(cores))
      }
      #bdist = parallel::mclapply(1:nsim, .compute_bstat, etahat0=etahat0, n=n, sighat=sighat,
      #                           dd=dd, ddwm=ddwm, qv=NULL, qvl=NULL, amat=amat, awmat=awmat,
      #                           bvec=bvec, imat=imat, dv=dv, lams=ps, w=NULL, arp=arp, p=p_optim, phi=phi1, 
      #                           family=family, mc.cores=(4L))
    } else {
      bdist = lapply(1:nsim, .compute_bstat, etahat0=etahat0, n=n, sighat=sighat, dd=dd, ddwm=ddwm, qv=NULL, qvl=NULL,
                     amat=amat, awmat=awmat, bvec=bvec, imat=imat, dv=dv, lams=ps, w=NULL, arp=arp, p=p_optim, phi=phi1, family=family)
    }
    bdist = simplify2array(bdist)
    pval = sum(bdist > bval) / nsim
  } else {
    pval = 1
  }
  #----------------------
  #for visualization
  #----------------------
  etacomps = etacomps0 = vector("list", length = capl)
  etahat0_surf = etahat_surf = NULL
  muhat0_surf = muhat_surf = NULL
  ahatc_noz = round(ahatc, 6) |> as.matrix()
  ahatl_noz = round(ahatl, 6) |> as.matrix()
  #print (ahatc_noz)
  #print (nz)
  if(nz > 0){
    #print (dim(ahatc_noz))
    ahatc_noz = ahatc_noz[-c(1:nz), ,drop=F] #|> as.vector()
    ahatl_noz = ahatl_noz[-c(1:nz), ,drop=F] #|> as.vector()
  }
  for(icomp in 1:capl){
    x = X[[icomp]]
    if(NCOL(x) == 1){
      #print (dim(dd_lst[[icomp]]))
      #print (dim(ahatc_noz[var_track == icomp]))
      etahat_icomp = dd_lst[[icomp]] %*% ahatc_noz[var_track == icomp,,drop=F]
      #print (dim(etahat_icomp))
      etacomps[[icomp]] = etahat_icomp
    }
    if(NCOL(x) == 2){
      kts = kts_lst[[icomp]]
      k1 = kts$k1
      k2 = kts$k2
      newd = expand.grid(k1, k2)
      newd = as.matrix(newd)
      
      #H0 surface
      x1p = newd[,1]
      x2p = newd[,2]
      if(attr(x, "parametric") == "warped_plane"){
        xm0p = cbind(1, x1p, x2p, x1p*x2p)
        if(icomp > 1){
          xm0p = xm0p[, -1]
        }
      }
      if(attr(x, "parametric") == "linear"){
        xm0p = cbind(1, x1p, x2p)
        if(icomp > 1){
          xm0p = xm0p[, -1]
        }
      }
      psurf0 = xm0p %*% ahatl_noz[var_track_param == icomp,,drop=F]
      etahat0_surf = matrix(psurf0, m1, m2)
      muhat0_surf = family$linkinv(etahat0_surf)
      
      #H1 surface
      pans = makedelta_wps(newd[,1], newd[,2], k1 = k1, k2 = k2, decreasing = c(FALSE, FALSE))
      ddp = pans$delta
      if(icomp > 1){
        ddp = ddp[, -1]
      }
      psurf = ddp %*% ahatc_noz[var_track == icomp,,drop=F]
      etahat_surf = matrix(psurf, m1, m2)
      muhat_surf = family$linkinv(etahat_surf)
      etacomps[[icomp]] = etahat_surf
      etacomps0[[icomp]] = etahat0_surf
    }
  }
  rslt = list(pval = pval, bval = bval, knots = kts, k1=k1, k2=k2, bmat = dd, wmat = wmat, dmat = dmat, ps = ps,
              xm0 = xm0, x1m = x1m, amat = amat, awmat = awmat, face = face, etahat = etahat, etahat0 = etahat0,
              etahat0_surf = etahat0_surf, etahat_surf = etahat_surf, muhat0_surf = muhat0_surf, muhat_surf = muhat_surf,
              sighat = sighat, edf = edf, edfu = edfu_use, lams = lams, gcvs = gcvs, edfs = edfs, ahatl = ahatl, ahatc = ahatc,
              phi = phi1, covmat = covmat, covmat_inv = covmat_inv, etacomps = etacomps, p_optim = p_optim, kts_lst = kts_lst, 
              var_track = var_track, var_track_param = var_track_param, etacomps0 = etacomps0, use_constreg=use_constreg)
  if(nz >= 1){
    #assume y is iid with common sig2
    if(!wt.iter){
      #should work for ar(p)?
      if(!arp || arp & p_optim == 0){
        covmat0 = sighat^2 * mat1 %*% t(mat1)
      } else {
        #print (sig2hat_z)
        #covmat0 = sig2hat_z * mat1 %*% t(mat1)
        covmat0 = mat1 %*% t(mat1) #seems working with unscaled covariance
        #covmat0 = sighat^2 * mat1 %*% t(mat1)
      }
    } else {
      wt = wt.fun(y, etahat, n, weights = rep(1, n), fml = family$family)
      covmat0 = mat1 %*% diag(wt) %*% t(mat1)
    }
    covmat = covmat0[(1):(nz), (1):(nz), drop=FALSE]
    sez = sqrt(diag(covmat))
    zcoefs = ahatc[(1):(nz)]
    tz = zcoefs / sez
    cpar = 1.2
    #test more!
    if(!wt.iter){
      if ((n - cpar * edf) <= 0) {
        pz = 2*(1 - pt(abs(tz), edf))
      } else {
        #why not pnorm?
        pz = 2*(1 - pt(abs(tz), n - cpar * edf))
      }
    } else {
      #why not pnorm?
      pz = 2*(1 - pnorm(abs(tz)))
    }
    rslt$sez = sez
    rslt$pz = pz
    rslt$tz = tz
    rslt$zcoefs = zcoefs
  } else {rslt$sez = NULL; rslt$pz = NULL; rslt$tz = NULL; rslt$zcoefs = NULL}
  return(rslt)
}

summary.testpar <- function(object,...) {
  if (!inherits(object, "testpar")) {
    stop("object must be of class 'testpar'")
  }
  # Significance stars
  pval <- object$pval
  if (is.na(pval)) {
    stars <- ""
  } else if (pval <= 0.001) {
    stars <- "***"
  } else if (pval <= 0.01) {
    stars <- "**"
  } else if (pval <= 0.05) {
    stars <- "*"
  } else if (pval <= 0.1) {
    stars <- "."
  } else {
    stars <- ""
  }
  
  cat("--------------------------------------------------------------\n")
  cat("|   Parametric Fit vs Smooth Constrained Fit Test             |\n")
  cat("--------------------------------------------------------------\n\n")
  
  # Composite H0 vs H1 with proper formatting
  if (!is.null(object$X) && is.list(object$X)) {
    h0_terms <- h1_terms <- character(length(object$X))
    
    for (i in seq_along(object$X)) {
      xi <- object$X[[i]]
      param <- attr(xi, "parametric")
      type  <- attr(xi, "type")
      varnames <- attr(xi, "nm")
      
      # Fallback name
      if (is.null(varnames)) {
        varnames <- paste0("x", i)
      }
      
      # Format variable name as "x1 * x2" if it's an interaction
      if (param == "warped_plane" && length(varnames) == 2) {
      #if (param == "warped_plane") {
        var <- paste(varnames, collapse = " * ")
      } else {
        var <- varnames[1]
      }
      
      h0_terms[i] <- paste0(var, ": ", param)
      h1_terms[i] <- paste0(var, ": ", type)
    }
    
    # Combine terms across X with "+"
    h0_str <- paste(h0_terms, collapse = " + ")
    h1_str <- paste(h1_terms, collapse = " + ")
    
    cat("  Null hypothesis (H0):", h0_str, "\n")
    cat("  Alternative (H1):    ", h1_str, "\n\n")
  }
  
  # Test statistic
  if (!is.null(object$bval)) {
    cat(sprintf("  %-32s %.4f\n", "Mixture-of-beta test statistic:", round(object$bval, 4)))
  }
  
  # Degrees of freedom
  if (!is.null(object$edf)) {
    cat(sprintf("  %-32s %s\n", "Effective degrees of freedom:", round(object$edf, 4)))
  }
  
  # p-value
  cat(sprintf("  %-32s %.4g %s\n", "p-value:", round(pval, 4), stars))
  cat("---------------------------------------------------------------\n")
  cat("Signif. codes:  0 '***'  0.001 '**'  0.01 '*'  0.05 '.'  0.1 ' ' 1\n")
  
  invisible(object)
}


#-----------------------------------------------------
#similar to glm.fit
#X is a list with attributes: shp, parametric, type
#-----------------------------------------------------
# testpar.fit = function(X, y, zmat = NULL, family = gaussian(link="identity"), 
#                        ps = NULL, edfu = NULL, nsim = 200, multicore = TRUE, 
#                        GCV = FALSE, lams = NULL, 
#                        arp = FALSE, p = NULL, space = "Q",...)
# {
#   #print (head(y))
#   #print (multicore)
#   wt.iter = ifelse(family$family == "gaussian", FALSE, TRUE)
#   extras = list(...)
#   #new:
#   #capl = NCOL(x)
#   capl = length(X)
#   n = length(y)
#   k1 = k2 = NULL
#   p_optim = 0
#   if(arp){
#     #if the user didn't provide p, then choose p from 0:3
#     if(is.null(p)){
#       p = 0:2 
#     } else if (p < 0 || p > 2) {
#       warning("The order in AR(p) must be an integer >= 0 and <= 2!")
#       p = 1
#     }
#   } else {
#     p = 0
#   }
#   #test more
#   # if(!wt.iter){
#   #   mu = mean(y)
#   #   y = y - mu
#   # }
#   #-------------------
#   #get xm, amat, dmat
#   #-------------------
#   mult = 2
#   nz = 0
#   if(!is.null(zmat)){
#     nz = NCOL(zmat)
#     #test more: center zmat to avoid solve() error
#     #zmat = scale(zmat, center = TRUE, scale = FALSE)
#   }
#   xm0 = NULL #combine dd columnwise
#   amat_lst = vector("list", length = capl)
#   awmat_lst = vector("list", length = capl)
#   dmat_lst = vector("list", length = capl)
#   dd_lst = kts_lst = vector("list", length = capl)
#   edfu_vec = numeric(capl)
#   var_track = var_track_param = NULL
#   iadd = ipr = 0
#   for(icomp in 1:capl){
#     x = X[[icomp]] #|> as.matrix()
#     if(NCOL(x) == 1){
#       iadd = iadd + 1
#       xu = unique(x)
#       n1 = length(xu)
#       type = attr(x, "type")
#       nkts = mult * switch(type, monotone = trunc(n1^(1/5)) + 6, convex = trunc(n1^(1/7)) + 6)
#       if(space == "Q"){
#         kts = quantile(xu, probs = seq(0, 1, length = nkts), names = FALSE)
#       }
#       if(space == "E"){
#         kts = 0:(nkts - 1) / (nkts - 1) * (max(x) - min(x)) + min(x)
#         #kts = seq.int(min(x), max(x), length = nkts)
#       }
#       kts_lst[[icomp]] = kts
#       #new: make delta here; combine it into xm
#       spl_degree = switch(type, monotone = 2L, convex = 3L)
#       spl_ord = spl_degree + 1
#       #dd_ans = bqspl(x, m=NULL, knots=kts)
#       #dd = dd_ans$bmat
#       dd = bSpline(x, knots = kts[-c(1, nkts)], degree = spl_degree, intercept = TRUE)
#       #dd = bs(x, knots = kts[-c(1, nkts)], degree = spl_degree, intercept = TRUE)
#       #xm = cbind(xm, dd)
#       #test!
#       if(icomp > 1){
#         dd = dd[, -1]
#       }
#       dd_lst[[icomp]] = dd
#       var_track = c(var_track, rep(icomp, NCOL(dd)))
#       #check more
#       #if(is.null(edfu)){
#       edfu = switch(type, monotone = (nkts/mult) , convex = (nkts/mult + 1))
#       edfu_vec[icomp] = edfu #+ nz
#       #}
#       shp = attr(x, "shape")
#       parametric = attr(x, "parametric")
#       #new: make amat here; add it to amat_lst
#       #print (kts)
#       # print (spl_ord)
#       # print (min(x))
#       # print (max(x))
#       amat = makeamat_1D(shp=shp, spl_ord=spl_ord, kts=kts, x=x, x1=min(x), xn=max(x))
#       #check more
#       #if(iadd > 1) {
#       if(icomp > 1){  
#         amat = amat[, -1]
#       }
#       amat_lst[[icomp]] = amat
#       #new: make awmat in case we have plane comps.
#       #new: make dmat here; add it to dmat_lst
#       nc = NCOL(dd)
#       dmat = makedmat_1D(nc, q = spl_ord)
#       dmat_lst[[icomp]] = dmat
#       #dv_lst[[icomp]] = crossprod(dmat)
#       #make xm0 here 
#       #xm0 = switch(attr(x, "parametric"), linear = cbind(1, x), quadratic = cbind(1, x, x^2))
#       xm0_icomp = switch(parametric, linear = cbind(1, x), quadratic = cbind(1, x, x^2))
#       #check more
#       #if(iadd > 1 |ipr > 1) {
#       if(icomp > 1){  
#         xm0_icomp = xm0_icomp[, -1]
#       }
#       xm0 = cbind(xm0, xm0_icomp) 
#       var_track_param = c(var_track_param, rep(icomp, NCOL(xm0_icomp)))
#       #ddwm_lst[[icomp]] = xm0_icomp 
#       #print (shp)
#       awmat = makeamat_1D_param(shp = shp, parametric = parametric, x1=min(x), xn=max(x))
#       #if(ipr > 1 | iadd > 1){
#       if(icomp > 1){  
#         awmat = awmat[, -1]
#       }
#       awmat_lst[[icomp]] = awmat
#     }
#     
#     use_constreg = FALSE
#     if(NCOL(x) == 2){
#       use_constreg = TRUE
#       shp = shp_pr = attr(x, "shape")
#       #print (shp_pr)
#       parametric = attr(x, "parametric")
#       ipr = ipr + 1
#       x1 = x[, 1]
#       x2 = x[, 2]
#       xu1 = unique(x1)
#       xu2 = unique(x2)
#       m1 = round(5*n^(1/6))
#       m2 = round(5*n^(1/6))
#       x1sc = (x1 - min(x1)) / (max(x1) - min(x1))
#       x2sc = (x2 - min(x2)) / (max(x2) - min(x2))
#       if(parametric == "linear"){
#         if(space == "Q"){
#           k1 = quantile(x1, probs = seq(0, 1, length = m1), names = FALSE)
#           k2 = quantile(x2, probs = seq(0, 1, length = m2), names = FALSE)
#         }
#         if(space == "E"){
#           k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) + min(x1)
#           k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2) - min(x2)) + min(x2)
#         }
#       } 
#       if(parametric == "warped_plane"){
#         if(space == "Q"){
#           k1 = quantile(x1sc, probs = seq(0, 1, length = m1), names = FALSE)
#           k2 = quantile(x2sc, probs = seq(0, 1, length = m2), names = FALSE)
#         }
#         if(space == "E"){
#           k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1sc) - min(x1sc)) + min(x1sc)
#           k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2sc) - min(x2sc)) + min(x2sc)
#         }
#       }
#       #if(is.null(edfu)){
#       edfu = 3*(n^(1/3))
#       edfu_vec[icomp] = edfu
#       #}
#       #new: make delta here; combine it into xm
#       if(parametric == "linear"){
#         sp = space
#         dd_ans = makedelta_wps(x1, x2, space = c(sp, sp), k1 = k1, k2 = k2, decreasing = c(FALSE, FALSE))
#       }
#       if(parametric == "warped_plane"){
#         sp = space
#         dd_ans = makedelta_wps(x1sc, x2sc, space = c(sp, sp), k1 = k1, k2 = k2, decreasing = c(FALSE, FALSE))
#       }
#       dd = dd_ans$delta
#       if(icomp > 1){
#         dd = dd[, -1]
#       }
#       dd_lst[[icomp]] = dd
#       var_track = c(var_track, rep(icomp, NCOL(dd)))
#       #new: need to re-define k1 and k2; some cells might be empty
#       k1 = dd_ans$k1
#       k2 = dd_ans$k2
#       kts = list(k1 = k1, k2 = k2)
#       kts_lst[[icomp]] = kts 
#       m1 = length(k1)
#       m2 = length(k2)
#       use_constreg = TRUE
#       #new: make amat here; add it to amat_lst
#       amat = makeamat_nonadd(k1, k2, shp_pr)
#       #check more
#       if(icomp > 1){
#         amat = amat[, -1]
#       }
#       amat_lst[[icomp]] = amat
#       #new: make dmat here; add it to dmat_lst
#       dmat = makedmat_2D(k1, k2)
#       #check more
#       if(icomp > 1){ 
#         dmat = dmat[, -1]
#       }
#       dmat_lst[[icomp]] = dmat
#       #dv_lst[[icomp]] = crossprod(dmat)
#       #make xm0 here 
#       xm0_icomp = switch(parametric, linear = cbind(1, x1, x2), warped_plane = cbind(1, x1sc, x2sc, x1sc*x2sc))
#       #check more
#       if(icomp > 1){
#         xm0_icomp = xm0_icomp[, -1]
#       }
#       xm0 = cbind(xm0, xm0_icomp)
#       var_track_param = c(var_track_param, rep(icomp, NCOL(xm0_icomp)))
#       #print (shp)
#       awmat = makeamat_1D_param(shp, parametric = parametric)
#       if(icomp > 1){ 
#         awmat = awmat[, -1]
#       }
#       awmat_lst[[icomp]] = awmat
#     } 
#     #make xm0 here 
#     #xm0 = switch(attr(x, "parametric"), linear = cbind(1, x), quadratic = cbind(1, x, x^2), linear_plane = cbind(1, x1, x2), warped_plane = cbind(1, x1sc, x2sc, x1sc*x2sc))
#   }
#   amat = Matrix::bdiag(amat_lst) |> as.matrix()
#   dmat = Matrix::bdiag(dmat_lst) |> as.matrix()
#   #dv = Matrix::bdiag(dv_lst) |> as.matrix()
#   
#   dd = do.call(cbind, dd_lst)
#   #dd = Matrix(dd, sparse = TRUE)
#   
#   nr_am = NROW(amat)
#   #----------------------
#   #create bvec for qprog
#   #----------------------
#   bvec = rep(0, nr_am) #will be different when using constreg
#   #----------------------------------------------------------------------------
#   #2.get H1 fit: cone \ L satisfying some shape constraint
#   #write a new function: fit.alt
#   #----------------------------------------------------------------------------
#   nc = nc_noz = NCOL(dd)
#   if(nz > 0){
#     # dd_1 = dd_lst[[1]]
#     # dd_1 = cbind(zmat, dd_1)
#     # dd_lst[[1]] = dd_1
#     
#     # dm_1 = dmat_lst[[1]]
#     # dm1_zero = matrix(0, nrow = nrow(dm_1), ncol = nz)
#     # dm_1 = cbind(dm1_zero, dm_1)
#     # dmat_lst[[1]] = dm_1
#     
#     dd = cbind(zmat, dd)
#     dmat_zero = matrix(0, nrow = nrow(dmat), ncol = nz)
#     dmat = cbind(dmat_zero, dmat)
#     amat_zero = matrix(0, nrow = nrow(amat), ncol = nz)
#     amat = cbind(amat_zero, amat)
#   }
#   dv = crossprod(dmat)
#   #new:
#   qv0 = crossprod(dd)
#   
#   #dv = Matrix(dv, sparse = TRUE)
#   #qv0 = Matrix(qv0, sparse = TRUE)
#   
#   nc_am = NCOL(amat)
#   imat = diag(nc_am)
#   if(GCV){
#     #temp:
#     edfu = sum(edfu_vec)
#     edfu_grid = c(edfu-1, edfu, edfu + (1:2) * 3)
#     if(is.null(lams)){
#       lams = sapply(edfu_grid, function(e){uniroot(f = .search_ps, qv0 = qv0, dv = dv, dd = dd, edfu = e, interval = c(1e-10, 2e+2), tol=1e-6)$root})
#       lams = c(1e-3, rev(lams))
#       if(arp){
#         #test!
#         lams = 2^(2:8)
#         #lams = 2^(1:7)
#         lams = lams/2^7/n^(2*spl_ord/(2*spl_ord+1))
#         # lams = 2^(1:8)
#         # lams = lams/2^8/n^(3/4)
#       }
#       gcvs_rslt = parallel::mclapply(lams, function(e) fit.hypo(p=p, dd=dd, y=y, amat=amat, bvec=bvec,
#                                                                 dv=dv, family=family, arp=arp, ps=e), mc.cores = (8L))
#       #gcvs_rslt = parallel::mclapply(lams, function(e) fit.hypo(p=p, dd=dd, y=y, amat=amat, bvec=bvec,
#       #                                                          dv=dv, family=family, arp=arp, ps=e), mc.cores = (8L))
#       gcvs = sapply(gcvs_rslt, function(rslti) rslti$gcv, simplify = TRUE)
#       edfs = sapply(gcvs_rslt, function(rslti) rslti$edf, simplify = TRUE)
#       choice_id = which(gcvs == min(gcvs))
#       ps = lams[choice_id] 
#       p_optim = p
#       c(sse1, etahat, ahatc, face, qv, edf, edfs, gcv, gcvs, sighat, edfu_use, covmat, covmat_inv, phi1, mat1) %<-% gcvs_rslt[[choice_id]][1:15]
#     }
#     #print (head(etahat))
#     etahat = as.matrix(etahat)
#     muhat = family$linkinv(etahat)
#     #only for H1 fit
#   } else {
#     #lams = ps
#     if(is.null(ps)){
#       #edfu = sum(edfu_vec) + nz
#       edfu = edfu_vec
#       #not right:...
#       #if(nz > 0){
#       #  edfu[1] = edfu_vec[1] + nz
#       #}
#       #if(!arp){
#       if(capl == 1){
#         #ps = uniroot(f = .search_ps, qv0 = qv0, dv = dv, dd = dd, edfu = edfu, interval = c(1e-10, 2e+2), tol = .Machine$double.eps^0.32)$root
#         ps = uniroot(f = .search_ps, qv0 = qv0, dv = dv, dd = dd, edfu = edfu, interval = c(1e-10, 2e+2), tol=1e-6)$root
#         #test!
#         #if(arp & p >= 1){
#         #print (ps)
#         #ps = ps * 0.6 #1/2 (2/3) will be unbiased for p=2 but inflated for p=1, 2/3 is inflated for p=1
#         #print (ps)
#         #}
#       } else if (capl > 1) {
#         ps = NULL
#         for(ic in 1:capl){
#           dd_ic = dd_lst[[ic]]
#           qv0_ic = crossprod(dd_ic)
#           dm_ic = dmat_lst[[ic]]
#           dv_ic = crossprod(dm_ic)
#           #psi = uniroot(f = .search_ps, qv0 = qv0_ic, dv = dv_ic, dd = dd_ic, edfu = edfu[ic], interval = c(1e-10, 2e+2), tol = .Machine$double.eps^0.32)$root
#           psi = uniroot(f = .search_ps, qv0 = qv0_ic, dv = dv_ic, dd = dd_ic, edfu = edfu[ic], interval = c(1e-10, 2e+2), tol=1e-6)$root
#           ps = c(ps, psi)
#         }
#         #ps = sapply(1:capl, function(ic)uniroot(f = .search_ps, qv0 = crossprod(dd_lst[[ic]]), dv = crossprod(dmat_lst[[ic]]), dd = dd_lst[[ic]], edfu = edfu[ic], interval = c(1e-10, 2e+2), tol = .Machine$double.eps^0.32)$root)
#       }
#       #}else{
#       #ps = 2/n^(2*spl_ord/(2*spl_ord+1)) #seems working
#       #  ps = 1/n^(2*spl_ord/(2*spl_ord+1)) #a little inflated
#       #}
#     } 
#     lams = ps
#     covmat = covmat_inv = phi = NULL
#     #new:
#     if(length(ps) >= 1){
#       for(ic in 1:capl){
#         dmat_lst[[ic]] = sqrt(ps[ic]) * dmat_lst[[ic]]
#       }
#       #already add zmat in the 1st element of dmat_lst
#       dmat = Matrix::bdiag(dmat_lst) |> as.matrix()
#       #recreate dmat
#       if(nz > 0){
#         dmat_zero = matrix(0, nrow = nrow(dmat), ncol = nz)
#         dmat = cbind(dmat_zero, dmat)
#       }
#       dv = crossprod(dmat)
#     }
#     if(length(p) == 1){ 
#       #arp with user-defined order, or !arp and p=0
#       #cat('call arp with user-defined order, ps=', ps, '\n')
#       ansc = fit.hypo(p=p, dd=dd, y=y, amat=amat, bvec=bvec, dv=dv, family=family, ps=ps, arp=arp)
#       c(sse1, etahat, ahatc, face, qv, edf, edfs, gcv, gcvs, sighat, edfu_use, covmat, covmat_inv, phi1, mat1) %<-% ansc[1:15]
#       p_optim = p
#     } else{
#       #p = 0:2
#       #new: test p=0 first
#       ansc = fit.hypo(p=0, dd=dd, y=y, amat=amat, bvec=bvec, dv=dv, family=family, ps=ps, arp=FALSE)
#       if(ansc$pval.ts > 0.05){
#         c(sse1, etahat, ahatc, face, qv, edf, edfs, gcv, gcvs, sighat, edfu_use, covmat, covmat_inv, phi1, mat1) %<-% ansc[1:15]
#         p_optim = 0
#       }else{
#         aic_rslt = parallel::mclapply(p[-1], fit.hypo, dd=dd, y=y, amat=amat, bvec=bvec, dv=dv, family=family, ps=ps, arp=arp, mc.cores=(8L))
#         aics = sapply(aic_rslt, function(rslti) rslti$aic, simplify = TRUE)
#         #print (aics)
#         choice_id = which(aics == min(aics))
#         p_optim = (p[-1])[choice_id] 
#         c(sse1, etahat, ahatc, face, qv, edf, edfs, gcv, gcvs, sighat, edfu_use, covmat, covmat_inv, phi1, mat1) %<-% (aic_rslt[[choice_id]])[1:15]
#       }
#     }
#     muhat = family$linkinv(etahat)
#   }
#   #----------------------------------------------------------------------------
#   #1.get H0 fit: linear space satisfying some shape constraint
#   #use sparse matrix later
#   #write a new function: fit.null
#   #----------------------------------------------------------------------------
#   wmat = x1m0 = x1m = NULL
#   if(!use_constreg){
#     pm0 = xm0 %*% solve(crossprod(xm0), t(xm0))
#     #x1m0 = dd - pm0 %*% dd
#     
#     # nc = ncol(dd)
#     # qr_dd = qr(dd)
#     # dd_r = qr_dd$rank
#     # if(dd_r < nc){
#     #   dd_2 = qr.Q(qr_dd, complete = TRUE)[, 1:dd_r]
#     #   rm_id = qr_dd$pivot[(dd_r+1):nc]
#     #   x1m0 = dd_2 - pm0 %*% dd_2
#     # }else{
#     x1m0 = dd - pm0 %*% dd
#     #}
#     
#     qr_x1m = qr(x1m0)
#     x1m = qr.Q(qr_x1m, complete = TRUE)[, 1:qr_x1m$rank]
#     # if(dd_r < nc){
#     #   xm1tb = crossprod(x1m, dd_2)
#     # } else {
#     xm1tb = crossprod(x1m, dd)
#     #}
#     qr_xm1tb = qr(t(xm1tb))
#     wmat = qr.Q(qr_xm1tb, complete = TRUE)[, -c(1:qr_xm1tb$rank), drop = FALSE]
#     # if(dd_r < nc){
#     #   ddwm = dd_2 %*% wmat
#     #   awmat = amat[,-rm_id] %*% wmat
#     # }else {
#     ddwm = dd %*% wmat
#     awmat = amat %*% wmat
#     #}
#   } else {
#     ddwm = xm0
#     awmat = Matrix::bdiag(awmat_lst) |> as.matrix()
#   }
#   if(nz > 0){
#     #add z before splines
#     ddwm = cbind(zmat, ddwm)
#     awmat_zero = matrix(0, nrow = nrow(awmat), ncol = nz)
#     awmat = cbind(awmat_zero, awmat)
#   }
#   ansl = fit.hypo(p=p_optim, dd=ddwm, y=y, amat=awmat, bvec=rep(0, nrow(awmat)), dv=NULL, family=family, ps=0, arp=arp)
#   sse0 = ansl$dev
#   etahat0 = ansl$etahat
#   muhat0 = family$linkinv(etahat0)
#   ahatl = ansl$ahat
#   qvl = ansl$qv
#   #cat(arp, '\n')
#   #----------------------------------------------------------------------------
#   #test: H0 vs H1
#   #----------------------------------------------------------------------------
#   if(wt.iter){
#     bval = (sse0 - sse1) / n #sse0 is llh0; sse1 is llh1
#   } else {
#     bval = (sse0 - sse1) / sse0
#   }
#   #temp
#   sm=1e-7
#   #sm = 1e-9
#   
#   #new
#   #Check if the user wants parallel computing
#   #assume user-default is FALSE; user needs to call option() to change this
#   use_parallel = isTRUE(getOption("cgam.parallel", FALSE))
#   #Decide how many cores to use
#   cores = 1
#   if(use_parallel){
#     #cores = getOption("cgam.cores", parallel::detectCores(logical = FALSE) - 1)
#     cores = min(8, getOption("cgam.cores", parallel::detectCores(logical = FALSE) - 1))
#     cores = max(1, cores)
#   }
#   
#   #Detect OS platform
#   is_windows = .Platform$OS.type == "windows"
#   #print (use_parallel)
#   #print (cores)
#   if(bval > sm){
#     #if(multicore){
#     if(use_parallel && cores > 1){
#       #tot_cores = parallel::detectCores()
#       #message(sprintf("Running in parallel using %d cores.", cores))
#       if(is_windows){
#         # Windows: use parLapply with a PSOCK cluster
#         cl = parallel::makeCluster(cores)
#         on.exit(parallel::stopCluster(cl))  # ensure cleanup
#         bdist = parallel::parLapply(cl, 1:nsim, .compute_bstat, etahat0=etahat0, n=n, sighat=sighat,
#                                     dd=dd, ddwm=ddwm, qv=NULL, qvl=NULL, amat=amat, awmat=awmat,
#                                     bvec=bvec, imat=imat, dv=dv, lams=ps, w=NULL, arp=arp, p=p_optim, phi=phi1, 
#                                     family=family)
#       } else {
#         # Unix/macOS: use mclapply
#         bdist = parallel::mclapply(1:nsim, .compute_bstat, etahat0=etahat0, n=n, sighat=sighat,
#                                    dd=dd, ddwm=ddwm, qv=NULL, qvl=NULL, amat=amat, awmat=awmat,
#                                    bvec=bvec, imat=imat, dv=dv, lams=ps, w=NULL, arp=arp, p=p_optim, phi=phi1, 
#                                    family=family, mc.cores=cores)
#       }
#     }else{
#       #message("Running in serial mode.")
#       bdist = lapply(1:nsim, .compute_bstat, etahat0=etahat0, n=n, sighat=sighat, dd=dd, ddwm=ddwm, qv=NULL, qvl=NULL,
#                      amat=amat, awmat=awmat, bvec=bvec, imat=imat, dv=dv, lams=ps, w=NULL, arp=arp, p=p_optim, phi=phi1, family=family)
#     }
#     bdist = simplify2array(bdist)
#     pval = sum(bdist > bval) / nsim
#   } else {
#     pval = 1
#   }
#   #----------------------
#   #for visualization
#   #----------------------
#   etacomps = etacomps0 = vector("list", length = capl)
#   etahat0_surf = etahat_surf = NULL
#   muhat0_surf = muhat_surf = NULL
#   ahatc_noz = round(ahatc, 6) |> as.matrix()
#   ahatl_noz = round(ahatl, 6) |> as.matrix()
#   #print (ahatc_noz)
#   #print (nz)
#   if(nz > 0){
#     #print (dim(ahatc_noz))
#     ahatc_noz = ahatc_noz[-c(1:nz), ,drop=F] #|> as.vector()
#     ahatl_noz = ahatl_noz[-c(1:nz), ,drop=F] #|> as.vector()
#     #new:
#     # dd_1 = dd_lst[[1]]
#     # dd_1 = dd_1[, -c(1:nz), drop=F]
#     # dd_lst[[1]] = dd_1
#   }
#   for(icomp in 1:capl){
#     x = X[[icomp]]
#     if(NCOL(x) == 1){
#       #print (dim(dd_lst[[icomp]]))
#       #print (dim(ahatc_noz[var_track == icomp]))
#       etahat_icomp = dd_lst[[icomp]] %*% ahatc_noz[var_track == icomp,,drop=F]
#       #print (dim(etahat_icomp))
#       etacomps[[icomp]] = etahat_icomp
#     }
#     if(NCOL(x) == 2){
#       kts = kts_lst[[icomp]]
#       k1 = kts$k1
#       k2 = kts$k2
#       newd = expand.grid(k1, k2)
#       newd = as.matrix(newd)
#       
#       #H0 surface
#       x1p = newd[,1]
#       x2p = newd[,2]
#       if(attr(x, "parametric") == "warped_plane"){
#         xm0p = cbind(1, x1p, x2p, x1p*x2p)
#         if(icomp > 1){
#           xm0p = xm0p[, -1]
#         }
#       }
#       if(attr(x, "parametric") == "linear"){
#         xm0p = cbind(1, x1p, x2p)
#         if(icomp > 1){
#           xm0p = xm0p[, -1]
#         }
#       }
#       #if(nz == 0){
#       #print (dim(xm0p))
#       #print (length(ahatl_noz))
#       #print (icomp)
#       psurf0 = xm0p %*% ahatl_noz[var_track_param == icomp,,drop=F]
#       # } else if (nz > 0){
#       #   psurf0 = xm0p %*% ahatl[-((ncol(xm0p) + 1):(ncol(xm0p) + nz))]
#       # }
#       etahat0_surf = matrix(psurf0, m1, m2)
#       muhat0_surf = family$linkinv(etahat0_surf)
#       
#       #H1 surface
#       pans = makedelta_wps(newd[,1], newd[,2], k1 = k1, k2 = k2, decreasing = c(FALSE, FALSE))
#       ddp = pans$delta
#       if(icomp > 1){
#         ddp = ddp[, -1]
#       }
#       #if(nz == 0){
#       psurf = ddp %*% ahatc_noz[var_track == icomp,,drop=F]
#       # } else {
#       #   psurf = ddp %*% ahatc[-((nc_noz + 1):(nc_noz + nz))]
#       # }
#       etahat_surf = matrix(psurf, m1, m2)
#       muhat_surf = family$linkinv(etahat_surf)
#       etacomps[[icomp]] = etahat_surf
#       etacomps0[[icomp]] = etahat0_surf
#     }
#   }
#   rslt = list(pval = pval, bval = bval, knots = kts, k1=k1, k2=k2, bmat = dd, wmat = wmat, dmat = dmat, ps = ps,
#               xm0 = xm0, x1m = x1m, amat = amat, awmat = awmat, face = face, etahat = etahat, etahat0 = etahat0,
#               etahat0_surf = etahat0_surf, etahat_surf = etahat_surf, muhat0_surf = muhat0_surf, muhat_surf = muhat_surf,
#               sighat = sighat, edf = edf, edfu = edfu_use, lams = lams, gcvs = gcvs, edfs = edfs, ahatl = ahatl, ahatc = ahatc,
#               phi = phi1, covmat = covmat, covmat_inv = covmat_inv, etacomps = etacomps, p_optim = p_optim, kts_lst = kts_lst, 
#               var_track = var_track, var_track_param = var_track_param, etacomps0 = etacomps0)
#   if(nz >= 1){
#     #assume y is iid with common sig2
#     if(!wt.iter){
#       #should work for ar(p)?
#       if(!arp || arp & p_optim == 0){
#         covmat0 = sighat^2 * mat1 %*% t(mat1)
#       } else {
#         #print (sig2hat_z)
#         #covmat0 = sig2hat_z * mat1 %*% t(mat1)
#         covmat0 = mat1 %*% t(mat1) #seems working with unscaled covariance
#         #covmat0 = sighat^2 * mat1 %*% t(mat1)
#       }
#     } else {
#       cicfamily = CicFamily(family)
#       wt.fun = cicfamily$wt.fun
#       wt = wt.fun(y, etahat, n, weights = rep(1, n), fml = family$family)
#       covmat0 = mat1 %*% diag(wt) %*% t(mat1)
#     }
#     covmat = covmat0[(1):(nz), (1):(nz), drop=FALSE]
#     sez = sqrt(diag(covmat))
#     zcoefs = ahatc[(1):(nz)]
#     tz = zcoefs / sez
#     cpar = 1.2
#     #test more!
#     if(!wt.iter){
#       if ((n - cpar * edf) <= 0) {
#         pz = 2*(1 - pt(abs(tz), edf))
#       } else {
#         #why not pnorm?
#         pz = 2*(1 - pt(abs(tz), n - cpar * edf))
#       }
#     } else {
#       #why not pnorm?
#       pz = 2*(1 - pnorm(abs(tz)))
#     }
#     rslt$sez = sez
#     rslt$pz = pz
#     rslt$tz = tz
#     rslt$zcoefs = zcoefs
#   } else {rslt$sez = NULL; rslt$pz = NULL; rslt$tz = NULL; rslt$zcoefs = NULL}
#   return(rslt)
# }

#####################################################################
#for parallel
#####################################################################
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

#####################################################################
#make the 1st derivative of the basis(edge) of the quadratic bspline at any point xi in the support of the predictor x
#xi is one point and its 1st derivative of a basis(edge) is not zero iff it falls between two knots 
#####################################################################
bqsplfirderi = function(x, xi, knots) {
  ############
  #make knots#   
  ############
  #if (is.null(knots)) {
  #knots are equally spaced on the support of min(x) to max(x)
  #	knots = 0:(m - 1) / (m - 1) * (max(x) - min(x)) + min(x)
  #t = 1:(m + 4) * 0
  #t[1] = t[2] = min(x)
  #t[m+3] = t[m+4] = max(x)
  #t[3:(m+2)] = knots
  #}
  m = length(knots)
  t = 1:(m+4) * 0
  t[1] = t[2] = min(x)
  t[m+3] = t[m+4] = max(x)
  t[3:(m+2)] = knots
  ################################################################################################
  #make bqsplines' 1st derivatives at any point in the support; row: some point xi; column: basis#
  ################################################################################################
  firder = matrix(0, nrow = 1, ncol = (m + 1))
  
  #1st edge
  bool = xi <= t[4] & xi >= t[3]
  if (bool) {	
    firder[,1] = 2 * (xi - t[4]) / (t[4] - t[2]) / (t[4] - t[3]) 
  }
  
  #2nd edge 1st part
  bool = xi <= t[4] & xi >= t[3]
  if (bool) {
    firder[,2] = 2 * ((t[4] - xi) / (t[4] - t[2]) / (t[4] - t[3]) - (xi - t[3]) / (t[5] - t[3]) / (t[4] - t[3]))
  }
  
  #2nd edge 2nd part
  bool = xi <= t[5] & xi >= t[4]
  if (bool) {
    firder[,2] = 2 * (xi - t[5]) / (t[5] - t[3]) / (t[5] - t[4])
  }
  
  #3rd edge to the (m-1)th edge
  for(i in 3:(m - 1)){
    #1st part
    bool = xi <= t[i+1] & xi >= t[i]
    if (bool) {
      firder[,i] = 2 * (xi - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i])
    }
    
    #2nd part
    bool = xi <= t[i+2] & xi >= t[i+1]
    if (bool) {
      firder[,i] = 2 * ((t[i+2] - xi) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) - (xi - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1]))
    }
    
    #3rd part
    bool = xi <= t[i+3] & xi >= t[i+2]
    if (bool) {
      firder[,i] = 2 * (xi - t[i+3]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
    }		
  }
  
  #the mth edge 1st part
  bool = xi <= t[m+1] & xi >= t[m]
  if (bool) {
    firder[,m] = 2 * (xi - t[m]) / (t[m+2] - t[m]) / (t[m+1] - t[m])
  }
  
  #the mth edge 2nd part 
  bool = xi <= t[m+2] & xi >= t[m+1]
  if (bool) {
    firder[,m] = 2 * ((t[m+2] - xi) / (t[m+2] - t[m+1]) / (t[m+2] - t[m]) - (xi - t[m+1]) / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1])) 
  }
  
  #the (m+1)th edge 
  bool = xi <= t[m+2] & xi >= t[m+1]
  if (bool) {
    firder[,m+1] = 2 * (xi - t[m+1]) / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1])
  }
  
  firder
}

##########################################################################################################################################
#make the 2nd derivative of the basis(edge) of the cubic bspline at any point xi in the support of the predictor x
#xi is one point and its 2nd derivative of a basis(edge) is not zero iff it falls between two knots 
##########################################################################################################################################
bcsplsecder = function(x, xi, knots) {
  ############
  #make knots#  
  ############
  #if (is.null(knots)) {
  #	knots = 0:(m - 1) / (m - 1) * (max(x) - min(x)) + min(x)
  #	t = 1:(m + 6)*0
  #	t[1:3] = min(x)
  #	t[(m+4):(m+6)] = max(x)
  #	t[4:(m+3)] = knots
  #}
  m = length(knots) 
  t = 1:(m + 6)*0
  t[1:3] = min(x)
  t[(m+4):(m+6)] = max(x)
  t[4:(m+3)] = knots
  #####################################################################
  #make bcsplines' 2nd derivatives at knots. row: knots; column: basis#
  #####################################################################
  secder =  matrix(0, nrow = 1, ncol = (m + 2)) 
  
  #1st edge: t4 and t5 #checked
  bool = xi < t[5] & xi >= t[4]
  #checked
  if (bool) {	
    secder[,1] = 6 * (t[5] - xi) / (t[5] - t[2]) / (t[5] - t[3]) / (t[5] - t[4])
  }
  
  #2nd edge: t4, t5 and t6 #checked
  bool = xi <= t[5] & xi >= t[4]
  #checked
  if (bool) {
    secder[,2] = 6 * ((xi - t[4]) / (t[6] - t[3]) / (t[6] - t[4]) / (t[5] - t[4]) - (t[5] - xi) / (t[6] - t[3]) / (t[5] - t[3]) / (t[5] - t[4]) - (t[5] - xi) / (t[5] - t[2]) / (t[5] - t[3]) / (t[5] - t[4])) 
  }	
  bool = xi <= t[6] & xi > t[5]
  #checked
  if (bool) {
    secder[,2] =  6 * (t[6] - xi) / (t[6] - t[3]) / (t[6] - t[4]) / (t[6] - t[5])  
  }
  
  #3rd edge: t4, t5, t6 and t7 #checked
  bool =  xi <= t[5] & xi >= t[4]
  #checked
  if (bool) {
    secder[,3] = -6 * ((xi - t[4]) / (t[7] - t[4]) / (t[6] - t[4]) / (t[5] - t[4]) + (xi - t[4]) / (t[6] - t[3]) / (t[6] - t[4]) / (t[5] - t[4]) + (xi - t[5]) / (t[6] - t[3]) / (t[5] - t[3]) / (t[5] - t[4]))
  }	
  bool = xi <= t[6] & xi > t[5]
  #checked
  if (bool) {
    secder[,3] = 6 * ((xi - t[5]) / (t[7] - t[4]) / (t[7] - t[5]) / (t[6] - t[5]) + (xi - t[6]) / (t[7] - t[4]) / (t[6] - t[4]) / (t[6] - t[5]) + (xi - t[6]) / (t[6] - t[3]) / (t[6] - t[4]) / (t[6] - t[5]))
  }
  bool = xi <= t[7] & xi > t[6]
  #checked 
  if (bool) {
    secder[,3] = 6 * (t[7] - xi) / (t[7] - t[4]) / (t[7] - t[5]) / (t[7] - t[6])
  }
  
  #4th  to (m-1)th edge: t[i], t[i+1], t[i+2], t[i+3] and t[i+4] #checked
  if (m > 4) {
    for (i in 4:(m-1)) {
      bool = xi <= t[i+1] & xi >= t[i]
      #checked
      if (bool) {
        secder[,i] = 6 * (xi - t[i]) / (t[i+3] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i])
      }	
      bool = xi <= t[i+2] & xi > t[i+1]
      #checked
      if (bool) {
        secder[,i] = -6 * ((xi - t[i+1]) / (t[i+4] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1]) + (xi - t[i+1]) / (t[i+3] - t[i]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1]) + (xi - t[i+2]) / (t[i+3] - t[i]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]))
      }	
      bool = xi <= t[i+3] & xi > t[i+2] 
      #checked
      if (bool) {
        secder[,i] = 6 * ((xi - t[i+2]) / (t[i+4] - t[i+1]) / (t[i+4] - t[i+2]) / (t[i+3] - t[i+2]) + (xi - t[i+3]) / (t[i+4] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2]) + (xi - t[i+3]) / (t[i+3] - t[i]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2]))
      }
      bool = xi <= t[i+4] & xi > t[i+3]
      #checked
      if (bool) {
        secder[,i] = 6 * (t[i+4] - xi) / (t[i+4] - t[i+1]) / (t[i+4] - t[i+2]) / (t[i+4] - t[i+3]) 
      }
    }
  }
  
  #mth edge: t[m], t[m+1], t[m+2] and t[m+3] #checked
  bool = xi <= t[m+1] & xi >= t[m]
  #checked
  if (bool) {
    secder[,m] = 6 * (xi - t[m]) / (t[m+3] - t[m]) / (t[m+2] - t[m]) / (t[m+1] - t[m])
  }	
  bool = xi <= t[m+2] & xi > t[m+1]
  #checked
  if (bool) {
    secder[,m] = -6 * ((xi - t[m+1]) / (t[m+4] - t[m+1]) / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1]) + (xi - t[m+1]) / (t[m+3] - t[m]) / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1]) + (xi - t[m+2]) / (t[m+3] - t[m]) / (t[m+2] - t[m]) / (t[m+2] - t[m+1]))
  }	
  bool = xi <= t[m+3] & xi > t[m+2] 
  #checked
  if (bool) {
    secder[,m] = 6 * ((xi - t[m+2]) / (t[m+4] - t[m+1]) / (t[m+4] - t[m+2]) / (t[m+3] - t[m+2]) + (xi - t[m+3]) / (t[m+4] - t[m+1]) / (t[m+3] - t[m+1]) / (t[m+3] - t[m+2]) + (xi - t[m+3]) / (t[m+3] - t[m]) / (t[m+3] - t[m+1]) / (t[m+3] - t[m+2]))
  }	
  
  #(m+1)th edge: t[m+1], t[m+2] and t[m+3] #checked
  bool = xi <= t[m+2] & xi >= t[m+1]
  #checked
  if (bool) {
    secder[,m+1] = 6 * (xi - t[m+1]) / (t[m+4] - t[m+1]) / (t[m+3] - t[m+1]) / (t[m+2] - t[m+1])
  }
  bool = xi <= t[m+3] & xi > t[m+2]
  #checked
  if (bool) {
    secder[,m+1] = -6 * ((xi - t[m+2]) / (t[m+5] - t[m+2]) / (t[m+4] - t[m+2]) / (t[m+3] - t[m+2]) + (xi - t[m+2]) / (t[m+4] - t[m+1]) / (t[m+4] - t[m+2]) / (t[m+3] - t[m+2]) + (xi - t[m+3]) / (t[m+4] - t[m+1]) / (t[m+3] - t[m+1]) / (t[m+3] - t[m+2]))
  }
  
  #(m+2)th edge: t[m+2] and t[m+3]
  bool = xi <= t[m+3] & xi >= t[m+2]
  #checked
  if (bool) {
    secder[,m+2] = 6 * (xi - t[m+2]) / (t[m+5] - t[m+2]) / (t[m+4] - t[m+2]) / (t[m+3] - t[m+2])
  }	
  secder
}

###########################
#bcspline first derivative#
###########################
bcsplfirderi = function(x, knots = NULL, xmin = 0, xmax = 1) {
  m = length(knots) 
  t = 1:(m + 6)*0
  t[1:3] = xmin #min(x)
  t[(m+4):(m+6)] = xmax #max(x)
  t[4:(m+3)] = knots
  
  nx = length(x)
  firder =  matrix(0, nrow = nx, ncol = (m + 2)) 
  
  #1st edge: t4 and t5 
  i = 1; j = i+1
  bool = x < t[5] & x >= t[4]
  firder[bool, 1] = -3 * (t[5] - x[bool])^2 / (t[5] - t[2]) / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
  
  #2nd edge: t4, t5 and t6 
  i = 2; j = i+1
  bool = x <= t[5] & x >= t[4]
  #if (bool) { 
  a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
  b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
  a1 = -2 * (t[i+3] - x[bool]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
  b1 = ((t[j+2] - x[bool])- (x[bool] - t[j])) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + ((t[j+3] - x[bool]) - (x[bool] - t[j+1])) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1]) 
  firder[bool, 2] = a / (t[i+3] - t[i]) + a1 * (x[bool] - t[i]) / (t[i+3] - t[i]) - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
  #}	
  bool = x <= t[6] & x > t[5]
  #if (bool) {
  b = (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
  b1 =  -2 * (t[j+3] - x[bool]) / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
  firder[bool, 2] = -b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
  #}
  
  #3rd edge: t4, t5, t6 and t7 
  i = 3; j = i+1
  bool =  x <= t[5] & x >= t[4]
  #if (bool) {
  a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
  b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j]) 
  a1 = ((t[i+2] - x[bool]) - (x[bool] - t[i])) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + ((t[i+3] - x[bool]) - (x[bool] - t[i+1])) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
  b1 = 2 * (x[bool] - t[j]) / (t[j+2] - t[j]) / (t[j+1] - t[j])  
  firder[bool, 3] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1 
  #}	
  bool = x <= t[6] & x > t[5]
  #if (bool) {
  a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
  b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
  a1 =  -2 * (t[i+3] - x[bool]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
  b1 = ((t[j+2] - x[bool]) - (x[bool] - t[j])) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + ((t[j+3] - x[bool]) - (x[bool] - t[j+1])) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])   
  firder[bool, 3] =  a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1 
  #}
  bool = x <= t[7] & x > t[6]
  #if (bool) {
  b = (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
  b1 = -2 *(t[j+3] - x[bool]) / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
  firder[bool, 3] = -b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
  #}
  
  #4th  to (m-1)th edge: t[i], t[i+1], t[i+2], t[i+3] and t[i+4] 
  if (m > 4) {
    for (i in 4:(m-1)) {
      j = i+1
      bool = x <= t[i+1] & x >= t[i]
      #if (bool) {
      a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
      a1 = 2 * (x[bool] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i])
      firder[bool, i] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1
      #}	
      bool = x <= t[i+2] & x > t[i+1]
      #if (bool) {
      a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
      b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
      a1 =  ((t[i+2] - x[bool]) - (x[bool] - t[i])) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + ((t[i+3] - x[bool]) - (x[bool] - t[i+1])) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
      b1 = 2 * (x[bool] - t[j]) / (t[j+2] - t[j]) / (t[j+1] - t[j])
      firder[bool, i] =  a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
      #}	
      bool = x <= t[i+3] & x > t[i+2] 
      #if (bool) {
      a =  (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
      b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
      a1 = -2 * (t[i+3] - x[bool]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
      b1 = ((t[j+2] - x[bool]) - (x[bool] - t[j])) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + ((t[j+3] - x[bool]) - (x[bool] - t[j+1])) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
      firder[bool, i] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
      #}
      bool = x <= t[i+4] & x > t[i+3]
      #if (bool) {
      b = (t[j+3] - x[bool])**2 / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
      b1 = -2 * (t[j+3] - x[bool]) / (t[j+3] - t[j+1]) / (t[j+3] - t[j+2])
      firder[bool, i] =  -b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
      #}
    }
  }
  #mth edge: t[m], t[m+1], t[m+2] and t[m+3] #checked
  i = m; j = i+1
  bool = x <= t[m+1] & x >= t[m]
  #if (bool) {
  a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
  a1 = 2 * (x[bool] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i])
  firder[bool, m] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1
  #}	
  bool = x <= t[m+2] & x > t[m+1]
  #if (bool) {
  a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
  b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
  a1 =  ((t[i+2] - x[bool]) - (x[bool] - t[i])) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + ((t[i+3] - x[bool]) - (x[bool] - t[i+1]) ) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1]) 
  b1 = 2 * (x[bool] - t[j]) / (t[j+2] - t[j]) / (t[j+1] - t[j])
  firder[bool, m] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1
  #}	
  bool = x <= t[m+3] & x > t[m+2] 
  #if (bool) {
  a = (t[i+3] - x[bool])**2 / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
  b = (x[bool] - t[j]) * (t[j+2] - x[bool]) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + (t[j+3] - x[bool]) * (x[bool] - t[j+1]) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
  a1 =  -2 * (t[i+3] - x[bool]) / (t[i+3] - t[i+1]) / (t[i+3] - t[i+2])
  b1 = ((t[j+2] - x[bool]) - (x[bool] - t[j])) / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) + ((t[j+3] - x[bool]) - (x[bool] - t[j+1])) / (t[j+3] - t[j+1]) / (t[j+2] - t[j+1])
  firder[bool, m] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1  - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1  
  #}	
  
  #(m+1)th edge: t[m+1], t[m+2] and t[m+3] 
  i = m+1; j = i+1
  bool = x <= t[m+2] & x >= t[m+1]
  #if (bool) {
  a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
  a1 = 2 * (x[bool] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i])
  firder[bool, m+1] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1
  #}
  bool = x <= t[m+3] & x > t[m+2]
  #if (bool) {
  a = (x[bool] - t[i]) * (t[i+2] - x[bool]) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + (t[i+3] - x[bool]) * (x[bool] - t[i+1]) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
  b = (x[bool] - t[j])**2 / (t[j+2] - t[j]) / (t[j+1] - t[j])
  a1 = ((t[i+2] - x[bool]) - (x[bool] - t[i])) / (t[i+2] - t[i]) / (t[i+2] - t[i+1]) + ((t[i+3] - x[bool]) - (x[bool] - t[i+1])) / (t[i+3] - t[i+1]) / (t[i+2] - t[i+1])
  b1 = 2 * (x[bool] - t[j]) / (t[j+2] - t[j]) / (t[j+1] - t[j])
  firder[bool, m+1] = a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1 - b / (t[i+4] - t[i+1]) + (t[i+4] - x[bool]) / (t[i+4] - t[i+1]) * b1  
  #}
  
  #(m+2)th edge: t[m+2] and t[m+3]
  i = m+2
  bool = x <= t[m+3] & x >= t[m+2]
  #if (bool) {
  a = (x[bool] - t[i])**2 / (t[i+2] - t[i]) / (t[i+1] - t[i])
  a1 = 2 * (x[bool] - t[i]) / (t[i+2] - t[i]) / (t[i+1] - t[i]) 
  firder[bool, m+2] =  a / (t[i+3] - t[i]) + (x[bool] - t[i]) / (t[i+3] - t[i]) * a1
  #}	
  firder
}

#------------------------------------------
#create penalty matrix
#------------------------------------------
# makedmat_1D = function(m1)
# {
#   #m1 = length(k1)
#   dmat = matrix(0, nrow = (m1-1), ncol = m1)
#   for (i in 2:m1) {
#     dmat[i-1, i-1] = 1; dmat[i-1, i] = -1
#   }
#   return (dmat)
# }

makedmat_1D = function(m, q = 2) {
  #dmat = NULL
  dmat = matrix(0, nrow = (m - q), ncol = m)
  
  # zero order
  if (q == 0) {
    # for (i in 1:m) {
    #   dmat[i, i] = 1
    # }
    diag(dmat) = 1
  }
  
  # first order
  if (q == 1) {
    for (i in 2:m) {
      #dmat[i-1, i-1] = 1; dmat[i-1, i] = -1
      dmat[i-1, i-1] = -1; dmat[i-1, i] = 1
    }
    #id1 = 1:(m-1)
    #id2 = 2:m
    #dmat[cbind(id1, id1)] = 1
    #dmat[cbind(id1, id2)] = -1
  }
  
  # second order
  if (q == 2) {
    for (i in 3:m) {
      dmat[i-2, i-2] = 1; dmat[i-2, i-1] = -2; dmat[i-2, i] = 1
    }
    #id1 = 1:(m-2)
    #id2 = 2:(m-1)
    #id3 = 3:m
    #dmat[cbind(id1, id1)] = 1
    #dmat[cbind(id1, id2)] = -2
    #dmat[cbind(id1, id3)] = 1
  }
  
  # third-order
  if (q == 3) {
    for (i in 4:m) {
      dmat[i-3, i-3] = 1; dmat[i-3, i-2] = -3; dmat[i-3, i-1] = 3; dmat[i-3, i] = -1
    }
    # id1 = 1:(m-3)
    # id2 = 2:(m-2)
    # id3 = 3:(m-1)
    # id4 = 4:m
    # dmat[cbind(id1, id1)] = 1
    # dmat[cbind(id1, id2)] = -3
    # dmat[cbind(id1, id3)] = 3
    # dmat[cbind(id1, id4)] = -1
  }
  
  # fourth-order (new)
  if (q == 4) {
    dm3 = matrix(0, nrow = (m - 3), ncol = m)
    #third-order
    for (i in 4:m) {
      dm3[i-3, i-3] = 1; dm3[i-3, i-2] = -3; dm3[i-3, i-1] = 3; dm3[i-3, i] = -1
    }
    
    # id1 = 1:(m-3)
    # id2 = 2:(m-2)
    # id3 = 3:(m-1)
    # id4 = 4:m
    # dm3[cbind(id1, id1)] = 1
    # dm3[cbind(id1, id2)] = -3
    # dm3[cbind(id1, id3)] = 3
    # dm3[cbind(id1, id4)] = -1
    
    dm1 = matrix(0, nrow = (m - 4), ncol = (m - 3))
    for(i in 2:(m - 3)){
      dm1[i-1, i-1] = 1; dm1[i-1, i] = -1
    }
    
    # id1 = 1:(m - 4)
    # id2 = 2:(m - 3)
    # dm1[cbind(id1, id1)] = 1
    # dm1[cbind(id1, id2)] = -1
    
    dmat = dm1 %*% dm3
  }
  
  return (dmat)
}

# 
# makedmat_1D = function(m, q = 2) {
#   #dmat = NULL
#   dmat = matrix(0, nrow = (m - q), ncol = m)
#   #  third-order
#   if (q == 3) {
#     for (i in 4:m) {
#       dmat[i-3, i-3] = 1; dmat[i-3, i-2] = -3; dmat[i-3, i-1] = 3; dmat[i-3, i] = -1
#     }
#   }
#   # second order
#   if (q == 2) {
#     for (i in 3:m) {
#       dmat[i-2, i-2] = 1; dmat[i-2, i-1] = -2; dmat[i-2, i] = 1
#     }
#   }
#   # first order
#   if (q == 1) {
#     for (i in 2:m) {
#       dmat[i-1, i-1] = 1; dmat[i-1, i] = -1
#     }
#   }
#   # zero order
#   if (q == 0) {
#     for (i in 1:m) {
#       dmat[i, i] = 1
#     }
#   }
#   
#   # fourth-order (new)
#   if (q == 4) {
#     dm3 = matrix(0, nrow = (m - 3), ncol = m)
#     #  third-order
#     for (i in 4:m) {
#       dm3[i-3, i-3] = 1; dm3[i-3, i-2] = -3; dm3[i-3, i-1] = 3; dm3[i-3, i] = -1
#     }
#     dm1 = matrix(0, nrow = (m - 4), ncol = (m - 3))
#     for(i in 2:(m - 3)){
#       dm1[i-1, i-1] = 1; dm1[i-1, i] = -1
#     }
#     dmat = dm1 %*% dm3
#   }
#   return (dmat)
# }

#------------------------------------------
#create penalty matrix for wps
#------------------------------------------
makedmat_2D = function(k1, k2)
{
  m1 = length(k1)
  m2 = length(k2)
  dmat = matrix(0, nrow = 2*(m1 * m2 - m1 - m2), ncol = m1*m2)
  irow = 0
  for(j in 1:m2){
    for(i in 1:(m1-2)){
      irow=irow+1
      dmat[irow,m2*(i+1)+j]=1/(k1[i+2]-k1[i+1])
      dmat[irow,m2*i+j]=-1/(k1[i+2]-k1[i+1])-1/(k1[i+1]-k1[i])
      dmat[irow,m2*(i-1)+j]=1/(k1[i+1]-k1[i])
    }
  }
  for(j in 1:(m2-2)){
    for(i in 1:m1){
      irow=irow+1
      dmat[irow,m2*(i-1)+j+2]=1/(k2[j+2]-k2[j+1])
      dmat[irow,m2*(i-1)+j+1]=-1/(k2[j+2]-k2[j+1])-1/(k2[j+1]-k2[j])
      dmat[irow,m2*(i-1)+j]=1/(k2[j+1]-k2[j])
    }
  }
  return (dmat)
}

#################################################################
#1D amat: cubic or quadratic splines; different from cgam bigdata
#used for splines
#################################################################
makeamat_1D = function(shp = 9, spl_ord = 3L, kts, x, x1, xn){
  sgn = 1
  #decr, conc
  if(shp %in% c(10, 12, 14, 16)){
    sgn = -1
  }
  if(spl_ord == 3L) {
    #cat('call makeamat_1D','\n')
    #print (kts)
    amat = Reduce(rbind, lapply(kts, function(e){bqsplfirderi(x, e, knots = kts)}))
    #print (amat)
    #nkts = length(kts)
    #amat = splines::splineDesign(knots = c(x1, x1, kts, xn, xn), x = kts, ord = spl_ord, derivs = rep(1, nkts), outer.ok = TRUE)
  }
  if(spl_ord == 4L) {
    #nkts = length(kts)
    #amat = splines::splineDesign(knots = c(x1, x1, kts, xn, xn), x = kts, ord = spl_ord, derivs = rep(2, nkts), outer.ok = TRUE)
    #x1_add = splines::splineDesign(knots = c(x1, x1, kts, xn, xn), x = x1, ord = spl_ord, derivs = 1, outer.ok = TRUE)
    #xn_add = splines::splineDesign(knots = c(x1, x1, kts, xn, xn), x = xn, ord = spl_ord, derivs = 1, outer.ok = TRUE)
    amat = Reduce(rbind, lapply(kts, function(e){bcsplsecder(x, e, knots = kts)}))
    # amat = NULL
    # for(kti in kts) {
    #   ami = bcsplsecder(x, kti, knots = kts)
    #   amat = rbind(amat, ami)
    # }
    x1_add = bcsplfirderi(x1, knots = kts, xmin = x1, xmax = xn)
    xn_add = bcsplfirderi(xn, knots = kts, xmin = x1, xmax = xn)
  }
  amat = sgn * amat
  #incr.conv
  if(shp == 13){
    amat = rbind(amat, x1_add)
  }
  #incr.conc
  if(shp == 14){
    amat = rbind(amat, xn_add)
  }
  #decr.conv
  if(shp == 15){
    amat = rbind(amat, -xn_add)
  }
  #decr.conc
  if(shp == 16){
    amat = rbind(amat, -x1_add)
  }
  return (amat)
}

#######################################################################################
#1D amat for xmat'beta
#2D amat when we have wps; singular problem; use_constreg=TRUE
#######################################################################################
makeamat_param = function(shp = c(1, 1), parametric = "linear", x1=NULL, xn=NULL){
  awmat = NULL
  if(length(shp) == 1){
    sgn = 1
    #use 1,2,3,4; don't want to add another attribute in s.incr; shp - 8
    shp = shp - 8
    if(shp %in% c(2, 4)){
      sgn = -1
    }
    
    if(parametric == "linear") {
      if(shp %in% c(1, 2)){
        awmat = matrix(0, nrow=1, ncol=2)
        awmat[1, 2] = 1
        awmat = sgn * awmat
      } #else {
      #awmat = matrix(0, nrow=1, ncol=2) #avoid error in test.fit when do awmat[,-1]
      #}
      
      #unconstrained: call lm
      if(shp %in% c(3, 4)){
        awmat = matrix(0, nrow=1, ncol=2)
      }
      
      #incr.conv 5 | incr.conc 6
      if(shp %in% c(5, 6)){
        awmat = matrix(0, nrow=1, ncol=2)
        awmat[1, 2] = 1
      }
      
      #decr.conv 7 | decr.conc 8
      if(shp %in% c(7, 8)){
        awmat = matrix(0, nrow=1, ncol=2)
        awmat[1, 2] = -1
      }
      
      #incr.conc 7
      # if(shp == 7){
      #   xn_add = c(0, 1, 2*xn)
      #   awmat = rbind(awmat, xn_add)
      # }
      
      #decr.conv 6
      # if(shp == 6){
      #   xn_add = c(0, 1, 2*xn)
      #   awmat = rbind(awmat, -xn_add)
      # }
      
      #decr.conc 8
      # if(shp == 8){
      #   x1_add = c(0, 1, 2*x1)
      #   awmat = rbind(awmat, -x1_add)
      # }
    }
    
    if(parametric == "quadratic") {
      if(shp %in% c(1, 2)){
        awmat = matrix(0, nrow=2, ncol=3)
        awmat[1, 2] = 1
        awmat[1, 3] = 2*x1
        awmat[2, 2] = 1
        awmat[2, 3] = 2*xn
        awmat = sgn * awmat
      }
      
      if(shp %in% c(3, 4)){
        awmat = matrix(0, nrow=1, ncol=3)
        awmat[1, 3] = 1
        awmat = sgn * awmat
      }
      
      #incr.conv 5
      if(shp == 5){
        awmat = matrix(0, nrow=1, ncol=3)
        awmat[1, 3] = 1
        x1_add = c(0, 1, 2*x1)
        awmat = rbind(awmat, x1_add)
      }
      #incr.conc 6
      if(shp == 6){
        awmat = matrix(0, nrow=1, ncol=3)
        awmat[1, 3] = 1
        awmat = -awmat
        xn_add = c(0, 1, 2*xn)
        awmat = rbind(awmat, xn_add)
      }
      #decr.conv 7
      if(shp == 7){
        awmat = matrix(0, nrow=1, ncol=3)
        awmat[1, 3] = 1
        xn_add = c(0, 1, 2*xn)
        awmat = rbind(awmat, -xn_add)
      }
      #decr.conc 8
      if(shp == 8){
        awmat = matrix(0, nrow=1, ncol=3)
        awmat[1, 3] = 1
        awmat = -awmat
        x1_add = c(0, 1, 2*x1)
        awmat = rbind(awmat, -x1_add)
      }
    }
    
    #new: cubic
    #no way to get increasing/decreasing when we have x^3
    #only shapes are 3 or 4
    #have to compute wmat and then get BW and AW
    if(parametric == "cubic"){
      if(shp %in% c(3, 5, 7)){
        awmat = matrix(0, nrow=2, ncol=4)
        awmat[1, 3] = 1
        awmat[1, 4] = 3*x1
        awmat[2, 3] = 1
        awmat[2, 4] = 3*xn
        
        #incr.conv 5
        if(shp == 5){
          x1_add = c(0, 1, 2*x1, 3*x1^2)
          awmat = rbind(awmat, x1_add)
        }
        
        #decr.conv 7
        if(shp == 7){
          xn_add = c(0, 1, 2*xn, 3*xn^2)
          awmat = rbind(awmat, -xn_add)
        }
      } 
      
      if(shp %in% c(4, 6, 8)){
        awmat = matrix(0, nrow=2, ncol=4)
        awmat[1, 3] = 1
        awmat[1, 4] = 3*x1
        awmat[2, 3] = 1
        awmat[2, 4] = 3*xn
        awmat = -awmat
        
        #incr.conc 6
        if(shp == 6){
          xn_add = c(0, 1, 2*xn, 3*xn^2)
          awmat = rbind(awmat, xn_add)
        }
        #decr.conc 8
        if(shp == 8){
          x1_add = c(0, 1, 2*x1, 3*x1^2)
          awmat = rbind(awmat, -x1_add)
        }
      }
      
      if (shp %in% c(1, 2)){
        awmat = matrix(0, nrow=1, ncol=4)
      }
    }
  }
  
  if(length(shp) == 2) {
    sgn = c(1, 1)
    if(shp[1] %in% c(2, 4)){
      sgn[1] = -1
    }
    if(shp[2] %in% c(2, 4)){
      sgn[2] = -1
    }
    #incr/decr
    if(shp[1] %in% c(1, 2) & shp[2] %in% c(1, 2)){
      if(parametric == "linear"){
        awmat = matrix(0, nrow = 2, ncol = 3)
        awmat[1, 2] = awmat[2, 3] = 1
        awmat[1, ] = awmat[1, ] * sgn[1]
        awmat[2, ] = awmat[2, ] * sgn[2]
      }
      
      if(parametric == "warped_plane"){ #x1 and x2 are scaled to be within [0,1]
        awmat = matrix(0, nrow = 4, ncol = 4)
        awmat[1, 2] = 1 * sgn[1]
        awmat[2, c(2, 4)] = 1 * sgn[1]
        awmat[3, 3] = 1 * sgn[2]
        awmat[4, 3:4] = 1 * sgn[2]
      }
    }
    #convex/concave: not sure what to do with s.conc.conv: no such a shape!
    #not implemented in testpar.fit 
    #b0+b1x1+b2x2+b3x1^2+b4x2^2: should be quadratic not linear
    #b0+b1x1+b2x2+b3x1^2+b4x2^2+b5x1x2
    if(shp[1] %in% c(3, 4) & shp[2] %in% c(3, 4)){
      # if(parametric == "linear"){
      #   awmat = matrix(0, nrow = 2, ncol = 5)
      #   awmat[1, 3] = awmat[2, 5] = 1
      #   awmat[1, ] = awmat[1, ] * sgn[1]
      #   awmat[2, ] = awmat[2, ] * sgn[2]
      # }
      
      #not test
      if(parametric == "quadratic"){
        awmat = matrix(0, nrow = 2, ncol = 5)
        awmat[1, 4] = awmat[2, 5] = 1
        awmat[1, ] = awmat[1, ] * sgn[1]
        awmat[2, ] = awmat[2, ] * sgn[2]
      }
      
      #should have another name....
      # if(parametric == "warped_plane"){
      #   awmat = matrix(0, nrow = 3, ncol = 6)
      #   awmat[1, 3] = awmat[2, 5] = awmat[3, 6] = 1
      #   awmat[1, ] = awmat[1, ] * sgn[1]
      #   awmat[2, ] = awmat[2, ] * sgn[2]
      #   awmat[3, ] = awmat[3, ] * sgn[3]
      # }
    } 
  }
  return (awmat)
}

###############################################################
#2D amat: wps splines used; all pairs of 1:8 are considered
###############################################################
makeamat_nonadd = function(k1, k2, shp_pr){
  m1 = length(k1)
  m2 = length(k2)
  
  #the code constrain the vertical/j dimension first
  shp_pr = rev(shp_pr)
  #amat = NULL
  if(shp_pr[1] %in% c(3, 4)) {
    nr1 = m2*(m1-2)  
  } else {
    nr1 = m2*(m1-1)
  }
  
  if(shp_pr[2] %in% c(3, 4)) {
    nr2 = m1*(m2-2)  
  } else {
    nr2 = m1*(m2-1)
  }
  
  amat = matrix(0, nrow = (nr1 + nr2), ncol = m1*m2)
  irow = 0
  
  #--------
  #dim 1
  #--------
  if(shp_pr[1] %in% c(1, 2)){
    for(j in 1:(m2-1)){
      for(i in 1:m1){
        irow = irow+1
        amat[irow, m2*(i-1)+j]=-1
        amat[irow, m2*(i-1)+j+1]=1
      }
    }
    if(shp_pr[1] == 2) {amat[1:irow, ] = -amat[1:irow, ]}
  }
  
  if(shp_pr[1] %in% c(3, 4)){
    for(j in 1:(m2-2)){
      for(i in 1:m1){
        irow = irow+1
        amat[irow,m2*(i-1)+j] = 1
        amat[irow,m2*(i-1)+j+1] = -2
        amat[irow,m2*(i-1)+j+2] = 1
      }
    }
    if(shp_pr[1] == 4) {amat[1:irow, ] = -amat[1:irow, ]}
  }
  
  if(shp_pr[1] %in% c(5, 8)){
    for(j in 1:(m2-2)){
      for(i in 1:m1){
        irow = irow+1
        amat[irow, m2*(i-1)+j] = 1
        amat[irow, m2*(i-1)+j+1] = -2
        amat[irow, m2*(i-1)+j+2] = 1
        if(j == 1) {
          irow = irow + 1
          amat[irow, (m2*(i-1)+j):(m2*(i-1)+j+1)] = c(-1, 1)
        }
      }
    }
    if(shp_pr[1] == 8) {amat[1:irow, ] = -amat[1:irow, ]}
  }
  
  if(shp_pr[1] %in% c(6, 7)){
    for(j in 1:(m2-2)){
      for(i in 1:m1){
        irow = irow+1
        amat[irow, m2*(i-1)+j] = 1
        amat[irow, m2*(i-1)+j+1] = -2
        amat[irow, m2*(i-1)+j+2] = 1
        if(j == (m2-2)) {
          irow = irow + 1
          amat[irow, (m2*(i-1)+j+1):(m2*(i-1)+j+2)] = c(1, -1)
        }
      }
    }
    if(shp_pr[1] == 7) {amat[1:irow, ] = -amat[1:irow, ]}
  }
  
  ed_dim1 = irow
  #--------
  #dim 2
  #--------
  if(shp_pr[2] %in% c(1, 2)){
    for(j in 1:m2){
      for(i in 1:(m1-1)){
        irow = irow+1
        amat[irow, m2*(i)+j] = 1
        amat[irow, m2*(i-1)+j] = -1
      }
    }
    if(shp_pr[2] == 2) {amat[(ed_dim1+1):irow, ] = -amat[(ed_dim1+1):irow, ]}
  }
  
  if(shp_pr[2] %in% c(3, 4)){
    for(j in 1:m2){
      for(i in 1:(m1-2)){
        irow = irow+1
        amat[irow, m2*(i-1)+j] = 1
        amat[irow, m2*(i)+j] = -2
        amat[irow, m2*(i+1)+j] = 1
      }
    }
    if(shp_pr[2] == 4) {amat[(ed_dim1+1):irow, ] = -amat[(ed_dim1+1):irow, ]}
  }
  
  if(shp_pr[2] %in% c(5, 8)){
    for(j in 1:m2){
      for(i in 1:(m1-2)){
        irow = irow+1
        amat[irow, m2*(i-1)+j] = 1
        amat[irow, m2*(i)+j] = -2
        amat[irow, m2*(i+1)+j] = 1
        if(i == 1) {
          irow = irow + 1
          amat[irow, m2*(i-1)+j] = -1
          amat[irow, m2*(i)+j] = 1
        }
      }
    }
    if(shp_pr[2] == 8) {amat[(ed_dim1+1):irow, ] = -amat[(ed_dim1+1):irow, ]}
  }
  
  if(shp_pr[2] %in% c(6, 7)){
    for(j in 1:m2){
      for(i in 1:(m1-2)){
        irow = irow+1
        amat[irow, m2*(i-1)+j] = 1
        amat[irow, m2*(i)+j] = -2
        amat[irow, m2*(i+1)+j] = 1
        if(i == (m1-2)) {
          irow = irow + 1
          amat[irow, m2*(i)+j] = 1
          amat[irow, m2*(i+1)+j] = -1
        }
      }
    }
    if(shp_pr[2] == 7) {amat[(ed_dim1+1):irow, ] = -amat[(ed_dim1+1):irow, ]}
  }
  
  return (amat)
}

#------------------------------------------
#one-sided hypothesis test
#not consider weights yet
#need to include the binomial case
#------------------------------------------
.compute_bstat <- function(iloop, etahat0, n, sighat, dd, ddwm, qv=NULL, qvl=NULL, amat, awmat, bvec, 
                           imat, dv, lams, w=NULL, family = gaussian(link = "identity"), 
                           arp = FALSE, p = 1, phi = NULL, covmat=NULL)
{
  muhat0 = etahat0
  if(family$family != "gaussian"){
    muhat0 = family$linkinv(etahat0)
  }
  ysim = ysim.fun(n, mu0 = muhat0, fml = family$family, phi = phi, sd = sighat)
  #----------
  #H0 fit
  #----------
  #bvec will be different if using constreg
  anslsim = fit.hypo(p=p, dd=ddwm, y=ysim, amat=awmat, bvec=rep(0, nrow(awmat)), dv=NULL, family=family, ps=0, arp=arp)
  ahatlsim = anslsim$ahat
  etah0 = ddwm %*% ahatlsim
  dev0sim = anslsim$dev
  
  #----------
  #H1 fit
  #----------
  anssim = fit.hypo(p=p, dd=dd, y=ysim, amat=amat, bvec=bvec, dv=dv, family=family, ps=lams[1], arp=arp)
  ahatsim = anssim$ahat
  etah = dd %*% ahatsim
  dev1sim = anssim$dev

  if(family$family == "gaussian"){
    #new:
    if(arp){
      bstati = (dev0sim - dev1sim) / dev0sim 
    }else{
      if (!is.null(w)) {
        sse0sim = sum(w * (ysim - etah0)^2)
        sse1sim = sum(w * (ysim - etah)^2)
      } else {
        sse0sim = sum((ysim - etah0)^2)
        sse1sim = sum((ysim - etah)^2)
      }
      bstati = (sse0sim - sse1sim) / sse0sim
    }
  } else {
    bstati = (dev0sim - dev1sim) / n
  }
  return (bstati)
}

#Yule-Walker
yw = function(e, p = 1,...){
  n = length(e)
  e = e - mean(e)
  gamma = NULL
  f = function(k, e){mean((e[(k+1):n]) * (e[1:(n-k)]))}
  #1. compute autocovariance
  kvec = 0:p
  gamma = sapply(kvec, f, e)
  #2. set up yw equation and solve for phi
  Gamma = toeplitz(gamma[1:p])
  phi = solve(Gamma, gamma[2:(p + 1)])
  #3. get white noise sigma2
  sigma2 = gamma[1] - drop(crossprod(gamma[2:(p + 1)], phi))
  #4. get covariance matrix for all n data points; ignore gamma 0 - gamma p
  #autocov = function(lag, phi, sigma2, p, n){
  #check more
  if (p == 1) {
    cov_vec = sapply(0:(n - 1), function(lag) sigma2 / (1 - phi^2) * phi^lag)
  } else {
    cov_vec = numeric(n)
    cov_vec[1:(p + 1)] = gamma
    for(i in (p + 2):n) {
      cov_vec[i] = crossprod(phi, cov_vec[(i-1):(i-1-p+1)])
    }
  }
  #}
  #cov_vec = sapply(0:(n - 1), autocov)
  covmat = toeplitz(cov_vec)
  #huan's paper uses correlation matrix, not covariance matrix:
  #no difference? makes a difference when check power
  covmat = covmat/gamma[1]
  #covmat = round(covmat, 10)
  #umat = t(chol(covi)) #L matrix in the paper
  #correction for bias?
  #sigma2 = sigma2/(1-sum(phi))
  #sigma2 = sigma2*n/(n-length(phi))
  rslt = list(gamma = gamma, phi = phi, sigma2 = sigma2, covmat = covmat)
}

#used for uinv of ar(1)
fuinv = function(n, phi, sig2hat = 1){
  #n = length(y)
  uinv = diag(rep(1/sig2hat^.5, n))
  uinv[1,1] = sqrt(1-phi^2)/sig2hat^.5
  uinv[cbind(2:n, 1:(n - 1))] = -phi/sig2hat^.5
  #for(irow in 2:n){uinv[irow,(irow-1)] = -phi/sig2hat^.5}
  return(uinv)
}

#-------------------------------------------------
#for testpar_tobe_in_cgam.R, ps is absorbed in dv
#-------------------------------------------------
fit.hypo = function(p = 0, dd, y, amat, bvec = NULL, dv = NULL, 
                    wmat=NULL, family = gaussian(link='identity'), 
                    ps = 0, arp = FALSE,...){
  wt.iter = ifelse(family$family == "gaussian", FALSE, TRUE)
  covmat = covmat_inv = phi = gamma = sigma2 = sighat = sig2hat_z = xtil = aic = pval.ts = NULL
  n = length(y)
  #print (n)
  if(!wt.iter){
    qv0 = crossprod(dd)
    #print (head(dd))
    if(!is.null(dv)){ #only for H1
      qv = qv0 + dv #ps is absorbed in dv
    } else {
      qv = qv0 # + 1e-4 * diag(nrow(qv0))
    }
    cv = crossprod(dd, y)
    #new:
    #cv = Matrix::crossprod(dd, y)
    #new:
    if(all(amat == 0)){
      ahat = solve(qv, cv)
      face = integer(0)
    } else {
      #ans = qprog(qv, cv, amat, bvec)
      #ahat = ans$thetahat
      #face = ans$face
      
      ans = solve.QP(qv, cv, t(amat), bvec)
      ahat = ans$solution
      face = ans$iact
    }
    
    #ans = solve.QP(qv, cv, t(amat), bvec)
    #ahat = ans$solution
    #face = ans$iact
    
    #new:
    #amat = Matrix::Matrix(amat, sparse = T)
    #ans = solve.QP(qv, cv, t(amat), bvec)
    #ahat = ans$solution
    #face = ans$iact
    etahat = dd %*% ahat
    etahat = as.numeric(etahat)
    sse = sum((y - etahat)^2)
    #print (sse)
    #new:
    pval.ts = Box.test(y-etahat, lag = 1)$p.value
    if(arp & p > 0){
      diff = sum((y-mean(y))^2)
      nrep = 0
      while(diff > 1e-5 & nrep < 20){
        nrep = nrep + 1
        oldeta = etahat
        #step 1: get the error term and estimate phi, gamma, sigma2 in arp
        e = y - oldeta
        ansyw = yw(e, p=p)
        covmat = ansyw$covmat
        phi = ansyw$phi
        sig2hat = ansyw$sigma2 #will have more bias than iid case
        gamma = ansyw$gamma
        
        if(p>1){
          #umat = t(chol(covmat)) #L matrix in the paper
          #uinv = solve(umat)
          
          umat = t(chol(covmat))
          iumat = diag(ncol(umat))
          uinv = forwardsolve(umat, iumat)
        } else if(p==1){
          uinv = fuinv(n = length(y), phi=phi, sig2hat = (sig2hat/gamma[1])) #scaled Sigma mat with gamma[1]
          #n = length(y)
          #uinv = diag(rep(1/sig2hat^.5, n))
          #uinv[1,1] = sqrt(1-phi^2)/sig2hat^.5
          #for(irow in 2:n){uinv[irow,(irow-1)] = -phi/sig2hat^.5}
        }
        
        #uinv = backsolve(umat, iumat)
        ytil = uinv %*% y
        xtil = uinv %*% dd
        cv = crossprod(xtil, ytil)
        qv0 = crossprod(xtil)
        #qv0 = Matrix::crossprod(xtil)
        
        if(!is.null(dv)){ #only for H1
          qv = qv0 + dv #ps is absorbed in dv
        } else {
          qv = qv0
        }
        
        #new:
        if(all(amat == 0)){
          ahat = solve(qv, cv)
          face = integer(0)
        } else {
          #ans = qprog(qv, cv, amat, bvec)
          #ahat = ans$thetahat
          #face = ans$face
          
          ans = solve.QP(qv, cv, t(amat), bvec)
          ahat = ans$solution
          face = ans$iact
        }
        
        #new:
        #ans = solve.QP(qv, cv, t(amat), bvec)
        #ahat = ans$solution
        #face = ans$iact
        
        etahat = dd %*% ahat
        etahat = as.numeric(etahat)
        diff = mean((etahat - oldeta)^2)
      }
      #print (nrep)
      etatil = xtil %*% ahat
      sse = sum((ytil - etatil)^2)
      covmat_inv = crossprod(uinv)
      #covmat_inv = round(covmat_inv, 10)
      
      #uinv_keep = uinv
      #dev0 = sum((y - etahat)^2)
      sig2hat = n/(n-p)*sig2hat
    }
  } else {
    n = length(y)
    muhat = rep.int(1/2, n)
    etahat = family$linkfun(muhat)
    w = 1:n*0 + 1
    gr = gr.fun(y, etahat, weights = w, fml = family$family)
    wt = wt.fun(y, etahat, n, weights = w, fml = family$family)
    cv = crossprod(dd, (wt * etahat - gr))
    qv0 = crossprod(dd)
    
    if(!is.null(dv)){ #only for H1
      qv = qv0 + dv #ps is absorbed in dv
    } else {
      qv = qv0
    }

    #new:
    if(all(amat == 0)){
      ahat = solve(qv, cv)
      face = integer(0)
    } else {
      #ans = qprog(qv, cv, amat, bvec)
      #ahat = ans$thetahat
      #face = ans$face
      
      ans = solve.QP(qv, cv, t(amat), bvec)
      ahat = ans$solution
      face = ans$iact
    }
    
    # ans = solve.QP(qv, cv, t(amat), bvec)
    # ahat = ans$solution
    # face = ans$iact
    
    etahat = dd %*% ahat
    #muhat0 = binomial(link = "logit")$linkinv(etahat0)
    etahat = as.matrix(etahat)
    muhat = family$linkinv(etahat)
    diff = 1
    nrep = 0
    sm = 1e-5
    #sm = 1e-8
    while (diff > sm & nrep < 20){
      oldmu = muhat
      nrep = nrep + 1
      gr = gr.fun(y, etahat, weights = w, fml = family$family)
      wt = wt.fun(y, etahat, n, weights = w, fml = family$family)
      cv = crossprod(dd, (wt * etahat - gr))
      #qv = t(dd) %*% diag(wt) %*% dd
      dd2 = apply(dd, 2, function(ddi) ddi * sqrt(wt))
      qv0 = crossprod(dd2)
      
      if(!is.null(dv)){ #only for H1
        qv = qv0 + dv #ps is absorbed in dv
      } else {
        qv = qv0
      }
      
      #new:
      if(all(amat == 0)){
        ahat = solve(qv, cv)
        face = integer(0)
      } else {
        #ans = qprog(qv, cv, amat, bvec)
        #ahat = ans$thetahat
        #face = ans$face
        
        ans = solve.QP(qv, cv, t(amat), bvec)
        ahat = ans$solution
        face = ans$iact
      }
      
      # ans = solve.QP(qv, cv, t(amat), bvec)
      # ahat = ans$solution
      # face = ans$iact
      
      etahat = dd %*% ahat
      #muhat0 = binomial(link = "logit")$linkinv(etahat0)
      etahat = as.matrix(etahat)
      muhat = family$linkinv(etahat)
      diff = mean((muhat - oldmu)^2)
    }
    #check!
    #want negative log-like
    llh = llh.fun(y, muhat, etahat, phihat=NULL, n, weights=NULL, fml = family$family) * n / 2
  }
  #print (sse)
  dev = ifelse(wt.iter, llh, sse)
  #----------------------------------------------------------------------------
  #get obs edf and sig2hat
  #----------------------------------------------------------------------------
  #umat = chol(qv)
  #uinv = solve(umat)
  #new: 2025/9
  umat = chol(qv)
  uinv = backsolve(umat, diag(ncol(qv)))
  atil = amat %*% uinv
  qvinv = tcrossprod(uinv)
  #new term
  dtd = crossprod(dd)
  bigp_base = qvinv %*% dtd
  edf_base = sum(diag(bigp_base))
  
  #umat = chol(qv)
  #iumat = diag(ncol(umat))
  #uinv = backsolve(umat, iumat)
  #atil = amat %*% uinv
  #xw_uinv = dd %*% uinv
  
  nc_am = NCOL(amat)
  imat = diag(nc_am)
  if (all(face == 0) | length(face) == 0) {
    #if(length(face) == 0) {
    pmat = imat
    edf = edf_base
  } else {
    dp = -t(atil)
    smat = dp[, face, drop = FALSE]
    #new: avoid solve error
    qr_smat = qr(smat)
    #already orthonormal here!
    smat = qr.Q(qr_smat)[, 1:(qr_smat$rank), drop = FALSE]
    # for(ic in 1:NCOL(smat)){
    #   ic_norm = crossprod(smat[,ic])
    #   smat[,ic] = smat[,ic]/sqrt(ic_norm[1,1])
    # }
    #pmat_polar = smat %*% solve(crossprod(smat), t(smat))
    pmat_polar = tcrossprod(smat)
    pmat = (imat - pmat_polar)
    bigp_trace_polor = uinv %*% tcrossprod(pmat_polar, uinv) %*% dtd
    edf = edf_base - sum(diag(bigp_trace_polor))
  }
  if(is.null(covmat_inv)){
    #bigp = xw_uinv %*% tcrossprod(pmat, xw_uinv)
    #edf = sum(diag(bigp))
  
    if(!wt.iter){
      sig2hat = dev / (n - 1 * edf)
    }
  } else {
    #arp will have an extra covmat_inv
    #bigp = xw_uinv %*% tcrossprod(pmat, xw_uinv) %*% covmat_inv
    if (all(face == 0) | length(face) == 0){
      bigp_trace = qvinv %*% crossprod(dd, covmat_inv) %*% dd
    } else {
      bigp_trace = (qvinv - uinv %*% tcrossprod(pmat_polar, uinv)) %*% crossprod(dd, covmat_inv) %*% dd
    }
    edf = sum(diag(bigp_trace))
    #edf = sum(diag(bigp))
    #no use:
    #sig2hat = dev / (n - 1 * edf)
    #sig2hat = dev0/(n-edf)
    sig2hat_z = dev / (n - 1 * edf)
  }
  
  edfs = edf
  gcv = dev / (1 - edf / n)^2
  gcvs = gcv
  
  #edfu: check more...
  #bigp_un = tcrossprod(xw_uinv)
  #if(ps > 0){
  #dv absorbs ps
  if(!is.null(dv)){
    #bigp_un = dd %*% solve(qv0 + ps * dv, t(dd))
    if(is.null(covmat_inv)){
      #bigp_un = dd %*% solve(qv, t(dd))
      bigp_un = qvinv %*% dtd #not the projection matrix, just to make it faster to get the trace
    } else {
      #for arp qv0 is t(dd) * Rinv * dd
      #bigp_un = dd %*% solve(qv, t(dd)) %*% covmat_inv
      bigp_un = qvinv %*% crossprod(dd, covmat_inv) %*% dd
    }
  } else {
    if(is.null(covmat_inv)){
      #bigp_un = dd %*% solve(qv0, t(dd))
      bigp_un = dd %*% MASS::ginv(qv0) %*% t(dd)
    } else {
      #bigp_un = dd %*% solve(qv0, t(dd)) %*% covmat_inv
      bigp_un = dd %*% MASS::ginv(qv0) %*% t(dd) %*% covmat_inv
    }
  }
  
  edfu_use = sum(diag(bigp_un))
  if(!wt.iter){
    sighat = sig2hat^.5
    # if(!is.null(covmat_inv)){
    #   aic = n*log(sig2hat*n)+2*(p+edf)
    # }else{
    #   aic = n*log(sig2hat*(n - 1 * edf))+2*(p+edf)
    # }
    #aic = n*log(sse) - n*log(n) + 2*(p + edf)
    #if(!is.null(covmat_inv)){
    aic = n*log(sig2hat)+2*(p+edf)
    #}
  }
  
  #mat1 is used for se of parametric covariates, same name as what is used in wps.fit
  #wps.fit doesn't have mat1 for binomial
  if(!arp || arp & p == 0){
    mat1 = uinv %*% pmat %*% t(uinv) %*% t(dd)
  }else if(arp & p > 0){
    #mat1 = uinv %*% pmat %*% t(uinv) %*% t(xtil)
    mat1 = uinv %*% pmat %*% t(uinv) %*% t(xtil*sqrt(gamma[1]))
  }
  rslt = list(dev=dev, etahat=etahat, ahat=ahat, face=face,
              qv=qv, edf=edf, edfs=edfs, gcv=gcv, gcvs=gcvs,
              sighat=sighat, edfu_use=edfu_use,
              covmat=covmat, covmat_inv=covmat_inv,
              phi=phi, mat1=mat1, pval.ts=pval.ts, sig2hat_z=sig2hat_z, aic=aic, gamma=gamma)
  #rslt = list(etahat = etahat, dev = dev, ahat = ahat, face = face, edf = edf, edfs = edfs,
  #            gcv = gcv, gcvs = gcv, edfu_use = edfu_use, sighat = sighat, covmat=covmat,
  #            covmat_inv = covmat_inv,
  #            qv = qv, phi = phi, gamma = gamma, sigma2 = sigma2,
  #            ddtil = xtil)
  return (rslt)
}

#------------------------------------------
#need a better way to find ps....
#------------------------------------------
.search_ps = function(pen, qv0, dv, dd, edfu, covmat_inv=NULL)
{
  #if(is.null(covmat_inv)){
  #qv = qv0 + pen * dv
  #val = sum(diag(dd %*% ginv(qv, t(dd)))) - edfu * (.85)
  #val = sum(diag(dd %*% solve(qv, t(dd)))) - edfu #edfu * (.8)
  #val = sum(diag(dd %*% solve(qv, t(dd)))) - edfu * (.85)
  
  qv = qv0 + pen * dv
  umat = chol(qv)
  qvinv = chol2inv(umat)
  val = sum(diag(qvinv %*% crossprod(dd))) - edfu * (.85)

  #cat('uniroot: ', val, '\n')
  #val = sum(diag(dd %*% solve(qv) %*% t(dd))) - edfu
  #}else{
  #  bigp_un = dd %*% solve(t(dd) %*% covmat_inv %*% dd + pen * dv, t(dd)) %*% covmat_inv
  #  val = sum(diag(bigp_un)) - edfu
  #}
  #val = .compute_edf(pen, qv0, dv, dd, y, amat, bvec, imat)$edfi - edfu
  return (val)
}

#---------------------------------------------
#gradient of hessian for binomial, poisson,etc
#---------------------------------------------
gr.fun <- function(y, etahat = NULL, weights = NULL, fml = object$family){
  n <- length(y)
  if (is.null(weights)) {
    weights <- 1:n*0 + 1
  }
  w <- weights
  if (fml == "poisson") {
    #cut_id = etahat > 5
    #if(sum(cut_id) == 0){
    gr <- w * (exp(etahat) -  y)
    #} else {
    #  gr <- 1:n*0
    #  gr[!cut_id] <- w[!cut_id] * (exp(etahat[!cut_id]) -  y[!cut_id])
    #  gr[cut_id] <- w[cut_id] * (exp(5) -  y[cut_id])
    #}
    #gr <- w * (exp(etahat) -  y)
  }
  if (fml == "binomial") {
    if (all(etahat == 0)) {
      gr <- w * (1/2 - y)
    } else {
      gr <- 1:n*0
      for (i in 1:n) {
        if (etahat[i] > 100) {
          gr[i] <- w[i] * (1 - y[i])
        } else {gr[i] <- w[i] * (exp(etahat[i]) / (1 + exp(etahat[i])) - y[i])}
      }
    }
  }
  if (fml == "Gamma") {
    gr <- w * (1 - y * exp(-etahat))
  }
  gr
}

#---------------------------------------------
#2nd deriv of hessian for binomial, poisson,etc
#---------------------------------------------
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
      for (i in 1:n) {
        if (etahat[i] > 100) {
          wt[i] <- 0
        } else {
          wt[i] <- w[i] * exp(etahat[i]) / ((1 + exp(etahat[i]))^2)
        }
      }
    }
  }
  if (fml == "gaussian") {
    wt <- w # (1:n*0 + 1) / w
  }
  if (fml == "Gamma") {
    wt <-  w * y * exp(-etahat)
  }
  wt <- as.vector(wt)
  wt
}

llh.fun <- function(y, muhat = NULL, etahat = NULL, phihat = NULL, n = NULL, weights = NULL, fml = object$family){
  sm <- 1e-7
  #sm <- 1e-5
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

#different from cgam, used to simulate y when computing mixture-beta dist
ysim.fun <- function(n, mu0 = NULL, fml = object$family, shp0 = NULL, sd = NULL, phi = NULL) {
  if (fml == "binomial") {
    #ysim <- 1:n*0
    #ysim[runif(n) < .5] <- 1
    ysim <- rbinom(n, size = 1, prob = mu0)
  }
  if (fml == "poisson") {
    if (!is.null(mu0)) {
      ysim <- rpois(n, mu0)
    }
  }
  if (fml == "gaussian") {
    if(!is.null(phi)){
      ysim <- mu0 + arima.sim(n = n, list(ar = phi), sd = sd)
    }else{
      ysim <- mu0 + rnorm(n, sd = sd)
    }
  }
  if (fml == "Gamma") {
    ysim <- rgamma(n, shape=1)
  }
  ysim
}