######
#cgam#
######
cgam <- function(formula, nsim = 1e+2, family = gaussian(), cpar = 1.2, data = NULL, weights = NULL, sc_x = FALSE, sc_y = FALSE)
{
  cl <- match.call()
  if (is.character(family)) 
     family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
     family <- family()
  if (is.null(family$family)) 
     stop("'family' not recognized!")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  if (family$family == "binomial") {
	if (class(y) == "factor") {
		y = ifelse(y == levels(y)[1], 0, 1)
	}
  }
  shapes1 <- NULL; shapes2 <- NULL
  xmat <- NULL; xnms <- NULL
  tr <- NULL; pl <- NULL; umb <- NULL
  tree.delta <- NULL; umbrella.delta <- NULL 
  tid1 <- NULL; tid2 <- NULL; tpos2 <- 0
  uid1 <- NULL; uid2 <- NULL; upos2 <- 0
  nums <- NULL; ks <- list(); sps <- NULL; xid <- 1
  zmat <- NULL; zid <- NULL; zid0 <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_param <- NULL; is_fac <- NULL; vals <- NULL; st <- 1; ed <- 1
  ztb <- list(); iztb <- 1
  for (i in 2:ncol(mf)) {
    if (is.numeric(attributes(mf[,i])$shape)) {
       shapes1 <- c(shapes1, attributes(mf[,i])$shape)
       xmat <- cbind(xmat, mf[,i])
       xnms <- c(xnms, attributes(mf[,i])$nm)
       nums <- c(nums, attributes(mf[,i])$numknots)
       sps <- c(sps, attributes(mf[,i])$space)
       ks[[xid]] <- attributes(mf[,i])$knots
       xid <- xid + 1
    }
    if (is.character(attributes(mf[,i])$shape)) {
       shapes2 <- c(shapes2, attributes(mf[,i])$shape)
       if (attributes(mf[,i])$shape == "tree") {
		pl <- c(pl, attributes(mf[,i])$pl)
		treei <- tree.fun(mf[,i], attributes(mf[,i])$pl)
		tree.delta <- rbind(tree.delta, treei)
		tpos1 <- tpos2 + 1
		tpos2 <- tpos2 + nrow(treei)
		tid1 <- c(tid1, tpos1)
		tid2 <- c(tid2, tpos2)
		tr <- cbind(tr, mf[,i])
       }
       if (attributes(mf[,i])$shape == "umbrella") {
		umbi <- umbrella.fun(mf[,i])
		umbrella.delta <- rbind(umbrella.delta, umbi)
		upos1 <- upos2 + 1
		upos2 <- upos2 + nrow(umbi)
		uid1 <- c(uid1, upos1)
		uid2 <- c(uid2, upos2)
		umb <- cbind(umb, mf[,i])	
	}
    }
    if (is.null(attributes(mf[,i])$shape)) {
	if (!is.null(names(mf)[i])) {
	  znms <- c(znms, names(mf)[i])
	}
        if (!is.matrix(mf[,i])) {
          zid <- c(zid, i)
	  is_param <- c(is_param, TRUE)
          if (is.factor(mf[,i])) {
	    is_fac <- c(is_fac, TRUE)
	    ch_char <- suppressWarnings(is.na(as.numeric(levels(mf[, i]))))
            if (any(ch_char)) {
	      vals <- c(vals, unique(levels(mf[, i]))[-1])
            } else {
	      vals <- c(vals, as.numeric(levels(mf[, i]))[-1])
	    }
            nlvs <- length(attributes(mf[,i])$levels)
	    ed <- st + nlvs - 2 
	    zid1 <- c(zid1, st)
	    zid2 <- c(zid2, ed)
	    st <- st + nlvs - 1
	    zmat0 <- model.matrix(~ mf[, i])[, -1, drop = FALSE]
	    zmat <- cbind(zmat, zmat0)
	    ztb[[iztb]] <- mf[,i]
	    iztb <- iztb + 1
          } else {
	    is_fac <- c(is_fac, FALSE)
            zmat <- cbind(zmat, mf[, i])
	    ztb[[iztb]] <- mf[,i]
            iztb <- iztb + 1
	    ed <- st
            zid1 <- c(zid1, st)
	    zid2 <- c(zid2, ed) 
	    st <- st + 1
	    vals <- c(vals, "")
         }
       } else {
	  is_param <- c(is_param, FALSE)
          is_fac <- c(is_fac, FALSE)
	  zmat0 <- mf[, i]
	  mat_cols <- ncol(zmat0)
	  mat_rm <- NULL
	  #rm_num <- 0
	  for (irm in 1:mat_cols) {
       	  	if (all(round(diff(zmat0[, irm]), 8) == 0)) {
                	mat_rm <- c(mat_rm, irm)
          	}
   	  }
	  if (!is.null(mat_rm)) {
	  	zmat0 <- zmat0[, -mat_rm, drop = FALSE]
		#rm_num <- rm_num + length(mat_rm)
	  }
	  zmat <- cbind(zmat, zmat0)
          ztb[[iztb]] <- mf[,i]
          iztb <- iztb + 1
	  vals <- c(vals, 1)
	  zid <- c(zid, i)
	  nlvs <- ncol(zmat0) + 1
	  ed <- st + nlvs - 2
	  zid1 <- c(zid1, st)
	  zid2 <- c(zid2, ed)
	  st <- st + nlvs - 1
      }
    }
  }
  dimnames(zmat)[[2]] <- NULL
  if (family$family == "binomial" | family$family == "poisson") {
     wt.iter = TRUE
  } else {wt.iter = FALSE}
  if (is.null(shapes1) & is.null(shapes2)) {
    nsim <- 0
  }
  #attr(xmat, "shape") <- shapes1
  xmat0 <- xmat; shapes0 <- shapes1; nums0 <- nums; ks0 <- ks; sps0 <- sps; xnms0 <- xnms; idx_s <- NULL; idx <- NULL
  if (any(shapes1 == 17)) {
    kshapes <- length(shapes1)
    obs <- 1:kshapes
    idx_s <- obs[which(shapes1 == 17)]; idx <- obs[which(shapes1 != 17)]
  
    xmat0[ ,1:length(idx_s)] <- xmat[ ,idx_s]
    shapes0[1:length(idx_s)] <- shapes1[idx_s]
    nums0[1:length(idx_s)] <- nums[idx_s]
    sps0[1:length(idx_s)] <- sps[idx_s]
    ks0[1:length(idx_s)] <- ks[idx_s]
    xnms0[1:length(idx_s)] <- xnms[idx_s]

    if (length(idx) > 0) {
      xmat0[ ,(1 + length(idx_s)):kshapes] <- xmat[ ,idx]
      shapes0[(1 + length(idx_s)):kshapes] <- shapes1[idx]
      nums0[(1 + length(idx_s)):kshapes] <- nums[idx]
      sps0[(1 + length(idx_s)):kshapes] <- sps[idx]
      ks0[(1 + length(idx_s)):kshapes] <- ks[idx]
      xnms0[(1 +length(idx_s)):kshapes] <- xnms[idx]
    }
    #xmat <- xmat0; nums <- nums0; ks <- ks0; sps <- sps0; xnms <- xnms0
  }
  shapes <- c(shapes1, shapes2)
  ans <- cgam.fit(y = y, xmat = xmat0, zmat = zmat, shapes = shapes0, numknots = nums0, knots = ks0, space = sps0, nsim = nsim, family = family, cpar = cpar, wt.iter = wt.iter, umbrella.delta = umbrella.delta, tree.delta = tree.delta, weights = weights, sc_x = sc_x, sc_y = sc_y, idx_s = idx_s, idx = idx)
  if (!is.null(uid1) & !is.null(uid2)) {
    uid1 <- uid1 + ans$d0 + ans$capm
    uid2 <- uid2 + ans$d0 + ans$capm
  }
  if (!is.null(tid1) & !is.null(tid2)) {
    tid1 <- tid1 + ans$d0 + ans$capm + ans$capu
    tid2 <- tid2 + ans$d0 + ans$capm + ans$capu 
  }
#new:
  knots <- ans$knots
  numknots <- ans$numknots
#new:
  if (length(knots) > 0) {
    names(knots) <- xnms
    #for (i in 1:length(knots)) {
    #  names(knots)[i] <- xnms[i]
    #  if (knots[[i]] == 0L) {
    #    knots[[i]] <- NULL
    #  }
    #}
  }
  rslt <- list(etahat = ans$etahat, muhat = ans$muhat, vcoefs = ans$vcoefs, xcoefs = ans$xcoefs, zcoefs = ans$zcoefs, ucoefs = ans$ucoefs, tcoefs = ans$tcoefs, coefs = ans$coefs, cic = ans$cic, d0 = ans$d0, edf0 = ans$edf0, etacomps = ans$etacomps, xmat = xmat, zmat = zmat, ztb = ztb, tr = tr, umb = umb, tree.delta = tree.delta, umbrella.delta = umbrella.delta, bigmat = ans$bigmat, shapes = shapes, shapesx = shapes1, wt = ans$wt, wt.iter = ans$wt.iter, family = ans$family, SSE0 = ans$sse0, SSE1 = ans$sse1, pvals.beta = ans$pvals.beta, se.beta = ans$se.beta, null_df = ans$df.null, df = ans$df, resid_df_obs = ans$resid_df_obs, null_deviance = ans$dev.null, deviance = ans$dev, tms = mt, capm = ans$capm, capms = ans$capms, capk = ans$capk, capt = ans$capt, capu = ans$capu, xid1 = ans$xid1, xid2 = ans$xid2, tid1 = tid1, tid2 = tid2, uid1 = uid1, uid2 = uid2, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2, nsim = nsim, xnms = xnms,  ynm = ynm, znms = znms, is_param = is_param, is_fac = is_fac, knots = knots, numknots = numknots, sps = sps, ms = ans$ms, cpar = ans$cpar, pl = pl, idx_s = idx_s, idx = idx)
  rslt$call <- cl
  class(rslt) <- "cgam"
  return (rslt) 
}

###############
#amat function#
###############
amat.fun <- function(x)
{
	obs <- 1:length(x)
	amat <- NULL
	xu <- unique(x)
	nx <- sort(xu)
	hd <- head(nx, 1)
	tl <- nx
	while (length(tl) > 1) {
		hd <- head(tl, 1)
		tl <- tl[-1]
		paired <- 0
		for (i in 1:length(tl)) {
			a1 <- 1:length(x)*0
			if (hd * tl[i] > 0) {
				if (hd < 0) {
		        		if (tl[i] > hd) {
						a1[min(obs[which(x == hd)])] <- -1; a1[min(obs[which(x == tl[i])])] <- 1 
					} else {
						a1[min(obs[which(x == hd)])] <- 1; a1[min(obs[which(x == tl[i])])] <- -1
					}
				}
				if (hd > 0) {
					if (tl[i] < hd) {
						a1[min(obs[which(x == hd)])] <- -1; a1[min(obs[which(x == tl[i])])] <- 1 
					} else {
						a1[min(obs[which(x == hd)])] <- 1; a1[min(obs[which(x == tl[i])])] <- -1
					}
				}
			#if (!all(a1 == 0) & paired == 0 ) {amat <- rbind(amat, a1); paired <- 1}
			}
			if (hd * tl[i] == 0) {
				if (hd == 0) {
					#if (tl[i] > 0) {
						a1[min(obs[which(x == hd)])] <- 1; a1[min(obs[which(x == tl[i])])] <- -1
					#}
				} else {
					a1[min(obs[which(x == hd)])] <- -1; a1[min(obs[which(x == tl[i])])] <- 1
				}
			}
			if (!all(a1 == 0) & paired == 0 ) {amat <- rbind(amat, a1); paired <- 1}
		}
	}
	dimnames(amat) <- NULL
	amat
}	

###############
#bmat function#
###############
bmat.fun <- function(x)
{
	obs <- 1:length(x)
	bmat <- NULL
	hd <- head(x, 1)
	tl <- x
	j <- 0
	while (length(tl) > 1) {
		hd <- head(tl, 1)
		tl <- tl[-1]
		paired <- 0
		j <- j + 1
		for (i in 1:length(tl)) {
			b1 <- 1:length(x)*0
			if (hd == tl[i]) {b1[j] <- -1; b1[j + i] <- 1}
			if (!all(b1 == 0) & paired == 0 ) {bmat <- rbind(bmat, b1); paired <- 1}
		}
	}
	dimnames(bmat) <- NULL
	bmat
}

##########
#cgam.fit#
##########
cgam.fit <- function(y, xmat, zmat, shapes, numknots, knots, space, nsim, family = gaussian(), cpar = 1.2, wt.iter = FALSE, umbrella.delta = NULL, tree.delta = NULL, weights = NULL, sc_x = FALSE, sc_y = FALSE, idx_s = NULL, idx = NULL) {
#print (weights)
        linkfun <- family$linkfun
	cicfamily <- CicFamily(family)
	llh.fun <- cicfamily$llh.fun
	etahat.fun <- cicfamily$etahat.fun
	gr.fun <- cicfamily$gr.fun
	wt.fun <- cicfamily$wt.fun
	zvec.fun <- cicfamily$zvec.fun
	muhat.fun <- cicfamily$muhat.fun
	ysim.fun <- cicfamily$ysim.fun
	deriv.fun <- cicfamily$deriv.fun
	dev.fun <- cicfamily$dev.fun 
	n <- length(y)
	sm <- 1e-7 
	#sm <- 1e-5
	capl <- length(xmat) / n
	if (capl < 1) {capl <- 0}
	if (round(capl, 8) != round(capl, 1)) {stop ("Incompatible dimensions for xmat!")}
#new:
if (capl > 0 & sc_x) {
	for (i in 1:capl) {xmat[,i] <- (xmat[,i] - min(xmat[,i])) / (max(xmat[,i]) - min(xmat[,i]))}
	#for (i in 1:capl) {xmat[,i] <- (xmat[,i] - mean(xmat[,i])) / sd(xmat[,i])}
	#for (i in 1:capl) {xmat[,i] <- xmat[,i] / sd(xmat[,i])}
}
#new:
#if (sd(y) > 1e+3) {
#	sc_y <- TRUE
#}
if (sc_y) {
	sc <- sd(y)
	y <- y / sc	
}
	capk <- length(zmat) / n
	if (capk < 1) {capk <- 0}
	if (round(capk, 8) != round(capk, 1)) {stop ("Incompatible dimensions for zmat!")}
#new:
	capls <- sum(shapes == 17)
####################################################
#get basis functions for the constrained components#
####################################################	
	delta <- NULL
	varlist <- NULL
	xid1 <- NULL; xid2 <- NULL; xpos2 <- 0
	knotsuse <- list(); numknotsuse <- NULL
	mslst <- list()
#new:
	capm <- 0
	capms <- 0
	if (capl - capls > 0) {
		del1_ans <- makedelta(xmat[, 1], shapes[1], numknots[1], knots[[1]], space = space[1])
		del1 <- del1_ans$amat
		knotsuse[[1]] <- del1_ans$knots 
		mslst[[1]] <- del1_ans$ms
		numknotsuse <- c(numknotsuse, length(del1_ans$knots))
        	m1 <- length(del1) / n
#new code: record the number of columns of del1 if shapes0[1] == 17:
		if (shapes[1] == 17) {capms <- capms + m1}
        	var1 <- 1:m1*0 + 1
		xpos1 <- xpos2 + 1
		xpos2 <- xpos2 + m1
		xid1 <- c(xid1, xpos1)
		xid2 <- c(xid2, xpos2)
		if (capl == 1) {
        		delta <- del1
         		varlist <- var1
          	} else {
	      		for (i in 2:capl) {
#new code:
	        		del2_ans <- makedelta(xmat[,i], shapes[i], numknots[i], knots[[i]], space = space[i])
				del2 <- del2_ans$amat
				knotsuse[[i]] <- del2_ans$knots
				mslst[[i]] <- del2_ans$ms
				numknotsuse <- c(numknotsuse, length(del2_ans$knots))
				m2 <- length(del2) / n
#new code: record the number of columns of del2 if shapes0[i] == 17:
				if (shapes[i] == 17) {capms <- capms + m2}
				xpos1 <- xpos2 + 1
				xpos2 <- xpos2 + m2
				xid1 <- c(xid1, xpos1)
				xid2 <- c(xid2, xpos2)
				delta <- rbind(del1, del2)
				varlist <- 1:(m1 + m2)*0
				varlist[1:m1] <- var1
				varlist[(m1 + 1):(m1 + m2)] <- (1:m2)*0 + i
				var1 <- varlist
				m1 <- m1 + m2
				del1 <- delta
	      		}
	    	}
		if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk > 0) {
			bigmat <- rbind(1:n*0 + 1, t(zmat), t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
			np <- 1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)  + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk == 0) {
			bigmat <- rbind(1:n*0 + 1, t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
			np <- 1 + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk > 0) {
			bigmat <- rbind(1:n*0 + 1, t(zmat), delta)
			np <- 1 + capk + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk == 0) {
			bigmat <- rbind(1:n*0 + 1, delta)
			np <- 1 + capms
		} else {
			print ("error in capk, shapes!")
		} 
#new:
	capm <- length(delta) / n - capms
	} else {
	  	if (capk + capls > 0) {
#new:
			if (capls  < 1 & capk > 0) {
          			bigmat <- rbind(1:n*0 + 1, t(zmat))
          			np <- 1 + capk
			} else if (capls > 0) {
				delta <- NULL; varlist <- NULL
				del1_ans <- makedelta(xmat[,1], 17, numknots[1], knots[[1]], space = space[1])
				del1 <- del1_ans$amat
				knotsuse[[1]] <- del1_ans$knots 
				mslst[[1]] <- del1_ans$ms
				numknotsuse <- c(numknotsuse, length(del1_ans$knots))
				m1 <- length(del1) / n
				var1 <- 1:m1*0 + 1
				xpos1 <- xpos2 + 1
				xpos2 <- xpos2 + m1
				xid1 <- c(xid1, xpos1)
				xid2 <- c(xid2, xpos2)
				if (capls == 1) {
        				delta <- del1
         				varlist <- var1
          			} else {
					for (i in 2:capls) {
	        				del2_ans <- makedelta(xmat[,i], 17, numknots[i], knots[[i]], space = space[i])
						del2 <- del2_ans$amat
						knotsuse[[i]] <- del2_ans$knots
						mslst[[i]] <- del2_ans$ms
						numknotsuse <- c(numknotsuse, length(del2_ans$knots))
						m2 <- length(del2) / n
						xpos1 <- xpos2 + 1
						xpos2 <- xpos2 + m2
						xid1 <- c(xid1, xpos1)
						xid2 <- c(xid2, xpos2)
						delta <- rbind(del1, del2)
						varlist <- 1:(m1 + m2)*0
						varlist[1:m1] <- var1
						varlist[(m1 + 1):(m1 + m2)] <- (1:m2)*0 + i
						var1 <- varlist
						m1 <- m1 + m2
						del1 <- delta
	      				}
				}
				if (capk < 1){
					bigmat <- rbind(1:n*0 + 1, delta)
					capms <- length(delta) / n
					np <- 1 + capms
				} else {
					bigmat <- rbind(1:n*0 + 1, t(zmat), delta)
					capms <- length(delta) / n
					np <- 1 + capk + capms 
				}			
			}
            	} else {bigmat <- matrix(1:n*0 + 1, nrow = 1); capm <- 0; capms <- 0; np <- 1}
	}
	if (!is.null(umbrella.delta)) {
		bigmat <- rbind(bigmat, umbrella.delta)
		capu <- length(umbrella.delta) / n
	} else {capu <- 0}
	if (!is.null(tree.delta)) {
		bigmat <- rbind(bigmat, tree.delta)
		capt <- length(tree.delta) / n
	} else {capt <- 0}
	if (!is.null(umbrella.delta) | !is.null(tree.delta)) 
		delta_ut <- rbind(umbrella.delta, tree.delta)
	if (capl + capk + capu + capt > 0) {
#		if (capl + capu + capt > 0) {
#new:
		if (capl - capls + capu + capt > 0) {
#new:initialize cvec
cvec <- NULL
			if (wt.iter) {
				etahat <- etahat.fun(n, y, fml = family$family)
				gr <- gr.fun(y, etahat, weights, fml = family$family)  
				wt <- wt.fun(etahat, n, weights, fml = family$family)     
				cvec <- wt * etahat - gr
			} else {wt <- wt.fun(etahat, n, weights, fml = family$family)}
			  	zvec <- zvec.fun(cvec, wt, y, fml = family$family)
#        		  	gmat <- t(bigmat %*% sqrt(diag(wt)))
#new: avoid memory allocation error
				gmat <- t(bigmat)
				for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
			  	dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
               		  	zsend <- gmat[, 1:np, drop = FALSE] 
				ans <- coneB(zvec, t(dsend), zsend)
			  	etahat <- t(bigmat) %*% ans$coefs
			  	if (wt.iter) {
					muhat <- muhat.fun(etahat, fml = family$family)
			  		diff <- 1
					if (family$family == "binomial") {
						mdiff <- abs(max(muhat) - 1) > sm	
					} else {mdiff <- TRUE}

					nrep <- 0
##########
#iterate!#	
##########
			  		while (diff > sm & mdiff & nrep < n^2){
						oldmu <- muhat	
						nrep <- nrep + 1
						gr <- gr.fun(y, etahat, weights, fml = family$family)	
						wt <- wt.fun(etahat, n, weights, fml = family$family) 
						cvec <- wt * etahat - gr
						#zvec <- cvec / sqrt(wt)
						zvec <- zvec.fun(cvec, wt, y, fml = family$family)						
#						gmat <- t(bigmat %*% sqrt(diag(wt)))
						gmat <- t(bigmat)
						for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
						dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
        	  				zsend <- gmat[, 1:np, drop = FALSE]
						ans <- coneB(zvec, t(dsend), zsend)
						etahat <- t(bigmat) %*% ans$coefs
						muhat <- muhat.fun(etahat, fml = family$family)
						diff <- mean((muhat - oldmu)^2)	
						mdiff <- abs(max(muhat) - 1)
						if (family$family == "binomial") {
							mdiff <- abs(max(muhat) - 1) > sm	
						} else {mdiff <- TRUE}
					}
			 	}
				yhat <- ans$yhat
				coefskeep <- ans$coefs 
########################
#if capk >= 0, we have:#
########################
				zcoefs <- coefskeep[1:(capk + 1)]
######################
#we will always have:#
######################
				vcoefs <- coefskeep[1:np]
#######################
#if capm > 0, we have:#
#######################
				xcoefs <- NULL
				#if (capm > 0) {
				#	xcoefs <- coefskeep[(np + 1):(np + capm)]
				#}
#new:
				if (capl > 0) {
					xcoefs <- coefskeep[(np - capms + 1):(np + capm)]
				}
#######################
#if capu > 0, we have:#
#######################

				ucoefs <- NULL
				if (capu > 0) {
					ucoefs <- coefskeep[(np + 1 + capm):(np + capm + capu)]
				}
#######################
#if capt > 0, we have:#
#######################
				tcoefs <- NULL
				if (capt > 0) {
					tcoefs <- coefskeep[(np + 1 + capm + capu):(np + capm + capu + capt)]
				}

#########################################################		  
#if we have at least one constrained predictor, we have:#
#########################################################

				thvecs <- NULL
				if (capl > 0) {	
#new code:
					dcoefs <- coefskeep[(np - capms + 1):(np + capm)]
					#dcoefs <- coefskeep[(np + 1):(np + capm)]	
	
#####################################################
#thvecs is f(x), where x has one of the eight shapes#
#####################################################

					thvecs <- matrix(nrow = capl, ncol = n)
	    				ncon <- 1
	    				for (i in 1:capl) {
	    	  				thvecs[i,] <- t(delta[varlist == i,]) %*% dcoefs[varlist == i]
#new:
	    	  				#if (shapes[i] > 2 & shapes[i] < 5) { 
						if (shapes[i] > 2 & shapes[i] < 5 | shapes[i] > 10 & shapes[i] < 13) { 
            		    				ncon <- ncon + 1
            		    				#thvecs[i,] <- thvecs[i,] + zcoefs[ncon] * xmat[,i]
							thvecs[i,] <- thvecs[i,] + vcoefs[capk + ncon] * xmat[,i]
            	  				}
	    				}
				}
#new:order thvecs back
if (!is.null(idx_s)) { 
	thvecs0 <- thvecs
	thvecs0[idx_s,] <- thvecs[1:length(idx_s), ]
	if (!is.null(idx)) {
		thvecs0[idx,] <- thvecs[(1+length(idx_s)):capl, ]
	}
	thvecs <- thvecs0
}
				thvecs_ut <- NULL
				if (capu + capt > 0) {
					thvecs_ut <- t(delta_ut) %*% coefskeep[(np + 1 + capm):(np + capm + capu + capt)]
				}
				if (!is.null(thvecs_ut)) {
					thvecs <- rbind(thvecs, t(thvecs_ut))
				}
#new: problem when not gaussian
if (sc_y) {
	y <- y*sc
	etahat <- etahat*sc
	for (i in 1:nrow(thvecs)) {
		thvecs[i,] <- thvecs[i,] * sc
	}
}
	  			etakeep <- etahat
				muhatkeep <- muhat.fun(etakeep, fml = family$family)
				wtkeep <- wt
				#llh <- llh.fun(y, muhat, etahat, n, weights, fml = family$family)
				llh <- llh.fun(y, muhatkeep, etakeep, n, weights, fml = family$family)
#print (sum(abs(coefskeep) > 0))
				df_obs <- sum(abs(coefskeep) > 0)
				if (family$family == "poisson") {		
					mu0 <- mean(y)
					eta0 <- log(mu0)
				} else {mu0 <- NULL}
		} else if (capk + capls > 0 & capl - capls + capt + capu == 0) {
			if (is.null(weights)) {
				weights <- 1:n*0 + 1
			}
			prior.w <- weights
			vmat <- t(bigmat[1:np, , drop = FALSE])
#new: avoid solve problem	
#tvmat <- check_irred(t(vmat), msg = FALSE)$edge
#vmat <- t(tvmat)		
			if (wt.iter) {
				nrep <- 0
				muhat <- mean(y) + 1:n*0
				etahat <- linkfun(muhat)
				diff <- 1
				if (family$family == "binomial") {
					mdiff <- abs(max(muhat) - 1) > sm	
				} else {mdiff <- TRUE}
				while (diff > sm & mdiff & nrep < n^2) {
					nrep <- nrep + 1
					oldmu <- muhat
					zhat <- etahat + (y - muhat) * deriv.fun(muhat, fml = family$family)				
					#w <- diag(as.vector(prior.w / deriv.fun(muhat)))		
					#w <- diag(as.vector(prior.w * (deriv.fun(muhat, fml = family$family))^(-1)))
					#b <- solve(t(vmat) %*% w %*% vmat) %*% t(vmat) %*% w %*% zhat
					w <- as.vector(prior.w * (deriv.fun(muhat, fml = family$family))^(-1))
 					tvmat <- t(vmat)
 					for (i in 1:n) {tvmat[,i] <- tvmat[,i] * w[i]}
					b <- solve(tvmat %*% vmat) %*% tvmat %*% zhat
					etahat <- vmat %*% b
					muhat <- muhat.fun(etahat, fml = family$family)		
					diff <- mean((muhat - oldmu)^2)	
					mdiff <- abs(max(muhat) - 1)
					if (family$family == "binomial") {
						mdiff <- abs(max(muhat) - 1) > sm	
					} else {mdiff <- TRUE}
				}
				zcoefs <- b[1:(capk + 1)]
				#se.beta <-  sqrt(diag(solve(t(vmat) %*% w %*% vmat)))[1:(capk + 1)]
				se.beta <- sqrt(diag(solve(tvmat %*% vmat)))[1:(capk + 1)]
				zstat <- zcoefs / se.beta
				pvals.beta <-  1 - pchisq(zstat^2, df = 1)
			} else {
				#w <- diag(prior.w)
				#b <- solve(t(vmat) %*% w %*% vmat) %*% t(vmat) %*% w %*% y
				w <- prior.w
				tvmat <- t(vmat)
				for (i in 1:n) {tvmat[,i] <- tvmat[,i] * w[i]}
				b <- solve(tvmat %*% vmat) %*% tvmat %*% y
				etahat <- vmat %*% b
				muhat <- muhat.fun(etahat, fml = family$family)	
				sdhat2 <- sum(prior.w * (y - muhat)^2) / (n - np)
				zcoefs <- b[1:(capk + 1)]
				se.beta <-  sqrt(diag(solve(tvmat %*% vmat) * sdhat2))[1:(capk + 1)]
				tstat <- zcoefs / se.beta
				pvals.beta <-  (1 - pt(abs(tstat), df = n - np)) * 2 
			}
#add thvecs if capls > 0:
			thvecs <- NULL
			if (capls > 0) {
				thvecs <- matrix(nrow = capls, ncol = n)
	    			dcoefs <- b[(capk + 2):np]
	    			for (i in 1:capls) {
					thvecs[i,] <- t(delta[varlist == i,]) %*% dcoefs[varlist == i]
	    			}
			}
#new:
if (sc_y) {
	y <- y*sc
	etahat <- etahat*sc
	muhat <- muhat.fun(etahat, fml = family$family)
	for (i in 1:nrow(thvecs)) {
		thvecs[i,] = thvecs[i,] * sc
	}
}
			llh <-  llh.fun(y, muhat, etahat, n, weights, fml = family$family)
			df_obs <- np
			dfmean <- np
			rslt <- new.env()
			rslt$family <- family 
			rslt$wt.iter <- wt.iter 
			#rslt$wt <- diag(w)
			rslt$wt <- w
			rslt$bigmat <- bigmat
			rslt$etahat <- etahat
			rslt$muhat <- muhat
			rslt$d0 <- np
			rslt$capm <- 0
			rslt$capms <- capms
			rslt$capk <- capk
			rslt$capu <- capu
			rslt$capt <- capt
			rslt$xid1 <- xid1 + np - capms
			rslt$xid2 <- xid2 + np - capms 
			rslt$dfmean <- dfmean	
			rslt$edf0 <- dfmean 
			#rslt$llh <- llh 
			#if (nsim > 0) {
			rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - 1.5 * (dfmean - np)))
			#}
			rslt$zcoefs <- zcoefs
			rslt$coefs <- b
			rslt$vcoefs <- b
			rslt$xcoefs <- b[(capk + 2):np]
			rslt$se.beta <- se.beta 
			rslt$pvals.beta <- pvals.beta 
			rslt$dev <- dev.fun(y, muhat, etahat, weights, fml = family$family)$dev
			rslt$dev.null <- dev.fun(y, muhat, etahat, weights, fml = family$family)$dev.null
			rslt$df <- n - np 
			rslt$df.null <- n - 1
			rslt$resid_df_obs <- n - np - 1.5 * (df_obs - np)
			rslt$vhat <- etahat 		
			rslt$vmat <- vmat	
			rslt$etacomps <- thvecs
			rslt$knots <- knotsuse
			rslt$numknots <- numknotsuse 
			rslt$ms <- mslst
			return (rslt) 
			
		}

##########
#get cic#
##########
		  if (capl - capls + capu + capt > 0 & nsim > 0) {
	  		dfs <- 1:nsim
	  		for (isim in 1:nsim) {
				#set.seed(123)
	    			ysim <- ysim.fun(n, mu0, fml = family$family)
		  		if (wt.iter) {
					etahat <- etahat.fun(n, ysim, fml = family$family)
					gr <- gr.fun(ysim, etahat, weights, fml = family$family)
					wt <- wt.fun(etahat, n, weights, fml = family$family)
					cvec <- wt * etahat - gr
				} else {wt <- wt.fun(etahat, n, weights, fml = family$family)}
					zvec <- zvec.fun(cvec, wt, ysim, fml = family$family)
#            				gmat <- t(bigmat %*% sqrt(diag(wt)))
					gmat <- t(bigmat)
					for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
           				dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
            				zsend <- gmat[, 1:np, drop = FALSE]
					ans <- try(coneB(zvec, t(dsend), zsend))
					if (class(ans) == "try-error") next 
				if (wt.iter) {
						etahat <- t(bigmat) %*% ans$coefs
						muhat <- muhat.fun(etahat, fml = family$family)
						diff <- 1
						if (family$family == "binomial") {
							mdiff <- abs(max(muhat) - 1) > sm	
						} else {mdiff <- TRUE}
##########
#iterate!#
##########
						nrep <- 0 
						while (diff > sm & nrep < n^2 & mdiff > sm) {
							nrep <- nrep + 1
							oldmu <- muhat	
							gr <- gr.fun(ysim, etahat, weights, fml = family$family)
							wt <- wt.fun(etahat, n, weights, fml = family$family)
							cvec <- wt * etahat - gr
							#zvec <- cvec / sqrt(wt)
							zvec <- zvec.fun(cvec, wt, y, fml = family$family)
#							gmat <- t(bigmat %*% sqrt(diag(wt)))
							gmat <- t(bigmat)
							for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
							dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
							zsend <- gmat[, 1:np, drop = FALSE]
							ans <- try(coneB(zvec, t(dsend), zsend))
							if (class(ans) == "try-error") next 
							etahat <- t(bigmat) %*% ans$coefs
							muhat <- muhat.fun(etahat, fml = family$family)
							diff <- mean((muhat - oldmu)^2)	
							if (family$family == "binomial") {
								mdiff <- abs(max(muhat) - 1) > sm	
							} else {mdiff <- TRUE}
						}
				}
	    			dfs[isim] <- sum(abs(ans$coefs) > 0)
	  		}
	  		dfmean <- mean(dfs)
#print (dfmean)
		   } else if (capl - capls + capu + capt > 0 & nsim == 0) {
			dfmean <- NULL
		   } 
###################################################
#if the user does not give any predictor, we have:#
###################################################
	} else {
		rslt <- new.env()
		rslt$family <- family 
		rslt$wt.iter <- wt.iter
		rslt$muhat <- 1:n*0 + mean(y)
		rslt$etahat <- linkfun(muhat)
		rslt$dfmean <- 1
		print ("No predictor is provided")
		return (rslt)	
	}
	if (capl - capls + capu + capt > 0) {
#new:
#exclude the case we only have z or unrestricted smooth
		xid1 <- xid1 + np - capms
		xid2 <- xid2 + np - capms 
		#xid1 <- xid1 + np
		#xid2 <- xid2 + np
		rslt <- new.env()		
		rslt$family <- family 
		rslt$wt.iter <- wt.iter
		rslt$wt <- wtkeep
		rslt$bigmat <- bigmat  
		rslt$etahat <- etakeep
		rslt$muhat <- muhatkeep 
		rslt$d0 <- np
		rslt$capm <- capm
		rslt$capms <- capms 
		rslt$capk <- capk
		rslt$capu <- capu
		rslt$capt <- capt
		rslt$xid1 <- xid1
		rslt$xid2 <- xid2 
		rslt$coefs <- coefskeep
		rslt$vcoefs <- vcoefs
		rslt$xcoefs <- xcoefs
		rslt$zcoefs <- zcoefs
		rslt$ucoefs <- ucoefs
		rslt$tcoefs <- tcoefs 
		rslt$dfmean <- dfmean
		#rslt$llh <- llh
		#if (nsim > 0) {
		if (!is.null(dfmean)) {
#check!
#new: use dfmean - np if n - np - 1.5 * (dfmean - np) < 0
			if ((n - np - 1.5 * (dfmean - np)) <= 0) {
				rslt$cic <- llh + log(1 + 2 * dfmean / (dfmean - np))
			} else {
        			rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - 1.5 * (dfmean - np)))
			}
		} else {rslt$cic <- NULL}
		#if (nsim > 0) {
			#rslt$edf0 <- dfmean - np
		#}
		rslt$edf0 <- dfmean 
		rslt$etacomps <- thvecs
		vmat <- t(bigmat[1:np, , drop = FALSE])
		rslt$vmat <- vmat
		if (is.null(weights)) {
			weights <- 1:n*0 + 1
		}
		prior.w <- weights
		#w <- diag(as.vector(prior.w / deriv.fun(muhatkeep, fml = family$family)))
		w <- as.vector(prior.w / deriv.fun(muhatkeep, fml = family$family))
###############################################
#the case capk = 0 and capk >= 0 are combined:#
###############################################
#debugged: vhat -> muhat.fun(vhat)
		vhat <- vmat %*% vcoefs	
		muvhat <- muhat.fun(vhat, fml = family$family)
#debugged: yhat -> muhatkeep
		sse1 <- sum(prior.w * (y - muhatkeep)^2)
		sse0 <- sum(prior.w * (y - muvhat)^2)
#new: use (df_obs - np) if (n - np - 1.5 * (df_obs - np)) <= 0
		#if ((n - np - 1.5 * (df_obs - np)) <= 0) {
		#	sdhat2 <- sse1 / (df_obs - np)
		#} else {
		#	sdhat2 <- sse1 / (n - np - 1.5 * (df_obs - np))
		#}
#new: 
		if ((n - np - cpar * df_obs) <= 0) {
			#sdhat2 <- sse1 / (df_obs - np)
			sdhat2 <- sse1 / df_obs
		} else {
			#sdhat2 <- sse1 / (n - np - 1.5 * (df_obs - np))
			sdhat2 <- sse1 / (n - np - cpar * df_obs)
		}
#debugged: vmat -> vmat and duse
#new: coefskeep include zcoefs; bigmat include vmat
		pmat <- vmat
		capbm <- length(bigmat) / n
		bigmat_nv <- bigmat[(np + 1):capbm, , drop = FALSE]	
		coefs_nv <- coefskeep[(np + 1):capbm]	
		duse <- coefs_nv > 1e-8
#print (coefskeep[1:np])
#print (dim(pmat))
		if (sum(duse) >= 1) {
			pmat = cbind(vmat, t(bigmat_nv[duse, , drop = FALSE]))
		}
#check!
#tpmat <- check_irred(t(pmat), msg = FALSE)$edge
#pmat <- t(tpmat)
		if (wt.iter) {
#			se2 <- solve(t(pmat) %*% w %*% pmat)
#new:
			tpmat <- t(pmat)
			for (i in 1:n) {tpmat[,i] <- tpmat[,i] * w[i]}
			se2 <- solve(tpmat %*% pmat)
		} else {
#			se2 <- solve(t(pmat) %*% diag(prior.w) %*% pmat) * sdhat2
			tpmat <- t(pmat)
			for (i in 1:n) {tpmat[,i] <- tpmat[,i] * prior.w[i]}
			se2 <- solve(tpmat %*% pmat) * sdhat2
		}		 		
		se.beta <- 1:(capk + 1)*0
		tstat <- 1:(capk + 1)*0
		pvals.beta <- 1:(capk + 1)*0
		rslt$zcoefs <- zcoefs 	
		for (i in 1:(capk + 1)) {
			se.beta[i] <- sqrt(se2[i,i])
			tstat[i] <- zcoefs[i] / se.beta[i]
#new code: n - np - 1.5 * (df_obs - np) must be positive
			if ((n - np - cpar * df_obs) <= 0) {
				pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  df_obs)) 
				warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
			} else {
				pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  n - np - cpar * df_obs)) 
			}
		}
		rslt$se.beta <- se.beta
		rslt$pvals.beta <- pvals.beta
		rslt$sse1 <- sse1
		rslt$sse0 <- sse0
		rslt$dev <- dev.fun(y, muhatkeep, etakeep, weights, fml = family$family)$dev
		rslt$dev.null <- dev.fun(y, muhatkeep, etakeep, weights, fml = family$family)$dev.null
		rslt$df <- n - np 
		rslt$df.null <- n - 1
		rslt$resid_df_obs <- n - np - cpar * df_obs
		rslt$df_obs <- df_obs
		rslt$vhat <- vhat
		if (length(knotsuse) == 0) {
			knotsuse <- NULL
		}
if (!is.null(idx_s)) {
	knotsuse0 <- knotsuse
	numknotsuse0 <- numknotsuse
	mslst0 <- mslst
	knotsuse0[idx_s] <- knotsuse[1:length(idx_s)]
	numknotsuse0[idx_s] <- numknotsuse[1:length(idx_s)]
	mslst0[idx_s] <- mslst[1:length(idx_s)]
	if (!is.null(idx)) { 
		knotsuse0[idx] <- knotsuse[(1+length(idx_s)):capl]
		numknotsuse0[idx] <- numknotsuse[(1+length(idx_s)):capl]
		mslst0[idx] <- mslst[(1+length(idx_s)):capl]
	}
	knotsuse <- knotsuse0
	numknotsuse <- numknotsuse0
	mslst <- mslst0
}
		rslt$knots <- knotsuse
		rslt$numknots <- numknotsuse 
		rslt$ms <- mslst
		rslt$capms <- capms
		rslt$cpar <- cpar 
		return (rslt)
	}
}


###########
#CicFamily#
###########
CicFamily <- function(object,...)UseMethod("CicFamily")
CicFamily <- function(object) {
  llh.fun <- function(y, muhat = NULL, etahat = NULL, n = NULL, weights = NULL, fml = object$family){
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
    llh 
  }

  etahat.fun <- function(n, y, fml = object$family){
    if (fml == "poisson") {
      etahat <- 1:n*0 + log(mean(y)) 
    } 
    if (fml == "binomial") {
      etahat <- 1:n*0 
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
	   for (i in 1:n) {
	     if (etahat[i] > 100) {
		 gr[i] <- w[i] * (1 - y[i]) 
	     } else {gr[i] <- w[i] * (exp(etahat[i]) / (1 + exp(etahat[i])) - y[i])}
           }         
         }
    }
    gr
  }

  wt.fun <- function(etahat = NULL, n = NULL, weights = NULL, fml = object$family){
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
    wt <- as.vector(wt)
    wt 
  }

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
    zvec 
  }

  muhat.fun <- function(etahat, wt = NULL, fml = object$family){
    n <- length(etahat)
    if (fml == "poisson") {
      muhat <- exp(etahat)
    }
    if (fml == "binomial") {
      muhat <- 1:n*0
      #muhat[wt == 0] <- 1
      for (i in 1:n) {
        if (etahat[i] > 100) {
	  muhat[i] <- 1
 	} else {
          muhat[i] <- exp(etahat[i]) / (1 + exp(etahat[i]))
        }
      }
    }
    if (fml == "gaussian") {
      muhat <- etahat
    }
   muhat 
  }

  ysim.fun <- function(n, mu0 = NULL, fml = object$family) {
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
    ysim 
  }

  deriv.fun <- function(muhat, fml = object$family) {
    if (fml == "binomial") {
	deriv <- 1 / (muhat * (1 - muhat))
    }
    if (fml == "poisson") {
	deriv <- 1 / muhat
    }
    if (fml == "gaussian") {
	deriv <- 1
    }
   deriv
  }

 dev.fun <- function(y, muhat, etahat, weights, fml = object$family){
  n <- length(y)
  sm <- 1e-7
  #sm <- 1e-5
  if (is.null(weights)) {
	weights <- 1:n*0 + 1
  }
  w <- weights
  vmat <- matrix(1:n*0 + 1, ncol = 1)
  if (fml == "poisson") {
        #dev <- 2 * sum(w * (y * log(y / muhat) - y + muhat))
	dev <- 0
	for (i in 1:n) {
	  if (y[i] == 0) {
            dev <- dev + 2 * w[i] * muhat[i]
          } else {
            dev <- dev + 2 * w[i] * (y[i] * log(y[i] / muhat[i]) - y[i] + muhat[i])
          }
	}
  }
  if (fml == "binomial") {
        dev <- 0
        for (i in 1:n) {
          if (y[i] == 0 & w[i] != 0) {
            dev <- dev + 2 * w[i] * log(w[i] / (w[i] - w[i] * muhat[i]))
          } else if (y[i] == 1 & w[i] != 0) {
              dev <- dev + 2 * w[i] * log(w[i] / (w[i] * muhat[i]))
          } else if (0 < y[i] & y[i] < 1 & w[i] != 0) {
              dev <- dev + 2 * w[i] * y[i] * log(w[i] * y[i] / (w[i] * muhat[i])) + 2 * (w[i] - w[i] * y[i]) * log((w[i] - w[i] * y[i]) / (w[i] - w[i] * muhat[i]))
          } else {
             stop ("y values must be 0 <= y <= 1!")
          }
       }
  }
  if (fml == "gaussian") {
        dev <- sum(w * (y - muhat)^2)
  }
###################
#get null deviance#
###################
  if (fml == "binomial" | fml == "poisson") {
      diff <- 1
      muhat0 <- mean(y) + 1:n*0
      if (fml == "poisson") {
         etahat0 <- log(muhat0)
      } 
      if (fml == "binomial") {
         etahat0 <- log(muhat0 / (1 - muhat0))
      } 		
      while (diff > sm) {
        oldmu <- muhat0
	zhat <- etahat0 + (y - muhat0) * deriv.fun(muhat0, fml = fml)		
#	wmat <- diag(as.vector(w / deriv.fun(muhat0, fml = fml)))			
#	b <- solve(t(vmat) %*% wmat %*% vmat) %*% t(vmat) %*% wmat %*% zhat
	wm <- as.vector(w / deriv.fun(muhat0, fml = fml))
        tvmat <- t(vmat)
        for (i in 1:n) {tvmat[,i] <- tvmat[,i] * wm[i]}
        b <- solve(tvmat %*% vmat) %*% tvmat %*% zhat
	etahat0 <- vmat %*% b
	muhat0 <- muhat.fun(etahat0, fml = fml)		
	diff <- mean((muhat0 - oldmu)^2)	
      }
      if (fml == "poisson") {
        #dev.null <- 2 * sum(w * (y * log(y / muhat0) - y + muhat0))
	dev.null <- 0
        for (i in 1:n) {
	  if (y[i] == 0) {
            dev.null <- dev.null + 2 * w[i] * muhat0[i]
          } else {
            dev.null <- dev.null + 2 * w[i] * (y[i] * log(y[i] / muhat0[i]) - y[i] + muhat0[i])
          }
	}
      }
      if (fml == "binomial") {
        dev.null <- 0
        for (i in 1:n) {
          if (y[i] == 0 & w[i] != 0) {
            dev.null <- dev.null + 2 * w[i] * log(w[i] / (w[i] - w[i] * muhat0[i]))
          } else if (y[i] == 1 & w[i] != 0) {
              dev.null <- dev.null + 2 * w[i] * log(w[i] / (w[i] * muhat0[i]))
          } else if (0 < y[i] & y[i] < 1 & w[i] != 0) {
              dev.null <- dev.null + 2 * w[i] * y[i] * log(w[i] * y[i] / (w[i] * muhat0[i])) + 2 * (w[i] - w[i] * y[i]) * log((w[i] - w[i] * y[i]) / (w[i] - w[i] * muhat0[i]))
          } else {
              stop ("y values must be 0 <= y <= 1!")
	  }
        }
      } 
  }
  if (fml == "gaussian") {
#     wmat <- diag(w)
#     b <- solve(t(vmat) %*% wmat %*% vmat) %*% t(vmat) %*% wmat %*% y
     tvmat <- t(vmat)
     for (i in 1:n) {tvmat[,i] <- tvmat[,i] * w[i]}
     b <- solve(tvmat %*% vmat) %*% tvmat %*% y
     etahat0 <- vmat %*% b
     muhat0 <- muhat.fun(etahat0, fml = fml)	
     dev.null <- sum(w * (y - muhat0)^2)
  }
  rslt <- new.env()
  rslt$dev <- dev
  rslt$dev.null <- dev.null 
  rslt
  }

  ans <- list(llh.fun = llh.fun, etahat.fun = etahat.fun, gr.fun = gr.fun, wt.fun = wt.fun, zvec.fun = zvec.fun, muhat.fun = muhat.fun, ysim.fun = ysim.fun, deriv.fun = deriv.fun, dev.fun = dev.fun)
  class(ans) <- "CicFamily"
  return (ans)
}

#######################
#eight shape functions#
####################### 
incr <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 1
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

decr <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 2
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x 
} 

conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 3
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 4
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

incr.conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 5
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

decr.conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 6
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

incr.conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 7
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

decr.conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 8
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

s.incr <- function(x, numknots = 0, knots = 0, space = "E")
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 9
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    x
}

s.decr <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 10
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 11
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 12
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.incr.conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 13
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.incr.conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 14
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.decr.conv <- function(x, numknots = 0, knots = 0, space = "E") 
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 15
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots	
    attr(x, "space") <- space
    x
}

s.decr.conc <- function(x, numknots = 0, knots = 0, space = "E") 
{
   cl <- match.call()
   pars <- match.call()[-1]
   attr(x, "nm") <- deparse(pars$x)
   attr(x, "shape") <- 16
   attr(x, "numknots") <- numknots
   attr(x, "knots") <- knots	
   attr(x, "space") <- space
   x
}

s <- function(x, numknots = 0, knots = 0, space = "E") 
{
   cl <- match.call()
   pars <- match.call()[-1]
   attr(x, "nm") <- deparse(pars$x)
   attr(x, "shape") <- 17
   attr(x, "numknots") <- numknots
   attr(x, "knots") <- knots	
   attr(x, "space") <- space
   x
}

######################################################
#tree function: give the tree shape to x and return x#
######################################################
tree <- function(x, pl = NULL)
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- "tree"
    if (is.null(pl)) {
	if (is.numeric(x)) {
		if (0%in%x) {
			pl <- 0
		} else {pl <- min(x)}
	} else {
		xu <- unique(x)
		pl <- xu[1]
	}
    } else {
	if (!(pl%in%x)) {
		stop ("placebo level is not a level of the tree variable!")
	}
    }
    attr(x, "pl") <- pl
    x

}

##################################################
#tree.fun: make delta to a tree ordering variable#
##################################################
#tree.fun <- function(x)
#{
#    if (min(x) != 0) {
#	stop ("A tree ordering variable must have its placebo equal to 0!")
#    }
#    if (!all(round(x, 0) == x)) {
#	stop ("All elements of a tree ordering variable must be integers!")
#    }
#    if (any(x < 0))
#	stop ("All elements of a tree ordering variable must be positive!")
#    nx <- x
#    obs <- 1:length(x)
#    delta <- matrix(0, nrow = length(attributes(factor(x))$levels) - 1, ncol = length(x))
#    pl <- min(nx)
#    for (i in 1:nrow(delta)) {
#       nx <- nx[which(nx != pl)]
#       pl <- min(nx)
#       index <- obs[x == pl]
#       delta[i, index] <- 1
#    }
#  attr(delta, "shape") <- "tree"
#  delta
#}

tree.fun <- function(x, pl = NULL)
{
  #if (pl %in% x) {
  if (is.null(pl)) {
	if (is.numeric(x)) {
		if (0%in%x) {
			pl <- 0
		} else {pl <- min(x)}
	} else {
		xu <- unique(x)
		pl <- xu[1]
	}
  } else {
	if (!(pl%in%x)) {
		stop ("placebo level is not a level of the tree variable!")
	}
  }
  pos <- x %in% pl
  xu <- unique(x)
  nx <- xu[xu != pl] 
  delta <- matrix(0, nrow = (length(xu) - 1), ncol = length(x))
  for (i in 1:nrow(delta)) {
  	xi = nx[i]
	delta[i, !pos & x %in% xi] <- 1
	pos <- pos | x %in% xi
  }
  #} else {stop ("placebo level is not a level of the tree variable!")}
  attr(delta, "shape") <- "tree"
  delta
}

##############################################################
#umbrella function: give the umbrella shape to x and return x#
##############################################################
umbrella <- function(x)
{
  cl <- match.call()
  pars <- match.call()[-1]
  attr(x, "nm") <- deparse(pars$x)
  attr(x, "shape") <- "umbrella"
  x
}

###########################################################
#umbrella.fun: make delta to an umbrella ordering variable#
###########################################################
umbrella.fun <- function(x)
{
	amat <- amat.fun(x)
	bmat <- bmat.fun(x)
	constr <- t(rbind(amat, bmat))
	vmat <- qr.Q(qr(constr), complete = TRUE)[, -(1:(qr(constr)$rank)), drop = FALSE]
	if (!is.null(bmat)) {
		wperp <- t(rbind(t(vmat), bmat))
		wmat <- qr.Q(qr(wperp), complete = TRUE)[, -(1:(qr(wperp)$rank)), drop = FALSE]	
		atil <- amat %*% wmat 
		delta <- t(wmat %*% t(atil) %*% solve(atil %*% t(atil)))
	} else {
		delta <- t(t(amat) %*% solve(amat %*% t(amat)))
	}
	attr(delta, "shape") <- "umbrella"
	delta
}

###################################################
#find delta for a specific predictor x and a shape#
################################################### 
makedelta = function(x, sh, numknots = 0, knots = 0, space = "E", suppre = FALSE, interp = FALSE) {
#if (!interp) {
#x = (x - min(x)) / (max(x) - min(x))
#}
	n = length(x)
# find unique x values
#round(x,8) will make 0 edge in amat!
	#xu = sort(unique(round(x, 8)))
#new: center and scale to avoid numerical instabillity
	xu = sort(unique(x))
	n1 = length(xu)
	sm = 1e-7
	ms = NULL
#  increasing or decreasing
	if (sh < 3) {
		amat = matrix(0, nrow = n1 - 1, ncol = n)
		for (i in 1: (n1 - 1)) {
			amat[i, x > xu[i]] = 1
		}
		if (sh == 2) {amat = -amat}
		for (i in 1:(n1 - 1)) {
#new: use ms in predict.cgam
			ms = c(ms, mean(amat[i, ]))
			amat[i, ] = amat[i, ] - mean(amat[i, ])
		}
	} else if (sh == 3 | sh == 4) {
#  convex or concave
		amat = matrix(0, nrow = n1 - 2 ,ncol = n)
		#for (i in 1: (n1 - 2)) {
		#	amat[i, x > xu[i]] = x[x > xu[i]] - xu[i]
		#}
		for (i in 1: (n1 - 2)) {
			amat[i, x > xu[i+1]] = x[x > xu[i+1]] - xu[i+1]
		}
		if (sh == 4) {amat = -amat}
		xm = cbind(1:n*0+1,x)
		xpx = solve(t(xm) %*% xm)
		pm = xm %*% xpx %*% t(xm)
#new: use ms in predict.cgam
		ms = amat %*% t(pm)
		#amat = amat - amat %*% t(pm)
		amat = amat - ms
	} else if (sh > 4 & sh < 9) {
		amat = matrix(0, nrow = n1 - 1, ncol = n)
		if (sh == 5) { ### increasing convex
			for (i in 1:(n1 - 1)) {
				amat[i, x > xu[i]] = (x[x > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			for (i in 1:(n1 - 1)) {
				ms = c(ms, mean(amat[i, ]))
				amat[i,] = amat[i,] - mean(amat[i,])
			}
		} else if (sh == 6) {  ## decreasing convex
			for (i in 1:(n1 - 1)) {
				amat[i, x < xu[i + 1]] = (x[x < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			for (i in 1:(n1 - 1)) {
				ms = c(ms, mean(amat[i, ]))			
				amat[i,] = amat[i,] - mean(amat[i, ])
			}
#print (ms)
		} else if (sh == 7) { ## increasing concave
			for (i in 1:(n1 - 1)) {
				amat[i, x < xu[i + 1]] = (x[x < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			for (i in 1:(n1 - 1)) {
				ms = c(ms, mean(amat[i, ]))
				amat[i,] = -amat[i,] + mean(amat[i,])
			}		
		} else if (sh == 8) {## decreasing concave
			for (i in 1:(n1 - 1)) {
				amat[i, x > xu[i]] = (x[x > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			for (i in 1:(n1 - 1)) {
				ms = c(ms, mean(amat[i, ]))
				amat[i,] = -amat[i,] + mean(amat[i,])
			}
		}
	} else if (sh > 8 & sh < 18) {
		#if (all(knots == 0) & numknots == 0) {
		if (length(knots) < 2 & numknots == 0) {
			if (sh == 9 | sh == 10) {#1 2
				k = trunc(n1^(1/5)) + 4
			} else {k = trunc(n1^(1/7) + 4)}
			if (space == "Q") {
				t = quantile(xu, probs = seq(0, 1, length = k), names = FALSE)
			}
			if (space == "E") {
				#t = 0:k / k * (max(x) - min(x)) + min(x)
				t = 0:(k-1) / (k-1) * (max(x) - min(x)) + min(x)
			} 
		#} else if (any(knots != 0) & numknots == 0) {
		} else if (length(knots) >= 2 & numknots == 0) {
			t = knots
		#} else if (all(knots == 0) & numknots != 0) {
		} else if (length(knots) < 2 & numknots != 0) {
			if (space == "Q") {
				t = quantile(xu, probs = seq(0, 1, length = numknots), names = FALSE)
			} 
			if (space == "E") {
				#k = numknots
#new: numknots should be the # of all knots
				k = numknots - 1
				#if (sh == 9 | sh == 10) {#1 2
				#	k = trunc(n1^(1/5)) + 4
				#} else {k = trunc(n1^(1/7) + 4)}
				t = 0:k / k * (max(x) - min(x)) + min(x)
			}
		#} else if (any(knots != 0) & numknots != 0) {
		} else if (length(knots) >= 2 & numknots != 0) {
			#t0 = quantile(xu, probs = seq(0, 1, length = numknots), names = FALSE)
			t = knots
			if (!suppre) {
				print("'knots' is used! 'numknots' is not used!")
			}
			#print ("'knots' is used!")
			#if (numknots != length(knots)) {
			#	if (!suppre) {
			#		print("length(knots) is not equal to 'numknots'! 'knots' is used!")
			#	}
			#} else if (any(t0 != knots)) {
			#	if (!suppre) {
			#		print("equal x-quantiles knots != 'knots'! 'knots' is used! ") 
			#	}
			#}
		}
		if (sh == 9) {#1			
			amat_ans = monincr(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 10) {#2
			amat_ans = mondecr(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 11) {#3
			amat_ans = convex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 12) {#4
			amat_ans = concave(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 13) {#5
			amat_ans = incconvex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 14) {#6
			amat_ans = incconcave(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
		} else if (sh == 15) {#7
			#amat_ans = -incconcave(x, t, interp)
			amat_ans = incconcave(x, t, interp)
			amat = -amat_ans$sigma
			if (!interp) {
				ms = -amat_ans$ms
			}
		} else if (sh == 16) {#8
			#amat_ans = -incconvex(x, t, interp)
			amat_ans = incconvex(x, t, interp)
			amat = -amat_ans$sigma
			if (!interp) {
				ms = -amat_ans$ms
			}
		} else if (sh == 17) {#unconstrained
			amat_ans = incconvex(x, t, interp)
			amat = amat_ans$sigma
			ms = amat_ans$ms
			#amat = -incconcave(x, t)
			#amat = rbind(x, t(bcspl(x, m = length(t), knots = t)$bmat)) 
			#amat = rbind(x, convex(x, t))
		}
	}
	#if (sh < 9) {
	#	rslt = list(amat = amat, knots = 0, ms = ms)
	#} else {
	#	rslt = list(amat = amat, knots = t, ms = ms)
	#}
	if (sh < 9) {t = 0}	
	rslt = list(amat = amat, knots = t, ms = ms)
	rslt
}


# Monotone increasing
monincr = function(xs, t, interp = FALSE) {
	n = length(xs)
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	knt = 1:m
	for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:(k-1)) {
		index = x >= t[1] & x <= t[j]
		sigma[j, index] = 0

		index = x > t[j] & x <= t[j+1]
		sigma[j, index] = (x[index] - t[j])^2 / (t[j+2] - t[j]) / (t[j+1] - t[j])

		index = x > t[j+1] & x <= t[j+2]
		sigma[j, index] = 1 - (x[index] - t[j+2])^2 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	    
		index = x > t[j+2] #& x <= t[m]
		sigma[j, index] = 1
	}
	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 0
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = (x[index] - t[k])^2 / (t[k+2] - t[k]) / (t[k+1] - t[k])
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = 1 - (x[index] - t[k+2])^2 / (t[k+2] - t[k+1]) / (t[k+2] - t[k])
	
	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = 1 - (t[2] - x[index])^2 / (t[2] - t[1])^2

	index = x > t[2] 
	sigma[k+1, index] = 1
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 0
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = (x[index] - t[k+1])^2 / (t[k+2] - t[k+1])^2
	
#new:
	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}	


########################################################
# Monotone decreasing
mondecr = function(xs, t, interp = FALSE) {
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k + 2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	#t = x[knt]
	for (j in 1:(k - 1)) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = 1

		index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = 1 - (x[index] - t[j])^2 / (t[j+2] - t[j]) / (t[j+1] - t[j])

	    	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = (x[index] - t[j+2])^2 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])

	    	index = x > t[j+2] 
	    	sigma[j, index] = 0
	}

	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 1
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = 1 - (x[index] - t[k])^2 / (t[k+2] - t[k]) / (t[k+1] - t[k])

	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = (x[index] - t[k+2])^2 / (t[k+2] - t[k+1]) / (t[k+2] - t[k])

	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = (t[2] - x[index])^2 / (t[2] - t[1])^2

	index = x > t[2] 
	sigma[k+1, index] = 0

	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 1
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = 1 - (x[index] - t[k+1])^2 / (t[k+2] - t[k+1])^2

	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}

########################################################
# Convex
convex = function(xs, t, interp = FALSE) {
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:(k-1)) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = 0
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = (x[index] - t[j])^3 / (t[j+2] - t[j]) / (t[j+1] - t[j]) / 3
	    
	   	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = x[index] - t[j+1] - (x[index] - t[j+2])^3 / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) / 3 + (t[j+1] - t[j])^2 / 3 /(t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])

	    	index = x > t[j+2]
	    	sigma[j, index] = (x[index] - t[j+1]) + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	}
	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 0
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = (x[index] - t[k])^3 / (t[k+2] - t[k]) / (t[k+1] - t[k]) / 3

	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = x[index] - t[k+1] - (x[index] - t[k+2])^3 / (t[k+2] - t[k]) / (t[k+2] - t[k+1]) / 3 + (t[k+1] - t[k])^2 / 3 / (t[k+2] -t[k]) - (t[k+2] - t[k+1])^2 / 3 / (t[k+2] - t[k])
	
	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = x[index] - t[1] + (t[2] - x[index])^3 / (t[2] - t[1])^2 / 3 - (t[2] - t[1]) / 3
	
	index = x > t[2] 
	sigma[k+1, index] = x[index] - t[1] - (t[2] - t[1]) / 3
		
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 0
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = (x[index] - t[k+1])^3 / (t[k+2] - t[k+1])^2 / 3
	
	ms = NULL
	if (!interp) {
		xm = cbind(1:n*0+1, x)
		pm = xm %*% solve(t(xm) %*% xm) %*% t(xm)
		ms = matrix(0, nrow = nrow(sigma), ncol = ncol(sigma))
		for (i in 1:m) {
			ms[i,] = pm %*% sigma[i,]
			ms[i,] = ms[i, rank(xs)]		
			sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}


########################################################
# Concave
concave = function(xs, t, interp = FALSE) {
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 2
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:m
	#for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	#t = x[knt]
	for (j in 1:k) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = x[index] - t[1]
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = t[j] - t[1] + ((t[j+1] - t[j])^3 - (t[j+1] - x[index])^3) / 3 / (t[j+1] - t[j]) / (t[j+2] - t[j]) + (x[index] - t[j]) * (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	    
	        index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + ((t[j+2] - t[j+1])^3 - (t[j+2] - x[index])^3) / 3 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])	
 	   
 	   	index = x > t[j+2]
 	   	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	}

	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = -(t[2] - x[index])^3 / 3 / (t[2] - t[1])^2
	
	index = x > t[2] 
	sigma[k+1, index] = 0
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = x[index] - t[1]
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = t[k+1] - t[1] + ((t[k+2] - t[k+1])^2 * (x[index] - t[k+1]) - (x[index] - t[k+1])^3 / 3) / (t[k+2] - t[k+1])^2
	
	ms = NULL
	if (!interp) {
		xm = cbind(1:n*0+1, x)
		pm = xm %*% solve(t(xm) %*% xm) %*% t(xm)
		ms = matrix(0, nrow = nrow(sigma), ncol = ncol(sigma))
		for (i in 1:m) {
			ms[i,] = pm %*% sigma[i,]
			ms[i,] = ms[i, rank(xs)]		
			sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - pm %*% sigma[i,]
			sigma[i,] = sigma[i, rank(xs)]
		}
	}	
	rslt = list(sigma = sigma, ms = ms)
	rslt
}

########################################################
# Increasing and Convex
incconvex = function(xs, t, interp = FALSE) {
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 3
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:(k+2)
	#for (i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:(k-1)) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = 0
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = (x[index] - t[j])^3 / (t[j+2] - t[j]) / (t[j+1] - t[j]) / 3
	    
	    	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = x[index] - t[j+1] - (x[index] - t[j+2])^3 / (t[j+2] - t[j]) / (t[j+2] - t[j+1]) / 3 + (t[j+1] - t[j])^2 / 3 /(t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	    
	  	index = x > t[j+2] 
	    	sigma[j, index] = (x[index] - t[j+1]) + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) - (t[j+2] - t[j+1])^2 / 3 / (t[j+2] - t[j])
	}
	index = x >= t[1] & x <= t[k]
	sigma[k, index] = 0
	
	index = x > t[k] & x <= t[k+1]
	sigma[k, index] = (x[index] - t[k])^3 / (t[k+2] - t[k]) / (t[k+1] - t[k]) / 3
	   
	index = x > t[k+1] & x <= t[k+2]
	sigma[k, index] = x[index] - t[k+1] - (x[index] - t[k+2])^3 / (t[k+2] - t[k]) / (t[k+2] - t[k+1]) / 3 + (t[k+1] - t[k])^2 / 3 / (t[k+2] - t[k]) - (t[k+2] - t[k+1])^2 / 3 / (t[k+2] - t[k])
	
	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = x[index] - t[1] + (t[2] - x[index])^3 / (t[2] - t[1])^2 / 3 - (t[2] - t[1]) / 3
	
	index = x > t[2] 
	sigma[k+1, index] = x[index] - t[1] - (t[2] - t[1]) / 3
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = 0
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = (x[index] - t[k+1])^3 / (t[k+2] - t[k+1])^2 / 3
	
	sigma[k+3,] = x

	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}

########################################################
# Increasing and Concave
incconcave = function(xs, t, interp = FALSE) {
#xs = (xs - min(xs)) / (max(xs) - min(xs))
	x = sort(xs)
	n = length(x)
	k = length(t) - 2
	m = k + 3
	sigma = matrix(0, nrow = m, ncol = n)
	obs = 1:n
	#knt = 1:(k+2)
	#for(i in 1:(k+2)) {knt[i] = min(obs[abs(x - t[i]) == min(abs(x - t[i]))])}
	for (j in 1:k) {
	 	index = x >= t[1] & x <= t[j]
	 	sigma[j, index] = x[index] - t[1]
	 	
	 	index = x > t[j] & x <= t[j+1]
	 	sigma[j, index] = t[j] - t[1] + ((t[j+1] - t[j])^3 - (t[j+1] - x[index])^3) / 3 / (t[j+1] - t[j]) / (t[j+2] - t[j]) + (x[index] - t[j]) * (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	    
	    	index = x > t[j+1] & x <= t[j+2]
	    	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + ((t[j+2] - t[j+1])^3 - (t[j+2] - x[index])^3) / 3 / (t[j+2] - t[j+1]) / (t[j+2] - t[j])
	   
	    	index = x > t[j+2] 
	    	sigma[j, index] = t[j] - t[1] + (t[j+1] - t[j])^2 / 3 / (t[j+2] - t[j]) + (t[j+2] - t[j+1]) * (t[j+1] - t[j]) / (t[j+2] - t[j]) + (t[j+2] - t[j+1])^2 / 3 /(t[j+2] - t[j])
	    
	}

	index = x >= t[1] & x <= t[2]
	sigma[k+1, index] = -(t[2] - x[index])^3 / 3 / (t[2] - t[1])^2
	
	index = x > t[2] 
	sigma[k+1, index] = 0
	
	index = x >= t[1] & x <= t[k+1]
	sigma[k+2, index] = x[index] - t[1]
	
	index = x > t[k+1] & x <= t[k+2]
	sigma[k+2, index] = t[k+1] - t[1] + ((t[k+2] - t[k+1])^2 * (x[index] - t[k+1]) - (x[index] - t[k+1])^3 / 3) / (t[k+2] - t[k+1])^2
	
	sigma[k+3, ] = x

	ms = NULL
	if (!interp) {
		ms = apply(sigma, 1, mean)
		for (i in 1:m) {
			sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	} else {
		for (i in 1:m) {
			#sigma[i,] = sigma[i,] - mean(sigma[i,])
			sigma[i,] = sigma[i, rank(xs)]
		} 
	}
	rslt = list(sigma = sigma, ms = ms)
	rslt
}


##############
#summary.cgam#
##############
summary.cgam <- function(object,...) {
	if (!is.null(object$zcoefs)) {
		family <- object$family 
		resid_df_obs <- object$resid_df_obs
		wt.iter <- object$wt.iter
		coefs <- object$zcoefs
		se <- object$se.beta
		tval <- coefs / se
		pvalbeta <- object$pvals.beta
		n <- length(coefs)
		sse0 <- object$SSE0
		sse1 <- object$SSE1
		cic <- object$cic
		deviance <- object$deviance
		null_deviance <- object$null_deviance
		df <- object$df
		null_df <- object$null_df
		zid <- object$zid
		#zid1 <- object$zid1 - 1 - length(shapes)
		#zid2 <- object$zid2 - 1 - length(shapes)
#new: zid1, zid2 just index zmat not bigmat
		zid1 <- object$zid1
		zid2 <- object$zid2
		tms <- object$tms
		is_param <- object$is_param
		is_fac <- object$is_fac
		vals <- object$vals
		if (wt.iter) {
			rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "z.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))			
			rownames(rslt1)[1] <- "(Intercept)"
			if (n > 1) {
				lzid <- length(zid1)
				for (i in 1:lzid) {
					pos1 <- zid1[i]; pos2 <- zid2[i]
					for (j in pos1:pos2) {
						if (!is_param[i]) {
							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")						
						} else {
							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")						
						}	
					}
				}
			}
			rslt1 <- as.matrix(rslt1)
		} else {
			rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
			rownames(rslt1)[1] <- "(Intercept)"
			if (n > 1) {
				lzid <- length(zid1)
				for (i in 1:lzid) {
					pos1 <- zid1[i]; pos2 <- zid2[i]; 
					for (j in pos1:pos2) {
						if (!is_param[i]) {
							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")						
						} else {
							rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")						
						}	
					}
				}
			}
			rslt1 <- as.matrix(rslt1)
		}
		#if (!is.null(sse0) & !is.null(sse1)) {
			#rslt2 <- cbind(SSE.Linear = sse0, SSE.Full = sse1)
#new:
		#	rslt2 <- data.frame("SSE.Linear" = sse0, "SSE.Full" = sse1)
		#	rownames(rslt2)[1] <- ""
		#	ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs, cic = cic, null_deviance = null_deviance, null_df = null_df, deviance = deviance, df = df, resid_df_obs = resid_df_obs, family = family) 
		#	class(ans) <- "summary.cgam"
		#	ans
		#} else {
			ans <- list(call = object$call, coefficients = rslt1, zcoefs = coefs, cic = cic, null_deviance = null_deviance, null_df = null_df, deviance = deviance, df = df, resid_df_obs = resid_df_obs, family = family)
			class(ans) <- "summary.cgam"
			ans
		#}
	} else {
		ans <- list(zcoefs = object$zcoefs)
		class(ans) <- "summary.cgam"
		ans
	}
}

#############
#summary.wps#
#############
summary.wps <- function(object,...) {
	if (!is.null(object$zcoefs)) {
		coefs <- object$zcoefs
		se <- object$se.beta
		#tval <- object$tz
		pvalbeta <- object$pvals.beta
		tval <- coefs / se
		n <- length(coefs)
		#sse0 <- object$SSE0
		#sse1 <- object$SSE1
		zid <- object$zid
#new: zid1, zid2 just index zmat not bigmat
		zid1 <- object$zid1
		zid2 <- object$zid2
		#tms <- object$tms
		#zmat <- object$zmat
		#is_mat <- object$is_mat
		is_param <- object$is_param
		is_fac <- object$is_fac
		vals <- object$vals
		tms <- object$tms
		rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
		rownames(rslt1)[1] <- "(Intercept)"
		if (n > 1) {
			lzid <- length(zid1)
			for (i in 1:lzid) {
				pos1 <- zid1[i]; pos2 <- zid2[i]
				for (j in pos1:pos2) {
					if (!is_param[i]) {
						rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")						
					} else {
						rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")					
					}	
				}
			}
		}
		rslt1 <- as.matrix(rslt1)
		#if (!is.null(sse0) & !is.null(sse1)) {
		#	rslt2 <- data.frame("SSE.Linear" = sse0, "SSE.Full" = sse1)
#new:
		#	rownames(rslt2)[1] <- ""
			#rslt2 <- as.matrix(rslt2)
		#	ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs) 
		#	class(ans) <- "summary.wps"
		#	ans
		#} else {
			ans <- list(call = object$call, coefficients = rslt1, zcoefs = coefs)
			class(ans) <- "summary.wps"
			ans
		#}
	} else {
		ans <- list(zcoefs = object$zcoefs)
		class(ans) <- "summary.wps"
		ans
	}
}


#####################
#print.summary.cgam #
#####################
print.summary.cgam <- function(x,...) {
	if (!is.null(x$zcoefs)) {
	#if (!is.null(x$se.beta)) {
		cat("Call:\n")
		print(x$call)
		cat("\n")
		cat("Coefficients:")
		cat("\n")
		printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
		cat("\n")
		if (x$family$family == "binomial") {
			message("(Dispersion parameter for binomial family taken to be 1)", "\n")
		}
		if (x$family$family == "poisson") {
			message("(Dispersion parameter for poisson family taken to be 1)", "\n")
		} 
		if (x$family$family == "gaussian") {
			message("(Dispersion parameter for gaussian family taken to be ", round(x$deviance/x$df,4), " )", "\n")
		} 
		message("Null deviance: ", round(x$null_deviance,4), " ", "on ", x$null_df, " ", "degrees of freedom")
		message("Residual deviance: ", round(x$deviance,4), " ", "on ", x$resid_df_obs, " ", "observed degrees of freedom")
		#if (is.null(x$cic)) {
		#	message("CIC value is not available when there is no shape-restricted predictor")
		#} else {message("CIC: ", round(x$cic,4))}
		if (!is.null(x$cic)) {
			message("CIC: ", round(x$cic,4))
		}
		#if (!is.null(x$residuals)) {
		#	cat("==============================================================", "\n")
		#	cat("Call:\n")
		#	print(x$call)
		#	cat("\n")
		#	printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
		#}
	} else {
		print ("No linear predictor is defined")	
		#print ("Residual degree of freedom is negative!")
	}
}

####################
#print.summary.wps #
####################
print.summary.wps <- function(x,...) {
	if (!is.null(x$zcoefs)) {
	#if (!is.null(x$se.beta)) {
		cat("Call:\n")
		print(x$call)
		cat("\n")
		cat("Coefficients:")
		cat("\n")
		printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
		#if (!is.null(x$residuals)) {
			#cat("==============================================================", "\n")
		#	cat("+--------------------------------------+\n")
  		#	cat("|         SSE.Linear vs SSE.Full       |\n")
#printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
  		#	cat("+--------------------------------------+\n")
			#cat("Call:\n")
			#print(x$call)
			#cat("\n")
		#	printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
		#}
	} else {
		print ("No linear predictor is defined")	
	}
}


##############
#predict.cgam#
##############
predict.cgam = function(object, newData,...) {
#print (is.data.frame(newData))
#print (newData)
#new: 
	if (!is.data.frame(newData)) {
		#newData = as.data.frame(newData)
		stop ("newData must be a data frame!")	
	}
	family = object$family
	cicfamily = CicFamily(family)
	muhat.fun = cicfamily$muhat.fun	
	#shapes = object$shapes
#new: only use shapes for x != umb or tree
	shapes = object$shapesx
	np = object$d0; capm = object$capm; capk = object$capk; capt = object$capt; capu = object$capu
#new:
	xid10 = object$xid1; xid20 = object$xid2; 
	uid1 = object$uid1; uid2 = object$uid2; tid1 = object$tid1; tid2 = object$tid2 
#new:
	xmat0 = object$xmat; knots0 = object$knots; numknots0 = object$numknots; sps0 = object$sps; ms0 = object$ms
#new: s(x)'s are re-ordered in cgam, we need re-order xmat0 etc into a order with 17's first to match bigmat
idx_s = object$idx_s; idx = object$idx
if (!is.null(idx_s)) {
	xmat00 = xmat0
	knots00 = knots0
	numknots00 = numknots0
	sps00 = sps0
	ms00 = ms0
	xmat00[,1:length(idx_s)] = xmat0[,idx_s] 
	knots00[1:length(idx_s)] = knots0[idx_s] 	
	numknots00[1:length(idx_s)] = numknots0[idx_s] 
	sps00[1:length(idx_s)] = sps0[idx_s] 
	ms00[1:length(idx_s)] = ms0[idx_s] 
	if (!is.null(idx)) {
		capl = ncol(xmat0)
		xmat00[,(1+length(idx_s)):capl] = xmat0[,idx] 
		knots00[(1+length(idx_s)):capl] = knots0[idx] 	
		numknots00[(1+length(idx_s)):capl] = numknots0[idx] 
		sps00[(1+length(idx_s)):capl] = sps0[idx] 
		ms00[(1+length(idx_s)):capl] = ms0[idx] 
	}
xmat0 = xmat00
knots0 = knots00
numknots0 = numknots00
sps0 = sps00
ms0 = ms00
}
	zmat = object$zmat; umb = object$umb; tr = object$tr
#new:
	ztb = object$ztb; zid1 = object$zid1; zid2 = object$zid2; iz = 1
	bigmat = object$bigmat; umbrella.delta = object$umbrella.delta; tree.delta = object$tree.delta
	coefs = object$coefs; zcoefs = object$zcoefs; vcoefs = object$vcoefs; xcoefs0 = object$xcoefs; ucoefs = object$ucoefs; tcoefs = object$tcoefs
	tt = object$tms
	if (!inherits(object, "cgam")) { 
	        warning("calling predict.cgam(<fake-cgam-object>) ...")
        }
	if (missing(newData) || is.null(newData)) {
		etahat = object$etahat
		muhat = muhat.fun(etahat, fml = family$family)
		ans = list(fit = muhat, etahat = etahat, newbigmat = object$bigmat)
		return (ans) 
	}
	Terms = delete.response(tt)
#model.frame will re-organize newData in the original order in formula
	m = model.frame(Terms, newData)
#print (m)
	newdata = m
#print (head(newdata))
#new:
	newx0 = NULL; newxv = NULL; newx = NULL; newx_s = NULL; newu = NULL; newt = NULL; newz = NULL; newv = NULL
	#newz = list(); iz = 1;
	rn = nrow(newdata)
#print (rn)
#new:
	newetahat = 0; newmuhat = 0
#new: newx_sbasis
	newxbasis = NULL; newx_sbasis = NULL; newubasis = NULL; newtbasis = NULL; newbigmat = NULL
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
#get shape attributes and elem out of newdata
	for (i in 1:ncol(newdata)) {
		if (is.null(attributes(newdata[,i])$shape)) {
			if (is.factor(newdata[,i])) {
				lvli = levels(newdata[,i])
				ztbi = levels(ztb[[iz]])
				newdatai = NULL
				if (!any(lvli %in% ztbi)) {
					stop ("new factor level must be among factor levels in the fit!")
				} else {
					id1 = which(ztbi %in% lvli)
					#if (any(id1 > 1)) {
#delete the base level
					klvls = length(ztbi)
					if (klvls > 1) {						
						newimat = matrix(0, nrow = rn, ncol = klvls-1)
					        for (i1 in 1:rn) {
							if (newdata[i1,i] != ztbi[1]) {
								id_col = which(ztbi %in% newdata[i1,i]) - 1
								newimat[i1,id_col] = 1
							}
						}
						#for (i1 in 1: klvls) {
						#	which(levels(newdatai) == ztbi[i1])
						#	for (i2 in 1:rn) {
						#		if (levels(newdatai[i2, i1]) == lvli[i2]) {
						#			if (lvli[i2] != ztbi[1]) {
						#				newimat[i2, i1] = 1
						#			}
						#		}
						#	}
						#}
						#ztb_use = ztbi[-1]
						#nc = length(id1[id1 > 1])
						#newdatai = matrix(0, nrow = rn, ncol = nc)
						#for (ic in 1:nc) {
						#	id2 = which(newdata[,i] == ztb_use[ic])
						#	newdatai[id2, ic] = 1
						#}
						newdatai = newimat
					}
				}
				#if (length(levels(newdata[,i])) > 2) {
				#	klvls = length(levels(newdata[,i]))
				#	vals = as.numeric(levels(newdata[,i]))
				#	newimat = matrix(0, nrow = rn, ncol = klvls - 1)
				#	for (i1 in 1: (klvls - 1)) {
				#		for (i2 in 1:rn) {
				#			if (newdatai[i2] == vals[i1 + 1]) {
				#				newimat[i2, i1] = 1
				#			}
				#		}
				#	}
				#	newdatai = newimat
				#}
			} else {
				newdatai = newdata[,i]
			}
			newz = cbind(newz, newdatai)
			iz = iz + 1
#print (head(newz))
		}
		if (is.numeric(attributes(newdata[,i])$shape)) {
			newx0 = cbind(newx0, newdata[,i])
			if ((attributes(newdata[,i])$shape > 2 & attributes(newdata[,i])$shape < 5) | (attributes(newdata[,i])$shape > 10 & attributes(newdata[,i])$shape < 13)) {
				newxv = cbind(newxv, newdata[,i])
			}
		} 
		if (is.character(attributes(newdata[,i])$shape)) {
       			if (attributes(newdata[,i])$shape == "tree") {
				newt = cbind(newt, newdata[,i])
			}
       			if (attributes(newdata[,i])$shape == "umbrella") {
				newu = cbind(newu, newdata[,i])	
			}
     		}
	}
#print (head(newt))
#print (head(newu))
#print (head(newx0))
#print (head(newxv))
#print (head(newz))
#print (head(newv))
#new: separate x and x_s, move shape 17 to the beginning
	if (!is.null(shapes)) {
		if (any(shapes == 17)) {
			kshapes <- length(shapes)
        		obs <- 1:kshapes
        		idx_s <- obs[which(shapes == 17)]; idx <- obs[which(shapes != 17)]	
			newx1 <- newx0
			shapes0 <- 1:kshapes*0
			newx1[,1:length(idx_s)] <- newx0[,idx_s]
			shapes0[1:length(idx_s)] <- shapes[idx_s]   
			if (length(idx) > 0) {
				newx1[,(1 + length(idx_s)):kshapes] <- newx0[,idx]
				shapes0[(1 + length(idx_s)):kshapes] <- shapes[idx]
    			}
			newx0 <- newx1; shapes <- shapes0
		}
#new:
		if (all(shapes < 9)) {
			newx = newx0
			xid1 = xid10; xid2 = xid20
			xmat = xmat0
#new:
			sh_x = shapes 
			ms_x = ms0
		} else if (all(shapes > 8)) {
			newx_s = newx0
			xid1_s = xid10; xid2_s = xid20
			xmat_s = xmat0
			numknots = numknots0
			knots = knots0
			sps = sps0
			sh = shapes
			ms = ms0
		} else if (any(shapes > 8) & any(shapes < 9)) {
			newx = newx0[, shapes < 9, drop = FALSE]
			xmat = xmat0[, shapes < 9, drop = FALSE]
#new:
			sh_x = shapes[shapes < 9]
			ms_x = ms0[shapes < 9]
			xid1 = xid10[shapes < 9]; xid2 = xid20[shapes < 9]

			newx_s = newx0[, shapes > 8, drop = FALSE]
			xmat_s = xmat0[, shapes > 8, drop = FALSE]
			sh = shapes[shapes > 8]
			ms = ms0[shapes > 8]
			xid1_s = xid10[shapes > 8]; xid2_s = xid20[shapes > 8]
			numknots = numknots0[shapes > 8]
			knots = knots0[shapes > 8]
			sps = sps0[shapes > 8]
		}
	}
	#if (!is.null(shapes) & any(shapes == 17)) {
	if (!is.null(shapes)) {
		vcoefs = vcoefs[1:(1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))]
	} else {vcoefs = vcoefs[1:(1 + capk)]}
	if (capk > 0) {
		vcoefs_nz = vcoefs[-c(2:(1 + capk))]
	} else {vcoefs_nz = vcoefs}
	newv = cbind(1:rn*0 + 1, newz, newxv)
	newv_nz = cbind(1:rn*0 + 1, newxv)
#print (newv_nz)
#print (vcoefs)
	#etahat_v = as.vector(newv %*% vcoefs)
	etahat_v_nz = as.vector(newv_nz %*% vcoefs_nz)
#print (etahat_v_nz)
#new:   
	etahat_s = 1:rn*0; newx_sbasis = NULL; xs_coefs = NULL; var_xs = NULL
	if (!is.null(newx_s)) {
 		ks = ncol(newx_s)
		del = NULL
		for (i in 1:ks) {
			xi = xmat_s[,i]
			nxi = newx_s[,i]
			if (any(nxi > max(xi)) | any(nxi < min(xi))) {
				stop ("No extrapolation is allowed in cgam prediction!")
			}
#new: scale accordingly
#nxi = (nxi - min(xi)) / (max(xi) - min(xi))
			pos1 = xid1_s[i] - (1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))
			pos2 = xid2_s[i] - (1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))
			xs_coefs = c(xs_coefs, xcoefs0[pos1:pos2])
			msi = ms[[i]]
			deli_ans = makedelta(nxi, sh[i], numknots[i], knots[[i]], space = sps[i], suppre = TRUE, interp = TRUE)
			deli = deli_ans$amat
			if (sh[i] > 10 & sh[i] < 13) {
				#x = xmat_s[,i]
				xs = sort(xi)
				ord = order(xi)
				nx = length(xi)
				obs = 1:nx
				nr = nrow(deli)
				nc = length(nxi)
				ms0 = matrix(0, nrow = nr, ncol = nc)
				for (i1 in 1:nc) {
					for (i2 in 1:nr) {
						ms0[i2, i1] = my_line(xp = nxi[i1], y = msi[i2, ][ord], x = xs, end = nx, start = 1)$yp 
					}
				}
				deli = deli - ms0
			} else {
				deli = deli - msi
			}
#new:
			var_xs = c(var_xs, 1:nrow(deli)*0 + i)
			del = rbind(del, deli)
		}
		newx_sbasis = t(del)
		etahat_s = as.vector(newx_sbasis %*% xs_coefs)
	}
	etahat_x = 1:rn*0; newxbasis = NULL; xcoefs = NULL; var_x = NULL
	if (!is.null(newx)) {
		kx = ncol(xmat)
		del = NULL
		for (i in 1:kx) {
			xi = xmat[,i]
#new:
			nxi = newx[,i]
			if (any(nxi > max(xi)) | any(nxi < min(xi))) {
				stop ("No extrapolation is allowed in cgam prediction!")
			}
#new: scale accordingly
#nxi = (nxi - min(xi)) / (max(xi) - min(xi))
			shi = sh_x[i]
			pos1 = xid1[i] - (1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))
			pos2 = xid2[i] - (1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13))
			xcoefs = c(xcoefs, xcoefs0[pos1:pos2])
			deli = pred_del(xi, shi, nxi, ms_x[[i]])
#new:
			var_x = c(var_x, 1:nrow(deli)*0 + i)
			del = rbind(del, deli)
		}
#new: 
		newxbasis = t(del)
		etahat_x = as.vector(newxbasis %*% xcoefs)
	}
#new:
	etahat_u = 1:rn*0; newuedge = NULL
	if (!is.null(newu)) {	
		for (j in 1:ncol(umb)) {
			u = umb[,j]
			nu = length(u)
			us = sort(u)
			ord = order(u)
			if (any(newu[,j] > max(u)) | any(newu[,j] < min(u))) {
				stop ("no extrapolation is allowed in cgam!")
			}
			#u_edges = t(umbrella.fun(u))
			pos1 = uid1[j]; pos2 = uid2[j]
			u_edges = t(bigmat[pos1:pos2, , drop = FALSE])
			nuedge = ncol(u_edges)
			newuedge0 = matrix(0, nrow = rn, ncol = nuedge)
			for (i in 1:nuedge) {
				ue = u_edges[,i]; udist = round(diff(ue), 4); s = sign(udist)
				ueord = ue[ord]; udistord = round(diff(ueord), 4); sord = sign(udistord)
				obs = 1:(nu-1)
				posin = obs[sord != 0] + 1
				pos = unique(c(1, posin, nu))
				uk = us[pos]
				uek = ueord[pos]
				npos = length(pos)
				nd = length(newu[,j])
				for (l in 1:(npos-1)) {
					if (any(newu[,j] == uk[npos])) {
						ids = which(newu[,j] == uk[npos])
						newuedge0[ids,i] = uek[npos]
					}
					if (any(newu[,j] >= uk[l]) & any(newu[,j] < uk[l+1])) {
						ids = which(newu[,j] >= uk[l] & newu[,j] < uk[l+1])
						newuedge0[ids,i] = uek[l]
					}
				}
			}
		newuedge = cbind(newuedge, newuedge0) 
		}
		newubasis = newuedge
		etahat_u = as.vector(newubasis %*% ucoefs)
	}
# we don't have a newtbasis
# estimate tree
	etahat_t = 1:rn*0 
#print (newt)
	if (!is.null(newt)) {
#each column of  tru is a "table" of all levels of a tree var
		tru = unique(tr)#; placebo = min(tru)
		#tr_etahat = 0
		for (i in 1: ncol(newt)) {
			pos1 = tid1[i] - np - capm - capu
			pos2 = tid2[i] - np - capm - capu
			newtu = unique(newt[,i])
			#tcoefi = tcoefs[pos1:pos2]
#check!	add 0 to be the coef for placebo	
			tcoefi = c(0, tcoefs[pos1:pos2]) 
			if (!all(newtu %in% tru[,i])) {
				stop ("new tree ordering factor must be among the old tree ordering factors!")
			}
			#placebo = min(tru[,i])
#new:
			placebo = object$pl[i]
			if (any(newtu != placebo)) {	
				#id_cf = which(tru[,i] %in% newtu) - 1 #no coef for placebo 
				#tr_etahat = tr_etahat + tcoefi[id_cf]
				#id_cf = which(newtu > placebo)
				id_cf = which(newtu != placebo)
				t_use = newtu[id_cf]
				nt = length(t_use)
				for (it in 1:nt) {
					etahat_t[newt[,i] == t_use[it]] = etahat_t[newt[,i] == t_use[it]] + tcoefi[tru[,i] == t_use[it]]
				}
			}
		}
		#etahat_t = tr_etahat
	} 
#print (etahat_v)
#print (etahat_s)
#print (etahat_x)
#print (etahat_u)
#print (etahat_t)
	etahat_z = 1:rn*0; zcf = NULL
#print (newt)
#print (newdata)
	if (!is.null(newz)) {
#delete the coef for 1
		zcoefs = zcoefs[-1]
		etahat_z = newz %*% zcoefs
		#ztbu = unique(ztb)#; placebo = min(ztbu)
		#lz = ncol(newz)
		#for (i in 1:lz) {
		#	pos1 = zid1[i]
		#	pos2 = zid2[i]
		#	zcoefi = zcoefs[pos1:pos2]
		#	zlvi = zcoefs[pos1:pos2]
		#	etahat_z[newz[,i] == 1] = etahat_z[newz[,i] == 1] + zcoefi
		#}
	} 
#print (etahat_z)
	etahat_v = etahat_v_nz + etahat_z	
	newetahat = etahat_v + etahat_s + etahat_x + etahat_u + etahat_t
	newmuhat = as.vector(muhat.fun(newetahat, fml = family$family))
	if (!is.null(newt)) {
		newtbasis = t(tree.delta[,1:nrow(newData),drop = FALSE])
	}
#print (newv)
	#newbigmat = t(cbind(newv, newx_sbasis, newxbasis, newubasis, newtbasis))
	#ans = list(v = newv, xbs = newxbasis, x_sbs = newx_sbasis, ubs = newubasis, tbs = newtbasis, bigmat = newbigmat, vcoefs = vcoefs, xcoefs = xcoefs, xs_coefs = xs_coefs, ucoefs = ucoefs, etahat = newetahat, muhat = newmuhat)
	#ans = list(muhat = newmuhat)
	muhat = newmuhat
	return (muhat) 
}

#############################
#predict delta for ordinal x#
#############################
pred_del = function(x, sh, xp, ms) {
	n = length(xp)
#x = (x - min(x)) / (max(x) - min(x))
	xu = sort(unique(x))
	n1 = length(xu)
	sigma = NULL
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
#  increasing or decreasing
	if (sh < 3) {
		sigma = matrix(0, nrow = n1 - 1, ncol = n)
		for (i in 1: (n1 - 1)) {
			sigma[i, xp > xu[i]] = 1
		}
		if (sh == 2) {sigma = -sigma}
		for (i in 1:(n1 - 1)) {sigma[i, ] = sigma[i, ] - ms[i]}
	} 
	if (sh == 3 | sh == 4) {
#  convex or concave
		sigma = matrix(0, nrow = n1 - 2, ncol = n)
		#for (i in 1: (n1 - 2)) {
		#	sigma[i, x > xu[i]] = x[x > xu[i]] - xu[i]
		#}
		for (i in 1: (n1 - 2)) {
			sigma[i, xp > xu[i+1]] = xp[xp > xu[i+1]] - xu[i+1]
		}
		if (sh == 4) {sigma = -sigma}
		#xm = cbind(1:n*0+1, xp)
		#xpx = solve(t(xm) %*% xm)
		#pm = xm %*% xpx %*% t(xm)
		#sigma = sigma - sigma %*% t(pm)
		xs = sort(x)
		ord = order(x)
		nx = length(x)
		obs = 1:nx
		m = nrow(ms)
		ms0 = matrix(0, nrow = m, ncol = n)
		for (i1 in 1:n) {
			for (i2 in 1:m) {
				ms0[i2, i1] = my_line(xp = xp[i1], y = ms[i2, ][ord], x = xs, end = nx, start = 1)$yp 
			}
		}
		sigma = sigma - ms0
	} 
	if (sh > 4 & sh < 9) {
		sigma = matrix(0, nrow = n1 - 1, ncol = n)
		if (sh == 5) { ### increasing convex
			for (i in 1:(n1 - 1)) {
				sigma[i, xp > xu[i]] = (xp[xp > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			for (i in 1:(n1 - 1)) {sigma[i,] = sigma[i,] - ms[i]}
		} else if (sh == 6) {  ## decreasing convex
			for (i in 1:(n1 - 1)) {
				sigma[i, xp < xu[i + 1]] = (xp[xp < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			for (i in 1:(n1 - 1)) {sigma[i,] = sigma[i,] - ms[i]}
#print (ms)
		} else if (sh == 7) { ## increasing concave
			for (i in 1:(n1 - 1)) {
				sigma[i, xp < xu[i + 1]] = (xp[xp < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			for (i in 1:(n1 - 1)) {sigma[i,] = -sigma[i,] + ms[i]}		
		} else if (sh == 8) {## decreasing concave
			for (i in 1:(n1 - 1)) {
				sigma[i, xp > xu[i]] = (xp[xp > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			for (i in 1:(n1 - 1)) {sigma[i,] = -sigma[i,] + ms[i]}
		}
	}
	return (sigma)
}
###########################################
#create a 3D plot for a cgam or wps object#
###########################################
plotpersp <- function(object,...)UseMethod("plotpersp")

################
#plotpersp.cgam#
################ 
plotpersp.cgam <- function(object, x1 = NULL, x2 = NULL, data = NULL, surface = "mu", categ = NULL, col = NULL, random = FALSE, x_grid = 10, y_grid = 10, xlim = range(x1), ylim = range(x2), zlim = NULL, xlab = NULL, ylab = NULL, zlab = NULL, th = NULL, ltheta = NULL, main = NULL, ticktype = "detailed",...) {
	if (!inherits(object, "cgam")) { 
	        warning("calling plotpersp(<fake-cgam-object>) ...")
        }
	#xmat <- object$xmat
	#cl <- match.call()
	#nms <- cl[-c(1, 2)]
	#lnms <- length(nms)
	#x1nm <- nms[1]$x
	#x1nm <- deparse(x1nm)
	#x2nm <- nms[2]$x
	#x2nm <- deparse(x2nm)
#new: default is plotpersp(object)	
	x1nm <- deparse(substitute(x1))
	x2nm <- deparse(substitute(x2))
#stop (print (x1nm))	
	xnms <- object$xnms
	xmat <- object$xmat
	if (x1nm == "NULL" | x2nm == "NULL") {
		if (length(xnms) >= 2) {
			x1nm <- xnms[1]
			x2nm <- xnms[2]	
			x1id <- 1
			x2id <- 2
			x1 <- xmat[, 1]
			x2 <- xmat[, 2]
		} else {stop ("Number of non-parametric predictors must >= 2!")}
	}
	ynm <- object$ynm
	#xmat <- object$xmat
#print (dim(xmat))
	is_fac <- object$is_fac
	is_param <- object$is_param
	family <- object$family
	fml <- family$family
	cicfamily <- CicFamily(family)
	muhat.fun <- cicfamily$muhat.fun
	znms <- object$znms
	kznms <- length(znms)
	if (!is.null(categ)) {
		if (!is.character(categ)) {
			warning("categ must be a character argument!")
		} else if (!any(znms == categ)) {
#in.or.out case
			#if (!is.null(attr(object, "sub"))) {
			#	if (!is.null(categ)) {
			#		categ = paste("factor(", categ, ")", sep = "")
			#	}
			#} else {
			#	warning(paste(categ, "is not an exact character name defined in the cgam fit!"))
			#	categ = NULL
			#}
			if (any(grepl(categ, znms))) {
				id = which(grepl(categ, znms))
				znmi = znms[id] 
				if (grepl("as.factor", znmi)) {
					categ = paste("as.factor(", categ, ")", sep = "")
				} else if (grepl("factor", znmi)) {
					categ = paste("factor(", categ, ")", sep = "")
				} else {print(paste(categ, "is not an exact character name defined in the cgam fit!"))}
			} else {print(paste(categ, "is not an exact character name defined in the cgam fit!"))}
		} else {
			obsz = 1:kznms
			zid = obsz[znms == categ]
#linear term:
			if (!(is_fac[zid])) {
				categ = NULL
			}
		}
	}
	shapes <- object$shapes
#new:	
	#zid1 = object$zid1 - 1 - length(shapes)
	#zid2 = object$zid2 - 1 - length(shapes)
	zid1 <- object$zid1
	zid2 <- object$zid2

	kznms <- length(znms)
	zmat <- object$zmat
	#zcoefs = object$zcoefs
#new: exclude the coef for the one vector
	zcoefs <- object$zcoefs[-1]
	#xmatnms <- object$xmatnms
	knms <- length(xnms)	
	obs <- 1:knms

	#if (!any(xmatnms == x1nm)) {
	#	warning(paste(x1nm, "is not an exact character name defined in the cgam fit!"))
	#}
	#if (!any(xmatnms == x2nm)) {
	#	warning(paste(x2nm, "is not an exact character name defined in the cgam fit!"))
	#}
	#x1id = obs[xmatnms == x1nm]
	#x2id = obs[xmatnms == x2nm]
	if (!is.null(data)) {
		if (!is.data.frame(data)) {
			stop ("User need to make the data argument a data frame with names for each variable!")
		}
		datnms <- names(data)
		if (!any(datnms == x1nm) | !any(datnms == x2nm)) {
			stop ("Check the accuracy of the names of x1 and x2!")
		}
		x1 <- data[ ,which(datnms == x1nm)]
		x2 <- data[ ,which(datnms == x2nm)]
		if (length(x1) != nrow(xmat)) {
			warning ("Number of observations in the data set is not the same as the number of elements in x1!")
		}
		#bool <- apply(xmat, 2, function(x) all(x1 == x))
		#if (any(bool)) {
			x1id <- obs[xnms == x1nm]
		#}
		if (length(x2) != nrow(xmat)) {
			warning ("Number of observations in the data set is not the same as the number of elements in x2!")
		}
		#bool <- apply(xmat, 2, function(x) all(x2 == x))
		#if (any(bool)) {
			x2id <- obs[xnms == x2nm]
		#}
	} else {
		#if (any(xmatnms == x1nm)) {
		#	x1id <- obs[xmatnms == x1nm]
		#} else {
		#	bool <- apply(xmat, 2, function(x) all(x1 == x))
		#	if (any(bool)) {
		#		x1id <- obs[bool]
		#	}
		#}
		#if (any(xmatnms == x2nm)) {
		#	x2id <- obs[xmatnms == x2nm]
		#} else {
		#	bool <- apply(xmat, 2, function(x) all(x2 == x))
		#	if (any(bool)) {
		#		x2id <- obs[bool]
		#	}
		#}
#new: x1 and x2 are in .Globe, not in formula
		if (all(xnms != x1nm)) {
			if (length(x1) != nrow(xmat)) {
				stop ("Number of observations in the data set is not the same as the number of elements in x1!")
			}
			bool <- apply(xmat, 2, function(x) all(x1 == x))
			if (any(bool)) {
				x1id <- obs[bool]
#change x1nm to be the one in formula
				x1nm <- xnms[bool] 
			} else {
				stop (paste(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the cgam fit!"))
			}
		} else {
			x1id <- obs[xnms == x1nm]
		}
		if (all(xnms != x2nm)) {
			if (length(x2) != nrow(xmat)) {
				stop ("Number of observations in the data set is not the same as the number of elements in x2!")
			}
			bool <- apply(xmat, 2, function(x) all(x2 == x))
			if (any(bool)) {
				x2id <- obs[bool]
				x2nm <- xnms[bool] 
			} else {
				stop (paste(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the cgam fit!"))
			}
		} else {
			x2id <- obs[xnms == x2nm]
		}
	}
#xmat is not the one in fit
#print (length(x1))
#print (length(x2))
	#xm <- cbind(x1, x2)
	xm <- xmat[, c(x1id, x2id)]
#print (all(xm == cbind(x1, x2)))
#print (head(cbind(x1, x2)))
	x1g <- 0:x_grid / x_grid * .95 * (max(xm[,1]) - min(xm[,1])) + min(xm[,1]) + .025 * (max(xm[,1]) - min(xm[,1]))
 	n1 <- length(x1g)
	x2g <- 0:y_grid / y_grid * .95 * (max(xm[,2]) - min(xm[,2])) + min(xm[,2]) + .025 * (max(xm[,2]) - min(xm[,2]))
	n2 <- length(x2g)	

	xgmat <- matrix(nrow = n1, ncol = n2)
	eta0 <- object$coefs[1]
	thvecs <- object$etacomps
	for (i2 in 1:n2) {
		for (i1 in 1:n1) {
			x1a <- max(xm[xm[,1] <= x1g[i1], 1])
			x1b <- min(xm[xm[,1] > x1g[i1], 1])
			v1a <- min(thvecs[x1id, xm[,1] == x1a])
			v1b <- min(thvecs[x1id, xm[,1] == x1b])
			alp <- (x1g[i1] - x1a) / (x1b - x1a)
			th1add <- (1 - alp) * v1a + alp * v1b
			x2a <- max(xm[xm[,2] <= x2g[i2],2])
			x2b <- min(xm[xm[,2] > x2g[i2],2])
			v2a <- min(thvecs[x2id, xm[,2] == x2a])
			v2b <- min(thvecs[x2id, xm[,2] == x2b])
			alp <- (x2g[i2] - x2a) / (x2b - x2a)
			th2add <- (1 - alp) * v2a + alp * v2b	
			xgmat[i1,i2] <- eta0 + th1add + th2add
		}
	}
	x3_add <- 0
	if (knms >= 3) {
		x3id <- obs[-c(x1id, x2id)]
		kx3 <- length(x3id)
		for (i in 1:kx3) {
			x3i <- xmat[, x3id[i]]
			x3i_use <- max(x3i[x3i <= median(x3i)])
			x3i_add <- min(thvecs[x3id[i], x3i == x3i_use])			
			x3_add <- x3_add + x3i_add
		}
	} 
	if (surface == "eta") {
		xgmat <- xgmat + as.numeric(x3_add)
	}	
	if (is.null(categ) & surface == "mu") {
		z_add <- 0
		if (!is.null(znms)) {
			kzids <- length(zid1)
			for (i in 1:kzids) {
				pos1 <- zid1[i]; pos2 <- zid2[i]
				zi <- zmat[, pos1:pos2, drop = FALSE]
				zcoefsi <- zcoefs[pos1:pos2]
				for (j in 1:ncol(zi)){
					uzij <- unique(zi[,j])
					kuzij <- length(uzij)
					nmodej <- sum(zi[,j] == uzij[1])
					zij_mode <- uzij[1]
					for (u in 2:kuzij) {
						if (sum(zi[,j] == uzij[u]) > nmodej) {
							zij_mode <- uzij[u]
							nmodej <- sum(zi[,j] == uzij[u])
						}
					}
					obsuzij <- 1:length(uzij)
					uzhatij <- uzij * zcoefsi[j] 
					zij_add <- uzhatij[obsuzij[uzij == zij_mode]]	
					z_add <- z_add + zij_add 
				}
			}
		}
		xgmat <- xgmat + as.numeric(x3_add) + as.numeric(z_add)
		for (i2 in 1:n2) {
			for (i1 in 1:n1) {
				xgmat[i1, i2] <- muhat.fun(xgmat[i1, i2], fml = fml)	
			}
		}
	} else if (!is.null(categ) & surface == "mu"){
		xgmats <- list()
		mins <- NULL; maxs <- NULL
		obsz <- 1:kznms
		zid <- obsz[znms == categ]
#print (class(znms[znms == categ]))
		pos1 <- zid1[zid]; pos2 <- zid2[zid]
#print (pos1)
#print (pos2)
		zi <- zmat[, pos1:pos2, drop = FALSE]
		z_add <- 1:nrow(zi)*0
		zcoefsi <- zcoefs[pos1:pos2]
#print (zcoefsi)
		for (j in 1:ncol(zi)) {
			zij <- zi[,j]
			zijhat <- zij * zcoefsi[j]
			z_add <- z_add + zijhat
		}
		z_add <- unique(z_add)
		kz_add <- length(z_add)
#new: plot the smallest one first:
		z_add <- z_add[order(z_add)]
#print (z_add)
		for (iz in 1:kz_add) {
			xgmats[[iz]] <- xgmat + as.numeric(x3_add) + z_add[iz]
			#mins <- c(mins, min(xgmats[[iz]]))
			#maxs <- c(maxs, max(xgmats[[iz]]))
			for (i2 in 1:n2) {
				for (i1 in 1:n1) {
					xgmats[[iz]][i1, i2] <- muhat.fun(xgmats[[iz]][i1, i2], fml = fml)	
				}
			}
			mins <- c(mins, min(xgmats[[iz]]))
			maxs <- c(maxs, max(xgmats[[iz]]))

		}
	}
	if (is.null(xlab)) {
		#xlab = deparse(x1nm)
		xlab <- x1nm
	}
	if (is.null(ylab)) {
		#ylab = deparse(x2nm)
		ylab <- x2nm
	}
	if (is.null(zlab)) {
		if (surface == "mu") {
			if (fml == "binomial") {
				zlab <- paste("Pr(", ynm, ")")
			} else if (fml == "poisson" | fml == "gaussian") {
				zlab <- paste("est mean of", ynm)
			}		
		}
		if (surface == "eta") {
			if (fml == "binomial") {
				zlab <- paste("est log odds ratio of", ynm)
			}  else if (fml == "poisson") {
				zlab <- paste("est log mean of", ynm)
			} else if (fml == "gaussian") {
				zlab <- paste("est mean of", ynm)
			}
		}
	}
	if (is.null(zlim)) {
		zlim <- range(xgmat, na.rm = TRUE)
	}
	palette <- c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
	if (!is.null(categ) & surface == "mu") {
		#palette = c("peachpuff", "lightblue", "grey", "wheat", "yellowgreen", "plum", "limegreen", "paleturqoise", "azure", "whitesmoke")
		kxgm <- length(xgmats)
		if (is.null(col)) {
			#if (kxgm == 2) {
			#	col = c("peachpuff", "lightblue")
			#} else if (kxgm == 3) { 
			#	col = c("peachpuff", "lightblue", "grey")
			#} else if (kxgm > 3 & kxgm < 11) {
			#	col = sample(palette, replace = FALSE)
			if (random) {
				col <- topo.colors(kxgm)
				#col <- sample(palette, size = kxgm, replace = FALSE)
#print (col)
			} else {
#print (kxgm)
				if (kxgm > 1 & kxgm < 11) {
					col <- palette[1:kxgm]
				} else {
					#integ <- floor(kxgm / 10)
					#rem <- kxgm %% 10
					#kint <- length(integ)
					#col = character(length = kxgm)
#print (col)
					#col <- NULL
					#for (i in 1:kint) {
#print (col[1 + (i - 1) * 10: i * 10])
#print (palette)
						#col[1 + (i - 1) * 10: i * 10] = palette
					#	col <- c(col, palette)
					#}
#print (col)
					#col[(kint * 10 + 1):kxgm] = palette[(kint * 10 + 1):kxgm]
					#col <- c(col, palette[1:rem])
#print ((kint * 10 + 1):kxgm)
#print (col)
#print (integ)
#new: use rainbow
					col <- topo.colors(kxgm)
				}
			} 
		} else {
			if (length(col) < kxgm) {
				#rem <- kxgm - length(col)
				#nrem <- length(rem)
				#rem_col <- palette[1:nrem]
				#col <- c(col, rem_col)
#new:
				col <- topo.colors(kxgm)
			} else if (length(col) > kxgm) {
				col <- col[1:kxgm]
				#print (paste("The first", kxgm, "colors are used!"))
			}
			#if (random) {
			#	print ("User defined colors are used!")
			#}
		}
#print (kxgm)
#new: set th for decr or incr
		decrs = shapes[c(x1id, x2id)] %in% c(2, 6, 8, 10, 15, 16)
		incrs = shapes[c(x1id, x2id)] %in% c(1, 5, 7, 9, 13, 14, 17)
		if (is.null(th) | !is.numeric(th)) {
			ang = NULL
			if (incrs[1] & incrs[2]) {
				if (is.null(ang)) {
					ang = -40
				}
			} else if (decrs[1] & incrs[2]) { 
				if (is.null(ang)) {
					ang = 40
				}
			} else if (incrs[1] & decrs[2]) {
				if (is.null(ang)) {
					ang = 240
				}
			} else if (decrs[1] & decrs[2]) {
				if (is.null(ang)) {
					ang = 140
				}
			} else {ang = -37}
		} else {ang = th}
		if (is.null(ltheta) | !is.numeric(ltheta)) {
			ltheta <- -135
		}
		for (i in 1:kxgm) {
#i = 1
#print (i)
			xgmat <- xgmats[[i]]
			#if (is.null(th) | !is.numeric(th)) {
			#	th <- -40
			#}
			#if (is.null(ltheta) | !is.numeric(ltheta)) {
			#	ltheta <- -135
			#}
#persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, theta = th)
#new: avoid thick labs
box = TRUE
axes = TRUE
if (i > 1) {
xlab = ylab = zlab = ""
box = FALSE
axes = FALSE
}
			persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, zlim = c(min(mins), max(maxs)), xlab = xlab, ylab = ylab, zlab = zlab, col = col[i], main = main, theta = ang, ltheta = ltheta, ticktype = ticktype, box = box, axes = axes,...)


			par(new = TRUE)
		}
	par(new = FALSE)
	} else {
		if (is.null(col)) {
			if (random) {
				col <- sample(palette, size = 1, replace = FALSE)
			} else {
				col <- "white"
			}
		} else {
			if (length(col) > 1) {
				col <- col[1]
				print ("The first color is used!")
			} 
			if (random) {
				print ("User defined color is used!")
			}
		}
		#if (is.null(th) | !is.numeric(th)) {
		#	th <- -40
		#}
		#if (is.null(ltheta) | !is.numeric(ltheta)) {
		#	ltheta <- -135
		#}
#new: set th for decr or incr
		decrs = shapes[c(x1id, x2id)] %in% c(2, 6, 8, 10, 15, 16)
		incrs = shapes[c(x1id, x2id)] %in% c(1, 5, 7, 9, 13, 14, 17)
		if (is.null(th) | !is.numeric(th)) {
			ang = NULL
			if (incrs[1] & incrs[2]) {
				if (is.null(ang)) {
					ang = -40
				}
			} else if (decrs[1] & incrs[2]) { 
				if (is.null(ang)) {
					ang = 40
				}
			} else if (incrs[1] & decrs[2]) {
				if (is.null(ang)) {
					ang = 240
				}
			} else if (decrs[1] & decrs[2]) {
				if (is.null(ang)) {
					ang = 140
				}
			} else {ang = -37}
		} else {ang = th}
		if (is.null(ltheta) | !is.numeric(ltheta)) {
			ltheta <- -135
		}
		persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, zlim = zlim, xlab = xlab, ylab = ylab, zlab = zlab, col = col, main = main, theta = ang, ltheta = ltheta, ticktype = ticktype,...)
	}
#print (col)
}

#################
#plotpersp.wps#
################ 
plotpersp.wps = function(object, x1 = NULL, x2 = NULL, data = NULL, surface = "C", categ = NULL, col = NULL, random = FALSE, xlim = range(x1), ylim = range(x2), zlim = NULL, xlab = NULL, ylab = NULL, zlab = NULL, th = NULL, ltheta = NULL, main = NULL, ticktype = "detailed",...) {
	if (!inherits(object, "wps")) { 
	        warning("calling plotpersp(<fake-wps-object>) ...")
        }
#new:
#if (is.null(x1) | is.null(x2)) {
#	stop ("x1 and x2 must be provided!")
#}
#print (is.null(x1))
	xnms = object$xnms
	xmat = object$xmat
#print (match.call())
#print (length(match.call()))
#print (class(match.call()[-c(1, 2)]))
#if () { 
#	cl = match.call()
#	nms = cl[-c(1, 2)]
	#lnms = length(nms)
#	x1nm0 = nms[1]$x
#	x1nm0 = deparse(x1nm0)
#	x2nm0 = nms[2]$x
#	x2nm0 = deparse(x2nm0)
#} else {
#new:
	x1nm0 = xnms[1]
	x2nm0 = xnms[2]
	x1 = xmat[, 1]
	x2 = xmat[, 2]
#}
	#knms = length(xnms)	
	#obs = 1:knms
	#x1nm0 = nms[1]$x
	#x1nm0 = deparse(x1nm0)
	#x2nm0 = nms[2]$x
	#x2nm0 = deparse(x2nm0)
	is_fac = object$is_fac
	#is_param = object$is_param
	ynm = object$ynm
#xmat is delta
	#delta = object$delta
	znms = object$znms
	decrs = object$decrs
	kznms = length(znms)
#zmat include 1 vector
	zmat = object$zmat
	zmat0 = zmat[, -1, drop = FALSE]
	zcoefs = object$zcoefs[-1]
	zid1 = object$zid1
	zid2 = object$zid2
	ah = object$coefs
	ahu = object$coefsu
	k1 = object$k1
	k2 = object$k2
	m1 = length(k1)
	m2 = length(k2)
	p = dim(zmat)[2]

	if (!is.null(categ)) {
		if (!is.character(categ)) {
			warning("categ must be a character argument!")
		} else if (!any(znms == categ)) {
#print ('TRUE')	
			warning(paste(categ, "is not an exact character name defined in the cgam fit!"))
			categ = NULL
		} else {
			obsz = 1:kznms
			zid = obsz[znms == categ]
			if (!(is_fac[zid])) {
				categ = NULL
			}
		}
	}
#new: switch xnms
	if (!is.null(data)) {
		if (!is.data.frame(data)) {
			stop ("User need to make the data argument a data frame with names for each variable!")
		}
		datnms = names(data)
		if (!any(datnms == x1nm0) | !any(datnms == x2nm0)) {
			stop ("Check the accuracy of the names of x1 and x2!")
		}
		x1 = data[ ,which(datnms == x1nm0)]
		x2 = data[ ,which(datnms == x2nm0)]
		#if (x1nm0 != xnms[1] & x2nm0 != xnms[2]) {
		#	x1nm = x2nm0
		#	x2nm = x1nm0
		#	tmp = x1
		#	x1 = x2
		#	x2 = tmp
		#} else {x1nm = x1nm0; x2nm = x2nm0}
	} else {
		if (all(xnms != x1nm0)) {
			#stop (paste(paste("'", x1nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
#new: in case of wrong data fame
			if (length(x1) != nrow(xmat)) {
				stop ("Number of observations in the data set is not the same as the number of elements in x1!")
			}
			bool = apply(xmat, 2, function(x) all.equal(x1, x))
			if (any(bool)) {
				#x1id = obs[bool]
				x1nm0 = xnms[bool] 
			} else {
				stop (paste(paste("'", x1nm0, "'", sep = ''), "is not a predictor defined in the wps fit!"))
			}
		} 
		if (all(xnms != x2nm0)) {
			#stop (paste(paste("'", x2nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
			if (length(x2) != nrow(xmat)) {
				stop ("Number of observations in the data set is not the same as the number of elements in x2!")
			}
			bool = apply(xmat, 2, function(x) all(x2 == x))
			if (any(bool)) {
				#x2id = obs[bool]
				x2nm0 = xnms[bool] 
			} else {
				stop (paste(paste("'", x2nm0, "'", sep = ''), "is not a predictor defined in the wps fit!"))
			}
		} 
	} 
	if (x1nm0 != xnms[1] & x2nm0 != xnms[2]) {
		x1nm = x2nm0
		x2nm = x1nm0
		tmp = x1
		x1 = x2
		x2 = tmp
	} else {x1nm = x1nm0; x2nm = x2nm0}
#print (dim(xmat))
	apl = 1:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
	aplu = 1:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
	apl[1] = ah[1]
	aplu[1] = ahu[1]
	apl[2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))] = ah[(p + 1):(m1 + m2 - 1 + (m1 - 1) * (m2 - 1) + p - 1)]
	aplu[2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))] = ahu[(p + 1):(m1 + m2 - 1 + (m1 - 1) * (m2 - 1) + p - 1)]
	mupl = matrix(apl[1], nrow = m1, ncol = m2)
	muplu = matrix(aplu[1], nrow = m1, ncol = m2)
	for (i1 in 2:m1) {
		mupl[i1, ] = mupl[i1, ] + apl[i1]
		muplu[i1, ] = muplu[i1, ] + aplu[i1]
	}
	for (i2 in 2:m2) {
		mupl[, i2] = mupl[, i2] + apl[m1 - 1 + i2]
		muplu[, i2] = muplu[, i2] + aplu[m1 - 1 + i2]
	}
	for (i1 in 2:m1) {
		for (i2 in 2:m2) {
			mupl[i1, i2] = mupl[i1, i2] + apl[m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2]
			muplu[i1, i2] = muplu[i1, i2] + aplu[m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2]
		}
	}
## reverse transform for decreasing
	if (is.null(th) | !is.numeric(th)) {
		ang = NULL
		if (!decrs[1] & !decrs[2]) {
			if (is.null(ang)) {
				ang = -40
			}
		} else if (decrs[1] & !decrs[2]) { 
			#k1 = -k1[m1:1]; 
			mupl = mupl[m1:1, 1:m2]; muplu = muplu[m1:1, 1:m2]
			if (is.null(ang)) {
				ang = 40
			}
		} else if (!decrs[1] & decrs[2]) {
			#k2 = -k2[m2:1]; 
			mupl = mupl[1:m1, m2:1]; muplu = muplu[1:m1, m2:1]
			if (is.null(ang)) {
				ang = 240
			}
		} else if (decrs[1] & decrs[2]) {
			#k1 = -k1[m1:1]; k2 = -k2[m2:1]; 
			mupl = mupl[m1:1, m2:1]; muplu = muplu[m1:1, m2:1]
			if (is.null(ang)) {
				ang = 140
			}
		}
	} else {ang = th}
	if (is.null(ltheta) | !is.numeric(ltheta)) {
		ltheta <- -135
	}
	if (is.null(categ)) {
		z_add = 0
		if (!is.null(znms)) {
			kzids = length(zid1)
			for (i in 1:kzids) {
				pos1 = zid1[i]; pos2 = zid2[i]
#zi is a factor 
				zi = zmat0[, pos1:pos2, drop = FALSE]
				zcoefsi = zcoefs[pos1:pos2]
				for (j in 1:ncol(zi)){
#find the 'mode' of the jth column of zi; add the coef corresponding to the 'mode'
					uzij = unique(zi[,j])
					kuzij = length(uzij)
					nmodej = sum(zi[,j] == uzij[1])
					zij_mode = uzij[1]
					for (u in 2:kuzij) {
						if (sum(zi[,j] == uzij[u]) > nmodej) {
							zij_mode = uzij[u]
							nmodej = sum(zi[,j] == uzij[u])
						}
					}
					obsuzij = 1:length(uzij)
					uzhatij = uzij * zcoefsi[j] 
					zij_add = uzhatij[obsuzij[uzij == zij_mode]]	
					z_add = z_add + zij_add 
				}
			}
		}
		mupl = mupl + as.numeric(z_add)
		muplu = muplu + as.numeric(z_add)
		mins = min(mupl); maxs = max(mupl)
		minsu = min(muplu); maxsu = max(muplu)
	} else {
		mupls = muplus = list()
		mins = maxs = NULL
		minsu = maxsu = NULL		
		obsz = 1:kznms
		zid = obsz[znms == categ]
		pos1 = zid1[zid]; pos2 = zid2[zid]
		zcoefsi = zcoefs[pos1:pos2]
#include the base level
		zcoefsi = c(0, zcoefsi)
		z_add = sort(zcoefsi)
		kz_add = length(z_add)
		for (iz in 1:kz_add) {
			mupls[[iz]] = mupl + z_add[iz]
			mins = c(mins, min(mupls[[iz]]))
			maxs = c(maxs, max(mupls[[iz]]))		
			muplus[[iz]] = muplu + z_add[iz]
			minsu = c(minsu, min(muplus[[iz]]))
			maxsu = c(maxsu, max(muplus[[iz]]))
		}

	}
	palette = c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
	if (is.null(categ)) {
		if (is.null(col)) {
			if (random) {
				col = sample(palette, size = 1, replace = FALSE)
			} else {
				col = "white"
			}
		} 
		if (surface == 'C') {
			musurf = mupl
			if (is.null(main)) {
				main = 'Constrained Warped-Plane Spline Surface'
			}
			zlim0 = c(min(mins), max(maxs))
		} else if (surface == 'U') {
			musurf = muplu
			if (is.null(main)) {
				main =  'Unconstrained Warped-Plane Spline Surface'
			}
			zlim0 = c(min(minsu), max(maxsu))
		}
		persp(k1, k2, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype,...)
	} else {
		kxgm = length(mupls)
		if (is.null(col)) {
			if (random) {
#new:
				col = topo.colors(kxgm)
				#col = sample(palette, size = kxgm, replace = FALSE)
			} else {
				if (kxgm > 1 & kxgm < 11) {
					col = palette[1:kxgm]
				} else {
					#integ = floor(kxgm / 10)
					#rem = kxgm %% 10
					#kint = length(integ)
					#col = NULL
					#for (i in 1:kint) {
					#	col = c(col, palette)
					#}
					#col = c(col, palette[1:rem])
#new:
					col = topo.colors(kxgm)
				}
			} 
		} else {
			if (length(col) < kxgm) {
				#rem = kxgm - length(col)
				#nrem = length(rem)
				#rem_col = palette[1:nrem]
				#col = c(col, rem_col)
#new:
				col = topo.colors(kxgm)
			} else if (length(col) > kxgm) {
				col = col[1:kxgm]
				#print (paste("The first", kxgm, "colors are used!"))
			}
			#if (random) {
				#print ("User defined colors are used!")
			#}
		}
		for (i in 1:kxgm) {
			mupli = mupls[[i]]
			muplui = muplus[[i]]
			if (surface == 'C') {
				musurf = mupli
				if (is.null(main)) {
					main = 'Constrained Warped-Plane Spline Surface'
				}
				zlim0 = c(min(mins), max(maxs))
			} else if (surface == 'U') {
				musurf = muplui
				if (is.null(main)) {
					main = 'Unconstrained Warped-Plane Spline Surface' 
				}				
				zlim0 = c(min(minsu), max(maxsu))
			}
			#par(mar = c(4, 2, 2, 2))
#print (sub)
#persp(k1, k2, musurf,  sub = sub)
			persp(k1, k2, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = col[i], cex.axis = .75, main = main, ticktype = ticktype,...)

			par(new = TRUE)
		}
		par(new = FALSE)
	}
}

#####
#wps#
#####
wps <- function(formula, family = gaussian(), data = NULL, weights = NULL, pnt = TRUE, pen = 0, cpar = 1.2)
{
  cl <- match.call()
  if (is.character(family)) 
     family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
     family <- family()
  if (is.null(family$family)) 
     stop("'family' not recognized!")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  if (family$family == "binomial") {
	if (class(y) == "factor") {
		y = ifelse(y == levels(y)[1], 0, 1)
	}
  }
  warp.delta <- NULL
  nums <- NULL; ks <- list(); sps <- NULL
  xmat <- NULL; x1 <- NULL; x2 <- NULL; xnms <- NULL; nks0 <- NULL; ks0 <- NULL; dcs <- NULL
  zmat <- NULL; zid <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_fac <- NULL; is_param <- NULL; vals <- NULL; st <- 1; ed <- 1
  for (i in 2: ncol(mf)) {
    if (is.character(attributes(mf[, i])$shape)) {
	sps <- attributes(mf[, i])$space
	dcs <- attributes(mf[, i])$decreasing
	ks0 <- attributes(mf[, i])$knots
	nks0 <- attributes(mf[, i])$numknots
	x1 <- (mf[, i])[, 1]
	x2 <- (mf[, i])[, 2]
	xmat <- cbind(x1, x2)
	xnms <- c(xnms, attributes(mf[, i])$name)
	ans_warp <- makedelta_wps(x1t = x1, x2t = x2, m1_0 = nks0[1], m2_0 = nks0[2], k1 = ks0$k1, k2 = ks0$k2, space = sps, decreasing = dcs)
	warp.delta0 <- ans_warp$delta
	k1 <- ans_warp$k1
	k2 <- ans_warp$k2
	ks[[1]] <- k1
	ks[[2]] <- k2
	warp.delta <- cbind(warp.delta, warp.delta0)
    }
    if (is.null(attributes(mf[, i])$shape)) {
	if (!is.null(names(mf)[i])) {
	  znms <- c(znms, names(mf)[i])
	}
        if (!is.matrix(mf[, i])) {
      		zid <- c(zid, i)
		is_param <- c(is_param, TRUE)
        	if (is.factor(mf[, i])) {
	  		is_fac <- c(is_fac, TRUE)
	 	 	ch_char <- suppressWarnings(is.na(as.numeric(levels(mf[, i]))))
          		if (any(ch_char)) {
	    			vals <- c(vals, unique(levels(mf[, i]))[-1])
          		} else {
	    			vals <- c(vals, as.numeric(levels(mf[, i]))[-1])
	  		}
          		nlvs <- length(attributes(mf[, i])$levels)
			ed <- st + nlvs - 2 
			zid1 <- c(zid1, st)
	  		zid2 <- c(zid2, ed)
	  		st <- st + nlvs - 1
	  		zmat0 <- model.matrix(~ mf[, i])[, -1, drop = FALSE]
	  		zmat <- cbind(zmat, zmat0)
    		} else {
			is_fac <- c(is_fac, FALSE)
			zmat <- cbind(zmat, mf[, i])
			ed <- st
			zid1 <- c(zid1, st)
			zid2 <- c(zid2, ed)	
			st <- st + 1		
			vals <- c(vals, "")
            	}
	} else {
	 	is_param <- c(is_param, FALSE)
		is_fac <- c(is_fac, FALSE)
	  	zmat0 <- mf[, i]
	  	zmat <- cbind(zmat, zmat0)
	  	#vals <- c(vals, 1)
	  	zid <- c(zid, i)
	  	nlvs <- ncol(zmat0) + 1
	  	ed <- st + nlvs - 2
	  	zid1 <- c(zid1, st)
	  	zid2 <- c(zid2, ed)
	  	st <- st + nlvs - 1
	}
    }
  }
  ans <- wps.fit(x1t = x1, x2t = x2, y = y, zmat = zmat, w = weights, pen = pen, pnt = pnt, cpar = cpar, decrs = dcs, delta = warp.delta, kts = ks)
  rslt <- list(k1 = ans$k1, k2 = ans$k2, muhat = ans$muhat, muhatu = ans$muhatu, SSE1 = ans$sse1, SSE0 = ans$sse0, edf = ans$edf, edfu = ans$edfu, delta = ans$delta, zmat = ans$zmat, xmat = xmat, coefs = ans$coefs, coefsu = ans$coefsu, zcoefs = ans$zcoefs, pvals.beta = ans$pz, se.beta = ans$sez, gcv = ans$gcv, gcvu = ans$gcvu, xnms = xnms, znms = znms, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2, ynm = ynm, decrs = dcs, tms = mt, is_param = is_param, is_fac = is_fac, d0 = ans$d0, pen = ans$pen, cpar = ans$cpar)
  rslt$call <- cl
  class(rslt) <- "wps"
  return (rslt) 
}

#####
wps.fit = function(x1t, x2t, y, zmat = NULL, w = NULL, pen = 0, pnt = TRUE, cpar = 1.2, decrs = c(FALSE, FALSE), delta = NULL, kts = NULL) {
	n = length(y)
	if (length(x1t) != n) {
		stop ("Error: length of x1 must equal length of y")
	}
	if (length(x2t) != n) {
		stop ("Error: length of x2 must equal length of y")
	}
# new
	k1 = kts[[1]]
	m1 = length(k1)
	k2 = kts[[2]]
	m2 = length(k2) 
	if (decrs[1]) {
		x1 = -x1t
		#if (!is.null(k1)) { 
			k1 = -k1[m1:1] 
		#}
	} else {x1 = x1t}
	if (decrs[2]) {
		x2 = -x2t
		#if (!is.null(k2)) {
			k2 = -k2[m2:1]
		#}
	} else {x2 = x2t}
	one = 1:n*0 + 1
	if (is.null(zmat)) {
		zmat = matrix(one, ncol = 1)
	} else {
		if (dim(zmat)[1] != n) {
			stop ("Error: # rows of zmat must equal length of y")
		}
		zproj = zmat %*% solve(crossprod(zmat), t(zmat))
		onep = one - zproj %*% one
# if the one vector is not in the space of zmat, then include the one vector in zmat 
		if (sum(onep^2) > 1e-12) {
			zmat = cbind(one, zmat)
		}
	}
	p = dim(zmat)[2]
	dimnames(zmat)[[2]] = NULL
# we already get warp.delta; the 1st col of xmat is one
	xmat0 = delta
#print (all(xmat0[, 1] == 1))
	if (dim(zmat)[2] > 1) {
		xmat = cbind(zmat, xmat0[, 2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))])
	} else {xmat = xmat0}
# constraint matrix
	amat = matrix(0, nrow = 2*m1*m2 - m1 - m2, ncol = m1*m2)
	amat[1, 2] = 1
	amat[2, m1 + 1] = 1
	irow = 2
	for (i2 in 2:m2) {
		irow = irow+1
		amat[irow, 2] = 1; amat[irow, m1 + m2 + i2 - 2] = 1
	}
	for (i1 in 2:m1) {
		irow = irow + 1
		amat[irow, m1 + 1] = 1; amat[irow, m1 + m2 + (i1 - 2) * (m2 - 1)] = 1
	}
	for (i1 in 2:(m1 - 1)) {
		irow = irow + 1
		amat[irow, i1] = -1; amat[irow, i1 + 1] = 1
	}
	for (i2 in 2:(m2 - 1)) {
		irow = irow + 1
		amat[irow, m1 + i2 - 1] = -1
		amat[irow, m1 + i2] = 1
	}
	for (i1 in 2:(m1 - 1)) {
		for (i2 in 2:m2) {
			irow = irow + 1
			amat[irow, i1] = -1
			amat[irow, i1 + 1] = 1
			amat[irow, m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2] = -1
			amat[irow, m1 + m2 - 2 + (i1 - 1) * (m2 - 1) + i2] = 1
		}
	}
	for (i2 in 2:(m2 - 1)) {
		for (i1 in 2:m1) {
			irow = irow + 1
			amat[irow, m1 + i2 - 1] = -1
			amat[irow, m1 + i2] = 1
			amat[irow, m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2] = -1
			amat[irow, m1 + m2 - 1 + (i1 - 2) * (m2 - 1) + i2] = 1
		}
	}
	amat0 = amat
	if (dim(zmat)[2] > 1) {
		amatz =  matrix(0, nrow = 2 * m1 * m2 - m1 - m2, ncol = dim(zmat)[2])
		amat = cbind(amatz, amat[, 2:(m1 * m2)])
	}
# penalty matrix
	kn1 = k1; kn2 = k2
	dmat = matrix(0, nrow = 2 * (m1 * m2 - m1 - m2), ncol = m1 * m2)
# by row
	irow = 1
	dmat[irow, 2] = -2; dmat[irow, 3] = 1
	for (i in 2:(m1 - 2)) {
		irow = irow + 1
		dmat[irow, i] = 1; dmat[irow, i + 1] = -2; dmat[irow, i + 2] = 1
	}
	for (ik2 in 2:m2) {
		irow = irow + 1
		dmat[irow, 2] = -2; dmat[irow, 3] = 1
		dmat[irow, m1 + m2 - 1 + ik2 - 1] = -2
		dmat[irow, m1 + 2 * m2 + ik2 - 3] = 1
		for (ik1 in 2:(m1 - 2)) {
			irow = irow + 1
			dmat[irow, ik1] = 1; dmat[irow, ik1 + 1] = -2; dmat[irow, ik1 + 2] = 1
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + ik2 - 1] = 1
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 1) + ik2 - 1] = -2
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * ik1 + ik2 - 1] = 1
		}
	}
# by col
	irow = irow + 1
	dmat[irow, m1 + 1] = -2; dmat[irow, m1 + 2] = 1
	for (i in 2:(m2 - 2)) {
		irow = irow + 1
		dmat[irow, m1 + i - 1] = 1; dmat[irow, m1 + i] = -2; dmat[irow, m1 + i + 1] = 1
	}
	for (ik1 in 2:m1) {
		irow = irow + 1
		dmat[irow, m1 + 1] = -2; dmat[irow, m1 + 2] = 1
		dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + 1] = -2
		dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + 2] = 1
		for (ik2 in 2:(m2 - 2)) {
			irow = irow + 1
			dmat[irow, m1 + ik2 - 1] = 1; dmat[irow, m1 + ik2] = -2; dmat[irow, m1 + ik2 + 1] = 1
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + ik2 - 1] = 1
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + ik2] = -2
			dmat[irow, m1 + m2 - 1 + (m2 - 1) * (ik1 - 2) + ik2 + 1] = 1
		}
	}
	dmat0 = dmat
	if (dim(zmat)[2] > 1) {
		dmatz = matrix(0, nrow = dim(dmat)[1], ncol = dim(zmat)[2])
		dmat = cbind(dmatz, dmat[, 2:(m1 * m2)])
	}
	k1 = kn1; k2 = kn2
#}
# weight
	if (is.null(w)) {
		xw0 = delta
		yw = y
		xw = xmat
		zw = zmat
	} else {
		if (min(w) > 1e-8) {
			yw = y * sqrt(w)
			xw = xmat
			zw = zmat
			xw0 = xmat0
			for (i in 1:n) {
				xw[i, ] = xmat[i, ] * sqrt(w[i])
				xw0[i, ] = xmat0[i, ] * sqrt(w[i])
			}
			for (i in 1:n) {
				zw[i, ] = zmat[i, ] * sqrt(w[i])
			}
		} else {
			xw0 = xmat0; xw = xmat; yw = y; zw = zmat
		}
	}
# transform to cone projection
sc = 1
if (round(pen, 6) > 1e-6) {
	ps = pen 
} else if (pnt & (round(pen, 6) == 0)) {
	mat = cbind(1, x1, x2, x1*x2) 
	mu_para = mat %*% solve(t(mat) %*% mat) %*% t(mat) %*% yw
	ssr = sum((yw - mu_para)^2)
	sc = ssr / (n-4)
	ps = 1 * sc
ps = max(1e-6, 1 * sc)
} else {ps = 1e-6} 
#	ps = max(1e-6, pen)
	qmat = t(xw) %*% xw + ps * t(dmat) %*% dmat	
	qmat0 = t(xw0) %*% xw0 + ps * t(dmat0) %*% dmat0
	umat = chol(qmat)
	#iumat = diag(ncol(umat))
        #uinv = backsolve(umat, iumat)
	#umat0 = chol(qmat0)
	uinv = solve(umat)
	#uinv0 = solve(umat0)
	atil = amat %*% uinv
	#atil0 = amat0 %*% uinv0
	cvec = t(uinv) %*% t(xw) %*% yw
	ans = coneA(cvec, atil)
	phihat = ans$thetahat
	ahat = uinv %*% phihat
	#cvec = crossprod(xw, yw)
	#cvec = t(xw) %*% yw	
#bv = rep(0, nrow(amat))
	#ans = qprog(qmat, cvec, amat, 1:nrow(amat)*0, msg = FALSE)
	#ahat = ans$thetahat
	muhat = xw %*% ahat	
# get trace of "proj" matrix
	sm = 1e-10
	nz = 1:dim(amat)[1] < sm	
	nz[amat %*% ahat > sm] = TRUE
	if (sum(nz) < dim(amat0)[1]) {
		gamj = -t(amat0[!nz, ])
		if (length(gamj) == dim(xw0)[2]) {
			gamj = matrix(gamj, ncol = 1)
		}
		aq = qr(gamj)
		if (aq$rank >= 1) {
			if (aq$rank == 1) {
				gamj = matrix(qr.Q(aq)[, 1:aq$rank, drop = FALSE], ncol = 1)
			} else {
				gamj = qr.Q(aq)[, 1:aq$rank, drop = FALSE]
			}	
			#pa = gamj %*% solve(crossprod(gamj), t(gamj))		
			pa = gamj %*% solve(t(gamj) %*% gamj) %*% t(gamj)	
			dj1 = dim(gamj)[1]; dj2 = dim(gamj)[2]
			imat = matrix(0, nrow = dj1, ncol = dj1) 
			for (i in 1:dj1) {imat[i, i] = 1}
			uvecs = matrix(rnorm(dj1 * (dj1 - dj2)), nrow = dj1)
			wmatt = (imat - pa) %*% uvecs
			aqt = qr(wmatt); wmat = qr.Q(aqt)
			b0mat = xw0 %*% wmat
		}
	} else {
		b0mat = xw0
		#inmat = matrix(0, nrow = dim(xmat0)[2], ncol = dim(xmat0)[2])
		#for (i in 1:dim(xmat0)[2]) {inmat[i, i] = 1}
		inmat = diag(dim(xmat0)[2])
		wmat = inmat		
	}
	p0mat = xmat0 %*% wmat %*% solve(t(wmat) %*% qmat0 %*% wmat) %*% t(wmat) %*% t(xmat0)
	#p0mat = xmat0 %*% wmat %*% solve(crossprod(wmat, qmat0 %*% wmat), t(wmat) %*% t(xmat0))	
	edf = sum(diag(p0mat)) + p - 1
	inmat = matrix(0, nrow = n, ncol = n)
	diag(inmat) = 1	
	sse1 = sum((yw - xw %*% ahat)^2)
#new:
zcoefs = ahat[1:p]
sse0 = sum((yw - zw %*% zcoefs)^2)
#new use edf instead if (n - cpar * edf) < 0
if ((n - p - cpar * edf) <= 0) {
	sig2hat = sse1 / edf
} else {
	sig2hat = sse1 / (n - p - cpar * edf)
}
	if (p >= 1) {
		pm = one %*% solve(crossprod(one)) %*% t(one)
		covmat = solve(t(zmat) %*% (inmat - p0mat + pm) %*% zmat)	
		sez = sqrt(diag(covmat) * sig2hat)	
		tz = zcoefs / sez
#new use edf instead if (n - cpar * edf) < 0
		if ((n - p - cpar * edf) <= 0) {
			pz = 2 * (1 - pt(abs(tz), edf))
			if (p > 1) {
				warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
			}
#print ('Check pz!')
		} else {
			pz = 2 * (1 - pt(abs(tz), n - p - cpar * edf))
		}
	}
# get unconstrained penalized estimator
	#prmatu = xw %*% solve(qmat, t(xw))
	prmatu = xw %*% solve(qmat) %*% t(xw)
	muhatu = prmatu %*% yw
	sseu = sum((yw - muhatu)^2)
	edfu = sum(diag(prmatu))
	ans = new.env()
	ans$k1 = k1
	ans$k2 = k2
	ans$muhat = muhat
	#ans$muplot = mupl
	ans$muhatu = muhatu
	#ans$muplotu = muplu
	ans$gcv = sse1 / (1 - edf / n)^2
	ans$gcvu = sseu / (1 - edfu/n)^2
	#ans$ssr = sse1
	ans$sse1 = sse1
	ans$sse0 = sse0
	ans$edf = edf
	#ans$nz = amat %*% ahat
	ans$edfu = edfu
	ans$coefs = ahat
	ans$coefsu = solve(qmat) %*% t(xw) %*% yw
#include coef for one vector
	#ans$zcoefs = ahat[1:p]
	ans$zcoefs = zcoefs
	ans$zmat = zmat 	
	ans$sighat = sig2hat
	ans$delta = xmat
	ans$pen = ps	
	ans$d0 = p
	ans$cpar = cpar 
	#ans$gcvus = gcvus
	#ans$lambdas_pen = lambdas_pen
	#if (p >= 2) {
	if (p >= 1) {
		ans$sez = sez
		ans$pz = pz
		ans$covmat = covmat
		ans$tz = tz
	} else {ans$sez = NULL; ans$pz = NULL; ans$covmat = NULL; ans$tz = NULL}
#print (sig2hat)
	return (ans)
}

####################################################################
#four monotonicity functions for warped-plane fit
####################################################################
ii <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "wps_ii"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "decreasing") <- c(FALSE, FALSE)
    return (xm)
}


di <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "wps_di"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "decreasing") <- c(TRUE, FALSE)
    return (xm)
}


dd <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "wps_dd"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "decreasing") <- c(TRUE, TRUE)
    return (xm)
}


###############################################################
#makedelta_wps: make delta to a pair of warped-plane variables#
###############################################################
makedelta_wps <- function(x1t, x2t, m1_0 = 0, m2_0 = 0, k1 = 0, k2 = 0, space = c("E", "E"), decreasing = c(FALSE, FALSE))
{
# x1 and x2 no need to sort
# if decreasing  (all calculations done for doubly-increasing case)
	n = length(x1t)
	if (decreasing[1]) {
#print (k1 != 0)
		x1 = -x1t
		#if (!is.null(k1)) {
		if (k1 != 0) { 
			m1 = length(k1); k1 = -k1[m1:1] 
		}
	} else {x1 = x1t}
	if (decreasing[2]) {
		x2 = -x2t
		#if (!is.null(k2)) {
		if (k2 != 0) {
			m2 = length(k2); k2 = -k2[m2:1]
		}
	} else {x2 = x2t}
# determine whether to use default knots or user-defined knots
# check if k1 and k2 > 1 elems
	make1 = FALSE   
	#if (is.null(k1)) {
	#if (k1 == 0) {
	if (all(k1 == 0)) {
		make1 = TRUE
	} else if (min(k1) > min(x1) | max(k1) < max(x1)) {
		warning ('Predictor should be within the range of the knots! User-defined knots not used!')
		make1 = TRUE
# new:	at least 4 kts?
	}# else if (length(k1) < 4) {
	#	make1 = TRUE
	#} else if (m1_0 >= 4) {
	#	make1 = TRUE
	#}
	make2 = FALSE
	#if (is.null(k2)) {
	#if (k2 == 0) {
	if (all(k2 == 0)) {
		make2 = TRUE
	} else if (min(k2) > min(x2) | max(k2) < max(x2)) {
		warning ('Predictor should be within the range of the knots! User-defined knots not used!')
		make2 = TRUE
	}# else if (length(k2) < 4) {
	#	make2 = TRUE
	#} else if (m2_0 >= 4) {
	#	make2 = TRUE
	#}
# add quantile part in it
	if (make1) {
# new:
#print (m1_0 == 12L)? !is.integer(m1_0) not work
#print (m1_0)
		if (m1_0  < 4 | round(m1_0, 0) != m1_0) {
		#if (m1_0  < 4) {
			if (m1_0 != 0) {
				warning ('At least four knots should be used! Number of knots is re-defined!')
			}
			m1 = 2 * round(n^(1/6)) + 4 
		} else {m1 = m1_0}
		if (space[1] == "Q") {
			k1 = quantile(unique(x1), probs = seq(0, 1, length = m1), names = FALSE)
		}
		if (space[1] == "E") {
			k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) + min(x1)
		} 
		#k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) + min(x1)
	} else { m1 = length(k1) }
	if (make2) {
#new:
		if (m2_0 < 4 | round(m2_0, 0) != m2_0) {
			if (m2_0 != 0) {
				warning ('At least four knots should be used! Number of knots is re-defined!')
			}
			m2 = 2 * round(n^(1/6)) + 4 
		} else {m2 = m2_0}
		if (space[2] == "Q") {
			k2 = quantile(unique(x2), probs = seq(0, 1, length = m2), names = FALSE)
		}
		if (space[2] == "E") {
			k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2) - min(x2)) + min(x2)
		} 
		#k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2) - min(x2)) + min(x2)	
	} else { m2 = length(k2) }
## check to see if empty knot intervals
	rm1 = k1[m1]; rm2 = k2[m2]
	keep = 1:m1 > 0
	for (i in 2:m1) {
		if (sum(x1 >= k1[i-1] & x1 < k1[i]) == 0) {
			keep[i] = FALSE
		}
	}
	k1 = k1[keep]; m1 = length(k1); k1[m1] = rm1
	keep = 1:m2 > 0
	for (i in 2:m2) {
		if (sum(x2 >= k2[i-1] & x2 < k2[i]) == 0) {
			keep[i] = FALSE
		}
	}
	k2 = k2[keep]; m2 = length(k2); k2[m2] = rm2
# make the basis functions
	b1 = matrix(0, nrow = n, ncol = m1)
	b2 = matrix(0, nrow = n, ncol = m2)
	for (i in 2:(m1 - 1)) {
		i1 = x1 >= k1[i-1] & x1 <= k1[i]
		b1[i1, i] = (x1[i1] - k1[i-1]) / (k1[i] - k1[i-1])
		i2 = x1 > k1[i] & x1 <= k1[i+1]
		b1[i2, i] = (k1[i+1] - x1[i2]) / (k1[i+1] - k1[i])
	}
	i1 = x1 >= k1[1] & x1 <= k1[2]
	b1[i1, 1] = (k1[2] - x1[i1]) / (k1[2] - k1[1])
	i2 = x1 > k1[m1-1] & x1 <= k1[m1]
	b1[i2, m1] = (x1[i2] - k1[m1-1]) / (k1[m1] - k1[m1-1])

	for (i in 2:(m2 - 1)) {
		i1 = x2 >= k2[i-1] & x2 <= k2[i]
		b2[i1, i] = (x2[i1] - k2[i-1]) / (k2[i] - k2[i-1])
		i2 = x2 > k2[i] & x2 <= k2[i+1]
		b2[i2, i] = (k2[i+1] - x2[i2]) / (k2[i+1] - k2[i])
	}
	i1 = x2 >= k2[1] & x2 <= k2[2]
	b2[i1, 1] = (k2[2] - x2[i1]) / (k2[2] - k2[1])
	i2 = x2 > k2[m2-1] & x2 <= k2[m2]
	b2[i2, m2] = (x2[i2] - k2[m2-1]) / (k2[m2] - k2[m2-1])
## design matrix
	xmat0 = matrix(nrow = n, ncol = m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
	xmat0[ ,1] = 1:n*0 + 1
	xmat0[ ,2:m1] = b1[ ,2:m1]
	xmat0[ ,(m1 + 1):(m1 + m2 - 1)] = b2[ ,2:m2]
	for (i in 1:(m1 - 1)) {
		xmat0[ ,(m1 + m2 + (i - 1) * (m2 - 1)):(m1 + m2 - 1 + i * (m2 - 1))] = b1[ ,i + 1] * b2[ ,2:m2]
	}
	#if (dim(zmat)[2] > 1) {
	#	xmat = cbind(zmat, xmat0[ ,2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))])
	#} else {xmat = xmat0}
	#bmat = xmat0[ ,2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))]
# ignore the zmat for now
	xmat = xmat0
# columns of delta are edges, different from other make_delta in cgam	
	delta = xmat
	ans = new.env()
	ans$delta = delta
	ans$k1 = k1
	ans$k2 = k2
	return (ans)	
	#attr(delta, "shape") = "warp"
	#delta
}


#################
#new coef method#
#################
coef.cgam <- function(object,...) {
  ans <- object$coefs
  ans	
}


###################
#new fitted method#
###################
fitted.cgam <- function(object,...) {
  ans <- object$muhat
  ans
}

fitted.wps <- function(object,...) {
  ans <- object$muhat
  ans
}

#############################
#shape selection part
#############################
#get(x = "s", pos = "package:cgam")
ShapeSelect <- function(formula, family = gaussian(), cpar = 1.2, data = NULL, weights = NULL, genetic = FALSE) {
	if (exists("s", parent.frame()) & class(get("s", envir = parent.frame())) == "function") {
		if (!identical(get("s", envir = parent.frame()), cgam::s)) {
			assign("s", cgam::s, envir = parent.frame())
		}
	}
	cl <- match.call()
	if (is.character(family)) 
		family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family)) 
		family <- family()
	if (is.null(family$family)) 
		stop("'family' not recognized!")
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
  	mf[[1L]] <- as.name("model.frame")
  	mf <- eval(mf, parent.frame())
  	ynm <- names(mf)[1]
  	mt <- attr(mf, "terms")
  	y <- model.response(mf, "any")
  	if (family$family == "binomial") {
		if (class(y) == "factor") {
			y <- ifelse(y == levels(y)[1], 0, 1)
		}
#new: test
		if (class(y) == "character") {
			y <- ifelse(factor(y) == levels(factor(y))[1], 0, 1)
		}
  	}
#print (head(y))
  	shpsx <- list(); shpsz <- list(); shpst <- list(); shpsvx <- NULL
	ix <- 1; iz <- 1; itr <- 1; sel <- FALSE
  	xmat <- NULL; xnms <- NULL
	vxmat <- NULL; vxnms <- NULL
  	zmat <- NULL; znms <- NULL; zfacs <- NULL
	vzmat <- NULL; vznms <- NULL; vzfacs <- NULL
	trmat <- NULL; trnms <- NULL; vtrmat <- NULL; vtrnms <- NULL
	for (i in 2:ncol(mf)) {
		if (!is.null(attributes(mf[, i])$type)) {
			if (attributes(mf[, i])$type == "fac" | attributes(mf[, i])$type == "lin") {
				if (attributes(mf[, i])$type == "fac") {
					zfacs <- c(zfacs, TRUE)
				} else {
					#if (is.character(mf[,i])) {
						#mf[,i] = ifelse(factor(mf[,i]) == levels(factor(mf[,i]))[1], 0, 1)
					#	nm <- attributes(mf[,i])$nm
					#	mf[,i] <- ifelse(factor(mf[,i]) == levels(factor(mf[,i]))[1], 0, 1)
					#	attr(mf[,i], "type") <- "lin"	
					#	attr(mf[,i], "shape") <- c(0, 1)
					#	attr(mf[,i], "nm") <- nm
					#	zfacs <- c(zfacs, TRUE)
					#} else {
						zfacs <- c(zfacs, FALSE)
					#}
				}
				zmat <- cbind(zmat, mf[, i])
       				znms <- c(znms, attributes(mf[, i])$nm)
				shpsz[[iz]] <- attributes(mf[, i])$shape
				iz <- iz + 1
				sel <- TRUE		
			}
			#if (attributes(mf[, i])$type == "lin") {
			#	zfacs <- c(zfacs, FALSE)
			#	zmat <- cbind(zmat, mf[, i])
#temp
       				#znms <- c(znms, names(mf)[i])
			#	znms <- c(znms, attributes(mf[, i])$nm)
			#	shpsz[[iz]] <- c(0, 1)
			#	iz <- iz + 1
			#}
			if (attributes(mf[, i])$type == "nparam") {
  				xmat <- cbind(xmat, mf[, i])
       				xnms <- c(xnms, attributes(mf[, i])$nm)
				shpsx[[ix]] <- attributes(mf[, i])$shape
				ix <- ix + 1
				sel <- TRUE
			}
			#if (attributes(mf[, i])$type == "ord.tree") {
			if (attributes(mf[, i])$type == "tree") {
				trmat <- cbind(trmat, mf[, i])
       				trnms <- c(trnms, attributes(mf[, i])$nm)
				shpst[[itr]] <- attributes(mf[, i])$shape
				itr <- itr + 1		
				sel <- TRUE
			}
		} else {
			if (is.numeric(attributes(mf[, i])$shape)) {
				shpsvx <- c(shpsvx, attributes(mf[, i])$shape)
				vxmat <- cbind(vxmat, mf[, i])
        			vxnms <- c(vxnms, attributes(mf[, i])$nm)
			} else if (is.character(attributes(mf[, i])$shape)) {
				if (attributes(mf[, i])$shape == "tree") {
					vtrmat <- cbind(vtrmat, mf[, i])
					vtrnms <- c(vtrnms, attributes(mf[, i])$nm)
				}
			} else {
				if (is.factor(mf[, i])) {
					vzfacs <- c(vzfacs, TRUE)
				} else {
					vzfacs <- c(vzfacs, FALSE)
				}
				vzmat <- cbind(vzmat, mf[, i])
#temp
       				vznms <- c(vznms, names(mf)[i])
			}
		}
	}
	if (!sel) {
		stop ("No variable to be selected! Use cgam instead!")
	}
	if (!is.null(xmat)) {
		capl <- ncol(xmat)
		bx <- sapply(shpsx, function(x) length(x))
		xnum <- rev(cumprod(bx))[1]
	} else {capl <- 0; xnum <- 1}
	if (!is.null(zmat)) {
		capk <- ncol(zmat)
		znum <- 2^capk
	} else {capk <- 0; znum <- 1}
	if (!is.null(trmat)) {
		capt <- ncol(trmat)
		trnum <- 3^capt
	} else {capt <- 0; trnum <- 1}
	nmod <- xnum * znum * trnum
	tm <- 0
	if (genetic) {
# call GA
		print ("Go genetic algorithm!")
		ptm <- proc.time()
		ans <- ConstrGA(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, nmod = nmod, zfacs = zfacs, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat)
		tm <- proc.time() - ptm
	} else {
		if (nmod <= 5e+2) {
			print ("Go through models one by one!")	
			Sys.sleep(1)
			ptm <- proc.time()
			ans <- ConstrALL(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, zfacs = zfacs, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat)
			tm <- tm + proc.time() - ptm
#could go one-by-one, but will suggest genetic
		} else if (nmod > 5e+2 & nmod <= 1e+6) {
			print ("Estimating the time of one fit!")
			ptm <- proc.time()	
			invisible(ConstrGA(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, nmod = nmod, zfacs = zfacs, time.est = TRUE, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat))
			tm <- proc.time() - ptm
			print (paste("Total time of one fit is roughly", round(tm[3], 1), "seconds"))
			if (round(tm[3] * nmod/60^2) > 1) {
				tm2 <- paste(round(tm[3] * nmod/60^2), "hours")
			} else {tm2 <- paste(round(tm[3] * nmod/60), "minutes")}
			msg <- paste("Number of all models is", nmod, "and the time for all fits is roughly", tm2, ". Go one by one?")
			res <- dlgMessage(msg, "yesno")$res
			if (res == "yes") {
				print ("Go one by one!")
				ptm <- proc.time()
				ans <- ConstrALL(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, zfacs = zfacs, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat) 
				tm <- proc.time() - ptm
			} else {
				print ("Go genetic!")
				ptm <- proc.time()
				ans <- ConstrGA(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, nmod = nmod, zfacs = zfacs, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat)
				tm <- proc.time() - ptm
			}
#too many models or memory problem
		} else {
			print ("Go genetic! Too many models!")	
			ptm <- proc.time()
			ans <- ConstrGA(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, nmod = nmod, zfacs = zfacs, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat)
			tm <- tm + proc.time() - ptm
		}
	}
#print (vznms)
	colnames(ans$pop2)[1:(capl+capk+capt)] <- c(xnms, znms, trnms)
#new:
top <- ans$pop2[1,]
if (!is.null(vxnms)) {
	lvx <- length(vxnms)
	nm0 <- colnames(top)
	for (i in 1:lvx) {
		#fi <- vxnms[i]
		pvi <- ShapeToChar(shpsvx[i], tag = "x")
		top <- cbind(top, pvi)
		#fi <- paste(pvi, "(" , vxnms[i], ")" , sep = "")
		#fm0 <- paste(fm0, fi, sep = "+")
	}
	colnames(top) <- c(nm0, vxnms)
}
if (!is.null(vtrnms)) {
	lvt <- length(vtrnms)
	nm0 <- colnames(top)
	for (i in 1:lvt) {
		#fi <- paste("tree", "(" , vtrnms[i], ")" , sep = "")
		#fm0 <- paste(fm0, fi, sep = "+")
		pvi <- "tree"
		top <- cbind(top, pvi)
	}
	colnames(top) <- c(nm0, vtrnms)
}
	rslt <- list(pop = ans$pop2, top = top, fitness = ans$fitness, tm = tm, xnms = xnms, znms = znms, trnms = trnms, zfacs = zfacs, mxf = ans$mxf, mnf = ans$mnf, GA = ans$GA, vzcat = ans$vzcat, vzmat = vzmat, shpsx = shpsx, vxnms = vxnms, vznms = vznms, vtrnms = vtrnms)
	form <- make_form(pvec = ans$pop2[1, 1:(capl+capk+capt)], pnms = colnames(ans$pop2)[1:(capl+capk+capt)], ynm = ynm, zfacs = zfacs, vznms = vznms, vxnms = vxnms, shpsvx = shpsvx, vtrnms = vtrnms)
	fm <- form$fm
#print (fm)
	zps <- form$zps
	if (is.null(fm)) {
		print ('No predictor is chosen')
		best.fit <- NULL	
	} else {
		best.fit <- cgam(formula = fm, nsim = 0, family = family, cpar = cpar, data = data, weights = weights)
		rslt$best.fit <- best.fit 
		rslt$fm <- fm
		id_noflat <- which(ans$pop[1, 1:capl, drop = FALSE] != 0)
		if (length(id_noflat) > 1) {
			msg <- paste("A perspective plot of the best fit?")
			res <- dlgMessage(msg, "yesno")$res
			if (res == "yes") {
				if (is.null(zps)) {
					plotpersp(best.fit)
				} else {plotpersp(best.fit, categ = zps[1])}
			} 
		}
	}
        rslt$best.fit <- best.fit
	rslt$call <- cl
	class(rslt) <- "shapeselect"
	return (rslt) 
}

######################
#extract the best fit#
######################
best.fit <- function(x) {
	if (!inherits(x, "shapeselect")) { 
	        stop("best.fit only works for an object of the ShapeSelect routine!")
        }
	object <- x$best.fit
	return (object)
}

##########################################
#new symbolic routine                   #
#s(): shape-restricted for genetic only #
#keep all 17 one shape routines for cgam#
#########################################
shapes <- function(x, set = "s.9")  
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    if (is.numeric(set)) {
	if (min(set) >= 0 & max(set) <= 16) {
		#shps <- c(0, set)
		shps <- set 
		shps <- unique(sort(shps))
		#attr(x, "char") <- NULL
	} else {stop ("Shape values for shape-restricted variables can only be between 0 and 16!")}
	attr(x, "type") <- "nparam"
    } else {
	attr(x, "type") <- "nparam"
	if (identical(set, "s.5")) {
		shps <- c(0, 9:12)
    	} else if (identical(set, "s.9")) {
		shps <- c(0, 9:16)
    	} else if (identical(set, "ord.5")) {
		shps <- c(0, 1:4)
    	} else if (identical(set, "ord.9")) {
		shps <- c(0, 1:8)
#print (shps)
    	} else if (identical(set, "tree")) {
		#shps <- c(0, "tree", "unordered")
		shps <- c(0, 1, 2)
		attr(x, "type") <- "tree"
    	} else {
		lshps <- length(set)
		shps <- unique(sort(unname(sapply(set, CharToShape))))
	}
    }
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- shps
    x
}


######################
#choose z as factors
######################
in.or.out <- function(z)  
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(z, "nm") <- deparse(pars$z)
if (is.factor(z)) {
    attr(z, "type") <- "fac"
} else {attr(z, "type") <- "lin"}
    #shps <- c(0, 1)
    attr(z, "shape") <- c(0, 1)
    z
}


#####################
#make a cgam formula#
#####################
make_form <- function(pvec = NULL, pnms = NULL, ynm = NULL, zfacs = NULL, vznms = NULL, vxnms = NULL, shpsvx = NULL, vtrnms = NULL) {
	fm <- NULL	
	zps <- NULL
	zi <- 1
	vs <- NULL
	rm <- c("flat", "out")
	pvec <- apply(as.matrix(pvec, nrow = 1), 2, as.character)
	if (all(pvec %in% rm)) {
		rslt <- list(fm = fm, zps = zps)
		#return (rslt)
	} else {
		if (any(pvec %in% rm)) {
			id_rm <- which(pvec %in% rm)
			pvec <- pvec[-id_rm]
			pnms <- pnms[-id_rm]
		}
		if (as.character(pvec[1]) == "in") {
			fm0 <- pnms[1]
			zps <- c(zps, fm0)
			zi <- zi + 1
#vs <- c(vs, pnms[1])
		} else if (as.character(pvec[1]) == "unordered") {
			fm0 <- paste("factor(", pnms[1], ")", sep = "")				
			zps <- c(zps, fm0)
#vs <- c(vs, pnms[1])
		} else {
			fm0 <- paste(as.character(pvec[1]), "(" , pnms[1], ")" , sep = "")
		}
		vs <- c(vs, pnms[1])
		if (length(pvec) > 1) {
#zi only works for zfacs, not tree unordered case
			#zi <- 1
			for (i in 2:length(pvec)) {
				#if (as.character(pvec[i]) == "flat") {
#not include the x
				#	next
				#} else 
				if (as.character(pvec[i]) == "in") {
					#if (zfacs[zi]) {				
					#	fi <- paste("factor(", pnms[i], ")", sep = "")
					#} else {fi <- pnms[i]}
					fi <- pnms[i]
					zps <- c(zps, fi)	
					zi <- zi + 1
				} else if (as.character(pvec[i]) == "unordered") {			
					fi <- paste("factor(", pnms[i], ")", sep = "")
					#fi <- pnms[i]
					zps <- c(zps, fi)
				} else {
#should work with tree?
					fi <- paste(as.character(pvec[i]), "(" , pnms[i], ")" , sep = "")
				}	
				fm0 <- paste(fm0, fi, sep = "+")	
			}
		}
		if (!is.null(vxnms)) {
			lvx <- length(vxnms)
			for (i in 1:lvx) {
				#fi <- vxnms[i]
				pvi <- ShapeToChar(shpsvx[i], tag = "x")
				fi <- paste(pvi, "(" , vxnms[i], ")" , sep = "")
				fm0 <- paste(fm0, fi, sep = "+")
			}
		}
		if (!is.null(vznms)) {
			lvz <- length(vznms)
			for (i in 1:lvz) {
				fi <- vznms[i]
				zps <- c(zps, fi)
				fm0 <- paste(fm0, fi, sep = "+")
			}
		}
		if (!is.null(vtrnms)) {
			lvt <- length(vtrnms)
			for (i in 1:lvt) {
				fi <- paste("tree", "(" , vtrnms[i], ")" , sep = "")
				fm0 <- paste(fm0, fi, sep = "+")
			}
		}
		fm <- as.formula(paste(ynm, "~", fm0, sep = "")) 
		rslt <- list(fm = fm, zps = zps)
	}
	return (rslt)
}


###################
#genetic algorithm#
###################
ConstrGA = function(y, xmat, zmat, trmat, family = gaussian(), shpsx = NULL, shpsvx = NULL, shpsz = NULL, shpst = NULL, cpar = 1.2, nmod = 2e+4, zfacs = NULL, time.est = FALSE, weights = NULL, vzmat = NULL, vzfacs = NULL, vxmat = NULL, vtrmat = NULL) {
	linkfun = family$linkfun
	cicfamily = CicFamily(family)
	llh.fun = cicfamily$llh.fun
	etahat.fun = cicfamily$etahat.fun
	gr.fun = cicfamily$gr.fun
	wt.fun = cicfamily$wt.fun
	zvec.fun = cicfamily$zvec.fun
	muhat.fun = cicfamily$muhat.fun
	ysim.fun = cicfamily$ysim.fun
	deriv.fun = cicfamily$deriv.fun
	dev.fun = cicfamily$dev.fun 
#print (head(xmat))
#print (head(zmat))
#print (head(y))
	n = length(y)
	#sm = 1e-7 
	capl = length(xmat) / n
	if (capl < 1) {capl = 0}
	if (round(capl, 8) != round(capl, 1)) {
		stop ("Incompatible dimensions for xmat!")
	}
#check!
	if (capl > 0) {
		for(i in 1:capl) {
			xmat[,i] = (xmat[,i] - min(xmat[,i])) / (max(xmat[,i]) - min(xmat[,i]))
			#xmat[,i] = (xmat[,i] - mean(xmat[,i])) / sd(xmat[,i])
			#xmat[,i] = xmat[,i] / sd(xmat[,i])
		}
	}
	caplv = length(vxmat) / n
	if (caplv < 1) {caplv = 0}
	if (round(caplv, 8) != round(caplv, 1)) {
		stop ("Incompatible dimensions for xmat!")
	}
#new:
	if (caplv > 0) {
		for(i in 1:caplv) {
			vxmat[,i] = (vxmat[,i] - min(vxmat[,i])) / (max(vxmat[,i]) - min(vxmat[,i]))
			#vxmat[,i] = (vxmat[,i] - mean(vxmat[,i])) / sd(vxmat[,i])
			#vxmat[,i] = vxmat[,i] / sd(vxmat[,i])
		}
	}
	capk = length(zmat) / n
	if (capk < 1) {capk = 0}
	if (round(capk, 8) != round(capk, 1)) {
		stop ("Incompatible dimensions for zmat!")
	}
	capkv = length(vzmat) / n
	if (capkv < 1) {capkv = 0}
	if (round(capkv, 8) != round(capkv, 1)) {
		stop ("Incompatible dimensions for vzmat!")
	}
	capt = length(trmat) / n
	if (capt < 1) {capt = 0}
	if (round(capt, 8) != round(capt, 1)) {
		stop ("Incompatible dimensions for trmat!")
	}
	captv = length(vtrmat) / n
	if (captv < 1) {captv = 0}
	if (round(captv, 8) != round(captv, 1)) {
		stop ("Incompatible dimensions for trmat!")
	}
#nrep=0
################################################################
##get basis functions for all allowed shapes for each component#
#not consider allowed shapes for now
################################################################
# get basis functions for the constrained components -- ordinal monotone
if (capl > 0) {
	delta = varlst = NULL
	if (1 %in% shpsx[[1]] | 2 %in% shpsx[[1]]) {
#print ('1')
		del1 = makedelta(xmat[, 1], 1)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
if (capl > 1) {
	for (i in 2:capl) {
		if (1 %in% shpsx[[i]] | 2 %in% shpsx[[i]]) {
#print ('1')
			del2 = makedelta(xmat[, i], 1)$amat
	 		m2 = length(del2) / n
		} else {del2 = NULL; m2 = 0}
	 	delta = rbind(del1, del2)
		if (m1 > 0 | m2 > 0) {
			varlist = 1:(m1+m2)*0
			if (m1 > 0) {		
				varlist[1:m1] = var1
			}
			if (m2 > 0) {
				varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
			}
#var1: keep track of x's in delta: 11112222...
			var1 = varlist
			m1 = m1 + m2
	 		del1 = delta
		}
	}
	#delta_om = delta
	#varlist_om = varlist
} #else {delta_om = del1; varlist_om = var1}
	delta_om = delta
	varlist_om = varlist
}
#print (shpsx[[1]])
# get basis functions for the constrained components -- smooth monotone
if (capl > 0) {
	delta = varlst = NULL
	if (9 %in% shpsx[[1]] | 10 %in% shpsx[[1]]) {
#print ('9')
		del1 = makedelta(xmat[, 1], 9)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
if (capl > 1) {
	for (i in 2:capl) {
		if (9 %in% shpsx[[i]] | 10 %in% shpsx[[i]]) {
#print ('9')
			del2 = makedelta(xmat[, i], 9)$amat
	 		m2 = length(del2) / n
		} else {del2 = NULL; m2 = 0}
	 	delta = rbind(del1, del2)
		if (m1 > 0 | m2 > 0) {
			varlist = 1:(m1+m2)*0
			if (m1 > 0) {	
				varlist[1:m1] = var1
			}
			if (m2 > 0) {	
				varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
			}			
			var1 = varlist
			m1 = m1 + m2
	 		del1 = delta
		}
	}
	#delta_sm = delta
	#varlist_sm = varlist
} #else {delta_sm = del1; varlist_sm = var1}
	delta_sm = delta
	varlist_sm = varlist
}
#print (delta_sm)
#print (varlist_sm)
# get basis functions for the constrained components -- ordinal convex
if (capl > 0) {
	delta = varlst = NULL
	if (3 %in% shpsx[[1]] | 4 %in% shpsx[[1]]) {
#print ('3')
		del1 = makedelta(xmat[, 1], 3)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
if (capl > 1) {
	for (i in 2:capl) {
		if (3 %in% shpsx[[i]] | 4 %in% shpsx[[i]]) {
#print ('3')
			del2 = makedelta(xmat[, i], 3)$amat
	 		m2 = length(del2) / n
		}
	 	delta = rbind(del1, del2)
		if (m1 > 0 | m2 > 0) {
			varlist = 1:(m1+m2)*0
			if (m1 > 0) {	
				varlist[1:m1] = var1
			}
			if (m2 > 0) {	
				varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
			}
			var1 = varlist
			m1 = m1 + m2
	 		del1 = delta
		}
	}
	#delta_ocv = delta
	#varlist_ocv = varlist
} #else {delta_ocv = del1; varlist_ocv = var1}
	delta_ocv = delta
	varlist_ocv = varlist
}
if (capl > 0) {
	delta = varlst = NULL
# get basis functions for the constrained components -- smooth convex
	if (11 %in% shpsx[[1]] | 12 %in% shpsx[[1]]) {
#print ('11')
		del1 = makedelta(xmat[, 1], 11)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
if (capl > 1) {
	for (i in 2:capl) {
		if (11 %in% shpsx[[i]] | 12 %in% shpsx[[i]]) {
#print ('11')
			del2 = makedelta(xmat[, i], 11)$amat
	 		m2 = length(del2) / n
		}
	 	delta = rbind(del1, del2)
		if (m1 > 0 | m2 > 0) {
			varlist = 1:(m1+m2)*0
			if (m1 > 0) {	
				varlist[1:m1] = var1
			}
			if (m2 > 0) {	
				varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
			}
			var1 = varlist
			m1 = m1 + m2
	 		del1 = delta
		}
	}
	#delta_scv = delta
	#varlist_scv = varlist
} #else {delta_scv = del1; varlist_scv = var1}
	delta_scv = delta
	varlist_scv = varlist
}
# get basis functions for the constrained components -- ordinal increasing convex
if (capl > 0) {
	delta = varlst = NULL
	if (5 %in% shpsx[[1]] | 8 %in% shpsx[[1]]) {
#print ('5')
		del1 = makedelta(xmat[, 1], 5)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
if (capl > 1) {
	for (i in 2:capl) {
		if (5 %in% shpsx[[i]] | 8 %in% shpsx[[i]]) {
#print ('5')
			del2 = makedelta(xmat[, i], 5)$amat
	 		m2 = length(del2) / n
		}
	 	delta = rbind(del1, del2)
		if (m1 > 0 | m2 > 0) {
			varlist = 1:(m1+m2)*0
			if (m1 > 0) {	
				varlist[1:m1] = var1
			}
			if (m2 > 0) {	
				varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
			}
			var1 = varlist
			m1 = m1 + m2
	 		del1 = delta
		}
	}
	#delta_oincv = delta
	#varlist_oincv = varlist
} #else {delta_oincv = del1; varlist_oincv = var1}
	delta_oincv = delta
	varlist_oincv = varlist
}
# get basis functions for the constrained components -- smooth increasing convex
if (capl > 0) {
	delta = varlst = NULL
	if (13 %in% shpsx[[1]] | 16 %in% shpsx[[1]]) {
#print ('13')
		del1 = makedelta(xmat[, 1], 13)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
if (capl > 1) {
	for (i in 2:capl) {
		if (13 %in% shpsx[[i]] | 16 %in% shpsx[[i]]) {
#print ('13')
			del2 = makedelta(xmat[, i], 13)$amat
	 		m2 = length(del2) / n
		}
		delta = rbind(del1, del2)
		if (m1 > 0 | m2 > 0) {
			varlist = 1:(m1+m2)*0
			if (m1 > 0) {
				varlist[1:m1] = var1
			}
			if (m2 > 0) {
				varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
			}
			var1 = varlist
			m1 = m1 + m2
		 	del1 = delta
		}
	}
	#delta_sincv = delta
	#varlist_sincv = varlist
} #else {delta_sincv = del1; varlist_sincv = var1}
	delta_sincv = delta
	varlist_sincv = varlist
}
# get basis functions for the constrained components -- ordinal decreasing convex
if (capl > 0) {
	delta = varlst = NULL
	if (6 %in% shpsx[[1]] | 7 %in% shpsx[[1]]) {
#print ('6')
		del1 = makedelta(xmat[, 1], 6)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
if (capl > 1) {
	for (i in 2:capl) {
		if (6 %in% shpsx[[i]] | 7 %in% shpsx[[i]]) {
#print ('6')
			del2 = makedelta(xmat[, i], 6)$amat
	 		m2 = length(del2) / n
		}
	 	delta = rbind(del1, del2)
		if (m1 > 0 | m2 > 0) {
			varlist = 1:(m1+m2)*0
			if (m1 > 0) {
				varlist[1:m1] = var1
			}
			if (m2 > 0) {
				varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
			}
			var1 = varlist
			m1 = m1 + m2
	 		del1 = delta
		}
	}
	#delta_odecv = delta
	#varlist_odecv = varlist
} #else {delta_odecv = del1; varlist_odecv = var1}
	delta_odecv = delta
	varlist_odecv = varlist
}
# get basis functions for the constrained components -- smooth increasing concave
if (capl > 0) {
	delta = varlst = NULL
	if (14 %in% shpsx[[1]] | 15 %in% shpsx[[1]]) {
#print ('14')
		del1 = makedelta(xmat[, 1], 14)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
if (capl > 1) {
	for (i in 2:capl) {
		if (14 %in% shpsx[[i]] | 15 %in% shpsx[[i]]) {
#print ('14')
			del2 = makedelta(xmat[, i], 14)$amat
	 		m2 = length(del2) / n
		}
	 	delta = rbind(del1, del2)
		if (m1 > 0 | m2 > 0) {
			varlist = 1:(m1+m2)*0
			if (m1 > 0) {
				varlist[1:m1] = var1
			}
			if (m2 > 0) {
				varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
			}
			var1 = varlist
			m1 = m1 + m2
	 		del1 = delta
		}
	}
	#delta_sincc = delta
	#varlist_sincc = varlist
} #else {delta_sincc = del1; varlist_sincc = var1}
	delta_sincc = delta
	varlist_sincc = varlist
}
#zmat: similar idea, zvars: columns for z1, z2,...,z's are not factorized for now
#print (zmat)
	if (capk > 0) {
#zcat: similar to the final delta of x's 
		zcat = NULL
		zvars = NULL	
		#st = 0
		for (k in 1:capk){
			zk = zmat[, k]
			is_fac = zfacs[k]
			if (is_fac) {
				zkmat = model.matrix(~ factor(zk))[, -1, drop = FALSE]
				zvars = c(zvars, rep(k, ncol(zkmat)))
			} else {
				zkmat = zk
				zvars = c(zvars, k)
			}
			zcat = cbind(zcat, zkmat)
			#if (is_fac) {
			#	zvars = c(zvars, rep(k, ncol(zkmat)))
			#} else {zvars = c(zvars, k)}
		}
	} 
#print (paste('capk:', capk))
	vzcat = NULL
	if (capkv > 0) {
		for (k in 1:capkv){
			vzk = vzmat[, k]
			is_fac = vzfacs[k]
			if (is_fac) {
				vzkmat = model.matrix(~ factor(vzk))[, -1, drop = FALSE]
			} else {
				vzkmat = vzk
			}
			vzcat = cbind(vzcat, vzkmat)
		}
	}
	vxcat = NULL
	if (caplv > 0) {
		for (l in 1:caplv) {
			xl = vxmat[, l]
			shpl = shpsvx[l]
			vxdd = t(makedelta(xl, shpl)$amat)
			if (shpl != 17) {
				#caplv = caplv - 1
				vxcat = cbind(vxcat, vxdd)
			} else if (shpl == 17) {
				capkv = capkv + 1
				vzcat = cbind(vzcat, vxdd)
			}
#print (dim(vxdd))
		}
	}
	if (capt > 0) {
#trcat: work the same way as zcat, include in bigmat as the final part
		trcat = NULL
		trvars = NULL	
		for (k in 1:capt){
			trk = trmat[, k]
			#trkmat = model.matrix(~ factor(trk))[, -1, drop = FALSE]
			trkmat = t(tree.fun(trk))
			trcat = cbind(trcat, trkmat)
			trvars = c(trvars, rep(k, ncol(trkmat)))
		}
	} 
	if (captv > 0) {
		vtrcat = NULL
		#vtrvars = NULL	
		for (k in 1:captv){
			vtrk = vtrmat[, k]
			#trkmat = model.matrix(~ factor(trk))[, -1, drop = FALSE]
			vtrkmat = t(tree.fun(vtrk))
			vtrcat = cbind(vtrcat, vtrkmat)
			#trvars = c(trvars, rep(k, ncol(trkmat)))
		}
	} 
# make the initial random population: capl: shapes for x's, capk: in or out for z's
	if (time.est) {
		npop = 1
	} else {
		#npop = 200 #* nmod / 2e+4
if (nmod < 2e+6) {
	npop = 200
} else {
	npop = 200 + round((nmod - 2e+6) / nmod * 20)
}
#print (npop)	
	}	
	popmat = matrix(0, nrow = npop, ncol = capl + capk + capt)
#print (capl)
	for (ipop in 1:npop) {
# trunc(runif(capl)*9) gives 0 ~ 8
# trunc(runif(capk)*2) gives 0 ~ 1
		#popmat[ipop, 1:capl] = trunc(runif(capl)*9)
		#if(capk>0){popmat[ipop,(capl+1):(capl+capk)]=trunc(runif(capk)*2)}
		if (capl > 0) {
			for (l in 1:capl) {
# include flat for each x
				popmat[ipop, l] = sample(shpsx[[l]], 1) 
			}	
		}			
		if (capk > 0) {
			popmat[ipop, (capl + 1):(capl + capk)] = trunc(runif(capk)*2)
		}
		if (capt > 0) {
			popmat[ipop, (capl + capk + 1):(capl + capk + capt)] = trunc(runif(capt)*3)
		}
	}
#print (popmat)
#print (class(popmat))
	fitness = 1:npop*(-1)
## keep track of fitnesses already calculated
	cicvals = matrix(0, ncol = 2, nrow = npop*100)
	kpop = matrix(0, nrow = 5000, ncol = capl + capk + capt)
	kfit = 1:5000
	robs = 1:(npop*100)
#check!
	#fits = vector("list", npop) 
	ivals = 0
## loop through the population
#base number's for x's
if (capl > 0) {
	bx = sapply(shpsx, function(x) length(x))
#number of x's for each base
	numx = unname(table(bx))
	ubx = unique(bx)
	lbx = length(ubx)
}
if (!time.est) {
cat(paste("Evaluating the fitness of the initial population....", "\n"))
}
	for (ipop in 1:npop) {
#if (ipop%%5 == 0) print (ipop)
#print (popmat[ipop,])
#print (capl)
# i1, i2: row numbers by base 10 instead of base 9 and 2
#imod is not right.... 		
		#i1 = 0
		#for (l in 1:capl) {i1 = i1 + 5^(capl - l) * popmat[ipop, l]}	
		i1 = 0 
		if (capl > 0) {
			for (i in 1:lbx) {
				capli = numx[i]
				bxi = ubx[i]
				popi = (popmat[ipop, 1:capl, drop = FALSE])[which(bx == bxi)]
				for (il in 1:capli) {i1 = i1 + bxi^(capli - il) * popi[il]}		
			}
		}	
  		i2 = 0
  		if (capk > 0) {
  			for (k in 1:capk) {i2 = i2 + 2^(capk - k) * popmat[ipop, capl + k]}
  		}
  		i3 = 0
  		if (capt > 0) {
  			for (tr in 1:capt) {i3 = i3 + 3^(capt - tr) * popmat[ipop, capl + capk + tr]}
  		}
  		#imod = i2 * 9^capl + i1 + 1
		mult = 1
if (capl > 0) {
		for (i in 1:lbx) {
			mult = mult * (ubx[i])^numx[i]
		}
}
		mult2 = 1
		if (capk > 0) {
			mult2 = 2^capk
		}
		#imod = i3 * 2^capk * mult + i2 * mult + i1 + 1
		imod = i3 * mult2 * mult + i2 * mult + i1 + 1
## if we've done this model before, read old CIC for fitness
  		if (!all(cicvals[, 1] != imod)) {
  			rownum = robs[cicvals[, 1] == imod]
  			fitness[ipop] = -cicvals[rownum, 2]
  		} else {	
#use 0 to represent all 'flat' 
			if (sum(popmat[ipop, ]) == 0) {
				#llh=-2*(nsuc*log(mu0)+(n-nsuc)*log(1-mu0))/n
#print ('test llh')
				llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), n = n, weights = weights, fml = family$family)				
				cic = llh + log(1 + 2 / (n - 1))
				fitness[ipop] = -cic
  			} else {
				delta = matrix(1:n*0+1, nrow = 1)
				if (capkv > 0) {
					delta = rbind(delta, t(vzcat))
				}
#print (head(t(delta)))
				if (capk > 0) {
  					if (sum(popmat[ipop, (capl+1):(capl+capk)]) > 0) {
  						#zuse = 1:st < 0
						zuse = 1:ncol(zcat) < 0
						for (i in 1:capk) {
							if (popmat[ipop, capl+i] == 1) {
								zuse[zvars == i] = TRUE
							}
						}
						delta = rbind(delta, t(zcat[, zuse]))
 						#delta = rbind(1:n*0+1, t(zcat[, zuse]))
 					} #else {delta = matrix(1:n*0+1, nrow = 1)}
 				}
#new: unordered part for tree = 2, equivalent to 'z'
				if (capt > 0) {
					if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
						truse = 1:ncol(trcat) < 0
						for (i in 1:capt) {
							if (popmat[ipop, capl+capk+i] == 2) {
								truse[trvars == i] = TRUE
							}
						}
				 		delta = rbind(delta, t(trcat[, truse]))
				 	} 
				}
				if (sum(popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13) > 0) {
 					usex = popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13
 					delta = rbind(delta, t(xmat[, usex]))
 				}
#new:
if (caplv > 0) {
	ch = shpsvx %in% c(3, 4, 11, 12)
	if (any(ch)) {
		vxin = vxmat[, which(ch)]
		delta = rbind(delta, t(vxin))
	}
}
				np = dim(delta)[1]
#print (popmat[ipop, ])
#new:
				if (capl > 0) {
					for (l in 1:capl) {
 						if (popmat[ipop, l] > 0) {
 							if (popmat[ipop, l] == 1) {#incr
#print ('1')
 								deladd = delta_om[varlist_om == l, ]
 							} else if (popmat[ipop, l] == 2) {#decr
#print ('2')
 								deladd = -delta_om[varlist_om == l, ]
 							} else if (popmat[ipop, l] == 3) {#conv
#print ('3')
 								deladd = delta_ocv[varlist_ocv == l, ]
 							} else if (popmat[ipop, l] == 4) {#conc
#print ('4')
 								deladd = -delta_ocv[varlist_ocv == l, ]
 							} else if (popmat[ipop, l] == 5) {#incr.conv
#print ('5')
 								deladd = delta_oincv[varlist_oincv == l, ]
 							} else if (popmat[ipop, l] == 8) {#decr.conc
#print ('8')
 								deladd = -delta_oincv[varlist_oincv == l, ]
 							} else if (popmat[ipop, l] == 6) {#decr.conv 
#print ('6')
 								deladd = delta_odecv[varlist_odecv == l, ]
 							} else if (popmat[ipop, l] == 7) {#incr.conc
#print ('7')
 								deladd = -delta_odecv[varlist_odecv == l, ]
 							} else if (popmat[ipop, l] == 9) {
#print ('9')
 								deladd = delta_sm[varlist_sm == l, ]
 							} else if (popmat[ipop, l] == 10) {
#print ('10')
 								deladd = -delta_sm[varlist_sm == l, ]
 							} else if (popmat[ipop, l] == 11) {
#print ('11')
 								deladd = delta_scv[varlist_scv == l, ]
 							} else if (popmat[ipop, l] == 12) {
#print ('12')
 								deladd = -delta_scv[varlist_scv == l, ]
 							} else if (popmat[ipop, l] == 13) {
#print ('13')	
 								deladd = delta_sincv[varlist_sincv == l, ]
 							} else if (popmat[ipop, l] == 16) {
#print ('16')
 								deladd = -delta_sincv[varlist_sincv == l, ]
 							} else if (popmat[ipop, l] == 14) {
#print ('14')
 								deladd = delta_sincc[varlist_sincc == l, ]
 							} else if (popmat[ipop, l] == 15) {
#print ('15')
 								deladd = -delta_sincc[varlist_sincc == l, ]
 							}
 							delta = rbind(delta, deladd)
#print (paste('dim: ', dim(delta)))
 						}	
 					} 
				}
				if (caplv > 0) {
					if (!is.null(vxcat)) {
						delta = rbind(delta, t(vxcat))
					}
				}
#tree part for tree var
				if (capt > 0) {
				  	if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
						truse = 1:ncol(trcat) < 0
						for (i in 1:capt) {
							if (popmat[ipop, capl+capk+i] == 1) {
								truse[trvars == i] = TRUE
							}
						}
				 		delta = rbind(delta, t(trcat[, truse]))
				 	} 
				}
				if (captv > 0) {
					delta = rbind(delta, t(vtrcat))
				}
#print (paste('dim: ', dim(t(trcat))))
#print (delta)
#print (cpar)
				ans = irls(y, delta, np, nsim = 100, family = family, weights = weights, cpar = cpar)
 				ivals = ivals + 1
 				cicvals[ivals, 2] = ans$cic
 				cicvals[ivals, 1] = imod
 				fitness[ipop] = -cicvals[ivals, 2]
 				kpop[ivals, ] = popmat[ipop, ]
				kfit[ivals] = -cicvals[ivals, 2]
				#fits[[ivals]] = ans$fhat
				
			}
		}				
 		#print (ipop - fitness[ipop])	
 	}
#iter = 0
mxf = NULL
mnf = NULL
#mdf = NULL
#mq3f = NULL
if (!time.est) {
##   Now loop through generations					
	nrep = 0
	obs = 1:npop
	q1 = trunc(npop/4)
	q2 = 2 * q1
	q3 = round(.98 * npop)
	nbig = 100
	nm = 10
	check = TRUE
	maxfit = -1000
	while (nrep < nbig & check) {
	#while (nrep < nbig) {
		ord = order(-fitness)
		popmat = popmat[ord, ,drop = FALSE]
		fitness = fitness[ord]
		nrep = nrep + 1
cat(paste("Iter =", nrep, " | Mean =", format(mean(fitness), digits = 6), " | Best =", format(max(fitness), digits = 6), "\n"))
# mutate!  randomly throughout population  ## don't mutate the most fit
# imut is the row # of popmat
# igene is the column # of x, z, and tree
		imut = trunc(runif(nm) * (npop - 1) + 2)  
		#igene = trunc(runif(nm) * (capk + capl) + 1)
#check!
igene = trunc(runif(nm) * (capk + capl + capt) + 1)
		for (im in 1:nm) {
			ipop = imut[im]
			if (igene[im] <= capl) {
				#popmat[imut[im], igene[im]] = trunc(runif(1)*9)
				#pos = imut[im]
				popmat[imut[im], igene[im]] = sample(shpsx[[igene[im]]], 1)
			} else if (igene[im] > capl & igene[im] <= (capl + capk)) {
				popmat[imut[im], igene[im]] = trunc(runif(1)*2)
			} else if (igene[im] > (capl + capk)) {
				popmat[imut[im], igene[im]] = trunc(runif(1)*3)
			}		
	  		#ipop = imut[im]
# find fitness of mutated phenotype
			#i1=0;for(l in 1:capl){i1=i1+9^(capl-l)*popmat[ipop,l]}		
	  		#i2=0
	  		#if(capk>0){
	  		#	for(k in 1:capk){i2=i2+2^(capk-k)*popmat[ipop,capl+k]}
	  		#}
	  		#imod=i2*9^capl+i1+1
			i1 = 0 
			if (capl > 0) {
				for (i in 1:lbx) {
					capli = numx[i]
					bxi = ubx[i]
					popi = (popmat[ipop, 1:capl, drop = FALSE])[which(bx == bxi)]
					for (il in 1:capli) {i1 = i1 + bxi^(capli - il) * popi[il]}		
				}
			}	
  			i2 = 0
  			if (capk > 0) {
  				for (k in 1:capk) {i2 = i2 + 2^(capk - k) * popmat[ipop, capl + k]}
  			}
  			i3 = 0
  			if (capt > 0) {
  				for (tr in 1:capt) {i3 = i3 + 3^(capt - tr) * popmat[ipop, capl + capk + tr]}
  			}
			mult = 1
			if (capl > 0) {
				for (i in 1:lbx) {
					mult = mult * (ubx[i])^numx[i]
				}
			}
			mult2 = 1
			if (capk > 0) {
				mult2 = 2^capk
			}
			imod = i3 * mult2 * mult + i2 * mult + i1 + 1
## if we've done this model before, read old CIC for fitness
  			if (!all(cicvals[, 1] != imod)) {
  				rownum = robs[cicvals[,1] == imod]
  				fitness[ipop] = -cicvals[rownum, 2]
  				#print("done")
## calculate fitness
  			} else {		
 	  			if (sum(popmat[ipop, ]) == 0) {
					#llh=-2*(nsuc*log(mu0)+(n-nsuc)*log(1-mu0))/n
					#cic=llh+log(1+2/(n-1))
	  				#fitness[ipop]=-cic
					llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), n = n, weights = weights, fml = family$family)				
					cic = llh + log(1 + 2 / (n - 1))
					fitness[ipop] = -cic
	  			} else {
	  				delta = matrix(1:n*0+1, nrow = 1)
					if (capkv > 0) {
						delta = rbind(delta, t(vzcat))
					}
	  				if (capk > 0) {
  						if (sum(popmat[ipop, (capl+1):(capl+capk)]) > 0) {
							zuse = 1:ncol(zcat) < 0
							for (i in 1:capk) {
								if (popmat[ipop, capl+i] == 1) {
									zuse[zvars == i] = TRUE
								}
							}
							delta = rbind(delta, t(zcat[, zuse]))
 							#delta = rbind(1:n*0+1, t(zcat[, zuse]))
 						} #else {
						#	delta = matrix(1:n*0+1, nrow = 1)
						#}
 					}
					if (capt > 0) {
						if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
							truse = 1:ncol(trcat) < 0
							for (i in 1:capt) {
								if (popmat[ipop, capl+capk+i] == 2) {
									truse[trvars == i] = TRUE
								}
							}
					 		delta = rbind(delta, t(trcat[, truse]))
					 	} 
					}
  					if (sum(popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13) > 0) {
 						usex = popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13
 						delta = rbind(delta, t(xmat[, usex]))
 					}
#new:
if (caplv > 0) {
	ch = shpsvx %in% c(3, 4, 11, 12)
	if (any(ch)) {
		vxin = vxmat[, which(ch)]
		delta = rbind(delta, t(vxin))
	}
}
					np = dim(delta)[1]
					if (capl > 0) {
 						for (l in 1:capl) {
 							if (popmat[ipop, l] > 0) {
 								if (popmat[ipop, l] == 1) {#incr
 									deladd = delta_om[varlist_om == l, ]
 								} else if (popmat[ipop, l] == 2) {#decr
 									deladd = -delta_om[varlist_om == l, ]
 								} else if (popmat[ipop, l] == 3) {#conv
 									deladd = delta_ocv[varlist_ocv == l, ]
 								} else if (popmat[ipop, l] == 4) {#conc
 									deladd = -delta_ocv[varlist_ocv == l, ]
 								} else if (popmat[ipop, l] == 5) {#incr.conv
 									deladd = delta_oincv[varlist_oincv == l, ]
 								} else if (popmat[ipop, l] == 8) {#decr.conc
 									deladd = -delta_oincv[varlist_oincv == l, ]
 								} else if (popmat[ipop, l] == 6) {#decr.conv 
 									deladd = delta_odecv[varlist_odecv == l, ]
 								} else if (popmat[ipop, l] == 7) {#incr.conc
 									deladd = -delta_odecv[varlist_odecv == l, ]
 								} else if (popmat[ipop, l] == 9) {
 									deladd = delta_sm[varlist_sm == l, ]
 								} else if (popmat[ipop, l] == 10) {
 									deladd = -delta_sm[varlist_sm == l, ]
 								} else if (popmat[ipop, l] == 11) {
 									deladd = delta_scv[varlist_scv == l, ]
 								} else if (popmat[ipop, l] == 12) {
 									deladd = -delta_scv[varlist_scv == l, ]
 								} else if (popmat[ipop, l] == 13) {
 									deladd = delta_sincv[varlist_sincv == l, ]
 								} else if (popmat[ipop, l] == 16) {
 									deladd = -delta_sincv[varlist_sincv == l, ]
 								} else if (popmat[ipop, l] == 14) {
 									deladd = delta_sincc[varlist_sincc == l, ]
 								} else if (popmat[ipop, l] == 15) {		
 									deladd = -delta_sincc[varlist_sincc == l, ]
 								}
 								delta = rbind(delta, deladd)
							}	
 						}
					}
					if (caplv > 0) {
						if (!is.null(vxcat)) {
							delta = rbind(delta, t(vxcat))
						}
					}
					if (capt > 0) {
					  	if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
							truse = 1:ncol(trcat) < 0
							for (i in 1:capt) {
								if (popmat[ipop, capl+capk+i] == 1) {
									truse[trvars == i] = TRUE
								}
							}
					 		delta = rbind(delta, t(trcat[, truse]))
					 	} 
					}
					if (captv > 0) {
						delta = rbind(delta, t(vtrcat))
					}
 					ans = irls(y, delta, np, nsim = 100, family = family, weights = weights, cpar = cpar)
 					ivals = ivals + 1
 					cicvals[ivals, 2] = ans$cic
 					cicvals[ivals, 1] = imod
 					fitness[ipop] = -cicvals[ivals, 2]
 					kpop[ivals, ] = popmat[ipop, ]
					kfit[ivals] = -cicvals[ivals, 2]
					#fits[[ivals]] = ans$fhat
 				
 				}
 			}
 		}
# reproduction: replace middle half with offspring: offspring combine elite with other  
		for (ipop in q1:q3) {
			mom = trunc(runif(1) * npop + 1)
			dad = trunc(runif(1) * q1 + 1)
			digits = runif(capl + capk + capt) > .5
			popmat[ipop, digits] = popmat[mom, digits]
			popmat[ipop, !digits] = popmat[dad, !digits]
# find fitness of baby
			#i1=0;for(l in 1:capl){i1=i1+9^(capl-l)*popmat[ipop,l]}		
	  		#i2=0
	  		#if(capk>0){
	  		#	for(k in 1:capk){i2=i2+2^(capk-k)*popmat[ipop,capl+k]}
	  		#}
	  		#imod=i2*9^capl+i1+1
			i1 = 0 
			if (capl > 0) {
				for (i in 1:lbx) {
					capli = numx[i]
					bxi = ubx[i]
					popi = (popmat[ipop, 1:capl, drop = FALSE])[which(bx == bxi)]
					for (il in 1:capli) {i1 = i1 + bxi^(capli - il) * popi[il]}		
				}
			}	
  			i2 = 0
  			if (capk > 0) {
  				for (k in 1:capk) {i2 = i2 + 2^(capk - k) * popmat[ipop, capl + k]}
  			}
  			i3 = 0
  			if (capt > 0) {
  				for (tr in 1:capt) {i3 = i3 + 3^(capt - tr) * popmat[ipop, capl + capk + tr]}
  			}
			mult = 1
if (capl > 0) {
			for (i in 1:lbx) {
				mult = mult * (ubx[i])^numx[i]
			}
}
			mult2 = 1
			if (capk > 0) {	
				mult2 = 2^capk		
			}
			imod = i3 * mult2 * mult + i2 * mult + i1 + 1
## if we've done this model before, read old CIC for fitness
  			if (!all(cicvals[, 1] != imod)) {
  				rownum = robs[cicvals[, 1] == imod]
  				fitness[ipop] = -cicvals[rownum, 2]
  				#print("done")
## calculate fitness
  			} else {
 	  			if (sum(popmat[ipop, ]) == 0) {
					#llh=-2*(nsuc*log(mu0)+(n-nsuc)*log(1-mu0))/n
					#cic=llh+log(1+2/(n-1))
	  				#fitness[ipop]=-cic
					llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), n = n, weights = weights, fml = family$family)				
					cic = llh + log(1 + 2 / (n - 1))
					fitness[ipop] = -cic
	  			} else {
	  				delta = matrix(1:n*0+1, nrow = 1)
					if (capkv > 0) {
						delta = rbind(delta, t(vzcat))
					}
	  				if (capk > 0) {
  						if (sum(popmat[ipop, (capl+1):(capl+capk)]) > 0) {
							zuse = 1:ncol(zcat) < 0
							for (i in 1:capk) {
								if (popmat[ipop, capl+i] == 1) {
									zuse[zvars == i] = TRUE
								}
							}
							delta = rbind(delta, t(zcat[, zuse]))
 							#delta = rbind(1:n*0+1, t(zcat[, zuse]))
 						} 
 					} 
					if (capt > 0) {
						if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
							truse = 1:ncol(trcat) < 0
							for (i in 1:capt) {
								if (popmat[ipop, capl+capk+i] == 2) {
									truse[trvars == i] = TRUE
								}
							}
					 		delta = rbind(delta, t(trcat[, truse]))
					 	} 
					}
  					if (sum(popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13) > 0) {
 						usex = popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13
 						delta = rbind(delta, t(xmat[, usex]))
 					}
#new:
if (caplv > 0) {
	ch = shpsvx %in% c(3, 4, 11, 12)
	if (any(ch)) {
		vxin = vxmat[, which(ch)]
		delta = rbind(delta, t(vxin))
	}
}
					np = dim(delta)[1]
					if (capl > 0) {
 						for (l in 1:capl) {
 							if (popmat[ipop, l] > 0) {
 								if (popmat[ipop, l] == 1) {#incr
 									deladd = delta_om[varlist_om == l, ]
 								} else if (popmat[ipop, l] == 2) {#decr
 									deladd = -delta_om[varlist_om == l, ]
 								} else if (popmat[ipop, l] == 3) {#conv
 									deladd = delta_ocv[varlist_ocv == l, ]
 								} else if (popmat[ipop, l] == 4) {#conc
 									deladd = -delta_ocv[varlist_ocv == l, ]
 								} else if (popmat[ipop, l] == 5) {#incr.conv
 									deladd = delta_oincv[varlist_oincv == l, ]
 								} else if (popmat[ipop, l] == 8) {#decr.conc
 									deladd = -delta_oincv[varlist_oincv == l, ]
 								} else if (popmat[ipop, l] == 6) {#decr.conv 
 									deladd = delta_odecv[varlist_odecv == l, ]
 								} else if (popmat[ipop, l] == 7) {#incr.conc
 									deladd = -delta_odecv[varlist_odecv == l, ]
 								} else if (popmat[ipop, l] == 9) {
 									deladd = delta_sm[varlist_sm == l, ]
 								} else if (popmat[ipop, l] == 10) {
 									deladd = -delta_sm[varlist_sm == l, ]
 								} else if (popmat[ipop, l] == 11) {
 									deladd = delta_scv[varlist_scv == l, ]
 								} else if (popmat[ipop, l] == 12) {
 									deladd = -delta_scv[varlist_scv == l, ]
 								} else if (popmat[ipop, l] == 13) {
 									deladd = delta_sincv[varlist_sincv == l, ]
 								} else if (popmat[ipop, l] == 16) {
 									deladd = -delta_sincv[varlist_sincv == l, ]
 								} else if (popmat[ipop, l] == 14) {
 									deladd = delta_sincc[varlist_sincc == l, ]
 								} else if (popmat[ipop, l] == 15) {
 									deladd = -delta_sincc[varlist_sincc == l, ]
 								}
 								delta = rbind(delta, deladd)
							}	
 						}
					}
					if (caplv > 0) {
						if (!is.null(vxcat)) {
							delta = rbind(delta, t(vxcat))
						}
					}
					if (capt > 0) {
					  	if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
							truse = 1:ncol(trcat) < 0
							for (i in 1:capt) {
								if (popmat[ipop, capl+capk+i] == 1) {
									truse[trvars == i] = TRUE
								}
							}
					 		delta = rbind(delta, t(trcat[, truse]))
					 	} 
					}
					if (captv > 0) {
						delta = rbind(delta, t(vtrcat))
					}
 					ans = irls(y, delta, np, nsim = 100, family = family, weights = weights, cpar = cpar)
 					ivals = ivals + 1
 					cicvals[ivals, 2] = ans$cic
 					cicvals[ivals, 1] = imod
 					fitness[ipop] = -cicvals[ivals, 2]
					kpop[ivals, ] = popmat[ipop, ]
					kfit[ivals] = -cicvals[ivals, 2]
					#fits[[ivals]] = ans$fhat
 				}
 			}
 		}	
## immigration
		for (ipop in (q3 + 1):npop) {
			if (capl > 0) {
				for (l in 1:capl) {
					popmat[ipop, l] = sample(shpsx[[l]], 1) 
				}	
			}	
			if (capk > 0) {
				popmat[ipop, (capl + 1):(capl + capk)] = trunc(runif(capk)*2)
			}
			if (capt > 0) {
				popmat[ipop, (capl + capk + 1):(capl + capk + capt)] = trunc(runif(capt)*3)
			}
# find fitness of new immigrant
			i1 = 0 
			if (capl > 0) {
				for (i in 1:lbx) {
					capli = numx[i]
					bxi = ubx[i]
					popi = (popmat[ipop, 1:capl, drop = FALSE])[which(bx == bxi)]
					for (il in 1:capli) {i1 = i1 + bxi^(capli - il) * popi[il]}		
				}
			}	
  			i2 = 0
  			if (capk > 0) {
  				for (k in 1:capk) {i2 = i2 + 2^(capk - k) * popmat[ipop, capl + k]}
  			}
  			i3 = 0
  			if (capt > 0) {
  				for (tr in 1:capt) {i3 = i3 + 3^(capt - tr) * popmat[ipop, capl + capk + tr]}
  			}
			mult = 1
if (capl > 0) {
			for (i in 1:lbx) {
				mult = mult * (ubx[i])^numx[i]
			}
}
			mult2 = 1
			if (capk > 0) {
				mult2 = 2^capk
			}
			imod = i3 * mult2 * mult + i2 * mult + i1 + 1
## if we've done this model before, read old CIC for fitness
  			if (!all(cicvals[, 1] != imod)) {
  				rownum = robs[cicvals[, 1] == imod]
  				fitness[ipop] = -cicvals[rownum, 2]
  			} else {		
	  			if (sum(popmat[ipop, ]) == 0) {
					llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), n = n, weights = weights, fml = family$family)				
					cic = llh + log(1 + 2 / (n - 1))
					fitness[ipop] = -cic
	  			} else {
	  				delta = matrix(1:n*0+1, nrow = 1)
					if (capkv > 0) {
						delta = rbind(delta, t(vzcat))
					}
	  				if (capk > 0) {
  						if (sum(popmat[ipop, (capl+1):(capl+capk)]) > 0) {
							zuse = 1:ncol(zcat) < 0
							for (i in 1:capk) {
								if (popmat[ipop, capl+i] == 1) {
									zuse[zvars == i] = TRUE
								}
							}
							delta = rbind(delta, t(zcat[, zuse]))
 							#delta = rbind(1:n*0+1, t(zcat[, zuse]))
 						} #else {delta = matrix(1:n*0+1, nrow = 1)}
 					}
					if (capt > 0) {
						if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
							truse = 1:ncol(trcat) < 0
							for (i in 1:capt) {
								if (popmat[ipop, capl+capk+i] == 2) {
									truse[trvars == i] = TRUE
								}
							}
					 		delta = rbind(delta, t(trcat[, truse]))
					 	} 
					}
  					if (sum(popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13) > 0) {
 						usex = popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13
 						delta = rbind(delta, t(xmat[, usex]))
 					}
#new:
if (caplv > 0) {
	ch = shpsvx %in% c(3, 4, 11, 12)
	if (any(ch)) {
		vxin = vxmat[, which(ch)]
		delta = rbind(delta, t(vxin))
	}
}
					np = dim(delta)[1]
					if (capl > 0) {
	 					for (l in 1:capl) {
	 						if (popmat[ipop, l] > 0) {
	 							if (popmat[ipop, l] == 1) {#incr
 									deladd = delta_om[varlist_om == l, ]
 								} else if (popmat[ipop, l] == 2) {#decr
 									deladd = -delta_om[varlist_om == l, ]
 								} else if (popmat[ipop, l] == 3) {#conv
 									deladd = delta_ocv[varlist_ocv == l, ]
 								} else if (popmat[ipop, l] == 4) {#conc
 									deladd = -delta_ocv[varlist_ocv == l, ]
 								} else if (popmat[ipop, l] == 5) {#incr.conv
 									deladd = delta_oincv[varlist_oincv == l, ]
 								} else if (popmat[ipop, l] == 8) {#decr.conc
 									deladd = -delta_oincv[varlist_oincv == l, ]
 								} else if (popmat[ipop, l] == 6) {#decr.conv 
 									deladd = delta_odecv[varlist_odecv == l, ]
 								} else if (popmat[ipop, l] == 7) {#incr.conc
 									deladd = -delta_odecv[varlist_odecv == l, ]
 								} else if (popmat[ipop, l] == 9) {
 									deladd = delta_sm[varlist_sm == l, ]
 								} else if (popmat[ipop, l] == 10) {
 									deladd = -delta_sm[varlist_sm == l, ]
 								} else if (popmat[ipop, l] == 11) {
 									deladd = delta_scv[varlist_scv == l, ]
 								} else if (popmat[ipop, l] == 12) {
 									deladd = -delta_scv[varlist_scv == l, ]
 								} else if (popmat[ipop, l] == 13) {
 									deladd = delta_sincv[varlist_sincv == l, ]
 								} else if (popmat[ipop, l] == 16) {
 									deladd = -delta_sincv[varlist_sincv == l, ]
 								} else if (popmat[ipop, l] == 14) {
 									deladd = delta_sincc[varlist_sincc == l, ]
 								} else if (popmat[ipop, l] == 15) {
 									deladd = -delta_sincc[varlist_sincc == l, ]
 								}
 								delta = rbind(delta, deladd)
							}	
 						}
					}
					if (caplv > 0) {
						if (!is.null(vxcat)) {
							delta = rbind(delta, t(vxcat))
						}
					}
					if (capt > 0) {
					  	if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
							truse = 1:ncol(trcat) < 0
							for (i in 1:capt) {
								if (popmat[ipop, capl+capk+i] == 1) {
									truse[trvars == i] = TRUE
								}
							}
					 		delta = rbind(delta, t(trcat[, truse]))
					 	} 
					}
					if (captv > 0) {
						delta = rbind(delta, t(vtrcat))
					}
  					ans = irls(y, delta, np, nsim = 100, family = family, weights = weights, cpar = cpar)
 					ivals = ivals + 1
 					cicvals[ivals, 2] = ans$cic
 					cicvals[ivals, 1] = imod
 					fitness[ipop] = -cicvals[ivals, 2]
					kpop[ivals, ] = popmat[ipop, ]
					kfit[ivals] = -cicvals[ivals, 2]
					#fits[[ivals]] = ans$fhat
				}
			}
		}		
## targeted mutation, add if capl > 0
		if (capl > 0) {
			#if (max(fitness) > maxfit + 1e-8) {
				ord = order(-fitness)
				popmat = popmat[ord, ,drop = FALSE]
				fitness = fitness[ord]
				repl = 0
				for (ic in 1:capl) {
				#for (ic in 1:(capl - 1)) {
					tpop = popmat[1, ]
#??
					#for (j in 1:9) {
					for (j in 1:length(shpsx[[ic]])) {
						#tpop[ic] = tpop[ic] + 1
						#if (tpop[ic] >= 9) {
						#if (tpop[ic] >= max(shpsx[[ic]])) {
							#tpop[ic] = 1
						#	tpop[ic] = min(shpsx[[ic]])
						#}
						tpop[ic] = sample(shpsx[[ic]], 1)
						#print(tpop)
# find fitness of mutated phenotype
#						i1=0;for(l in 1:capl){i1=i1+9^(capl-l)*tpop[l]}		
#	 					i2=0
#	  					if(capk>0){
#	  						for(k in 1:capk){i2=i2+2^(capk-k)*tpop[capl+k]}
#	  					}
#	  					imod=i2*9^capl+i1+1
						i1 = 0 
						if (capl > 0) {
							for (i in 1:lbx) {
								capli = numx[i]
								bxi = ubx[i]
								popi = (popmat[ipop, 1:capl, drop = FALSE])[which(bx == bxi)]
								for (il in 1:capli) {i1 = i1 + bxi^(capli - il) * popi[il]}		
							}
						}	
		  				i2 = 0
		  				if (capk > 0) {
		  					for (k in 1:capk) {i2 = i2 + 2^(capk - k) * popmat[ipop, capl + k]}
		  				}
		  				i3 = 0
		  				if (capt > 0) {
		  					for (tr in 1:capt) {i3 = i3 + 3^(capt - tr) * popmat[ipop, capl + capk + tr]}
		  				}
						mult = 1
						for (i in 1:lbx) {
							mult = mult * (ubx[i])^numx[i]
						}
						mult2 = 1
						if (capk > 0) {
							mult2 = 2^capk
						}
						imod = i3 * mult2 * mult + i2 * mult + i1 + 1
## if we've done this model before, read old CIC for fitness
	  					if (!all(cicvals[, 1] != imod)) {
	  						rownum = robs[cicvals[, 1] == imod]
	  						tfit = -cicvals[rownum, 2]
		  				} else {	
 		  					if (sum(tpop) == 0) {
								#llh=-2*(nsuc*log(mu0)+(n-nsuc)*log(1-mu0))/n
								llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), n = n, weights = weights, fml = family$family)				
								cic = llh + log(1 + 2 / (n - 1))
 								tfit = -cic
 							} else {
  								delta = matrix(1:n*0+1, nrow = 1)
								if (capkv > 0) {
									delta = rbind(delta, t(vzcat))
								}
  								if (capk > 0) {
									if (sum(tpop[(capl+1):(capl+capk)]) > 0) {
										zuse = 1:ncol(zcat) < 0
										for (i in 1:capk) {
											if (tpop[capl+i] == 1) {
												zuse[zvars == i] = TRUE
											}
										}
										delta = rbind(delta, t(zcat[, zuse]))
 										#delta = rbind(1:n*0+1, t(zcat[, zuse]))
 									} #else {delta = matrix(1:n*0+1, nrow = 1)}
 								} 
								if (capt > 0) {
									if (sum(tpop[(capl+capk+1):(capl+capk+capt)]) > 0) {
										truse = 1:ncol(trcat) < 0
										for (i in 1:capt) {
											if (tpop[capl+capk+i] == 2) {
												truse[trvars == i] = TRUE
											}
										}
 										delta = rbind(delta, t(trcat[, truse]))
 									}
								}
								if (sum(tpop[1:capl] > 2 & tpop[1:capl] < 5 | tpop[1:capl] > 10 & tpop[1:capl] < 13) > 0) {
 									usex = tpop[1:capl] > 2 & tpop[1:capl] < 5 | tpop[1:capl] > 10 & tpop[1:capl] < 13
 									delta = rbind(delta, t(xmat[, usex]))
 								}
#new:
if (caplv > 0) {
	ch = shpsvx %in% c(3, 4, 11, 12)
	if (any(ch)) {
		vxin = vxmat[, which(ch)]
		delta = rbind(delta, t(vxin))
	}
}
								np = dim(delta)[1]
 								for (l in 1:capl) {
 									if (tpop[l] > 0) {
 										if (tpop[l] == 1) {#incr
 											deladd = delta_om[varlist_om == l, ]
 										} else if (tpop[l] == 2) {#decr
 											deladd = -delta_om[varlist_om == l, ]
 										} else if (tpop[l] == 3) {#conv
 											deladd = delta_ocv[varlist_ocv == l, ]
 										} else if (tpop[l] == 4) {#conc
 											deladd = -delta_ocv[varlist_ocv == l, ]
 										} else if (tpop[l] == 5) {#incr.conv
 											deladd = delta_oincv[varlist_oincv == l, ]
 										} else if (tpop[l] == 8) {#decr.conc
 											deladd = -delta_oincv[varlist_oincv == l, ]
 										} else if (tpop[l] == 6) {#decr.conv 
 											deladd = delta_odecv[varlist_odecv == l, ]
 										} else if (tpop[l] == 7) {#incr.conc
 											deladd = -delta_odecv[varlist_odecv == l, ]
 										} else if (tpop[l] == 9) {
 											deladd = delta_sm[varlist_sm == l, ]
 										} else if (tpop[l] == 10) {
 											deladd = -delta_sm[varlist_sm == l, ]
 										} else if (tpop[l] == 11) {
 											deladd = delta_scv[varlist_scv == l, ]
 										} else if (tpop[l] == 12) {
 											deladd = -delta_scv[varlist_scv == l, ]
 										} else if (tpop[l] == 13) {
 											deladd = delta_sincv[varlist_sincv == l, ]
 										} else if (tpop[l] == 16) {
 											deladd = -delta_sincv[varlist_sincv == l, ]
 										} else if (tpop[l] == 14) {
 											deladd = delta_sincc[varlist_sincc == l, ]
 										} else if (tpop[l] == 15) {
 											deladd = -delta_sincc[varlist_sincc == l, ]
 										}
 										delta = rbind(delta, deladd)
									}	
 								}
								if (caplv > 0) {
									if (!is.null(vxcat)) {
										delta = rbind(delta, t(vxcat))
									}
								}
#nsim2 = 100 for target part
								if (capt > 0) {
									if (sum(tpop[(capl+capk+1):(capl+capk+capt)]) > 0) {
										truse = 1:ncol(trcat) < 0
										for (i in 1:capt) {
											if (tpop[capl+capk+i] == 1) {
												truse[trvars == i] = TRUE
											}
										}
 										delta = rbind(delta, t(trcat[, truse]))
 									}
								}
								if (captv > 0) {
									delta = rbind(delta, t(vtrcat))
								}
 								ans = irls(y, delta, np, nsim = 100, family = family, weights = weights, cpar = cpar)
 								ivals = ivals + 1
 								cicvals[ivals, 2] = ans$cic
 								cicvals[ivals, 1] = imod
								kpop[ivals, ] = tpop
								kfit[ivals] = -cicvals[ivals, 2]
 								tfit = -cicvals[ivals, 2]
								#fits[[ivals]] = ans$fhat
								#print(tfit)
 								if (tfit > fitness[1]) {
  									#points(npop - repl, tfit, pch = 5, col = "darkgreen")
 									popmat[npop - repl, ] = tpop
 									fitness[npop - repl] = tfit
 									repl = repl + 1
 								}
							}						
						}	
					}	
 				}
		}
		maxfit = max(fitness)
		if (fitness[q1] > maxfit - 1e-8){check = FALSE}
		#if (fitness[q1] > maxfit - 1e-6){check = FALSE}
mxf = c(mxf, max(fitness))
mnf = c(mnf, mean(fitness))
#mdf = c(mdf, quantile(fitness, probs = .5))
#mq3f = c(mdf, quantile(fitness, probs = .75))
	}
}
ord = order(-fitness)
# sort the pop by fitness
popmat = popmat[ord, ,drop = FALSE]
fitness = fitness[ord]
rslt = new.env()
rslt$fitness = fitness
#rslt$fhat = fits
pop2 = matrix(0, nrow = npop, ncol = (capl+capk+capt))
if (capl > 0) {
	for (i in 1:npop) {
		pop2[i, 1:capl] = apply(popmat[i, 1:capl, drop = FALSE], 2, ShapeToChar)
	}
	#names(pop2)[1:capl] = xnms
}
if (capk > 0) {
	for (i in 1:npop) {
		pop2[i, (capl+1):(capl+capk)] = apply(popmat[i, (capl+1):(capl+capk), drop = FALSE], 2, function(num) ShapeToChar(num, tag = "z"))
	}
	#names(pop2)[(capl+1):(capl+capk)] = znms
}
if (capt > 0) {
	for (i in 1:npop) {
		pop2[i, (capl+capk+1):(capl+capk+capt)] = apply(popmat[i, (capl+capk+1):(capl+capk+capt), drop = FALSE], 2, function(num) ShapeToChar(num, tag = "tree"))
	}
	#names(pop2)[(capl+capk+1):(capl+capk+capt)] = trnms
}
#colnames(pop2) = c(xnms, znms, trnms)
pop2 = as.data.frame(pop2, stringsAsFactors = FALSE)
rslt$pop = cbind(popmat, fitness)
pop2 = cbind(pop2, fitness)
rslt$pop2 = pop2
#rslt$top = pop2[1, stringsAsFactors = FALSE]
rslt$mxf = mxf
rslt$mnf = mnf
#rslt$mdf = mdf
#rslt$mq3f = mq3f
rslt$GA = TRUE
rslt$vzcat = vzcat
#print (paste('nrep: ', nrep))
#class(rslt) = "cgam"
return (rslt)
}


##############################
#go through models one-by-one#
##############################
ConstrALL = function(y, xmat, zmat, trmat, family = gaussian(), shpsx = NULL, shpsvx = NULL, shpsz = NULL, shpst = NULL, cpar = 1.2, zfacs = NULL, weights = NULL, vzmat = NULL, vzfacs = NULL, vxmat = NULL, vtrmat = NULL) {
	linkfun = family$linkfun
	cicfamily = CicFamily(family)
	llh.fun = cicfamily$llh.fun
	etahat.fun = cicfamily$etahat.fun
	gr.fun = cicfamily$gr.fun
	wt.fun = cicfamily$wt.fun
	zvec.fun = cicfamily$zvec.fun
	muhat.fun = cicfamily$muhat.fun
	ysim.fun = cicfamily$ysim.fun
	deriv.fun = cicfamily$deriv.fun
	dev.fun = cicfamily$dev.fun 
	n = length(y)
	sm = 1e-7 
	capl = length(xmat) / n
	if (capl < 1) {capl = 0}
	if (round(capl, 8) != round(capl, 1)) {
		stop ("Incompatible dimensions for xmat!")
	}
#check!
	if (capl > 0) {
		for(i in 1:capl) {
			xmat[,i] = (xmat[,i] - min(xmat[,i])) / (max(xmat[,i]) - min(xmat[,i]))
			#xmat[,i] = (xmat[,i] - mean(xmat[,i])) / sd(xmat[,i])
			#xmat[,i] = xmat[,i] / sd(xmat[,i])
		}
	}
	caplv = length(vxmat) / n
	if (caplv < 1) {caplv = 0}
	if (round(caplv, 8) != round(caplv, 1)) {
		stop ("Incompatible dimensions for zmat!")
	}
#new:
	if (caplv > 0) {
		for(i in 1:caplv) {
			vxmat[,i] = (vxmat[,i] - min(vxmat[,i])) / (max(vxmat[,i]) - min(vxmat[,i]))
			#vxmat[,i] = (vxmat[,i] - mean(vxmat[,i])) / sd(vxmat[,i])
			#vxmat[,i] = vxmat[,i] / sd(vxmat[,i])
		}
	}
	capk = length(zmat) / n
	if (capk < 1) {capk = 0}
	if (round(capk, 8) != round(capk, 1)) {
		stop ("Incompatible dimensions for zmat!")
	}
	capkv = length(vzmat) / n
	if (capkv < 1) {capkv = 0}
	if (round(capkv, 8) != round(capkv, 1)) {
		stop ("Incompatible dimensions for zmat!")
	}
	capt = length(trmat) / n
	if (capt < 1) {capt = 0}
	if (round(capt, 8) != round(capt, 1)) {
		stop ("Incompatible dimensions for trmat!")
	}
	captv = length(vtrmat) / n
	if (captv < 1) {captv = 0}
	if (round(captv, 8) != round(captv, 1)) {
		stop ("Incompatible dimensions for trmat!")
	}
################################################################
##get basis functions for all allowed shapes for each component#
#not consider allowed shapes for now
################################################################
# get basis functions for the constrained components -- ordinal monotone
if (capl > 0) {
	delta = varlist = NULL
	if (1 %in% shpsx[[1]] | 2 %in% shpsx[[1]]) {
#print ('1')
		del1 = makedelta(xmat[, 1], 1)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
	if (capl > 1) {
		for (i in 2:capl) {
			if (1 %in% shpsx[[i]] | 2 %in% shpsx[[i]]) {
				del2 = makedelta(xmat[, i], 1)$amat
		 		m2 = length(del2) / n
			} else {del2 = NULL; m2 = 0}
		 	delta = rbind(del1, del2)
			if (m1 > 0 | m2 > 0) {
				varlist = 1:(m1+m2)*0
				if (m1 > 0) {		
					varlist[1:m1] = var1
				}
				if (m2 > 0) {
					varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
				}
				var1 = varlist
				m1 = m1 + m2
		 		del1 = delta
			}
		}	
	}
	delta_om = delta
	varlist_om = varlist
}
# get basis functions for the constrained components -- smooth monotone
if (capl > 0) {
	delta = varlist = NULL
	if (9 %in% shpsx[[1]] | 10 %in% shpsx[[1]]) {
		del1 = makedelta(xmat[, 1], 9)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
	if (capl > 1) {
		for (i in 2:capl) {
			if (9 %in% shpsx[[i]] | 10 %in% shpsx[[i]]) {
				del2 = makedelta(xmat[, i], 9)$amat
		 		m2 = length(del2) / n
			} else {del2 = NULL; m2 = 0}
		 	delta = rbind(del1, del2)
			if (m1 > 0 | m2 > 0) {
				varlist = 1:(m1+m2)*0
				if (m1 > 0) {	
					varlist[1:m1] = var1
				}
				if (m2 > 0) {	
					varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
				}			
				var1 = varlist
				m1 = m1 + m2
		 		del1 = delta
			}
		}
	} 
	delta_sm = delta
	varlist_sm = varlist
}
# get basis functions for the constrained components -- ordinal convex
if (capl > 0) {
	delta = varlist = NULL
	if (3 %in% shpsx[[1]] | 4 %in% shpsx[[1]]) {
		del1 = makedelta(xmat[, 1], 3)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
	if (capl > 1) {
		for (i in 2:capl) {
			if (3 %in% shpsx[[i]] | 4 %in% shpsx[[i]]) {
				del2 = makedelta(xmat[, i], 3)$amat
		 		m2 = length(del2) / n
			}
		 	delta = rbind(del1, del2)
			if (m1 > 0 | m2 > 0) {
				varlist = 1:(m1+m2)*0
				if (m1 > 0) {	
					varlist[1:m1] = var1
				}
				if (m2 > 0) {	
					varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
				}
				var1 = varlist
				m1 = m1 + m2
		 		del1 = delta
			}
		}
	} 
	delta_ocv = delta
	varlist_ocv = varlist
}
# get basis functions for the constrained components -- smooth convex
if (capl > 0) {
	delta = varlist = NULL
	if (11 %in% shpsx[[1]] | 12 %in% shpsx[[1]]) {
		del1 = makedelta(xmat[, 1], 11)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1 
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
	if (capl > 1) {
		for (i in 2:capl) {
			if (11 %in% shpsx[[i]] | 12 %in% shpsx[[i]]) {
				del2 = makedelta(xmat[, i], 11)$amat
		 		m2 = length(del2) / n
			}
		 	delta = rbind(del1, del2)
			if (m1 > 0 | m2 > 0) {
				varlist = 1:(m1+m2)*0
				if (m1 > 0) {	
					varlist[1:m1] = var1
				}
				if (m2 > 0) {	
					varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
				}
				var1 = varlist
				m1 = m1 + m2
		 		del1 = delta
			}
		}
	}
	delta_scv = delta
	varlist_scv = varlist
}
if (capl > 0) {
	delta = varlist = NULL
# get basis functions for the constrained components -- ordinal increasing convex
	if (5 %in% shpsx[[1]] | 8 %in% shpsx[[1]]) {
		del1 = makedelta(xmat[, 1], 5)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
	if (capl > 1) {
		for (i in 2:capl) {
			if (5 %in% shpsx[[i]] | 8 %in% shpsx[[i]]) {
				del2 = makedelta(xmat[, i], 5)$amat
		 		m2 = length(del2) / n
			}
		 	delta = rbind(del1, del2)
			if (m1 > 0 | m2 > 0) {
				varlist = 1:(m1+m2)*0
				if (m1 > 0) {	
					varlist[1:m1] = var1
				}
				if (m2 > 0) {	
					varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
				}
				var1 = varlist
				m1 = m1 + m2
		 		del1 = delta
			}
		}
	} 
	delta_oincv = delta
	varlist_oincv = varlist
}
if (capl > 0) {
	delta = varlist = NULL
# get basis functions for the constrained components -- smooth increasing convex
	if (13 %in% shpsx[[1]] | 16 %in% shpsx[[1]]) {
		del1 = makedelta(xmat[, 1], 13)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
	if (capl > 1) {
		for (i in 2:capl) {
			if (13 %in% shpsx[[i]] | 16 %in% shpsx[[i]]) {
				del2 = makedelta(xmat[, i], 13)$amat
		 		m2 = length(del2) / n
			}
			delta = rbind(del1, del2)
			if (m1 > 0 | m2 > 0) {
				varlist = 1:(m1+m2)*0
				if (m1 > 0) {
					varlist[1:m1] = var1
				}
				if (m2 > 0) {
					varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
				}
				var1 = varlist
				m1 = m1 + m2
			 	del1 = delta
			}
		}
	} 	
	delta_sincv = delta
	varlist_sincv = varlist
}
if (capl > 0) {
	delta = varlist = NULL
# get basis functions for the constrained components -- ordinal decreasing convex
	if (6 %in% shpsx[[1]] | 7 %in% shpsx[[1]]) {
		del1 = makedelta(xmat[, 1], 6)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
	if (capl > 1) {
		for (i in 2:capl) {
			if (6 %in% shpsx[[i]] | 7 %in% shpsx[[i]]) {
				del2 = makedelta(xmat[, i], 6)$amat
		 		m2 = length(del2) / n
			}
		 	delta = rbind(del1, del2)
			if (m1 > 0 | m2 > 0) {
				varlist = 1:(m1+m2)*0
				if (m1 > 0) {
					varlist[1:m1] = var1
				}
				if (m2 > 0) {
					varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
				}
				var1 = varlist
				m1 = m1 + m2
		 		del1 = delta
			}
		}
	} 
	delta_odecv = delta
	varlist_odecv = varlist
}
if (capl > 0) {
	delta = varlist = NULL
# get basis functions for the constrained components -- smooth increasing concave
	if (14 %in% shpsx[[1]] | 15 %in% shpsx[[1]]) {
		del1 = makedelta(xmat[, 1], 14)$amat
		m1 = length(del1) / n
		var1 = 1:m1*0 + 1
	} else {del1 = NULL; m1 = 0; var1 = 0}
	delta = del1
	varlist = var1
	if (capl > 1) {
		for (i in 2:capl) {
			if (14 %in% shpsx[[i]] | 15 %in% shpsx[[i]]) {
				del2 = makedelta(xmat[, i], 14)$amat
		 		m2 = length(del2) / n
			}
		 	delta = rbind(del1, del2)
			if (m1 > 0 | m2 > 0) {
				varlist = 1:(m1+m2)*0
				if (m1 > 0) {
					varlist[1:m1] = var1
				}
				if (m2 > 0) {
					varlist[(m1+1):(m1+m2)] = (1:m2)*0+i
				}
				var1 = varlist
				m1 = m1 + m2
		 		del1 = delta
			}
		}
	}
	delta_sincc = delta
	varlist_sincc = varlist
}
	if (capk > 0) {
		zcat = NULL
		zvars = NULL	
		for (k in 1:capk){
			zk = zmat[, k]
			is_fac = zfacs[k]
			if (is_fac) {
				zkmat = model.matrix(~ factor(zk))[, -1, drop = FALSE]
				zvars = c(zvars, rep(k, ncol(zkmat)))
			} else {
				zkmat = zk
				zvars = c(zvars, k)
			}
			zcat = cbind(zcat, zkmat)
		}
	} 
	vzcat = NULL
	if (capkv > 0) {
		for (k in 1:capkv){
			vzk = vzmat[, k]
			is_fac = vzfacs[k]
			if (is_fac) {
				vzkmat = model.matrix(~ factor(vzk))[, -1, drop = FALSE]
			} else {
				vzkmat = vzk
			}
			vzcat = cbind(vzcat, vzkmat)
		}
	}
	vxcat = NULL
	if (caplv > 0) {
		for (l in 1:caplv) {
			xl = vxmat[, l]
			shpl = shpsvx[l]
			vxdd = t(makedelta(xl, shpl)$amat)
			if (shpl != 17) {
				#caplv = caplv - 1
				vxcat = cbind(vxcat, vxdd)
			} else if (shpl == 17) {
				capkv = capkv + 1
				vzcat = cbind(vzcat, vxdd)
			}
#print (dim(vxdd))
		}
	}
	if (capt > 0) {
#trcat: work the same way as zcat, include in bigmat as the final part
		trcat = NULL
		trvars = NULL	
		for (k in 1:capt){
			trk = trmat[, k]
			#trkmat = model.matrix(~ factor(trk))[, -1, drop = FALSE]
			trkmat = t(tree.fun(trk))
			trcat = cbind(trcat, trkmat)
			trvars = c(trvars, rep(k, ncol(trkmat)))
		}
	} 
	if (captv > 0) {
		vtrcat = NULL
		#vtrvars = NULL	
		for (k in 1:captv){
			vtrk = vtrmat[, k]
			#trkmat = model.matrix(~ factor(trk))[, -1, drop = FALSE]
			vtrkmat = t(tree.fun(vtrk))
			vtrcat = cbind(vtrcat, vtrkmat)
			#trvars = c(trvars, rep(k, ncol(trkmat)))
		}
	} 
# make the population of all models
	listAll = vector("list", capl + capk + capt)
	if (capl > 0) {
		listAll[1:capl] = shpsx
	}
	if (capk > 0) {
		listAll[(capl+1):(capl+capk)] = shpsz
	}
	if (capt > 0) {
		listAll[(capl+capk+1):(capl+capk+capt)] = shpst
	}
	popmat = as.matrix(expand.grid(listAll, stringsAsFactors = FALSE))
#ord = ncol(popmat):1
#pop0 = unname(popmat[, ord])
	npop = nrow(popmat)
	popmat = cbind(popmat, 1:npop*0)
	fitness = 1:npop*(-1)
## keep track of fitnesses already calculated
cat(paste("Evaluating the fitness of all models!", "\n"))
for (ipop in 1:npop) { 
#print (ipop)
	if (sum(popmat[ipop, ]) == 0) {
		llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), n = n, weights = weights, fml = family$family)				
		cic = llh + log(1 + 2 / (n - 1))	
		fitness[ipop] = -cic
		popmat[ipop,  (capl+capk+capt+1)] = -cic
		
  	} 
	delta = matrix(1:n*0+1, nrow = 1)
	if (capkv > 0 ) {
		delta = rbind(delta, t(vzcat))
	}
	if (capk > 0) {
  		if (sum(popmat[ipop, (capl+1):(capl+capk)]) > 0) {
  			#zuse = 1:st < 0
			zuse = 1:ncol(zcat) < 0
			for (i in 1:capk) {
				if (popmat[ipop, capl+i] == 1) {
					zuse[zvars == i] = TRUE
				}
			}
			delta = rbind(delta, t(zcat[, zuse]))
 			# delta = rbind(1:n*0+1, t(zcat[, zuse]))
 		}# else {delta = matrix(1:n*0+1, nrow = 1)}
 	} #else {delta = matrix(1:n*0+1, nrow = 1)}
	if (capt > 0) {
		if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
			truse = 1:ncol(trcat) < 0
			for (i in 1:capt) {
				if (popmat[ipop, capl+capk+i] == 2) {
					truse[trvars == i] = TRUE
				}
			}
			delta = rbind(delta, t(trcat[, truse]))
		} 
	}
	if (sum(popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13) > 0) {
 		usex = popmat[ipop, 1:capl] > 2 & popmat[ipop, 1:capl] < 5 | popmat[ipop, 1:capl] > 10 & popmat[ipop, 1:capl] < 13
 		delta = rbind(delta, t(xmat[, usex]))
 	}
#new:
if (caplv > 0) {
	ch = shpsvx %in% c(3, 4, 11, 12)
	if (any(ch)) {
		vxin = vxmat[, which(ch)]
		delta = rbind(delta, t(vxin))
	}
}
	np = dim(delta)[1]
 	if (capl > 0) {
		for (l in 1:capl) {
 			if (popmat[ipop, l] > 0) {
 				if (popmat[ipop, l] == 1) {#incr
#print ('1')
 					deladd = delta_om[varlist_om == l, ]
 				} else if (popmat[ipop, l] == 2) {#decr
#print ('2')
 					deladd = -delta_om[varlist_om == l, ]
 				} else if (popmat[ipop, l] == 3) {#conv
#print ('3')
 					deladd = delta_ocv[varlist_ocv == l, ]
 				} else if (popmat[ipop, l] == 4) {#conc
#print ('4')
 					deladd = -delta_ocv[varlist_ocv == l, ]
 				} else if (popmat[ipop, l] == 5) {#incr.conv
#print ('5')
 					deladd = delta_oincv[varlist_oincv == l, ]
 				} else if (popmat[ipop, l] == 8) {#decr.conc
#print ('8')
 					deladd = -delta_oincv[varlist_oincv == l, ]
 				} else if (popmat[ipop, l] == 6) {#decr.conv 
#print ('6')
 					deladd = delta_odecv[varlist_odecv == l, ]
 				} else if (popmat[ipop, l] == 7) {#incr.conc
#print ('7')
 					deladd = -delta_odecv[varlist_odecv == l, ]
 				} else if (popmat[ipop, l] == 9) {
#print ('9')
 					deladd = delta_sm[varlist_sm == l, ]
 				} else if (popmat[ipop, l] == 10) {
#print ('10')
 					deladd = -delta_sm[varlist_sm == l, ]
 				} else if (popmat[ipop, l] == 11) {
#print ('11')
 					deladd = delta_scv[varlist_scv == l, ]
 				} else if (popmat[ipop, l] == 12) {
#print ('12')
 					deladd = -delta_scv[varlist_scv == l, ]
 				} else if (popmat[ipop, l] == 13) {
#print ('13')	
 					deladd = delta_sincv[varlist_sincv == l, ]
 				} else if (popmat[ipop, l] == 16) {
#print ('16')
 					deladd = -delta_sincv[varlist_sincv == l, ]
 				} else if (popmat[ipop, l] == 14) {
#print ('14')
 					deladd = delta_sincc[varlist_sincc == l, ]
 				} else if (popmat[ipop, l] == 15) {
#print ('15')
 					deladd = -delta_sincc[varlist_sincc == l, ]
 				}
 				delta = rbind(delta, deladd)
#print (paste('dim: ', dim(delta)))
 			}	
 		} 
	}
	if (caplv > 0) {
		if (!is.null(vxcat)) {
			delta = rbind(delta, t(vxcat))
		}
	}
	if (capt > 0) {
		if (sum(popmat[ipop, (capl+capk+1):(capl+capk+capt)]) > 0) {
			truse = 1:ncol(trcat) < 0
			for (i in 1:capt) {
				if (popmat[ipop, capl+capk+i] == 1) {
					truse[trvars == i] = TRUE
				}
			}
			delta = rbind(delta, t(trcat[, truse]))
		} 
	}
	if (captv > 0) {
		delta = rbind(delta, t(vtrcat))
	}
	ans = irls(y, delta, np, nsim = 100, family = family, weights = weights, cpar = cpar) 	
	fitness[ipop] = -ans$cic
	popmat[ipop,  (capl+capk+capt+1)] = -ans$cic
}
ord = order(-fitness)
popmat = popmat[ord, ,drop = FALSE]
fitness = fitness[ord]
rslt = new.env()
rslt$pop = unname(popmat)
rslt$fitness = fitness
pop2 = matrix(0, nrow = npop, ncol = (capl+capk+capt))
for (i in 1:npop) {
	pop2[i, 1:capl] = apply(popmat[i, 1:capl, drop = FALSE], 2, ShapeToChar)
}
if (capk > 0) {
	for (i in 1:npop) {
		pop2[i, (capl+1):(capl+capk)] = apply(popmat[i, (capl+1):(capl+capk), drop = FALSE], 2, function(num) ShapeToChar(num, tag = "z"))
	}
}
if (capt > 0) {
	for (i in 1:npop) {
		pop2[i, (capl+capk+1):(capl+capk+capt)] = apply(popmat[i, (capl+capk+1):(capl+capk+capt), drop = FALSE], 2, function(num) ShapeToChar(num, tag = "tree"))
	}
}
#colnames(pop2) = c(xnms, znms, trnms)
#rslt$pop0 = pop0
pop2 = as.data.frame(pop2, stringsAsFactors = FALSE)
pop2 = cbind(pop2, fitness)
rslt$pop2 = pop2
#rslt$top = pop2[1,]
rslt$GA = FALSE
rslt$vzcat = vzcat
#class(rslt) = "cgam"
return (rslt)
}

##############################
#tranform shapes back to char
##############################
ShapeToChar = function(shp, tag = "x") {
	#if (max(shp) > 16 | min(shp) < 0) {
	#	stop ('No such a shape! A shape value must be between 0 and 16.')
	#}
	if (tag == "x") {
		if (shp == 0) {
			shp = 18
		}
		switch(shp, 
			cs1 = {ch = 'incr'},
			cs2 = {ch = 'decr'},
			cs3 = {ch = 'conv'},
			cs4 = {ch = 'conc'},
			cs5 = {ch = 'incr.conv'},
			cs6 = {ch = 'decr.conv'},
			cs7 = {ch = 'incr.conc'},
			cs8 = {ch = 'decr.conc'},
			cs9 = {ch = 's.incr'},
			cs10 = {ch = 's.decr'},
			cs11 = {ch = 's.conv'},
			cs12 = {ch = 's.conc'},
			cs13 = {ch = 's.incr.conv'},
			cs14 = {ch = 's.incr.conc'},
			cs15 = {ch = 's.decr.conv'},
			cs16 = {ch = 's.decr.conc'},
			cs17 = {ch = 's'},
			cs18 = {ch = 'flat'}
			#{print ('No such a shape')}
		)
	} else if (tag == "z") {
		if (shp == 0) {
			shp = 2
		}
		switch(shp, 
			cs1 = {ch = 'in'},
			cs2 = {ch = 'out'}
		)
	} else if (tag == "tree") {
		if (shp == 0) {
			shp = 3
		}
		switch(shp, 
			cs1 = {ch = 'tree'},
			cs2 = {ch = 'unordered'},
			cs3 = {ch = 'out'}
		)
	}
	return (ch)
}

#
CharToShape = function(ch) {
	shp = NULL	
	if (ch == 'flat') {
		shp = 0
	} else if (ch == 'incr') {
		shp = 1
	} else if (ch == 'decr') {
		shp = 2
	} else if (ch == 'conv') {
		shp = 3
	} else if (ch == 'conc') {
		shp = 4
	} else if (ch == 'incr.conv') {
		shp = 5
	} else if (ch == 'decr.conv') {
		shp = 6
	} else if (ch == 'incr.conc') {
		shp = 7
	} else if (ch == 'decr.conc') {
		shp = 8
	} else if (ch == 's.incr') {
		shp = 9
	} else if (ch == 's.decr') {
		shp = 10
	} else if (ch == 's.conv') {
		shp = 11
	} else if (ch == 's.conc') {
		shp = 12
	} else if (ch == 's.incr.conv') {
		shp = 13
	} else if (ch == 's.incr.conc') {
		shp = 14
	} else if (ch == 's.decr.conv') {
		shp = 15
	} else if (ch == 's.decr.conc') {
		shp = 16
	} else if (ch == 's') {
		shp = 17
	} else {
		stop ('shape not defined!')	
	}
	return (shp)
}

############################################################
#plot fitness vs iterations                                #
############################################################
plot.shapeselect = function(x,...)
{
  object = x
  if (!object$GA) {
  	stop("plot.shapeselect won't work for objects not fitted by genetic algorithm ...")
  }
  par(mfrow = c(1, 1))
  #iters = 1:100
  cex.points = 0.7
  col = c("slategrey","mediumorchid4")
  pch = c(16, 1)
  lty = c(1, 2)
  mx = object$mxf
  mn = object$mnf
  iters = 1:length(mx)
  ylim = c(min(mn), max(mx) + (max(mx) - min(mn)) * .01) 
  plot(iters, mx, type = "n", ylim = ylim, xlab = "Generation", ylab = "Fitness value")
  grid() 
  legend("bottomright", bty = "n", legend = c("Best", "Mean"), col = col, pch = pch, lty = lty, pt.cex = cex.points, inset = 0.01)  
  points(iters, mx, type = "o", pch = pch[1], lty = lty[1], col = col[1], cex = cex.points)
  for (i in 1:length(mx)) {	
	if (i < length(mx)) {
  		points(iters[i], mn[i], type = "o", pch = pch[2], lty = lty[2], col = col[2], cex = cex.points)
		segments(iters[i], mn[i], iters[i+1], mn[i+1], col = "mediumorchid4", lty = lty[2])
	} else {
		points(iters[i], mn[i], type = "o", pch = "*", col = col[2], cex = 2.1)
		Sys.sleep(.5)
		points(iters[i], mn[i], type = "o", pch = "*", col = "white", cex = 2.1)
		Sys.sleep(.5)
		points(iters[i], mn[i], type = "o", pch = "*", col = col[2], cex = 2.5)
		Sys.sleep(.5)
		points(iters[i], mn[i], type = "o", pch = "*", col = "white", cex = 2.5)
		Sys.sleep(.5)
		points(iters[i], mn[i], type = "o", pch = "*", col =  col[2], cex = 2.9)
		#Sys.sleep(.5)
		#points(iters[i], mn[i], type = "o", pch = "*", col = "white", cex = 2.9)
		#Sys.sleep(.5)
		#points(iters[i], mn[i], type = "o", pch = "*", col =  col[2], cex = 3.3)		
	}
	Sys.sleep(.03)
  }
}

########################################################################
# iteratively re-weighted least squares -- binary likelihood
########################################################################
irls = function(y, bigmat, np, nsim = 100, family = gaussian(), weights = NULL, cpar = 1.2) {
	linkfun = family$linkfun
	cicfamily = CicFamily(family)
	llh.fun = cicfamily$llh.fun
	etahat.fun = cicfamily$etahat.fun
	gr.fun = cicfamily$gr.fun
	wt.fun = cicfamily$wt.fun
	zvec.fun = cicfamily$zvec.fun
	muhat.fun = cicfamily$muhat.fun
	ysim.fun = cicfamily$ysim.fun
	deriv.fun = cicfamily$deriv.fun
	dev.fun = cicfamily$dev.fun 
	if (family$family == "binomial" | family$family == "poisson") {
		wt.iter = TRUE
	} else {wt.iter = FALSE}
	n = length(y)
	if (is.null(weights)) {
		weights = 1:n*0 + 1
	}
	sm = 1e-6
	m = dim(bigmat)[1] - np
# new: initialize cvec
	cvec = NULL
	if (wt.iter) {
		etahat = etahat.fun(n, y, fml = family$family)
		gr = gr.fun(y, etahat, weights, fml = family$family)  
		wt = wt.fun(etahat, n, weights, fml = family$family)     
		cvec = wt * etahat - gr
	} else {wt = wt.fun(etahat, n, weights, fml = family$family)}
	zvec = zvec.fun(cvec, wt, y, fml = family$family)
        gmat = t(bigmat %*% sqrt(diag(wt)))
	if (m > 0) {
#np always >= 1	
		dsend = gmat[, (np + 1):(np + m), drop = FALSE]
        	zsend = gmat[, 1:np, drop = FALSE] 
		ans = coneB(zvec, t(dsend), zsend, msg = FALSE)
#ans = hingep(zvec, t(dsend), zsend)
		coef = ans$coefs
		etahat = t(bigmat) %*% coef 
		muhat = muhat.fun(etahat, fml = family$family)
		if (wt.iter) {
			#muhat = muhat.fun(etahat, fml = family$family)
			diff = 1
			if (family$family == "binomial") {
				mdiff = abs(max(muhat) - 1) > sm	
			} else {mdiff = TRUE}
			nrep = 0
##########
#iterate!#	
##########
			while (diff > sm & mdiff & nrep < n^2) {
				oldmu = muhat	
				nrep = nrep + 1
				gr = gr.fun(y, etahat, weights, fml = family$family)	
				wt = wt.fun(etahat, n, weights, fml = family$family) 
				cvec = wt * etahat - gr
				zvec = zvec.fun(cvec, wt, y, fml = family$family)						
				gmat = t(bigmat %*% sqrt(diag(wt)))
				dsend = gmat[, (np + 1):(np + m), drop = FALSE]
        			zsend = gmat[, 1:np, drop = FALSE] 
				ans = coneB(zvec, t(dsend), zsend, msg = FALSE)
#ans = hingep(zvec, t(dsend), zsend)
				coef = ans$coefs
				etahat = t(bigmat) %*% coef
				muhat = muhat.fun(etahat, fml = family$family)
				diff = mean((muhat - oldmu)^2)	
				mdiff = abs(max(muhat) - 1)
				if (family$family == "binomial") {
					mdiff = abs(max(muhat) - 1) > sm	
				} else {mdiff = TRUE}
			}
		}	
	} else {
		prior.w = weights
		vmat = t(bigmat[1:np, , drop = FALSE])	
		w = diag(prior.w)
		coef = solve(t(vmat) %*% w %*% vmat) %*% t(vmat) %*% w %*% y
		etahat = vmat %*% coef
		muhat = muhat.fun(etahat, fml = family$family)			
		if (wt.iter) {
			nrep = 0
			muhat = mean(y) + 1:n*0
			etahat = linkfun(muhat)
			diff = 1
			if (family$family == "binomial") {
				mdiff = abs(max(muhat) - 1) > sm	
			} else {mdiff = TRUE}
			while (diff > sm & mdiff & nrep < n^2) {
				nrep = nrep + 1
				oldmu = muhat
				zhat = etahat + (y - muhat) * deriv.fun(muhat, fml = family$family)				
				#w <- diag(as.vector(prior.w / deriv.fun(muhat)))		
				w = diag(as.vector(prior.w * (deriv.fun(muhat, fml = family$family))^(-1)))
				coef = solve(t(vmat) %*% w %*% vmat) %*% t(vmat) %*% w %*% zhat
				etahat = vmat %*% coef
				muhat = muhat.fun(etahat, fml = family$family)		
				diff = mean((muhat - oldmu)^2)	
				mdiff = abs(max(muhat) - 1)
				if (family$family == "binomial") {
					mdiff = abs(max(muhat) - 1) > sm	
				} else {mdiff = TRUE}
			}
		}  
	}
	#muhat = muhat.fun(etahat, fml = family$family)
	llh = llh.fun(y, muhat, etahat, n, weights, fml = family$family)
	mukeep = muhat
	coefkeep = coef
	if (np > 0) {zcoefs = coefkeep[1:np]}
### get edf0
	if (m > 0) {
		dfs = 1:nsim*0
		sm = 1e-5
		if (family$family == "poisson") {		
			mu0 = mean(y)
			eta0 = log(mu0)
		} else {mu0 = NULL}
		for (isim in 1:nsim) {
#print (isim)
			ysim = ysim.fun(n, mu0 = mu0, fml = family$family)
		  	if (wt.iter) {
				etahat = etahat.fun(n, ysim, fml = family$family)
				gr = gr.fun(ysim, etahat, weights, fml = family$family)
				wt = wt.fun(etahat, n, weights, fml = family$family)
				cvec = wt * etahat - gr
			} else {wt = wt.fun(etahat, n, weights, fml = family$family)}
			zvec = zvec.fun(cvec, wt, ysim, fml = family$family)
            		gmat = t(bigmat %*% sqrt(diag(wt)))
           		dsend = gmat[, (np + 1):(np + m), drop = FALSE]
            		zsend = gmat[, 1:np, drop = FALSE]
			ans = try(coneB(zvec, t(dsend), zsend, msg = FALSE))
#ans = coneB(zvec, t(dsend), zsend)
#save(zvec,file='zv1.Rda')
#save(dsend,file='ds1.Rda')
#save(zsend,file='zs1.Rda')
#ans = hingep(zvec, t(dsend), zsend)
			if (class(ans) == "try-error") next 
			if (wt.iter) {
				etahat = t(bigmat) %*% ans$coefs
				muhat = muhat.fun(etahat, fml = family$family)
				diff = 1
				if (family$family == "binomial") {
					mdiff = abs(max(muhat) - 1) > sm	
				} else {mdiff = TRUE}
##########
#iterate!#
##########
				nrep = 0 
				while (diff > sm & nrep < n^2 & mdiff > sm) {
					nrep = nrep + 1
					oldmu = muhat	
					gr = gr.fun(ysim, etahat, weights, fml = family$family)
					wt = wt.fun(etahat, n, weights, fml = family$family)
					cvec = wt * etahat - gr
					#zvec <- cvec / sqrt(wt)
					zvec = zvec.fun(cvec, wt, y, fml = family$family)
					gmat = t(bigmat %*% sqrt(diag(wt)))
					dsend = gmat[, (np + 1):(np + m), drop = FALSE]
					zsend = gmat[, 1:np, drop = FALSE]
					ans = try(coneB(zvec, t(dsend), zsend, msg = FALSE))
#ans = coneB(zvec, t(dsend), zsend)
#save(zvec,file='zv2.Rda')
#save(dsend,file='ds2.Rda')
#save(zsend,file='zs2.Rda')
#ans = hingep(zvec, t(dsend), zsend)
					if (class(ans) == "try-error") next 
					etahat = t(bigmat) %*% ans$coefs
					muhat = muhat.fun(etahat, fml = family$family)
					diff = mean((muhat - oldmu)^2)	
					if (family$family == "binomial") {
						mdiff = abs(max(muhat) - 1) > sm	
					} else {mdiff = TRUE}
				}
			}
	    		dfs[isim] = sum(abs(ans$coefs) > 0)
		}
		dfmean = mean(dfs)
	} else {dfmean = np}
	rslt = new.env()
	rslt$edf0 = dfmean
	rslt$llh = llh
	rslt$fhat = mukeep
	rslt$coefs = coefkeep
	#cdose$cic = llh+log(1+2*dfmean/(n-np-1.5*(dfmean-np)))
	if ((n - np - cpar * (dfmean - np)) <= 0) {
		rslt$cic <- llh + log(1 + 2 * dfmean / (dfmean - np))
	} else {
        	rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - cpar * (dfmean - np)))
	}	
	rslt
}


