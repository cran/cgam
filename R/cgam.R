cgam <- function(formula, nsim = 1e+2, family = gaussian(), data = NULL, weights = NULL)
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
  mf <- mf[c(1L,m)]
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  if (family$family == "binomial") {
	if (class(y) == "factor") {
		y = ifelse(y == levels(y)[1], 0, 1)
	}
  }
  shapes1 <- NULL
  shapes2 <- NULL
  xmat <- NULL
  tr <- NULL
  umb <- NULL
  tree.delta <- NULL
  umbrella.delta <- NULL 
  tid1 <- NULL; tid2 <- NULL; tpos2 <- 0
  uid1 <- NULL; uid2 <- NULL; upos2 <- 0
  for (i in 2: ncol(mf)) {
    if (is.numeric(attributes(mf[,i])$shape)) {
       shapes1 <- c(shapes1, attributes(mf[,i])$shape)
       xmat <- cbind(xmat, mf[,i])
     }
    if (is.character(attributes(mf[,i])$shape)) {
       shapes2 <- c(shapes2, attributes(mf[,i])$shape)
       if (attributes(mf[,i])$shape == "tree") {
		tree.delta <- rbind(tree.delta, tree.fun(mf[,i]))
		tpos1 <- tpos2 + 1
		tpos2 <- tpos2 + nrow(tree.fun(mf[,i]))
		tid1 <- c(tid1, tpos1)
		tid2 <- c(tid2, tpos2)
		#tree.delta <- tree.fun(mf[,i])
		tr <- cbind(tr, mf[,i])
	}
       if (attributes(mf[,i])$shape == "umbrella") {
		#umbrella.delta <- umbrella.fun(mf[,i])
		umbrella.delta <- rbind(umbrella.delta, umbrella.fun(mf[,i]))
		upos1 <- upos2 + 1
		upos2 <- upos2 + nrow(umbrella.fun(mf[,i]))
		uid1 <- c(uid1, upos1)
		uid2 <- c(uid2, upos2)
		umb <- cbind(umb, mf[,i])	
	}
     }
  }
  shapes <- c(shapes1, shapes2) #
  attr(xmat, "shape") <- shapes1
  zmat <- NULL
  zid <- NULL; zid0 <- NULL; nzid <- NULL; dist <- 0; nzid0 <- 0; znames <- NULL
  for (i in 2: ncol(mf)) {
    if (is.null(attributes(mf[,i])$shape)) {
    #if (!all(round(diff(mf[, i]), 8) == 0)) {
      if (!is.matrix(mf[,i])) {
        #new code to make zmat
        if (is.factor(mf[,i])) {
          nlvs <- length(attributes(mf[,i])$levels)
          zid0 <- i + 0:(nlvs-2) + dist
          zid <- c(zid, zid0)
          dist <- nlvs - 2
          nzid0 <- length(zid0)
          nzid <- c(nzid, nzid0)
        } else {
	    zid <- c(zid, i)
        }
      zmat <- model.matrix(mt, mf)[,zid]
#print (zmat)
      } else {
          zmat <- cbind(zmat, mf[,i])
	  zid <- c(zid, 1:ncol(mf[,i]))
          nzid <- length(zid)
	  znames <- c(znames, paste(as.vector(attributes(mt)$term.labels)[i-1], 1:ncol(mf[,i]), sep = ""))
	  colnames(zmat) <- znames
#print (zmat)
      }
    }
    #}
  }
  if (is.matrix(zmat) & !is.null(zmat)){
    nzmat <- zmat
    if (qr(nzmat)$rank != ncol(nzmat)) {
       stop("zmat should be full column rank!")
    }
    mat_cols <- ncol(nzmat)
    mat_rows <- nrow(nzmat)
    mat_rm <- NULL
    for (i in 1:mat_cols) {
  #     #if (all.equal(diff(nzmat[, i]), rep(0, mat_rows-1)) == TRUE) {
       if (all(round(diff(nzmat[, i]), 8) == 0)) {
              mat_rm <- c(mat_rm, i)
       }
    }
    if (!is.null(mat_rm)) {
            nzmat <- nzmat[, -mat_rm, drop = FALSE]
    }
    zmat <- nzmat
  }
#print (zmat)
  if (family$family == "binomial"|family$family == "poisson") {
     wt.iter = TRUE
  } else {wt.iter = FALSE}
  if (is.null(shapes1) & is.null(shapes2)) {
    nsim <- 0
  }
  ans <- cgam.fit(y = y, xmat = xmat, zmat = zmat, shapes = shapes1, nsim = nsim, family = family, wt.iter = wt.iter, umbrella.delta = umbrella.delta, tree.delta = tree.delta, weights = weights)
  if (!is.null(uid1) & !is.null(uid2)) {
    uid1 <- uid1 + ans$d0 + ans$capm
    uid2 <- uid2 + ans$d0 + ans$capm
  }
  if (!is.null(tid1) & !is.null(tid2)) {
    tid1 <- tid1 + ans$d0 + ans$capm + ans$capu
    tid2 <- tid2 + ans$d0 + ans$capm + ans$capu 
  }
  rslt <- list(vhat = ans$vhat, etahat = ans$etahat, muhat = ans$muhat, vcoefs = ans$vcoefs, xcoefs = ans$xcoefs, zcoefs = ans$zcoefs, ucoefs = ans$ucoefs, tcoefs = ans$tcoefs, coefs = ans$coefs, cic = ans$cic, d0 = ans$d0, edf0 = ans$edf0, etacomps = ans$etacomps, xmat = xmat, zmat = zmat, tr = tr, umb = umb, tree.delta = tree.delta, umbrella.delta = umbrella.delta, bigmat = ans$bigmat, shapes = shapes, wt = ans$wt, wt.iter = ans$wt.iter, family = ans$family, SSE0 = ans$sse0, SSE1 = ans$sse1, pvals.beta = ans$pvals.beta, se.beta = ans$se.beta, null_df = ans$df.null, df = ans$df, resid_df_obs = ans$resid_df_obs, null_deviance = ans$dev.null, deviance = ans$dev, tms = mt, capm = ans$capm, capk = ans$capk, capt = ans$capt, capu = ans$capu, xid1 = ans$xid1, xid2 = ans$xid2, tid1 = tid1, tid2 = tid2, uid1 = uid1, uid2 = uid2, zid = zid, nzid = nzid, nsim = nsim)
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


cgam.fit <- function(y, xmat, zmat, shapes, nsim, family = gaussian(), wt.iter = FALSE, umbrella.delta = NULL, tree.delta = NULL, weights = NULL) {
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
	sm <- 5e-6 
	capl <- length(xmat) / n
	if (capl < 1) {capl <- 0}
	if (round(capl,8) != round(capl,1)) {stop ("Incompatible dimensions for xmat!")}
	capk <- length(zmat) / n
	if (capk < 1) {capk <- 0}
	if (round(capk,8) != round(capk,1)) {stop ("Incompatible dimensions for zmat!")}

####################################################
#get basis functions for the constrained components#
####################################################

	delta = NULL
	varlist = NULL
	xid1 = NULL; xid2 = NULL; xpos2 <- 0
	if (capl > 0) {
		del1 <- makedelta(xmat[,1], shapes[1])
        	m1 <- length(del1) / n
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
	        		del2 <- makedelta(xmat[,i], shapes[i])
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

#########################
# make constraint matrix#
#########################

		if (sum(shapes > 2 & shapes < 5) > 0 & capk > 0) {
        		bigmat <- rbind(1:n*0 + 1, t(zmat), t(xmat[,shapes > 2 & shapes < 5]), delta)
        		np <- 1 + capk + sum(shapes > 2 & shapes < 5)
		} else if (sum(shapes > 2 & shapes < 5) > 0 & capk == 0) {
			bigmat <- rbind(1:n*0 + 1, t(xmat[, shapes > 2 & shapes < 5]), delta)
			np <- 1 + sum(shapes > 2 & shapes < 5)
		} else if (sum(shapes > 2 & shapes < 5) == 0 & capk > 0) {
			bigmat <- rbind(1:n*0 + 1, t(zmat), delta)
			np <- 1 + capk
		} else if (sum(shapes > 2 & shapes < 5) == 0 & capk == 0) {
			bigmat <- rbind(1:n*0 + 1, delta)
			np <- 1
		} else {
        	  	print ("error in capk,shapes!")
        	  }
	capm <- length(delta) / n
	} else {
	  	if (capk > 0) {
          		bigmat <- rbind(1:n*0 + 1, t(zmat))
          		capm <- 0
          		np <- 1 + capk
            	} else { bigmat <- matrix(1:n*0 + 1, nrow = 1); capm <- 0; np <- 1}
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
		if (capl + capu + capt > 0) { 
			if (wt.iter) {
				etahat <- etahat.fun(n, y)
				gr <- gr.fun(y, etahat)  
				wt <- wt.fun(etahat, n, weights)     
				cvec <- wt * etahat - gr
			} else {wt <- 1:n*0 + 1}
			  	zvec <- zvec.fun(cvec, wt, y)
        		  	gmat <- t(bigmat %*% sqrt(diag(wt)))
			  	dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
               		  	zsend <- gmat[, 1:np, drop = FALSE] 
				ans <- coneB(zvec, t(dsend), zsend)
			  	etahat <- t(bigmat) %*% ans$coefs
			  	if (wt.iter) {
					muhat <- muhat.fun(etahat)
			  		diff <- 1
					nrep <- 0
##########
#iterate!#	
##########
			  		while (diff > sm & nrep < n^2){
						oldmu <- muhat	
						nrep <- nrep + 1
						gr <- gr.fun(y, etahat)	
						wt <- wt.fun(etahat, n, weights) 
						cvec <- wt * etahat - gr
						zvec <- cvec / sqrt(wt)
						gmat <- t(bigmat %*% sqrt(diag(wt)))
						dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
        	  				zsend <- gmat[, 1:np, drop = FALSE]
						ans <- coneB(zvec, t(dsend), zsend)
						etahat <- t(bigmat) %*% ans$coefs
						muhat <- muhat.fun(etahat)
						diff <- mean((muhat - oldmu)^2)	
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
				if (capm > 0) {
					xcoefs <- coefskeep[(np + 1):(np + capm)]
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
					dcoefs <- coefskeep[(np + 1):(capm + np)]	
	
#####################################################
#thvecs is f(x), where x has one of the eight shapes#
#####################################################

					thvecs <- matrix(nrow = capl, ncol = n)
	    				ncon <- 1
	    				for (i in 1:capl) {
	    	  				thvecs[i,] <- t(delta[varlist == i,]) %*% dcoefs[varlist == i]
	    	  				if (shapes[i] > 2 & shapes[i] < 5) { 
            		    				ncon <- ncon + 1
            		    				#thvecs[i,] <- thvecs[i,] + zcoefs[ncon] * xmat[,i]
							thvecs[i,] <- thvecs[i,] + vcoefs[capk + ncon] * xmat[,i]
            	  				}
	    				}
				}
				thvecs_ut <- NULL
				if (capu + capt > 0) {
					thvecs_ut <- t(delta_ut) %*% coefskeep[(np + 1 + capm):(np + capm + capu + capt)]
				}
				if (!is.null(thvecs_ut)) {
					thvecs <- rbind(thvecs, t(thvecs_ut))
				}
				llh <-  llh.fun(y, muhat, etahat, n, weights)
	  			etakeep <- etahat
				muhatkeep <- muhat.fun(etakeep)
				wtkeep <- wt
				df_obs <- sum(abs(coefskeep) > 0)
				if (family$family == "poisson") {		
					mu0 <- mean(y)
					eta0 <- log(mu0)
				}
		} else if (capk > 0 & (capl + capt + capu) == 0) {
			if (is.null(weights)) {
				weights <- 1:n*0 + 1
			}
			prior.w <- weights
			vmat <- t(bigmat[1:np, , drop = FALSE])
			if (wt.iter) {
				diff <- 1
				muhat <- mean(y) + 1:n*0
				etahat <- linkfun(muhat)		
				while (diff > sm) {
					oldmu <- muhat
					zhat <- etahat + (y - muhat) * deriv.fun(muhat)				
					w <- diag(as.vector(prior.w / deriv.fun(muhat)))			
					b <- solve(t(vmat) %*% w %*% vmat) %*% t(vmat) %*% w %*% zhat
					etahat <- vmat %*% b
					muhat <- muhat.fun(etahat)		
					diff <- mean((muhat - oldmu)^2)	
				}
				se.beta <-  sqrt(diag(solve(t(vmat) %*% w %*% vmat)))
				zstat <- b / se.beta
				pvals.beta <-  1 - pchisq(zstat^2, df = 1)
			} else {
				w <- diag(prior.w)
				b <- solve(t(vmat) %*% w %*% vmat) %*% t(vmat) %*% w %*% y
				etahat <- vmat %*% b
				muhat <- muhat.fun(etahat)	
				sdhat2 <- sum(prior.w * (y - muhat)^2) / (n - np)
				se.beta <-  sqrt(diag(solve(t(vmat) %*% w %*% vmat) * sdhat2))
				tstat <- b / se.beta
				pvals.beta <-  (1 - pt(abs(tstat), df = n - np)) * 2 
			}
			llh <-  llh.fun(y, muhat, etahat, n, weights)
			df_obs <- np
			dfmean <- np
			rslt <- new.env()
			rslt$family <- family 
			rslt$wt.iter <- wt.iter 
			rslt$wt <- w
			rslt$bigmat <- bigmat
			rslt$etahat <- etahat
			rslt$muhat <- muhat
			rslt$d0 <- np
			rslt$capm <- capm
			rslt$capk <- capk
			rslt$capu <- capu
			rslt$capt <- capt
			rslt$xid1 <- xid1
			rslt$xid2 <- xid2 
			rslt$dfmean <- dfmean	
			rslt$edf0 <- dfmean 
			#rslt$llh <- llh 
			#if (nsim > 0) {
				rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - 1.5 * (dfmean - np)))
			#}
			rslt$zcoefs <- b
			rslt$se.beta <- se.beta 
			rslt$pvals.beta <- pvals.beta 
			rslt$dev <- dev.fun(y, muhat, etahat, weights)$dev
			rslt$dev.null <- dev.fun(y, muhat, etahat, weights)$dev.null
			rslt$df <- n - np 
			rslt$df.null <- n - 1
			rslt$resid_df_obs <- n - np - 1.5 * (df_obs - np)
			rslt$vhat <- etahat 			
			return (rslt) 
			
		}

##########
#get edf0#
##########
		  if (capl + capu + capt > 0 & nsim > 0) {
	  		dfs <- 1:nsim
	  		for (isim in 1:nsim) {
				#set.seed(123)
	    			ysim <- ysim.fun(n, mu0)
		  		if (wt.iter) {
					etahat <- etahat.fun(n, ysim)
					gr <- gr.fun(ysim, etahat)
					wt <- wt.fun(etahat, n, weights)
					cvec <- wt * etahat - gr
				} else {wt <- 1:n*0 + 1}
					zvec <- zvec.fun(cvec, wt, ysim)
            				gmat <- t(bigmat %*% sqrt(diag(wt)))
           				dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
            				zsend <- gmat[, 1:np, drop = FALSE]
					ans <- try(coneB(zvec, t(dsend), zsend))
					if (class(ans) == "try-error") next 
				if (wt.iter) {
						etahat <- t(bigmat) %*% ans$coefs
						muhat <- muhat.fun(etahat)
						diff <- 1
##########
#iterate!#
##########
						nrep <- 0 
						while (diff > sm & nrep < n^2) {	
							nrep <- nrep + 1
							oldmu <- muhat	
							gr <- gr.fun(ysim, etahat)
							wt <- wt.fun(etahat, n, weights)
							cvec <- wt * etahat - gr
							zvec <- cvec / sqrt(wt)
							gmat <- t(bigmat %*% sqrt(diag(wt)))
							dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
							zsend <- gmat[, 1:np, drop = FALSE]
							ans <- try(coneB(zvec, t(dsend), zsend))
							if (class(ans) == "try-error") next 
							etahat <- t(bigmat) %*% ans$coefs
							muhat <- muhat.fun(etahat)
							diff <- mean((muhat - oldmu)^2)	
						}
				}
	    			dfs[isim] <- sum(abs(ans$coefs) > 0)
	  		}
	  		dfmean <- mean(dfs)
		   } else if (capl + capu + capt > 0 & nsim == 0) {
			dfmean <- NULL
		   } 
			#else {dfmean <- np} 
		   

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
	if(capl + capu + capt > 0){
		xid1 <- xid1 + np
		xid2 <- xid2 + np
		rslt <- new.env()		
		rslt$family <- family 
		rslt$wt.iter <- wt.iter
		rslt$wt <- wtkeep
		rslt$bigmat <- bigmat  
		rslt$etahat <- etakeep
		rslt$muhat <- muhatkeep 
		rslt$d0 <- np
		rslt$capm <- capm
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
        		rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - 1.5 * (dfmean - np)))
		} else {rslt$cic <- NULL}
		#if (nsim > 0) {
			#rslt$edf0 <- dfmean - np
		#}
		rslt$edf0 <- dfmean 
		rslt$etacomps <- thvecs
		vmat <- t(bigmat[1:np, , drop = FALSE])
		if (is.null(weights)) {
			weights <- 1:n*0 + 1
		}
		prior.w <- weights
		w <- diag(as.vector(prior.w / deriv.fun(muhatkeep)))

###############################################
#the case capk = 0 and capk >= 0 are combined:#
###############################################

		vhat <- vmat %*% vcoefs	
		sse1 <- sum(prior.w * (y - yhat)^2)
		sse0 <- sum(prior.w * (y - vhat)^2)
		if ((n - np - 1.5 * (df_obs - np)) <= 0) {
			sdhat2 <- sse1
		} else {
			sdhat2 <- sse1 / (n - np - 1.5 * (df_obs - np))
		}
		if (wt.iter) {
			se2 <- solve(t(vmat) %*% w %*% vmat)
		} else {
			se2 <-  solve(t(vmat) %*% diag(prior.w) %*% vmat) * sdhat2
		}		 			
		se.beta <- 1:(capk + 1)*0
		tstat <- 1:(capk + 1)*0
		pvals.beta <- 1:(capk + 1)*0
		rslt$zcoefs <- zcoefs 
#if ((n - np - 1.5 * (df_obs - np)) < 0) {
#	rslt$se.beta <- NULL
#	rslt$pvals.beta <- NULL
#} else {
		for (i in 1:(capk + 1)) {
			se.beta[i] <- sqrt(se2[i,i])
			tstat[i] <- zcoefs[i] / se.beta[i]
			pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  n - np - 1.5 * (df_obs - np))) 
		}
		rslt$se.beta <- se.beta
		rslt$pvals.beta <- pvals.beta
#}	
		rslt$sse1 <- sse1
		rslt$sse0 <- sse0
		rslt$dev <- dev.fun(y, muhatkeep, etakeep, weights)$dev
		rslt$dev.null <- dev.fun(y, muhatkeep, etakeep, weights)$dev.null
		rslt$df <- n - np 
		rslt$df.null <- n - 1
		rslt$resid_df_obs <- n - np - 1.5 * (df_obs - np)
		rslt$df_obs <- df_obs
		rslt$vhat <- vhat
		return (rslt)
	}
}


CicFamily <- function(object,...)UseMethod("CicFamily")
CicFamily <- function(object) {
  llh.fun <- function(y, muhat = NULL, etahat = NULL, n = NULL, weights = NULL, fml = object$family){
    sm <- 5e-6
    if (is.null(weights)) {
	weights <- 1:n*0 + 1
    }
    w <- weights
    if (fml == "poisson") {
      llh <- 2 * sum(w * (muhat - y * etahat)) / n
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
          llh <- log(sum(w * (y - etahat)^2)) - sum(log(w)) / n
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

  gr.fun <- function(y, etahat = NULL, fml = object$family){
    n <- length(y)
    if (fml == "poisson") {
       gr <- exp(etahat) -  y 
    }
    if (fml == "binomial") {
       if (all(etahat == 0)) { 
         gr <- 1/2 - y
       } else {
	   gr <- 1:n*0
	   for (i in 1:n) {
	     if (etahat[i] > 100) {
		 gr[i] <- 1 - y[i] 
	     } else {gr[i] <- exp(etahat[i]) / (1 + exp(etahat[i])) - y[i]}
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
      wt <- w * exp(etahat)
    }
    if (fml == "binomial") {
      if (all(etahat == 0)){
        wt <- 1:n*0 + 1/4
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
      wt <- (1:n*0 + 1) / w 
    }
    wt <- as.vector(wt)
    wt 
  }

  zvec.fun <- function(cvec = NULL, wt = NULL, y, sm = 5e-6, fml = object$family) {
    n <- length(y)
    if (fml == "gaussian") {
      zvec <- y
    }
    if (fml == "poisson") {
      zvec <- cvec / wt 
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
  sm <- 5e-6
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
          if (y[i] == 0) {
            dev <- dev + 2 * w[i] * log(w[i] / (w[i] - w[i] * muhat[i]))
          } else if (y[i] == 1) {
              dev <- dev + 2 * w[i] * log(w[i] / (w[i] * muhat[i]))
          } else if (0 < y[i] & y[i] < 1) {
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
	zhat <- etahat0 + (y - muhat0) * deriv.fun(muhat0)		
	wmat <- diag(as.vector(w / deriv.fun(muhat0)))			
	b <- solve(t(vmat) %*% wmat %*% vmat) %*% t(vmat) %*% wmat %*% zhat
	etahat0 <- vmat %*% b
	muhat0 <- muhat.fun(etahat0)		
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
          if (y[i] == 0) {
            dev.null <- dev.null + 2 * w[i] * log(w[i] / (w[i] - w[i] * muhat0[i]))
          } else if (y[i] == 1) {
              dev.null <- dev.null + 2 * w[i] * log(w[i] / (w[i] * muhat0[i]))
          } else if (0 < y[i] & y[i] < 1) {
              dev.null <- dev.null + 2 * w[i] * y[i] * log(w[i] * y[i] / (w[i] * muhat0[i])) + 2 * (w[i] - w[i] * y[i]) * log((w[i] - w[i] * y[i]) / (w[i] - w[i] * muhat0[i]))
          } else {
              stop ("y values must be 0 <= y <= 1!")
	  }
        }
      } 
  }
  if (fml == "gaussian") {
     wmat <- diag(w)
     b <- solve(t(vmat) %*% wmat %*% vmat) %*% t(vmat) %*% wmat %*% y
     etahat0 <- vmat %*% b
     muhat0 <- muhat.fun(etahat0)	
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
incr <- function(x) 
{
    attr(x, "shape") <- 1
    x 
}

decr <- function(x)
{
    attr(x, "shape") <- 2
    x 
} 

conv <- function(x)
{
    attr(x, "shape") <- 3
    x
}

conc <- function(x)
{
    attr(x, "shape") <- 4
    x
}

incr.conv <- function(x)
{
    attr(x, "shape") <- 5
    x
}

decr.conv <- function(x)
{
    attr(x, "shape") <- 6
    x
}

incr.conc <- function(x)
{
    attr(x, "shape") <- 7
    x
}

decr.conc <- function(x)
{
    attr(x, "shape") <- 8
    x
}

######################################################
#tree function: give the tree shape to x and return x#
######################################################
tree <- function(x)
{
   attr(x, "shape") <- "tree"
   x

}

##################################################
#tree.fun: make delta to a tree ordering variable#
##################################################
tree.fun <- function(x)
{
    if (min(x) != 0) {
	stop ("A tree ordering variable must have its placebo equal to 0!")
    }
    if (!all(round(x,0) == x)) {
	stop ("All elements of a tree ordering variable must be integers!")
    }
    #if (any(x < 0))
	#stop ("All elements of a tree ordering variable must be positive!")
    nx <- x
    obs <- 1:length(x)
    delta <- matrix(0, nrow = length(attributes(factor(x))$levels) - 1, ncol = length(x))
    pl <- min(nx)
    for (i in 1:nrow(delta)) {
       nx <- nx[which(nx != pl)]
       pl <- min(nx)
       index <- obs[x == pl]
       delta[i, index] <- 1
    }
  attr(delta, "shape") <- "tree"
  delta
}

##############################################################
#umbrella function: give the umbrella shape to x and return x#
##############################################################
umbrella <- function(x)
{
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
makedelta=function(x,sh){
	n=length(x)
# find unique x values
	xu=sort(unique(x))
	n1=length(xu)
	sm=1e-8
#  increasing or decreasing
	if(sh<3){
		amat=matrix(0,nrow=n1-1,ncol=n)
		for(i in 1:(n1-1)){
			amat[i,x>xu[i]]=1
		}
		if(sh==2){amat=-amat}
		for(i in 1:(n1-1)){amat[i,]=amat[i,]-mean(amat[i,])}
	}else if(sh==3|sh==4){
#  convex or concave
		amat=matrix(0,nrow=n1-2,ncol=n)
		for(i in 1:(n1-2)){
			amat[i,x>xu[i]]=x[x>xu[i]]-xu[i]
		}
		if(sh==4){amat=-amat}
		xm=cbind(1:n*0+1,x)
		xpx=solve(t(xm)%*%xm)
		pm=xm%*%xpx%*%t(xm)
		dmat=amat-amat%*%t(pm)
	}else if(sh>4){
		amat=matrix(0,nrow=n1-1,ncol=n)
		if(sh==5){ ### increasing convex
			for(i in 1:(n1-1)){
				amat[i,x>xu[i]]=(x[x>xu[i]]-xu[i])/(max(x)-xu[i])
			}
			for(i in 1:(n1-1)){amat[i,]=amat[i,]-mean(amat[i,])}
		}else if(sh==6){  ## decreasing convex
			for(i in 1:(n1-1)){
				amat[i,x<xu[i+1]]=(x[x<xu[i+1]]-xu[i+1])/(min(x)-xu[i+1])
			}
			for(i in 1:(n1-1)){amat[i,]=amat[i,]-mean(amat[i,])}
		}else if(sh==7){ ## increasing concave
			for(i in 1:(n1-1)){
				amat[i,x<xu[i+1]]=(x[x<xu[i+1]]-xu[i+1])/(min(x)-xu[i+1])
			}
			for(i in 1:(n1-1)){amat[i,]=-amat[i,]+mean(amat[i,])}		
		}else if(sh==8){## decreasing concave
			for(i in 1:(n1-1)){
				amat[i,x>xu[i]]=(x[x>xu[i]]-xu[i])/(max(x)-xu[i])
			}
			for(i in 1:(n1-1)){amat[i,]=-amat[i,]+mean(amat[i,])}
		}
	}
	amat
}

###################
#summary functions#
###################
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
		if (!is.null(x$residuals)) {
			cat("==============================================================", "\n")
			cat("Call:\n")
			print(x$call)
			cat("\n")
			printCoefmat(x$residuals, P.values = TRUE, has.Pvalue = TRUE)
		}
	} else {
		print ("No predictor is defined")	
		#print ("Residual degree of freedom is negative!")
	}
}

summary.cgam <- function(object,...) {
	if (!is.null(object$zcoefs)) {
		family <- object$family 
		resid_df_obs <- object$resid_df_obs
		wt.iter <- object$wt.iter
		coefs <- object$zcoefs
		se <- object$se.beta
#if (is.null(se)) {
#	print ("Residual degree of freedom is negative!")
#}
		tval <- coefs / se
		pvalbeta <- object$pvals.beta
		n <- length(coefs)
		sse0 <- object$sse0
		sse1 <- object$sse1
		#llh <- object$llh
		cic <- object$cic
		deviance <- object$deviance
		null_deviance <- object$null_deviance
		df <- object$df
		null_df <- object$null_df
		zid <- object$zid
		nzid <- object$nzid
		tms <- object$tms
		nsim <- object$nsim
		if (wt.iter) {
			rslt1 <- data.frame("Estimate" = round(coefs,4), "StdErr" = round(se,4), "z.value" = round(tval,4), "p.value" = round(pvalbeta,4))			
			rownames(rslt1)[1] <- "(Intercept)"
			if (n > 1){
				#num <- 2:n
				lzid <- length(zid)
				if (is.null(nzid)) {
					for (i in 1:lzid) {
						rownames(rslt1)[i+1] <- attributes(tms)$term.labels[zid[i]-1]
					}
				} else if (all(round(nzid, 8) == 1)) {
					for (i in 1:lzid) {
						rownames(rslt1)[i+1] <- paste(attributes(tms)$term.labels[zid[i]-1], 1, sep = "")
					}				 
				} 
			}
			rslt1 <- as.matrix(rslt1)
		} else {
			rslt1 <- data.frame("Estimate" = round(coefs,4), "StdErr" = round(se,4), "t.value" = round(tval,4), "p.value" = round(pvalbeta,4))
			rownames(rslt1)[1] <- "(Intercept)"
			if (n > 1){
				lzid <- length(zid)
				if (is.null(nzid)) {
					for (i in 1:lzid) {
						rownames(rslt1)[i+1] <- attributes(tms)$term.labels[zid[i]-1]
					}
				} else if (all(round(nzid, 8) == 1)) {
					for (i in 1:lzid) {
						rownames(rslt1)[i+1] <- paste(attributes(tms)$term.labels[zid[i]-1], 1, sep = "")
					}				 
				} 
			}
		rslt1 <- as.matrix(rslt1)
	}
		if (!is.null(sse0) & !is.null(sse1)) {
			rslt2 <- cbind(SSE.Linear = sse0, SSE.Full = sse1)
			ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs, cic = cic, null_deviance = null_deviance, null_df = null_df, deviance = deviance, df = df, resid_df_obs = resid_df_obs, family = family) 
			class(ans) <- "summary.cgam"
			ans
		} else {
			ans <- list(call = object$call, coefficients = rslt1, zcoefs = coefs, cic = cic, null_deviance = null_deviance, null_df = null_df, deviance = deviance, df = df, resid_df_obs = resid_df_obs, family = family)
			class(ans) <- "summary.cgam"
			ans
		}
	} else {
		ans <- list(zcoefs = object$zcoefs)
		class(ans) <- "summary.cgam"
		ans
	}
}

##############
#predict.cgam#
##############
predict.cgam = function(object, newData,...) {
	family = object$family
	cicfamily = CicFamily(family)
	muhat.fun = cicfamily$muhat.fun	
	shapes = object$shapes
	np = object$d0; capm = object$capm; capk = object$capk; capt = object$capt; capu = object$capu
	xid1 = object$xid1; xid2 = object$xid2; uid1 = object$uid1; uid2 = object$uid2; tid1 = object$tid1; tid2 = object$tid2 
	xmat = object$xmat; zmat = object$zmat; umb = object$umb; tr = object$tr
	bigmat = object$bigmat; umbrella.delta = object$umbrella.delta; tree.delta = object$tree.delta
	coefs = object$coefs; zcoefs = object$zcoefs; vcoefs = object$vcoefs; xcoefs = object$xcoefs; ucoefs = object$ucoefs; tcoefs = object$tcoefs
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
	m = model.frame(Terms, newData)
	newdata = m
	newx = NULL; newu = NULL; newt = NULL; newz = NULL; newv = NULL
	rn = nrow(newdata)
	newetahat = NULL; newmuhat = NULL
	newxbasis = matrix(nrow = nrow(newdata), ncol = capm); newubasis = matrix(nrow = nrow(newdata), ncol = capu); newtbasis = NULL; newbigmat = NULL
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
	for (i in 1:ncol(newdata)) {
		if (is.null(attributes(newdata[,i])$shape)) {
			newz = cbind(newz, newdata[,i])
			newv = cbind(newv, newz)
		}
		if (is.numeric(attributes(newdata[,i])$shape)) {
			newx = cbind(newx, newdata[,i])
			if (attributes(newdata[,i])$shape > 2 & attributes(newdata[,i])$shape < 5) {
				newv = cbind(newv, newdata[,i])
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
	newv = cbind(1:rn*0 + 1, newv)
	etahat = 1:rn*0
#######################################
#make newdata into same column as xmat#
#######################################
	#if (ncol(xmat) < 1) {
	#	stop ("There should be at least one constrained predictor to use predict.cgam!")
	#}
	if (!is.null(newx)) {
		newedge = NULL
		for (j in 1:ncol(xmat)) {
			x = xmat[,j]; nx = length(x); xs = sort(x)
			xu = unique(x); nxu = length(xu); xus = sort(xu)
			pos1 = xid1[j]; pos2 = xid2[j]
			x_edges = t(bigmat[pos1:pos2, , drop = FALSE])
			nedge = ncol(x_edges)
			ord = order(x)
			ordu = order(xu)
			if (any(newx[,j] > max(x)) | any(newx[,j] < min(x))) {
				stop ("No extrapolation is allowed in cgam prediction!")
			}
			newedge0 = matrix(0, nrow = rn, ncol = nedge)
#######################################
#store the 2nd edge to check incr/decr#
#######################################
			if (nedge >= 2) {
				e2 = x_edges[, 2, drop = FALSE]; dist2 = round(diff(e2), 4); s2 = sign(dist2)
				e2ord = e2[ord]; dist2ord = round(diff(e2ord), 4); s2ord = sign(dist2ord)
			}
			for (i in 1:nedge) {	
				e = x_edges[, i, drop = FALSE]; dist = round(diff(e), 4); s = sign(dist)
				eord = e[ord]; distord = round(diff(eord), 4); sord = sign(distord)
				eu = unique(e); distu = round(diff(eu), 4); su = sign(distu)
				euord = eu[ordu]; distuord = round(diff(euord), 4); suord = sign(distuord)
				if (nx == 2) {
					newedge0[,i] = my_line(newx[,j], eord, xs, nx, 1)$yp
				} else {
###############################
#concave and convex: all edges#
###############################
					if (any(sord == 1) & any(sord == -1)) {
						for (l in 1:(nxu-2)) {
							if (suord[l] != suord[l+1]) {
								pos = l+1
								break
							}
						}
						if (any(newx[,j] >= xus[pos])) {
							ids = which(newx[,j] >= xus[pos])
							newedge0[ids,i] = my_line(newx[ids,j], y = euord, x = xus, end = nxu, start = pos)$yp
						} 
						if (any(newx[,j] < xus[pos])) {
							ids = which(newx[,j] < xus[pos])
							newedge0[ids,i] = my_line(newx[ids,j], y = euord, x = xus, end = pos, start = 1)$yp
						}  
#########################################################
#incr and decr:all edges except for the 1st and the last#
#########################################################
					} else if (sord[1] == 0 & sord[nx-1] == 0 & (sum(sord != 0) == 1)) {
						for (l in 1:(nx-2)) {
							if (sord[l] != sord[l+1]) {
								pos0 = l + 1
								pos = pos0 + 1
								break 
							} 
						}
						if (any(newx[,j] >= xs[pos])) {
							ids = which(newx[,j] >= xs[pos])
							newedge0[ids,i] = eord[pos]
						} 
						if (any(newx[,j] < xs[pos])) {
							ids = which(newx[,j] < xs[pos])
							newedge0[ids,i] = eord[pos0]
						}
#########################################################
#incr.conc,incr.conv,decr.conc, and decr.conv: all edges#
#########################################################
					} else { 
						pos0 = NULL
						for (l in 2:(nx-1)) {
							if (all(sord[l:(nx-1)] == 0) && any(sord[1:(l-1)] != 0)) {
								pos0 = l 
								break
							}
						}
						if (is.null(pos0)) {
							pos0 = NULL
							for (l in (nx-1):2) {
								if (all(sord[1:(l-1)] == 0) && any(sord[l:(nx-1)] != 0)) {
									pos0 = l
									break
								}
							}		
							if (is.null(pos0) & (all(suord > 0) | all(suord < 0))) {
								newedge0[,i] = my_line(newx[,j], y = eord, x = xs, end = nx, start = 1)$yp
							} else {
								obs1 = pos0:(nx-1)
								del = NULL
								if (any(sord[pos0:(nx-1)] == 0)) {
									id = obs1[sord[pos0:(nx-1)] == 0]
									del = c(del, id+1)	
								}
								if (is.null(del)) {
									pos = pos0
									if (any(newx[,j] >= xs[pos])) {
										ids = which(newx[,j] >= xs[pos])
										newedge0[ids,i] = my_line(newx[ids,j], y = eord, x = xs, end = nx, start = pos)$yp
									} 
									if (any(newx[,j] < xs[pos])) {
										ids = which(newx[,j] < xs[pos])
										newedge0[ids,i] = my_line(newx[ids,j], y = eord, x = xs, end = pos, start = 1)$yp
									} 
								} else {
									xc = xs[-del]
									nxc = length(xc)
									ec = eord[-del]
									distc = round(diff(ec), 4)
									sc = sign(distc)
									for (l in 1:(nxc-2)) {
										if (sc[l] != sc[l+1]) {
											pos = l+1
											break
										}
									}
									if (any(newx[,j] >= xc[pos])) {
										ids = which(newx[,j] >= xc[pos])
										newedge0[ids,i] = my_line(newx[ids,j], y = ec, x = xc, end = nxc, start = pos)$yp
									} 
									if (any(newx[,j] < xc[pos])) {
										ids = which(newx[,j] < xc[pos])
										newedge0[ids,i] = my_line(newx[ids,j], y = ec, x = xc, end = pos, start = 1)$yp
									}
								}
							}
						} else {
							obs2 = 1:(pos0-1)
							del = NULL
							if (any(sord[1:(pos0-1)] == 0)) {
								id = obs2[sord[1:(pos0-1)] == 0]
								del = c(del, id+1)
							}
							if (is.null(del)) {
								pos = pos0
								if (any(newx[,j] >= xs[pos])) {
									ids = which(newx[,j] >= xs[pos])
									newedge0[ids,i] = my_line(newx[ids,j], y = eord, x = xs, end = nx, start = pos)$yp
								} 
								if (any(newx[,j] < xs[pos])) {
									ids = which(newx[,j] < xs[pos])
									newedge0[ids,i] = my_line(newx[ids,j], y = eord, x = xs, end = pos, start = 1)$yp
								} 
							} else {
								xc = xs[-del]
								nxc = length(xc) 
								ec = eord[-del]
								distc = round(diff(ec), 4)
								sc = sign(distc)
								for (l in 1:(nxc-2)) {
									if (sc[l] != sc[l+1]) {
										pos = l+1
										break
									}
								}
								if (any(newx[,j] >= xc[pos])) {
									ids = which(newx[,j] >= xc[pos])
									newedge0[ids,i] = my_line(newx[ids,j], y = ec, x = xc, end = nxc, start = pos)$yp
								} 
								if (any(newx[,j] < xc[pos]))  {
									ids =  which(newx[,j] < xc[pos])
									newedge0[ids,i] = my_line(newx[ids,j], y = ec, x = xc, end = pos, start = 1)$yp
								}
							}
						}
					}
				}
			}		
##################################################################
#make corrections to the 1st and the last element of the newedge0# 
#for the 1st and the last edge of the incr and the decr case     #
##################################################################
			if (nx > 2 & nedge >= 2) {
				e1 = x_edges[, 1, drop = FALSE]; e1ord = e1[ord]
				ek = x_edges[, nedge, drop = FALSE]; ekord = ek[ord]
				pos01 = 1; pos1 = 2; pos0k = nx-1; posk = nx
				if (nx == 3 & nedge >= 2) { 
					if ((s2ord[1] != s2ord[nx-1]) & any(c(s2ord[1], s2ord[nx-1]) == 0)) {
						if (any(newx[,j] < xs[pos1] & newx[,j] >= xs[pos01])) {
							ids = which(newx[,j] < xs[pos1] & newx[,j] >= xs[pos01])
							newedge0[ids, 1] = e1ord[pos01]
						} 
						if (any(newx[,j] >= xs[pos0k] & newx[,j] < xs[posk])) {
							ids = which(newx[,j] >= xs[pos0k] & newx[,j] < xs[posk])
							newedge0[ids,nedge] = ekord[pos0k]
						}
					}
				} else if (nx > 3 & nedge >= 2) {
					if (s2ord[1] == 0 & s2ord[nx-1] == 0) {
						if (any(newx[,j] < xs[pos1] & newx[,j] >= xs[pos01])) {
							ids = which(newx[,j] < xs[pos1] & newx[,j] >= xs[pos01])
							newedge0[ids,1] = e1ord[pos01]
						} 
				 		if (any(newx[,j] >= xs[pos0k] & newx[,j] < xs[posk])) {
							ids = which(newx[,j] >= xs[pos0k] & newx[,j] < xs[posk])
							newedge0[ids,nedge] = ekord[pos0k]
						}
					}
				}
			}	
		newedge = cbind(newedge, newedge0)
		newxbasis[, (xid1[j]-np):(xid2[j]-np)] = newedge0
		}
		etahat = cbind(newv, newedge) %*% c(vcoefs, xcoefs) 
	}
	if (!is.null(newu)) {
		newuedge = NULL
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
				ue = u_edges[, i, drop = FALSE]; udist = round(diff(ue), 4); s = sign(udist)
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
		newubasis[, (uid1[j]-np-capm):(uid2[j]-np-capm)] = newuedge0 
		}
		etahat = etahat + newuedge %*% ucoefs
	}
	if (!is.null(newt)) {
		tru = unique(tr); placebo = min(tru)
		pos_end = 0
		tr_etahat = 0
		for (i in 1: ncol(newt)) {
			#te = t(tree.fun(tr[,i]))
			pos1 = tid1[i]; pos2 = tid2[i]
			te = t(bigmat[pos1:pos2, , drop = FALSE])
			pos_start = pos_end + 1
			pos_end = pos_end + ncol(te)
			newtu = unique(newt[,i])
			for (inewt in 1: length(newtu)) {
				if (!any(tru[,i] == newtu[inewt])) {
					stop ("new tree ordering factor must be among the old tree ordering facors!")
				} 
			} 
			if (any(newt[,i] != placebo)) {	
				tr_k = sum(newtu != placebo)				
				tr_etahat = tr_etahat + te[, 1:tr_k , drop = FALSE] %*% c(tcoefs[pos_start:pos_end][1:tr_k])
			}
		}
		etahat = etahat + tr_etahat
	} 
	#newetahat = c(newetahat, etahat)
	newetahat = etahat 	
	if (is.null(newx)) {
		newetahat = newv %*% vcoefs + newetahat
	}
	newmuhat = muhat.fun(newetahat, fml = family$family)
	if (!is.null(newt)) {
		newtbasis = t(tree.delta[ ,1:nrow(newData), drop = FALSE])
	}
	#ans = new.env()
	if (!is.null(newt)) {
		newbigmat = t(cbind(newv, newxbasis, newubasis, newtbasis))
	} else {
		newbigmat = t(cbind(newv, newxbasis, newubasis))
	}
	ans = list(newv = newv, newxbasis = newxbasis, newubasis = newubasis, newtbasis = newtbasis, newbigmat = newbigmat, etahat = newetahat, fit = newmuhat)
	return (ans) 
}
