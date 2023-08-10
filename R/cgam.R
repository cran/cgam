######
#cgam#
######
cgam <- function(formula, cic = FALSE, nsim = 100, family = gaussian, cpar = 1.5, data = NULL, weights = NULL, sc_x = FALSE, sc_y = FALSE, pnt = TRUE, pen = 0, var.est = NULL, gcv = FALSE, pvf = TRUE)
{
	cl <- match.call()
    if (is.character(family))
    	family <- get(family, mode = "function", envir = parent.frame())
  	if (is.function(family))
     	family <- family()
  	if (is.null(family$family))
    	stop("'family' not recognized!")
    if (family$family == "ordered") {
		rslt <- cgam.polr(formula, data = data, weights = weights, family = family, nsim = nsim, cpar = cpar)
		rslt$call <- cl
		return (rslt)
	} else {
		labels <- NULL
  		mf <- match.call(expand.dots = FALSE)
		m <- match(c("formula", "data"), names(mf), 0L)
  		mf <- mf[c(1L, m)]
  		mf[[1L]] <- as.name("model.frame")
  		mf <- eval(mf, parent.frame())
 		ynm <- names(mf)[1]
  		mt <- attr(mf, "terms")
  		y <- model.response(mf, "any")
  		if (family$family == "binomial") {
			#if (class(y) == "factor") {
  		if(inherits(y, "factor")){
				y <- ifelse(y == levels(y)[1], 0, 1)
			}
  		}
#additive
  		shapes1_add <- NULL; shapes2_add <- NULL; shapes_add <- NULL
  		xmat_add <- NULL; xnms_add <- NULL
  		tr <- NULL; pl <- NULL; umb <- NULL
  		tree.delta <- NULL; umbrella.delta <- NULL
  		tid1 <- NULL; tid2 <- NULL; tpos2 <- 0
  		uid1 <- NULL; uid2 <- NULL; upos2 <- 0
  		nums_add <- NULL; ks_add <- list(); sps_add <- NULL; xid_add <- 1
  		zmat <- NULL; zid <- NULL; zid0 <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_param <- NULL; is_fac <- NULL; vals <- NULL; st <- 1; ed <- 1
  		ztb <- list(); iztb <- 1
		iadd <- 0
		varlist_add <- NULL
        xmat0_add <- NULL
#warp
		warp.delta <- NULL
  		nums_wp <- NULL; ks_wp <- list(); ks_wps <- list(); sps_wp <- NULL
  		xmat_wp <- NULL; x1_wp <- NULL; x2_wp <- NULL; xnms_wp <- NULL; nks0_wp <- NULL; ks0_wp <- NULL; dc <- NULL
  		x1_wps <- NULL; x2_wps <- NULL; dcss <- list()
  		#zmat <- NULL; zid <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_fac <- NULL; is_param <- NULL; vals <- NULL; st <- 1; ed <- 1
		iwps <- 0
		varlist_wps <- NULL
#tri
		tri.delta <- NULL
  		nums_tri <- NULL; ks_tri <- list(); ks_tris <- list(); nks_tris <- list(); sps_tri <- NULL
  		xmat_tri <- NULL; x1_tri <- NULL; x2_tri <- NULL; xnms_tri <- NULL; nks0_tri <- NULL; ks0_tri <- NULL; cvs <- NULL
  		x1_tris	<- NULL; x2_tris <- NULL; cvss <- list()
  		#zmat <- NULL; zid <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_fac <- NULL; is_param <- NULL; vals <- NULL; st <- 1; ed <- 1
  		#print (dim(mf))
		itri <- 0
		#print (head(mf))
  		for (i in 2:ncol(mf)) {
#additive
#print (attributes(mf[,i])$class)
		if (!is.null(attributes(mf[,i])$categ)) {
    		if (is.numeric(attributes(mf[,i])$shape) & (attributes(mf[,i])$categ == "additive")) {
       			labels <- c(labels, "additive")
       			shapes1_add <- c(shapes1_add, attributes(mf[,i])$shape)
       			xmat_add <- cbind(xmat_add, mf[,i])
       			xnms_add <- c(xnms_add, attributes(mf[,i])$nm)
       			nums_add <- c(nums_add, attributes(mf[,i])$numknots)
       			sps_add <- c(sps_add, attributes(mf[,i])$space)
       			ks_add[[xid_add]] <- attributes(mf[,i])$knots
       			xid_add <- xid_add + 1
    		}
    		if (is.character(attributes(mf[,i])$shape) & (attributes(mf[,i])$categ == "additive")) {
    			labels <- c(labels, "additive")
       			shapes2_add <- c(shapes2_add, attributes(mf[,i])$shape)
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
#warp
			if (is.character(attributes(mf[, i])$shape) & (attributes(mf[,i])$categ == "warp")) {
				iwps <- iwps + 1
				labels <- c(labels, rep(paste("warp", iwps, sep = "_"), 2))
    			sps_wp <- attributes(mf[, i])$space
				dcs <- attributes(mf[, i])$decreasing
				dcss[[iwps]] <- dcs
				ks0_wp <- attributes(mf[, i])$knots
				nks0_wp <- attributes(mf[, i])$numknots
				x1_wp <- (mf[, i])[, 1]
				x1_wps <- cbind(x1_wps, x1_wp)
				x2_wp <- (mf[, i])[, 2]
				x2_wps <- cbind(x2_wps, x2_wp)
				xmat_wp <- cbind(xmat_wp, mf[, i])
				xnms_wp <- c(xnms_wp, attributes(mf[, i])$name)
				#print (xnms_wp)
				#print (dim(xmat_wp))
				#print (nks0_wp)
				#print (ks0_wp)
				#print (sps_wp)
				ans_warp <- makedelta_wps(x1t = x1_wp, x2t = x2_wp, m1_0 = nks0_wp[1], m2_0 = nks0_wp[2], k1 = ks0_wp$k1, k2 = ks0_wp$k2, space = sps_wp, decreasing = dcs)
				warp.delta0 <- ans_warp$delta
				k1 <- ans_warp$k1
				k2 <- ans_warp$k2
				#print (k1)
				#print (k2)
				ks_wp[[1]] <- k1
				ks_wp[[2]] <- k2
				ks_wps[[iwps]] <- ks_wp
				if (iwps > 1) {
					warp.delta0 <- warp.delta0[,-1]
				}
				varlist_wps <- c(varlist_wps, 1:ncol(warp.delta0)*0 + iwps)
				warp.delta <- cbind(warp.delta, warp.delta0)
				#save(warp.delta, file='warpd.Rda')
    			#print (dim(warp.delta))
    			#print (ks_wps)
    		}
#tri
			if (is.character(attributes(mf[, i])$shape) & (attributes(mf[,i])$categ == "tri")) {
				itri <- itri + 1
				labels <- c(labels, rep(paste("tri", itri, sep = "_"), 2))
				sps_tri <- attributes(mf[, i])$space
				cvs <- attributes(mf[, i])$cvs
				cvss[[itri]] <- cvs
				ks0_tri <- attributes(mf[, i])$knots
				nks0_tri <- attributes(mf[, i])$numknots
				ks_tris[[itri]] <- ks0_tri
				nks_tris[[itri]] <- nks0_tri
				x1_tri <- (mf[, i])[, 1]
				x1_tris <- cbind(x1_tris, x1_tri)
				x2_tri <- (mf[, i])[, 2]
				x2_tris <- cbind(x2_tris, x2_tri)
				#xmat_tri <- cbind(x1_tri, x2_tri)
				xmat_tri <- cbind(xmat_tri, mf[, i])
				xnms_tri <- c(xnms_tri, attributes(mf[, i])$name)
    		}
    		#print (class((mf[, i])))
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
  		if (family$family == "binomial" | family$family == "poisson" | family$family == "Gamma") {
    		wt.iter = TRUE
  		} else {wt.iter = FALSE}
  	  	#attr(xmat, "shape") <- shapes1
  	  	#print (labels)
  		if (any(labels == "additive") | !is.null(zmat)) {
  			xmat0_add <- xmat_add; shapes0_add <- shapes1_add; nums0_add <- nums_add; ks0_add <- ks_add; sps0_add <- sps_add; xnms0_add <- xnms_add; idx_s <- NULL; idx <- NULL
  			if (any(shapes1_add == 17)) {
    			kshapes <- length(shapes1_add)
    			obs <- 1:kshapes
    			idx_s <- obs[which(shapes1_add == 17)]; idx <- obs[which(shapes1_add != 17)]
    			xmat0_add[ ,1:length(idx_s)] <- xmat_add[ ,idx_s]
    			shapes0_add[1:length(idx_s)] <- shapes1_add[idx_s]
    			nums0_add[1:length(idx_s)] <- nums_add[idx_s]
    			sps0_add[1:length(idx_s)] <- sps_add[idx_s]
    			ks0_add[1:length(idx_s)] <- ks_add[idx_s]
    			xnms0_add[1:length(idx_s)] <- xnms_add[idx_s]
    			if (length(idx) > 0) {
      				xmat0_add[ ,(1 + length(idx_s)):kshapes] <- xmat_add[ ,idx]
      				shapes0_add[(1 + length(idx_s)):kshapes] <- shapes1_add[idx]
      				nums0_add[(1 + length(idx_s)):kshapes] <- nums_add[idx]
      				sps0_add[(1 + length(idx_s)):kshapes] <- sps_add[idx]
      				ks0_add[(1 + length(idx_s)):kshapes] <- ks_add[idx]
      				xnms0_add[(1 +length(idx_s)):kshapes] <- xnms_add[idx]
    			}
    			#xmat <- xmat0; nums <- nums0; ks <- ks0; sps <- sps0; xnms <- xnms0
  			}
  			shapes_add <- c(shapes1_add, shapes2_add)
  		}
  		#print (labels)
  		xnms <- c(xnms_wp, xnms_add, xnms_tri)
  		xmat <- cbind(xmat_wp, xmat_add, xmat_tri)
  		colnames(xmat) <- xnms
  		#boolz <- is.null(labels) & !is.null(zmat)
  		if (all(labels == "additive")) {
  			if (is.null(shapes_add)) {
    			nsim <- 0
  			}
  			ans <- cgam.fit(y = y, xmat = xmat0_add, zmat = zmat, shapes = shapes0_add, numknots = nums0_add, knots = ks0_add, space = sps0_add, nsim = nsim, family = family, cpar = cpar, wt.iter = wt.iter, umbrella.delta = umbrella.delta, tree.delta = tree.delta, weights = weights, zid = zid, zid1 = zid1, zid2 = zid2, sc_x = sc_x, sc_y = sc_y, idx_s = idx_s, idx = idx, var.est = var.est, data = data)
  			if (!is.null(uid1) & !is.null(uid2)) {
    			uid1 <- uid1 + ans$d0 + ans$capm
    			uid2 <- uid2 + ans$d0 + ans$capm
  			}
  			if (!is.null(tid1) & !is.null(tid2)) {
    			tid1 <- tid1 + ans$d0 + ans$capm + ans$capu
    			tid2 <- tid2 + ans$d0 + ans$capm + ans$capu
  			}
#new:
  			knots_add <- ans$knots
  			numknots_add <- ans$numknots
#new:
  			if (length(knots_add) > 0) {
    			names(knots_add) <- xnms_add
  			}
  			rslt <- list(etahat = ans$etahat, muhat = ans$muhat, vcoefs = ans$vcoefs,
  			             xcoefs = ans$xcoefs, zcoefs = ans$zcoefs, ucoefs = ans$ucoefs,
  			             tcoefs = ans$tcoefs, coefs = ans$coefs, cic = ans$cic, d0 = ans$d0,
  			             edf0 = ans$edf0, etacomps = ans$etacomps, y = y,
  			             xmat_add = xmat_add, zmat = zmat, ztb = ztb, tr = tr, umb = umb,
  			             tree.delta = tree.delta, umbrella.delta = umbrella.delta, bigmat = ans$bigmat,
  			             shapes = shapes_add, shapesx = shapes1_add, shapesx0 = shapes0_add,
  			             prior.w = weights, wt = ans$wt, wt.iter = ans$wt.iter,
  			             family = family, SSE0 = ans$sse0, SSE1 = ans$sse1,
  			             pvals.beta = ans$pvals.beta, se.beta = ans$se.beta,
  			             df.null = ans$df.null, df = ans$df, df.residual = ans$df.residual,
  			             edf = ans$df_obs, null.deviance = ans$dev.null, deviance = ans$dev,
  			             tms = mt, capm = ans$capm, capms = ans$capms, capk = ans$capk, capt = ans$capt,
  			             capu = ans$capu, xid1 = ans$xid1, xid2 = ans$xid2, tid1 = tid1, tid2 = tid2,
  			             uid1 = uid1, uid2 = uid2, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2,
  			             nsim = nsim, xnms_add = xnms_add, xnms = xnms, ynm = ynm, znms = znms,
  			             is_param = is_param, is_fac = is_fac, knots = knots_add, numknots = numknots_add,
  			             sps = sps_add, ms = ans$ms, cpar = ans$cpar, pl = pl, idx_s = idx_s, idx = idx,
  			             xmat0 = ans$xmat2, knots0 = ans$knots2, numknots0 = ans$numknots2, sps0 = ans$sps2,
  			             ms0 = ans$ms2, phihat = ans$phihat, pvs = ans$pvs,
  			             s.edf = ans$s.edf, bstats = ans$bstats, pvsz = ans$pvsz,
  			             z.edf = ans$z.edf, fstats = ans$fstats, vh = ans$vh, kts.var = ans$kts.var,
  			             sc = ans$sc, sc_y = sc_y)
  			rslt$call <- cl
  			class(rslt) <- "cgam"
		} else if (any(grepl("warp", labels, fixed = TRUE)) & !any(grepl("tri", labels, fixed = TRUE))) {
            nprs <- sum(grepl("warp", labels, fixed = TRUE)) / 2
            delta_add <- NULL
            ms <- NULL
            np_add <- 0
            pb <- 0
            knotsuse_add <- NULL
            if (any(labels == "additive")) {
                #delta_add <- NULL
                if (length(shapes0_add) > 0) {
                    ans_add <- make_delta_add(xmat = xmat0_add, shapes = shapes0_add, numknots = nums0_add, knots = ks0_add, space = sps0_add)
                    delta_add <- ans_add$bigmat
                    ms <- ans_add$mslst
                    #np_add: smooth only, conv and conc
                    np_add <- ans_add$np
                    pb <- nrow(delta_add)
                    varlist_add <- ans_add$varlist
                    knotsuse_add <- ans_add$knotsuse
                }
            }
            delta_ut <- NULL
            if (!is.null(umbrella.delta)) {
                delta_add <- rbind(delta_add, umbrella.delta)
                delta_ut <- rbind(delta_ut, umbrella.delta)
            }
            if (!is.null(tree.delta)) {
                delta_add <- rbind(delta_add, tree.delta)
                delta_ut <- rbind(delta_ut, tree.delta)
            }
            ans <- wps.fit(x1t = x1_wps, x2t = x2_wps, y = y, zmat = zmat, xmat_add = xmat_add, delta_add = delta_add, delta_ut = delta_ut, varlist_add = varlist_add, shapes_add = shapes_add, np_add = np_add, w = weights, pen = pen, pnt = pnt, cpar = cpar, decrs = dcss, delta = warp.delta, kts = ks_wps, wt.iter = wt.iter, family = family, cic = cic, nsim = nsim, nprs = nprs, idx_s = idx_s, idx = idx, gcv = gcv, pvf = pvf)
            rslt <- list(ks_wps = ans$kts, muhat = ans$muhat, SSE1 = ans$sse1, SSE0 = ans$sse0, edfc = ans$edf, edf0 = ans$edf0, delta = ans$delta, y = y, zmat = ans$zmat, ztb = ztb, zmat_0 = ans$zmat_0, xmat_wp = xmat_wp, xmat_add = xmat_add, xmat0_add = xmat0_add, xmat = xmat, shapes = shapes_add, coefs = ans$coefs, zcoefs = ans$zcoefs, pvals.beta = ans$pz,pval=ans$pval,bval=ans$bval, se.beta = ans$sez, gcv = ans$gcv, xnms_wp = xnms_wp, xnms_add = xnms_add, xnms = xnms, xid_add = xid_add, znms = znms, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2, ynm = ynm, decrs = dcss, tms = mt, mf = mf, is_param = is_param, is_fac = is_fac, d0 = ans$d0, pb = pb, pen = ans$pen, cpar = ans$cpar, wt.iter = wt.iter, cic = ans$cic, family = family, nprs_wps = nprs, varlist_wps = varlist_wps, coef_wp = ans$coef_wp, varlist_add = varlist_add, np_add = np_add, coef_add = ans$coef_add, coef_ut = ans$coef_ut, etacomps = ans$etacomps, amat = ans$amat, dmat = ans$dmat, prior.w = weights, sig2hat = ans$sig2hat, ms = ms, knotsuse_add = knotsuse_add, vmat = ans$vmat)
            class(rslt) <- "wps"
  			if (!is.null(delta_add)) {
  				class(rslt) <- c("wps", "cgam")
  				if (min(which(labels == "additive")) < min(which(grepl("warp", labels, fixed = TRUE)))) {
  					class(rslt) <- c("cgam", "wps")
  				}
  			}
  		} else if (any(grepl("tri", labels, fixed = TRUE))) {
  			nprs <- sum(grepl("tri", labels, fixed = TRUE)) / 2
  			nprs_wps <- sum(grepl("warp", labels, fixed = TRUE)) / 2
  			amat_wp <- NULL
  			dmat_wp <- NULL
			delta_add <- NULL
            knotsuse_add <- NULL
            ms <- NULL
  			np_add <- 0
			pb <- 0
			if (any(labels == "additive")) {
				if (length(shapes0_add) > 0) {
					ans_add <- make_delta_add(xmat = xmat0_add, shapes = shapes0_add, numknots = nums0_add, knots = ks0_add, space = sps0_add)
					delta_add <- ans_add$bigmat
					np_add <- ans_add$np
					pb <- nrow(delta_add)
					varlist_add <- ans_add$varlist
                    ms <- ans_add$mslst
                    knotsuse_add <- ans_add$knotsuse
				}
			}
			if (!is.null(umbrella.delta)) {
				delta_add <- rbind(delta_add, umbrella.delta)
			}
			if (!is.null(tree.delta)) {
				delta_add <- rbind(delta_add, tree.delta)
			}
			if (any(grepl("warp", labels, fixed = TRUE))) {
				ans_wps <- makeamat_wps(ks_wps, nprs_wps)
				amat_wp <- ans_wps$amat
				dmat_wp <- ans_wps$dmat
				valist_wps <- ans_wps$valist_wps
			}
  			ans <- trispl.fit(x1t = x1_tris, x2t = x2_tris, y = y, zmat = zmat, xmat_add = xmat_add, delta_add = delta_add, varlist_add = varlist_add, shapes_add = shapes_add, np_add = np_add, xmat_wp = xmat_wp, delta_wp = warp.delta, varlist_wps = varlist_wps, amat_wp = amat_wp, dmat_wp = dmat_wp, w = weights, lambda = pen, pnt = TRUE, cpar = cpar, cvss = cvss, delta = tri.delta, kts = ks_tris, nkts = nks_tris, wt.iter = wt.iter, family = family, nsim = nsim, nprs = nprs)
  			rslt <- list(muhat = ans$muhat, etahat = ans$etahat, trimat = ans$trimat, y = y, zmat = ans$zmat, ztb = ztb, capk = ans$capk, xmat_tri = xmat_tri, xmat_add = xmat_add, xmat0_add = xmat0_add, xmat_wp = xmat_wp, xmat = xmat, shapes = shapes_add, coefs = ans$thhat, zcoefs = ans$zcoefs, se.beta = ans$se.beta, pvals.beta = ans$pvals.beta, gcv = ans$gcv, edf = ans$edf, edf0 = ans$edf0, cic = ans$cic, sse0 = ans$sse0, sse1 = ans$sse1, family = family, nprs = nprs, nprs_wps = nprs_wps, varlist_wps = varlist_wps, varlist_add = varlist_add, varlist_tri = ans$varlist, xnms_tri = xnms_tri, xnms_add = xnms_add, xnms_wp = xnms_wp, xnms = xnms, znms = znms, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2, ynm = ynm, decrs = dcss, cvss = cvss, tms = mt, is_param = is_param, is_fac = is_fac, d0 = ans$d0, pen = ans$pen, cpar = cpar, wt.iter = wt.iter, np_add = np_add, pb = pb, coef_add = ans$coef_add, coef_tri = ans$coef_tri, coef_wp = ans$coef_wp, etacomps = ans$etacomps, ks_wps = ks_wps, knots_lst = ans$knots_lst, m12_lst = ans$m12_lst, trimat_lst = ans$trimat_lst, bmat_lst = ans$bmat_lst, capk_lst = ans$capk_lst, prior.w = weights, dmatc = ans$dmatc, amatc = ans$amatc, pmatc = ans$pmatc, sig2hat = ans$sig2hat, knotsuse_add = knotsuse_add, ms = ms)
  			class(rslt) <- "trispl"
  			if (!is.null(delta_add) & is.null(warp.delta)) {
  				if (labels[1] == "additive") {
  					class(rslt) <- c("cgam", "trispl")
  				} else {class(rslt) <- c("trispl", "cgam")}
  			}
  			if (!is.null(warp.delta) & is.null(delta_add)) {
  				#if (labels[1] == "warp") {
  				if (grepl("warp", labels[1], fixed = TRUE)) {
  					class(rslt) <- c("wps", "trispl")
  				} else {class(rslt) <- c("trispl",  "wps")}
  			}
  			if (!is.null(warp.delta) & !is.null(delta_add)) {
  				if (labels[1] == "additive") {
  					class(rslt) <- c("cgam", "trispl", "wps")
  				} else if (grepl("warp", labels[1], fixed = TRUE)) {
  					class(rslt) <- c( "wps", "trispl", "cgam")
  				} else {class(rslt) <- c("trispl", "cgam", "wps")}
  			}
		}
	}
	rslt$call <- cl
	rslt$labels <- labels
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
cgam.fit <- function(y, xmat, zmat, shapes, numknots, knots, space, nsim, family = gaussian(), cpar = 1.2, wt.iter = FALSE, umbrella.delta = NULL, tree.delta = NULL, weights = NULL, zid = zid, zid1 = zid1, zid2 = zid2, sc_x = FALSE, sc_y = FALSE, idx_s = NULL, idx = NULL, var.est = NULL, data = NULL) {
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
    sc <- 1
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
		if(shapes[1] >= 9 & shapes[1] <= 17) {
		  numknotsuse <- c(numknotsuse, length(del1_ans$knots))
		} else {
		  numknotsuse <- c(numknotsuse, nrow(del1))
		}

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
				if(shapes[i] >= 9 & shapes[i] <= 17) {
				  numknotsuse <- c(numknotsuse, length(del2_ans$knots))
				} else {
				  numknotsuse <- c(numknotsuse, nrow(del2))
				}
				#numknotsuse <- c(numknotsuse, length(del2_ans$knots))
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
#new: initialize face
            face <- NULL
			if (wt.iter) {
				etahat <- etahat.fun(n, y, fml = family$family)
				gr <- gr.fun(y, etahat, weights, fml = family$family)
				wt <- wt.fun(y, etahat, n, weights, fml = family$family)
				cvec <- wt * etahat - gr
			} else {wt <- wt.fun(y, etahat, n, weights, fml = family$family)}
			  	zvec <- zvec.fun(cvec, wt, y, fml = family$family)
#        		gmat <- t(bigmat %*% sqrt(diag(wt)))
#new: avoid memory allocation error
				gmat <- t(bigmat)
				for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
			  	dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
                zsend <- gmat[, 1:np, drop = FALSE]
                #ans <- coneB(zvec, t(dsend), zsend)
                ans <- coneB(zvec, dsend, zsend)
                face <- ans$face
			  	etahat <- t(bigmat) %*% ans$coefs
#new: for monotinic variance estimation
                vh <- NULL
                kts.var <- NULL
                if (!is.null(var.est)) {
                    #x.var <- var.est
                    #test more:
                    if (!is.null(data)) {
                        nms <- names(data)
                        nm <- attributes(var.est)$nm
                        if (nm %in% nms) {
                            x.var <- data[[nm]]
                            attributes(x.var) <- attributes(var.est)
                        }
                    } else {x.var <- var.est}
                    db.exp <- attributes(x.var)$db.exp
                    kts.var <- attributes(x.var)$var.knots
                    shape.var <- attributes(x.var)$shape
                    iter <- 0
                    diff.var <- 100
                    oldeta <- etahat
                    evec0 <- y - etahat
                    #bigwt <- 10*max(evec0)
                    bigwt <- 1000
                    while(diff.var > 1e-8 & iter < 10) {
                        iter <- iter + 1
                        evec <- y - etahat
                        fit.var <- varest(evec, x.var, shape=shape.var, var.knots=kts.var, db.exp=db.exp)
                        v1 <- fit.var$vhat
                        wt0 <- 1/v1
                        wt0[wt0 > bigwt] <- bigwt
                        wt <- wt.fun(y, etahat, n, weights = wt0, fml = family$family)
                        zvec <- zvec.fun(cvec, wt, y, fml = family$family)
                        gmat <- t(bigmat)
                        for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
                        dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
                        zsend <- gmat[, 1:np, drop = FALSE]
                        ans <- coneB(zvec, dsend, zsend)
                        etahat <- t(bigmat) %*% ans$coefs
                        diff.var <- mean((etahat - oldeta)^2)
                        oldeta <- etahat
                    }
                    kts.var <- fit.var$var.knots
                    vh <- v1
                }
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
						wt <- wt.fun(y, etahat, n, weights, fml = family$family)
						cvec <- wt * etahat - gr
						#zvec <- cvec / sqrt(wt)
						zvec <- zvec.fun(cvec, wt, y, fml = family$family)
#						gmat <- t(bigmat %*% sqrt(diag(wt)))
						gmat <- t(bigmat)
						for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
						dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
        	  			zsend <- gmat[, 1:np, drop = FALSE]
                        #ans <- coneB(zvec, t(dsend), zsend)
                        ans <- coneB(zvec, dsend, zsend, face = face)
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
#print (idx_s)
#print (idx)
				#if (!is.null(idx_s)) {
				if (length(idx_s) > 0) {
					thvecs0 <- thvecs
					thvecs0[idx_s,] <- thvecs[1:length(idx_s), ]
					#if (!is.null(idx)) {
					if (length(idx) > 0) {
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
				df_obs <- sum(abs(coefskeep) > 0)
				#llh <- llh.fun(y, muhat, etahat, n, weights, fml = family$family)
				if (family$family == "Gamma") {
					if ((n - cpar * df_obs) <= 0) {
						phihat <- sum(((y - muhatkeep) / muhatkeep)^2) / df_obs
					} else {
						phihat <- sum(((y - muhatkeep) / muhatkeep)^2) / (n - cpar * df_obs)
					}
				} else {phihat <- NULL}
				#print (phihat)
				llh <- llh.fun(y, muhatkeep, etakeep, phihat, n, weights, fml = family$family)
				if (family$family == "poisson") {
					mu0 <- mean(y)
					eta0 <- log(mu0)
				} else {mu0 <- NULL}
                #new: get p-values for smooth components
                    pvs <- NULL
                    s.edf <- NULL
                    bstats <- NULL
                    if (is.numeric(shapes)) {
                        #changed to include shape==17
                        #if (all(shapes >= 9 & shapes <= 17)) {
                        if (any(shapes >= 9 & shapes <= 17)) {
                            for (i in 1:capl) {
                              if (shapes[i] >= 1 & shapes[i] <= 17 ){
                                  ansi <- cgam.pv(y=y, xmat=xmat, zmat=zmat, shapes=shapes, delta=delta, np=np, capms=capms, numknotsuse=numknotsuse, varlist=varlist, family=family, weights=weights, test_id=i, nsims=100)
                                  pvi <- ansi$pv
                                  edfi <- 1.5*sum(xcoefs[varlist == i] > 1e-8)
                                  bstati <- ansi$bstat
                                  pvs <- c(pvs, pvi)
                                  s.edf <- c(s.edf, edfi)
                                  bstats <- c(bstats, bstati)
                              }
                              #print (pvs)
                            }
                        }
                    }
                    #new: get p-values for categorical predictors, not for each level
                    pvsz <- NULL
                    z.edf <- NULL
                    fstats <- NULL
                    #temp: use cgam.pvz only when there's smooth component
                    if (is.numeric(shapes)) {
                        #changed to include shape==17
                        if (all(shapes >= 9 & shapes <= 17)) {
                            if (is.null(weights)) {
                                weights <- 1:n*0 + 1
                            }
                            prior.w <- weights
                            if (capk > 0) {
                                sse1 <- sum(prior.w * (y - muhatkeep)^2)
                                ansi <- cgam.pvz(y=y, bigmat=bigmat, df_obs=df_obs, sse1=sse1, np=np, zid=zid, zid1=zid1, zid2=zid2, muhat = muhatkeep, etahat = etakeep, coefskeep = coefskeep, wt.iter=wt.iter, family=family, weights=weights)
                                pvsz <- ansi$pvs
                                z.edf <- ansi$edfs
                                fstats <- ansi$fstats
                            }
                        }
                    }
                } else if (capk + capls > 0 & capl - capls + capt + capu == 0) {
                    if (is.null(weights)) {
                        weights <- 1:n*0 + 1
                    }
                    prior.w <- weights
                    vmat <- t(bigmat[1:np, , drop = FALSE])
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
                            #print (tvmat)
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
#new: gamma only
                        if (family$family == "Gamma") {
                            phihat <- sum(((y - muhatkeep) / muhatkeep)^2) / (n - np)
                        } else {phihat <- NULL}
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
                        se.beta <- sqrt(diag(solve(tvmat %*% vmat) * sdhat2))[1:(capk + 1)]
                        tstat <- zcoefs / se.beta
                        pvals.beta <-  (1 - pt(abs(tstat), df = n - np)) * 2
#new: for var estimation
                        vh <- NULL
                        kts.var <- NULL
                        if (!is.null(var.est)) {
                            x.var <- var.est
                            db.exp <- attributes(x.var)$db.exp
                            kts.var <- attributes(x.var)$var.knots
                            shape.var <- attributes(x.var)$shape
                            iter <- 0
                            diff.var <- 100
                            oldeta <- etahat
                            evec0 <- y - oldeta
                            #bigwt <- 10*max(evec0)
                            bigwt <- 1000
                            while(diff.var > 1e-8 & iter < 10) {
                                iter <- iter + 1
                                evec <- y - etahat
                                fit.var <- varest(evec, x.var, shape=shape.var, var.knots=kts.var, db.exp=db.exp)
                                v1 <- fit.var$vhat
                                #w <- 1/v1
                                w0 <- 1/v1
                                w0[w0 > bigwt] <- bigwt
                                w <- w0
                                tvmat <- t(vmat)
                                for (i in 1:n) {tvmat[,i] <- tvmat[,i] * w[i]}
                                b <- solve(tvmat %*% vmat) %*% tvmat %*% y
                                etahat <- vmat %*% b
                                diff.var <- mean((etahat - oldeta)^2)
                                oldeta <- etahat
                                #print (diff.var)
                            }
                            kts.var <- fit.var$var.knots
                            vh <- v1
                        }
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
                llh <- llh.fun(y, muhat, etahat, phihat, n, weights, fml = family$family)
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
                if (family$family == "binomial") {
                    nobs <- sum(prior.w)
                } else {nobs <- n}
                if ((nobs - np - cpar * (dfmean - np)) <= 0) {
                    rslt$cic <- llh + log(1 + 2 * dfmean / (dfmean - np))
                } else {
                    #new:
                    if (n<=200){cpar=1.5}
                    rslt$cic <- llh + log(1 + 2 * dfmean / (nobs - np - cpar * (dfmean - np)))
                }
                #rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - 1.5 * (dfmean - np)))
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
                #rslt$df.residual <- n - np - 1.5 * (df_obs - np)
                rslt$df.residual <- n - cpar * df_obs
                rslt$vhat <- etahat
                rslt$vmat <- vmat
                rslt$etacomps <- thvecs
                rslt$knots <- knotsuse
                rslt$numknots <- numknotsuse
                rslt$ms <- mslst
                #new: used for predict when there are only shape == 17 x's
                #rslt$xmat <- xmat
                rslt$xmat2 <- xmat
                rslt$knots2 <- knotsuse
                rslt$numknots2 <- numknotsuse
                rslt$ms2 <- mslst
                rslt$vh <- vh
                rslt$kts.var <- kts.var
                return (rslt)
            }
##########
#get cic#
##########
		  if (capl - capls + capu + capt > 0 & nsim > 0) {
	  		dfs <- 1:nsim
#new: initialize face
            face <- NULL
	  		for (isim in 1:nsim) {
				#set.seed(123)
	    		ysim <- ysim.fun(n, mu0, fml = family$family)
		  		if (wt.iter) {
					etahat <- etahat.fun(n, ysim, fml = family$family)
					gr <- gr.fun(ysim, etahat, weights, fml = family$family)
					wt <- wt.fun(ysim, etahat, n, weights, fml = family$family)
					cvec <- wt * etahat - gr
				} else {wt <- wt.fun(ysim, etahat, n, weights, fml = family$family)}
					zvec <- zvec.fun(cvec, wt, ysim, fml = family$family)
#            		gmat <- t(bigmat %*% sqrt(diag(wt)))
					gmat <- t(bigmat)
					for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
           			dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
            		zsend <- gmat[, 1:np, drop = FALSE]
                    #ans <- try(coneB(zvec, t(dsend), zsend))
                    ans <- try(coneB(zvec, dsend, zsend))
                    face <- ans$face
					#if (class(ans) == "try-error") next
          if(inherits(ans, "try-error")) next
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
							wt <- wt.fun(ysim, etahat, n, weights, fml = family$family)
							cvec <- wt * etahat - gr
							#zvec <- cvec / sqrt(wt)
							zvec <- zvec.fun(cvec, wt, y, fml = family$family)
#							gmat <- t(bigmat %*% sqrt(diag(wt)))
							gmat <- t(bigmat)
							for (i in 1:n) {gmat[i,] <- bigmat[,i] * sqrt(wt[i])}
							dsend <- gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
							zsend <- gmat[, 1:np, drop = FALSE]
                            #ans <- try(coneB(zvec, t(dsend), zsend))
                            ans <- try(coneB(zvec, dsend, zsend, face=face))
							#if (class(ans) == "try-error") next
              if(inherits(ans, "try-error")) next
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
		rslt$phihat <- phihat
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
			if ((n - np - cpar * (dfmean - np)) <= 0) {
				rslt$cic <- llh + log(1 + 2 * dfmean / (dfmean - np))
			} else {
                #new:
                if (n<=200){cpar=1.5}
        		rslt$cic <- llh + log(1 + 2 * dfmean / (n - np - cpar * (dfmean - np)))
			}
		} else {rslt$cic <- NULL}
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
#new:
		#if ((n - np - cpar * df_obs) <= 0) {
		if ((n - cpar * df_obs) <= 0) {
			#sdhat2 <- sse1 / (df_obs - np)
			sdhat2 <- sse1 / df_obs
		} else {
			#sdhat2 <- sse1 / (n - np - 1.5 * (df_obs - np))
			#sdhat2 <- sse1 / (n - np - cpar * df_obs)
			sdhat2 <- sse1 / (n - cpar * df_obs)
		}
#debugged: vmat -> vmat and duse
#new: coefskeep include zcoefs; bigmat include vmat
		pmat <- vmat
		capbm <- length(bigmat) / n
		bigmat_nv <- bigmat[(np + 1):capbm, , drop = FALSE]
		coefs_nv <- coefskeep[(np + 1):capbm]
		duse <- coefs_nv > 1e-8
		if (sum(duse) >= 1) {
			pmat <- cbind(vmat, t(bigmat_nv[duse, , drop = FALSE]))
		}
		if (wt.iter) {
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
			#if ((n - np - cpar * df_obs) <= 0) {
			if ((n - cpar * df_obs) <= 0) {
				pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  df_obs))
				warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
			} else {
				#pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  n - np - cpar * df_obs))
				pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  n - cpar * df_obs))
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
		#rslt$df.residual <- n - np - cpar * df_obs
		rslt$df.residual <- n - cpar * df_obs
		rslt$df_obs <- df_obs
		rslt$vhat <- vhat
		if (length(knotsuse) == 0) {
			knotsuse <- NULL
		}
#new: order back
		knotsuse2 <- knotsuse
		numknotsuse2 <- numknotsuse
		mslst2 <- mslst
		xmat2 <- xmat
		#if (!is.null(idx_s)) {
		if (length(idx_s) > 0) {
			knotsuse0 <- knotsuse
			numknotsuse0 <- numknotsuse
			mslst0 <- mslst
			knotsuse0[idx_s] <- knotsuse[1:length(idx_s)]
			numknotsuse0[idx_s] <- numknotsuse[1:length(idx_s)]
			mslst0[idx_s] <- mslst[1:length(idx_s)]
			#if (!is.null(idx)) {
			if (length(idx) > 0) {
				knotsuse0[idx] <- knotsuse[(1+length(idx_s)):capl]
				numknotsuse0[idx] <- numknotsuse[(1+length(idx_s)):capl]
				mslst0[idx] <- mslst[(1+length(idx_s)):capl]
			}
			knotsuse <- knotsuse0
			numknotsuse <- numknotsuse0
			mslst <- mslst0
		}
# the following has shape = 17 at the beginning: knots2 etc
		rslt$knots2 <- knotsuse2
		rslt$numknots2 <- numknotsuse2
		rslt$ms2 <- mslst2
		rslt$xmat2 <- xmat2
		rslt$knots <- knotsuse
		rslt$numknots <- numknotsuse
		rslt$ms <- mslst
		rslt$capms <- capms
		rslt$cpar <- cpar
        rslt$pvs <- pvs
        rslt$s.edf <- s.edf
        rslt$bstats <- bstats
        #new:
        rslt$pvsz <- pvsz
        rslt$z.edf <- z.edf
        rslt$fstats <- fstats
        rslt$vh <- vh
        rslt$kts.var <- kts.var
        rslt$sc <- sc
		#rslt$sdhat2 <- sdhat2
		return (rslt)
	}
}

###########
#CicFamily#
###########
CicFamily <- function(object,...)UseMethod("CicFamily")
CicFamily <- function(object) {
#temp:
#if (object$family == "ordered") {
#	object$family = "gaussian"
#}
#new:
  if (object$family == "gaussian") {
  	linkfun <- function (mu) mu
  }
  if (object$family == "poisson") {
  	linkfun <- function (mu) log(mu)
  }
  if (object$family == "binomial") {
  	linkfun <- binomial()$linkfun
  }
  if (object$family == "Gamma") {
  	#linkfun <- function (mu) log(mu)
     linkfun <- Gamma(link="log")$linkfun
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

  muhat.fun <- function(etahat, wt = NULL, fml = object$family){
    n <- length(etahat)
    if (fml == "poisson") {
      muhat <- exp(etahat)
    }
    if (fml == "binomial") {
      muhat <- 1:n*0
      id1 = which(etahat > 100)
      id2 = which(etahat <= 100)
      muhat[id1] = 1
      muhat[id2] = exp(etahat[id2]) / (1 + exp(etahat[id2]))
      #muhat[wt == 0] <- 1
      #for (i in 1:n) {
      #  if (etahat[i] > 100) {
      #      muhat[i] <- 1
      #  } else {
      #    muhat[i] <- exp(etahat[i]) / (1 + exp(etahat[i]))
      #  }
      #}
    }
    if (fml == "gaussian") {
      muhat <- etahat
    }
    if (fml == "Gamma") {
      muhat <- exp(etahat)
    }
   muhat
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
    if (fml == "Gamma") {
    	deriv <- 1 / muhat
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
  if (fml == "Gamma") {
    dev <- 0
	for (i in 1:n) {
	  if (y[i] == 0) {
	  		sm <- 1e-4
            dev <- dev + 2 * w[i] * ((sm - muhat[i]) / muhat[i] - log(sm / muhat[i]))
       } else {
            dev <- dev + 2 * w[i] * ((y[i] - muhat[i]) / muhat[i] - log(y[i] / muhat[i]))
       }
	}
  }
###################
#get null deviance#
###################
  if (fml == "binomial" | fml == "poisson" | fml == "Gamma") {
      diff <- 1
      muhat0 <- mean(y) + 1:n*0
      if (fml == "poisson") {
         etahat0 <- log(muhat0)
      }
      if (fml == "binomial") {
         etahat0 <- log(muhat0 / (1 - muhat0))
      }
      if (fml == "Gamma") {
    	etahat0 <- log(muhat0)
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
     if (fml == "Gamma") {
      	dev.null <- 0
		for (i in 1:n) {
	  		if (y[i] == 0) {
	  			sm <- 1e-4
            	dev.null <- dev.null + 2 * w[i] * ((sm - muhat0[i]) / muhat0[i] - log(sm / muhat0[i]))
       		} else {
            	dev.null <- dev.null + 2 * w[i] * ((y[i] - muhat0[i]) / muhat0[i] - log(y[i] / muhat0[i]))
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
  ans <- list(llh.fun = llh.fun, etahat.fun = etahat.fun, gr.fun = gr.fun, wt.fun = wt.fun, zvec.fun = zvec.fun, muhat.fun = muhat.fun, ysim.fun = ysim.fun, deriv.fun = deriv.fun, dev.fun = dev.fun, linkfun = linkfun)
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
    x
}

#s.incr <- function(x, numknots = 0, knots = 0, space = "E")
#{
#    cl <- match.call()
#    pars <- match.call()[-1]
#    attr(x, "nm") <- deparse(pars$x)
#    attr(x, "shape") <- 9
#    attr(x, "numknots") <- numknots
#    attr(x, "knots") <- knots
#    attr(x, "space") <- space
#    attr(x, "categ") <- "additive"
#    #class(x) <- "additive"
#    x
#}

#s.decr <- function(x, numknots = 0, knots = 0, space = "E")
#{
#    cl <- match.call()
#    pars <- match.call()[-1]
#    attr(x, "nm") <- deparse(pars$x)
#    attr(x, "shape") <- 10
#    attr(x, "numknots") <- numknots
#    attr(x, "knots") <- knots
#    attr(x, "space") <- space
#    attr(x, "categ") <- "additive"
#    #class(x) <- "additive"
#    x
#}


s.incr <- function(x, numknots = 0, knots = 0, var.knots = 0, space = "Q", db.exp = FALSE)
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 9
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    attr(x, "categ") <- "additive"
    attr(x, "db.exp") <- db.exp
    attr(x, "var.knots") <- var.knots
    #class(x) <- "additive"
    x
}

s.decr <- function(x, numknots = 0, knots = 0, var.knots = 0, space = "Q", db.exp = FALSE)
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 10
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    attr(x, "categ") <- "additive"
    attr(x, "db.exp") <- db.exp
    attr(x, "var.knots") <- var.knots
    #class(x) <- "additive"
    x
}


s.conv <- function(x, numknots = 0, knots = 0, space = "Q")
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 11
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
    x
}

s.conc <- function(x, numknots = 0, knots = 0, space = "Q")
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 12
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
    x
}

s.incr.conv <- function(x, numknots = 0, knots = 0, space = "Q")
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 13
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
    x
}

s.incr.conc <- function(x, numknots = 0, knots = 0, space = "Q")
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 14
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
    x
}

s.decr.conv <- function(x, numknots = 0, knots = 0, space = "Q")
{
    cl <- match.call()
    pars <- match.call()[-1]
    attr(x, "nm") <- deparse(pars$x)
    attr(x, "shape") <- 15
    attr(x, "numknots") <- numknots
    attr(x, "knots") <- knots
    attr(x, "space") <- space
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
    x
}

s.decr.conc <- function(x, numknots = 0, knots = 0, space = "Q")
{
   cl <- match.call()
   pars <- match.call()[-1]
   attr(x, "nm") <- deparse(pars$x)
   attr(x, "shape") <- 16
   attr(x, "numknots") <- numknots
   attr(x, "knots") <- knots
   attr(x, "space") <- space
   attr(x, "categ") <- "additive"
   #class(x) <- "additive"
   x
}

s <- function(x, numknots = 0, knots = 0, space = "Q")
{
   cl <- match.call()
   pars <- match.call()[-1]
   attr(x, "nm") <- deparse(pars$x)
   attr(x, "shape") <- 17
   attr(x, "numknots") <- numknots
   attr(x, "knots") <- knots
   attr(x, "space") <- space
   attr(x, "categ") <- "additive"
   #class(x) <- "additive"
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
    #stop (print (x))
    #print (class(x))
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
    attr(x, "categ") <- "additive"
    #class(x) <- "additive"
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
  attr(x, "categ") <- "additive"
  #class(x) <- "additive"
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
		if (!interp) {
			for (i in 1:(n1 - 1)) {
#new: use ms in predict.cgam
				ms = c(ms, mean(amat[i, ]))
				amat[i, ] = amat[i, ] - mean(amat[i, ])
			}
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
		if (!interp) {
			ms = amat %*% t(pm)
			#amat = amat - amat %*% t(pm)
			amat = amat - ms
		}
	} else if (sh > 4 & sh < 9) {
		amat = matrix(0, nrow = n1 - 1, ncol = n)
		if (sh == 5) { ### increasing convex
			for (i in 1:(n1 - 1)) {
				amat[i, x > xu[i]] = (x[x > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			if (!interp) {
				for (i in 1:(n1 - 1)) {
					ms = c(ms, mean(amat[i, ]))
					amat[i,] = amat[i,] - mean(amat[i,])
				}
			}
		} else if (sh == 6) {  ## decreasing convex
			for (i in 1:(n1 - 1)) {
				amat[i, x < xu[i + 1]] = (x[x < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			if (!interp) {
				for (i in 1:(n1 - 1)) {
					ms = c(ms, mean(amat[i, ]))
					amat[i,] = amat[i,] - mean(amat[i, ])
				}
			}
#print (ms)
		} else if (sh == 7) { ## increasing concave
			for (i in 1:(n1 - 1)) {
				amat[i, x < xu[i + 1]] = (x[x < xu[i + 1]] - xu[i + 1]) / (min(x) - xu[i + 1])
			}
			if (!interp) {
				for (i in 1:(n1 - 1)) {
					ms = c(ms, mean(amat[i, ]))
					amat[i,] = -amat[i,] + mean(amat[i,])
				}
			}
		} else if (sh == 8) {## decreasing concave
			for (i in 1:(n1 - 1)) {
				amat[i, x > xu[i]] = (x[x > xu[i]] - xu[i]) / (max(x) - xu[i])
			}
			if (!interp) {
				for (i in 1:(n1 - 1)) {
					ms = c(ms, mean(amat[i, ]))
					amat[i,] = -amat[i,] + mean(amat[i,])
				}
			}
		}
	} else if (sh > 8 & sh < 18) {
        #new: add two knots
		#if (all(knots == 0) & numknots == 0) {
		if (length(knots) < 2 & numknots == 0) {
			if (sh == 9 | sh == 10) {#1 2
                #k = trunc(n1^(1/5)) + 4
                if (n1 <= 50) {
                    k = 5
                } else if (n1>50 && n1<100) {
                    k = 6
                } else if (n1>= 100 && n1<200) {
                    k = 7
                } else {
                    k = trunc(n1^(1/5)) + 6
                }
			} else {
                #k = trunc(n1^(1/7) + 4)
                if (n1 <= 50) {
                    k = 5
                } else if (n1>50 && n1<100) {
                    k = 6
                } else if (n1>= 100 && n1<200) {
                    k = 7
                } else {
                    k = trunc(n1^(1/7)) + 6
                }
            }
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

#####
#make basis for additive components
make_delta_add = function(xmat, shapes, numknots, knots, space) {
	capl = ncol(xmat)
	if (capl < 1) {capl = 0}
	if (round(capl, 8) != round(capl, 1)) {stop ("Incompatible dimensions for xmat!")}
	capls = sum(shapes == 17)
	delta = NULL
	varlist = NULL
	knotsuse = list(); numknotsuse = NULL
	mslst = list()
	capm = 0
	capms = 0
	del1_ans = makedelta(xmat[, 1], shapes[1], numknots[1], knots[[1]], space = space[1], interp = FALSE)
	del1 = del1_ans$amat
	knotsuse[[1]] = del1_ans$knots
	mslst[[1]] = del1_ans$ms
	numknotsuse = c(numknotsuse, length(del1_ans$knots))
    m1 = nrow(del1)
#new code: record the number of columns of del1 if shapes0[1] == 17:
	if (shapes[1] == 17) {capms = capms + m1}
    var1 = 1:m1*0 + 1
	if (capl == 1) {
        delta = del1
        varlist = var1
    } else {
	    for (i in 2:capl) {
	        del2_ans = makedelta(xmat[,i], shapes[i], numknots[i], knots[[i]], space = space[i], interp = FALSE)
			del2 = del2_ans$amat
			knotsuse[[i]] = del2_ans$knots
			mslst[[i]] = del2_ans$ms
			numknotsuse = c(numknotsuse, length(del2_ans$knots))
			m2 = nrow(del2)
#new code: record the number of columns of del2 if shapes0[i] == 17:
			if (shapes[i] == 17) {capms = capms + m2}
			delta = rbind(del1, del2)
			varlist = 1:(m1 + m2)*0
			varlist[1:m1] = var1
			varlist[(m1 + 1):(m1 + m2)] = (1:m2)*0 + i
			var1 = varlist
			m1 = m1 + m2
			del1 = delta
	    }
	}
	if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0) {
		bigmat = rbind(t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
		np = sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)  + capms
	} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0) {
		bigmat = delta
		np = capms
	}
	rslt = list(bigmat=bigmat, np=np, varlist=varlist, knotsuse=knotsuse, mslst=mslst)
	return(rslt)
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
	sigma[k+1, index] = x[index] - t[1] + (t[2] - x[index])^3 / (t[2] - t[1])^2 / 3 - (t[2] - t[1]) / 3 #-(t[2]-t[1])^3/(t[2]-t[1])^2/3 #

	index = x > t[2]
	sigma[k+1, index] = x[index] - t[1] - (t[2] - t[1]) / 3 #-(t[2]-t[1])^3/(t[2]-t[1])^2/3#

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
			#rng=max(sigma[i,])
			#sigma[i,]=sigma[i,]/rng
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
        df.residual <- object$df.residual
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
        null.deviance <- object$null.deviance
        df <- object$df
        df.null <- object$df.null
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
        pvs <- object$pvs
        s.edf <- object$s.edf
        bstats <- object$bstats
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
            #new:
            rslt2 <- NULL
            if (!is.null(pvs)) {
                rslt2 <- data.frame("edf" = round(s.edf, 4), "mixture of Beta" = round(bstats, 4), "p.value" = round(pvs, 4))
                #rownames(rslt2) <- attributes(tms)$term.labels
                #debugged: check more
                if (!is.null(zid)) {
                    rownames(rslt2) <- (attributes(tms)$term.labels)[-(zid-1)]
                } else {
                    rownames(rslt2) <- (attributes(tms)$term.labels)
                    #if (any(object$shapes < 9) & length(pvs) < length(object$shapes)) {
                    #if (any(object$shapes < 9)) {
                    #  rm_id = which(object$shapes < 9)
                    #  rownames(rslt2) <- (attributes(tms)$term.labels)[-rm_id]
                    #} else {
                    #  rownames(rslt2) <- (attributes(tms)$term.labels)
                    #}
                }
            }

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
            rslt2 <- NULL
            if (!is.null(pvs)) {
                rslt2 <- data.frame("edf" = round(s.edf, 4), "mixture of Beta" = round(bstats, 4), "p.value" = round(pvs, 4))
                #debugged: check more
                if (!is.null(zid)) {
                    rownames(rslt2) <- (attributes(tms)$term.labels)[-(zid-1)]
                } else {
                    #new: remove the group name for cgamm
                    if(inherits(object, "cgamm")){
                      rownames(rslt2) <- rev(rev((attributes(tms)$term.labels))[-1])
                    } else {
                      rownames(rslt2) <- (attributes(tms)$term.labels)
                    }
                    #if (any(object$shapes < 9) & length(pvs) < length(object$shapes)) {
                    #if (any(object$shapes < 9)) {
                    #  rm_id = which(object$shapes < 9)
                    #  rownames(rslt2) <- (attributes(tms)$term.labels)[-rm_id]
                    #} else {
                    #  rownames(rslt2) <- (attributes(tms)$term.labels)
                    #}
                }
            }
        }
        #if (!is.null(sse0) & !is.null(sse1)) {
        #rslt2 <- cbind(SSE.Linear = sse0, SSE.Full = sse1)
        #new:
        #	rslt2 <- data.frame("SSE.Linear" = sse0, "SSE.Full" = sse1)
        #	rownames(rslt2)[1] <- ""
        #	ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs, cic = cic, null.deviance = null.deviance, df.null = df.null, deviance = deviance, df = df, df.residual = df.residual, family = family)
        #	class(ans) <- "summary.cgam"
        #	ans
        #} else {
        ans <- list(call = object$call, coefficients = rslt1, coefficients2 = rslt2,
                    zcoefs = coefs, cic = cic, null.deviance = null.deviance, df.null = df.null,
                    deviance = deviance, df = df, df.residual = df.residual, family = family)
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
  rslt1 <- NULL
  rslt2 <- NULL
  cic <- NULL
  is_param <- object$is_param
  is_fac <- object$is_fac
  vals <- object$vals
  tms <- object$tms
  if (!is.null(object$zmat)) {
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
    #new:
    cic <- object$cic
    rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
    #rownames(rslt1)[1] <- "(Intercept)"
    if (n >= 1) {
      lzid <- length(zid1)
      for (i in 1:lzid) {
        pos1 <- zid1[i]; pos2 <- zid2[i]
        for (j in pos1:pos2) {
          if (!is_param[i]) {
            rownames(rslt1)[j] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j], sep = "")
          } else {
            rownames(rslt1)[j] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")
          }
        }
      }
    }
    rslt1 <- as.matrix(rslt1)
    #ans <- list(call = object$call, coefficients = rslt1, coefficients2 = rslt2, zcoefs = coefs, cic = cic)
    #class(ans) <- "summary.wps"
    #ans
    #}
  }
  if (!is.null(object$pval)) {
    #new:
    pval <- object$pval; bval <- object$bval; edf <- object$edfc
    if (!is.null(pval)) {
      rslt2 <- data.frame("edf" = round(edf, 4), "Beta" = round(bval, 4), "p.value" = round(pval, 4))
      #debugged: check more
      zid <- object$zid
      if (!is.null(zid)) {
        rownames(rslt2) <- (attributes(tms)$term.labels)[-(zid-1)]
      } else {
        rownames(rslt2) <- (attributes(tms)$term.labels)
      }
    }
  }
  ans <- list(call = object$call, coefficients = rslt1, coefficients2 = rslt2, cic = cic)
  class(ans) <- "summary.wps"
  ans
}

#####################
#print.cgam #
#####################
print.cgam = function (x, ...)
{
    print(x$family)
    cat("Call:\n")
    print(x$call)
    cat("\n")
    if(!is.null(x$deviance)) {
        cat("Null deviance: ",round(x$null.deviance,4), "", "on", x$df.null, "", "degrees of freedom", "\n")
        cat("Residual deviance: ",round(x$deviance,4), "", "on", x$df.residual, "", "observed degrees of freedom", "\n")
    }
    if (!is.null(x$cic)) {
        cat("CIC: ", round(x$cic,4))
    }
    invisible(x)
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
            cat("(Dispersion parameter for binomial family taken to be 1)", "\n")
        }
        if (x$family$family == "poisson") {
            cat("(Dispersion parameter for poisson family taken to be 1)", "\n")
        }
        if (x$family$family == "gaussian") {
            cat("(Dispersion parameter for gaussian family taken to be ", round(x$deviance/x$df,4),")","\n", sep="")
        }
        cat("\n")
        cat("Null deviance: ",round(x$null.deviance,4), "", "on", x$df.null, "", "degrees of freedom", "\n")
        cat("Residual deviance: ",round(x$deviance,4), "", "on", x$df.residual, "", "observed degrees of freedom", "\n")
        #if (is.null(x$cic)) {
        #	message("CIC value is not available when there is no shape-restricted predictor")
        #} else {message("CIC: ", round(x$cic,4))}
        if (!is.null(x$coefficients2)) {
            cat("\n")
            #cat("Approximate significance of smooth terms: \n")
            cat("Approximate significance of constrained components: \n")
            printCoefmat(x$coefficients2, P.values = TRUE, has.Pvalue = TRUE)
        }
        if (!is.null(x$cic)) {
            cat("CIC: ", round(x$cic,4))
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
  if (!is.null(x$coefficients)) {
    #if (!is.null(x$se.beta)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Coefficients:")
    cat("\n")
    printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
    if (!is.null(x$cic)) {
      cat("CIC: ", round(x$cic,4), "\n", sep = "")
    }
  }
  if (!is.null(x$coefficients2)) {
    cat("\n")
    cat("Approximate significance of constrained surface component: \n")
    printCoefmat(x$coefficients2, P.values = TRUE, has.Pvalue = TRUE)
  }
}


##############
#predict.cgam#
##############
predict.cgam = function(object, newData, interval = c("none", "confidence", "prediction"), type = c("response", "link"), level = 0.95, n.mix = 500,...) {
  #print (is.data.frame(newData))
  #print (newData)
  #new:
  family = object$family
  cicfamily = CicFamily(family)
  muhat.fun = cicfamily$muhat.fun
  if (!inherits(object, "cgam")) {
    warning("calling predict.cgam(<fake-cgam-object>) ...")
  }
  if (missing(newData) || is.null(newData)) {
    #if (missing(newData)) {
    etahat = object$etahat
    muhat = muhat.fun(etahat, fml = family$family)
    #ans = list(fit = muhat, etahat = etahat, newbigmat = object$bigmat)
    ans = list(fit = muhat)
    return (ans)
  }
  if (!is.data.frame(newData)) {
    #newData = as.data.frame(newData)
    stop ("newData must be a data frame!")
  }
  #shapes = object$shapes
  #new: used for ci
  prior.w = object$prior.w
  vh = object$vh
  y = object$y
  #new: only use shapes for x != umb or tree
  shapes = object$shapesx
  np = object$d0; capm = object$capm; capms = object$capms; capk = object$capk; capt = object$capt; capu = object$capu
  #new:
  xid10 = object$xid1; xid20 = object$xid2;
  uid1 = object$uid1; uid2 = object$uid2; tid1 = object$tid1; tid2 = object$tid2
  #new:
  xmat0 = object$xmat0; knots0 = object$knots0; numknots0 = object$numknots0; sps0 = object$sps0; ms0 = object$ms0
  zmat = object$zmat; umb = object$umb; tr = object$tr
  #new:
  ztb = object$ztb; zid1 = object$zid1; zid2 = object$zid2; iz = 1
  bigmat = object$bigmat; umbrella.delta = object$umbrella.delta; tree.delta = object$tree.delta
  coefs = object$coefs; zcoefs = object$zcoefs; vcoefs = object$vcoefs; xcoefs0 = object$xcoefs; ucoefs = object$ucoefs; tcoefs = object$tcoefs
  tt = object$tms
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
  idx_s <- NULL
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
    #new:xmat is ordinal, xmat_s is smooth
    xmat = xmat_s = NULL
    newx = newx_s = NULL
    knots = NULL
    ms_x = ms = NULL
    sh_x = sh = NULL
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
        ms2 = matrix(0, nrow = nr, ncol = nc)
        for (i1 in 1:nc) {
          for (i2 in 1:nr) {
            ms2[i2, i1] = my_line(xp = nxi[i1], y = msi[i2, ][ord], x = xs, end = nx, start = 1)$yp
          }
        }
        deli = deli - ms2
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
  #muhat = newmuhat
  if ("none" %in% interval) {
    ans = list(fit = newmuhat, object = object)
    #new: add prediction interval
  } else if (interval == "confidence" | interval == "prediction") {
    n = ncol(bigmat)
    nc = ncol(xmat0)
    spl = splpr = NULL
    m_acc = 0
    object$ms0 -> ms0
    #sh = 17 first
    object$shapesx0 -> shapes
    #new:
    varlist = NULL
    for (i in 1:nc) {
      msi = ms0[[i]]
      shpi = shapes[i]
      ki = knots0[[i]]
      xi = xmat0[,i]
      xipr = newx0[,i]
      deli = makedelta(xi, shpi, knots = ki)
      #print (i)
      #if (i == 2) {
      #	stop (print (msi))
      #}
      spli = deli$amat
      varlist = c(varlist, rep(i, nrow(spli)))
      if (shpi >= 9) {
        dpri = makedelta(xipr, shpi, knots = ki, suppre = TRUE, interp = TRUE)
        splpri = dpri$amat
        if (shpi > 10 & shpi < 13) {
          xs = sort(xi)
          ord = order(xi)
          #nx = length(xi)
          #nx is n
          obs = 1:n
          nr = nrow(splpri)
          #nc = length(xipr)
          #rn is the length of a new vector
          ms2 = matrix(0, nrow = nr, ncol = rn)
          for (i1 in 1:rn) {
            for (i2 in 1:nr) {
              ms2[i2, i1] = my_line(xp = xipr[i1], y = msi[i2, ][ord], x = xs, end = n, start = 1)$yp
            }
          }
          splpri = splpri - ms2
        } else {
          splpri = splpri - msi
        }
      } else {splpri = pred_del(xi, shpi, xipr, msi)}
      mi = dim(spli)[1]
      m_acc = m_acc + mi
      spl = rbind(spl, spli)
      splpr = rbind(splpr, splpri)
    }
    capk = object$capk
    p = 1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)
    zmat = t(bigmat[1:p, ,drop = FALSE])
    if (family$family == "gaussian") {
      #nloop = 1000
      nloop = n.mix
      #nsec will be too big for ordinal, ignore for now
      nsec = 2^m_acc
      #print (dim(zmat))
      if (is.null(prior.w)) {
        prior.w = rep(1, n)
      }
      if (!all(prior.w == 1)) {
        for (i in 1:n) {
          spl[,i] = spl[,i] * sqrt(prior.w[i])
          zmat[i,] = zmat[i,] * sqrt(prior.w[i])
        }
      }
      if (!is.null(vh)) {
        wvec = 1/vh
      } else if (is.null(vh) & !is.null(prior.w)) {
        wvec = prior.w
      } else if (is.null(vh) & is.null(prior.w)) {
        wvec = rep(1, n)
      }
      if (!all(wvec == 1)) {
        for (i in 1:n) {
          spl[,i] = spl[,i] * sqrt(wvec[i])
          zmat[i,] = zmat[i,] * sqrt(wvec[i])
        }
      }
      yw = y * sqrt(wvec)
      #ans1 = coneB(yw, spl, zmat)
      ans1 = coneB(yw, t(spl), zmat)
      face = ans1$face
      muhat = ans1$yhat
      #muhat = object$muhat
      muhat0 = muhat / sqrt(wvec)
      #muhat0 = object$muhat
      #muhat = muhat0 * sqrt(prior.w)
      #print (ans1$df)
      #print (ans1$coef)
      sighat = sqrt(sum((yw - muhat)^2) / (n - 1.5*ans1$df))  ## varest is too large but "conservative"
      #sighat = sqrt(sum((yw - muhat)^2) / (n - 1.5*object$df_obs))
      nv = p
      ## estimate sector probabilties
      #set.seed(1)
      #sector = Matrix(1:nsec*0,nrow=1)
      #sector = big.matrix(1:nsec*0,nrow=1)
      #sector = 1:nsec*0
      #for (iloop in 1:nloop) {
      #    ysim = muhat0 + rnorm(n)*sighat
      #    ysim = ysim * sqrt(wvec)
      #    #ans = coneB(ysim, spl, zmat)
      #    ans = coneB(ysim, t(spl), zmat, face = face)
      #    cf = round(ans$coefs[(nv+1):(m_acc+nv)], 10)
      #    sec = 1:m_acc*0
      #    sec[cf > 0] = 1
      #    r = makebin(sec) + 1
      #    sector[r] = sector[r] + 1
      #}
      ### calculate the mixture hat matrix:
      #bsec = matrix(0, nrow=nsec,ncol=2)
      #bsec[,1] = 1:nsec
      #bsec[,2] = sector / nsec
      #keep = sector > 0
      #bsec = bsec[keep,]
      #ns = dim(bsec)[1]
      #bsec[,2] = bsec[,2] / sum(bsec[,2])


      #new: not use sector = 1:nsec*0; simulate for 1,000 times and record the faces used more than once
      # if times / nloop < 1e-3, then delete
      sector = NULL
      times = NULL
      df = NULL
      #first column shows sector; second column shows times
      for (iloop in 1:nloop) {
        ysim = muhat0 + rnorm(n)*sighat
        ysim = ysim * sqrt(prior.w)
        ans = coneB(ysim, t(spl), zmat, face = face)
        cf = round(ans$coefs[(nv+1):(m_acc+nv)], 10)
        sec = 1:m_acc*0
        sec[cf > 0] = 1
        sector = rbind(sector, sec)
        r = makebin(sec) + 1
        #sector[r] = sector[r] + 1
        if (iloop == 1) {
          df = rbind(df, c(r, 1))
        } else {
          if (r %in% df[,1]) {
            ps = which(df[,1] %in% r)
            df[ps,2] = df[ps,2] + 1
          } else {
            df = rbind(df, c(r, 1))
          }
        }
      }

      #remove times/sum(times) < 1e-3??
      sm_id = which((df[,2]/nloop) < 1e-3)
      if (any(sm_id)) {
        df = df[-sm_id, ,drop=FALSE]
      }

      #new:
      ns = nrow(df)
      bsec = df
      bsec[,2] = bsec[,2] / sum(bsec[,2])
      ord = order(bsec[,1])
      bsec = bsec[ord, ,drop=FALSE]

      ### calculate the mixture cov(alpha) matrix:
      obs = 1:m_acc;oobs = 1:(m_acc+nv)
      acov = matrix(0, nrow = m_acc+nv, ncol = m_acc+nv)
      for (is in 1:ns) {
        if (bsec[is,2] > 0) {
          jvec = getbin(bsec[is,1], m_acc)
          if (sum(jvec) == 1) {
            smat = cbind(zmat, spl[which(jvec==1),])
          } else if (sum(jvec) == 0) {
            smat = zmat
          } else {
            smat = cbind(zmat, t(spl[which(jvec==1),]))
          }
          acov1 = bsec[is,2]*solve(t(smat)%*%smat)
          acov2 = matrix(0,nrow=m_acc+nv,ncol=m_acc+nv)
          jobs = 1:(m_acc+nv)>0
          jm = 1:m_acc>0
          jm[obs[jvec==0]] = FALSE
          jobs[(nv+1):(m_acc+nv)] = jm
          nobs = oobs[jobs==TRUE]
          for (i in 1:sum(jobs)) {
            acov2[nobs[i],jobs] = acov1[i,]
          }
          acov = acov+acov2
        }
      }
      acov = acov*sighat^2
      #only xmatpr and splpr have interpolation points
      xmatpr = cbind(newv, t(splpr))
      #muhatpr = xmatpr %*% ans1$coefs
      #temp:
      muhatpr = xmatpr %*% object$coefs[1:ncol(xmatpr), ,drop=FALSE]
      #new: C.I. level
      mult = qnorm((1 - level)/2, lower.tail=FALSE)
      #hl = 2*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
      if (interval == "confidence") {
        hl = mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
      }
      if (interval == "prediction") {
        hl = mult*sqrt(sighat^2+diag(xmatpr%*%acov%*%t(xmatpr)))
      }
      #ans = list(fit = muhatpr, lower = muhatpr - hl, upper = muhatpr + hl)
    } else {
      ## initialize
      m = nrow(bigmat)
      #amat is different from the original code, because we put zmat at the first p columns in cgam
      #amat = matrix(0, nrow = (m-p), ncol = m)
      #for(i in (p+1):m){amat[i,i]=1}
      #the constraint is amat%*%bh >= 0; not constrain vmat
      amat = diag(m-p)
      zerom = matrix(0, nrow=nrow(amat), ncol=p)
      amat = cbind(zerom, amat)
      #alp0 = 1:m*0
      #eta0 = del%*%alp0
      #eta0 = 1:n*0
      #mu0 = 1:n*0 + 0.5
      z = 1:n*0
      ## get the dimension of the face
      deltil = t(bigmat)
      w = object$wt
      #print (any(w < 0))
      for(i in 1:m){deltil[,i]=deltil[,i]*sqrt(w)}
      umat = chol(crossprod(deltil))
      uinv = solve(umat)
      atil = amat%*%uinv
      bh = coef(object)
      #check
      etahat = object$etahat
      muhat = fitted(object)
      #z = etahat + (y - muhat)/muhat/(1-muhat)
      #ch = muhat > 1e-5 & muhat < 1-(1e-5)
      #z[ch] = etahat[ch] + (y[ch]-muhat[ch])/muhat[ch]/(1-muhat[ch])
      #z[!ch] = etahat[!ch]
      #w[ch] = muhat[ch]*(1-muhat[ch])
      #w[!ch] = 1e-5
      #z=eta0+(y-mu0)/mu0/(1-mu0)
      ch = muhat > 1e-5 & muhat < 1-(1e-5)
      z[ch] = etahat[ch] + (y[ch]-muhat[ch])/muhat[ch]/(1-muhat[ch])
      z[!ch] = etahat[!ch]
      w[ch] = muhat[ch]*(1-muhat[ch])
      w[!ch] = 1e-5
      ztil = z*sqrt(w)

      #z = etahat + (y - muhat) / family$variance(muhat)
      #print (family$variance(muhat)): sometimes give a value that is almost zero and it creates NA in thhat
      #ztil = z*sqrt(w)

      rw = round(amat%*% bh, 6) == 0
      #print (rw)
      #print (sum(rw))
      if(sum(rw) == 0){
        raj = 0
        edf = m
      } else {
        ajc = atil[rw,]
        raj = rankMatrix(t(ajc))[1]
        edf = min(m, 1.2*(m-raj))
      }
      ## thhat is final cone projection
      thhat = deltil %*% bh
      #print (thhat)
      #print (ztil)
      shat = sum((ztil - thhat)^2)/(n - edf)
      #print (shat)
      nloop = 100
      cmat = matrix(0, nrow = m, ncol = m)
      for(iloop in 1:nloop){
        zsim = thhat + rnorm(n)*shat
        zsimtil = t(uinv) %*% t(deltil) %*% zsim
        #if (iloop == 1) {
        asim = coneA(zsimtil, atil)
        #    face = asim$face
        #} else {
        #    if (length(face) >= 1) {
        #        asim = coneA(zsimtil, atil, face = face)
        #    } else{asim = coneA(zsimtil, atil)}
        #}
        ## get the matrix wm: cols span row space orthogonal to rows of A_J^c
        rw = round(atil %*% asim$thetahat, 6) == 0
        #print (shat)
        #print (rw)
        if(sum(rw) == 0){
          #pjmat = t(atil) %*% solve(atil %*% t(atil)) %*% atil
          pjmat = t(atil) %*% solve(tcrossprod(atil), atil)
        } else {
          ajc = atil[rw,]
          if(length(ajc) == m){
            wm = ajc/sqrt(sum(ajc^2))
            wm = matrix(wm, ncol=1)
          }else{
            aqr = qr(t(ajc))
            wm = qr.Q(aqr)[,1:aqr$rank]
          }
          pjmat = -wm%*%t(wm)
          for(i in 1:m){pjmat[i,i] = pjmat[i,i]+1}
        }
        cmat = cmat + uinv %*% pjmat %*% t(uinv)/nloop
      }
      cmat = cmat*shat
      #coveta = diag(del%*%cmat%*%t(del))
      xmatpr = cbind(newv, t(splpr))
      #coveta = diag(del%*%cmat%*%t(del))
      coveta = diag(xmatpr%*%cmat%*%t(xmatpr))
      etapr = xmatpr%*%bh
      #muhatpr = 1 - 1/(1+exp(etapr))
      muhatpr = family$linkinv(etapr)
      #new: C.I. level
      mult = qnorm((1 - level)/2, lower.tail=FALSE)
      #if (interval == "confidence") {
      uppeta = etapr + mult*sqrt(coveta)
      loweta = etapr - mult*sqrt(coveta)
      uppmu = family$linkinv(uppeta)
      lowmu = family$linkinv(loweta)
      #uppmu = 1 - 1/(1+exp(uppeta))
      #lowmu = 1 - 1/(1+exp(loweta))
      #}
    }
    #new: thvs are for individual c.i. bands and surfaces
    #smooth additive terms only, no umbrella or tree
    if (family$family == "gaussian") {
      dcoefs = object$coefs[(np - capms + 1):(np + capm),,drop=FALSE]
      etacomps = object$etacomps
      capl = nrow(etacomps)
      thvs_upp = thvs_lwr = matrix(0, nrow = capl, ncol = rn)
      ncon = 1
      # temp: constrained smooth terms covariance
      acov2 = acov[-c(1:np),-c(1:np)]
      for (i in 1:capl) {
        ddi = t(splpr[varlist == i, ,drop=FALSE])
        acovi = acov2[varlist == i, varlist == i]
        hli = mult*sqrt(diag(ddi %*% acovi %*% t(ddi)))
        ei = ddi%*%dcoefs[varlist == i, ,drop=FALSE]
        thvs_upp[i,] = ei + hli
        thvs_lwr[i,] = ei - hli
      }

      # order thvs back
      if (length(idx_s) > 0) {
        thvs0_upp = thvs_upp
        thvs0_upp[idx_s,] = thvs_upp[1:length(idx_s), ]
        thvs0_lwr = thvs_lwr
        thvs0_lwr[idx_s,] = thvs_lwr[1:length(idx_s), ]
        if (length(idx) > 0) {
          thvs0_upp[idx,] = thvs_upp[(1+length(idx_s)):capl, ]
          thvs0_lwr[idx,] = thvs_lwr[(1+length(idx_s)):capl, ]
        }
        thvs_upp = thvs0_upp
        thvs_lwr = thvs0_lwr
      }
      #new: problem when not gaussian
      if (object$sc_y) {
        for (i in 1:nrow(thvs_upp)) {
          thvs_upp[i,] = thvs_upp[i,] * object$sc
          thvs_lwr[i,] = thvs_lwr[i,] * object$sc
        }
      }
    }
    #new: NASA's idea for monotonic CI
    if (family$family == "gaussian") {
      lwr = muhatpr - hl
      upp = muhatpr + hl
    } else {
      lwr = loweta
      upp = uppeta
    }
    #for now, only handles one predictor
    if (length(object$shapes) == 1) {
      ord = order(newx0) #need a better way for > 1 predictor case

      lwr_tmp = lwr[ord]
      upp_tmp = upp[ord]
      lwr_tmp_u = unique(lwr_tmp)
      upp_tmp_u = unique(upp_tmp)
      check_ps_lwr = diff(lwr_tmp_u)
      check_ps_upp = diff(upp_tmp_u)
      if (object$shapes == 9) {
        ps_lwr = which(check_ps_lwr < 0)
        n_ps_lwr = length(ps_lwr)
        if(n_ps_lwr > 0){
          #NASA's for loop, the easier way is not working now....
          n = length(lwr_tmp)
          for(i in 1:(n - 1)){
            if(lwr_tmp[i + 1] < lwr_tmp[i]){
              lwr_tmp[i + 1] = lwr_tmp[i]
            }
          }
        }
        #lwr_tmp[ps_lwr + 1] = lwr_tmp[ps_lwr]

        ps_upp = which(check_ps_upp < 0)
        n_ps_upp = length(ps_upp)
        if(n_ps_upp > 0) {
          n = length(upp_tmp)
          for(i in (n - 1):1){
            if(upp_tmp[i] > upp_tmp[i + 1]){
              upp_tmp[i] = upp_tmp[i + 1]
            }
          }
        }
        #upp_tmp[ps_upp] = upp_tmp[ps_upp + 1]
      }
      if (object$shapes == 10) {
        ps_lwr = which(check_ps_lwr > 0)
        n_ps_lwr = length(ps_lwr)
        if(n_ps_lwr > 0){
          n = length(lwr_tmp)
          for(i in (n - 1):1){
            if(lwr_tmp[i] < lwr_tmp[i + 1]){
              lwr_tmp[i] = lwr_tmp[i + 1]
            }
          }
        }
        #lwr_tmp[ps_lwr] = lwr_tmp[ps_lwr + 1]

        ps_upp = which(check_ps_upp > 0)
        n_ps_upp = length(ps_upp)
        if(n_ps_upp > 0) {
          n = length(upp_tmp)
          for(i in 1:(n - 1)){
            if(upp_tmp[i + 1] > upp_tmp[i]){
              upp_tmp[i + 1] = upp_tmp[i]
            }
          }
        }
        #upp_tmp[ps_upp + 1] = upp_tmp[ps_upp]
      }

      lwr = lwr_tmp[order(ord)]
      upp = upp_tmp[order(ord)]
      #loweta = lwr
      #uppeta = upp

      #lwr = lwr_tmp
      #upp = upp_tmp
    }

    if ("response" %in% type) {
      lwr = family$linkinv(lwr)
      upp = family$linkinv(upp)
    }

    if (family$family == "gaussian") {
      ans = list(fit = muhatpr, lower = lwr, upper = upp, acov = acov, object = object, mult = mult, thvs_upp = thvs_upp, thvs_lwr = thvs_lwr)
      #ans = list(fit = muhatpr, lower = muhatpr - hl, upper = muhatpr + hl, acov = acov, object = object, mult = mult, thvs_upp = thvs_upp, thvs_lwr = thvs_lwr)
    } else { #else if (family$family == "binomial") {
      if ("response" %in% type) {
        ans = list(fit = muhatpr, lower = lwr, upper = upp, object = object, mult = mult)
        #ans = list(fit = muhatpr, lower = lowmu, upper = uppmu, object = object, mult = mult)
      } else if ("link" %in% type) {
        ans = list(fit = etapr, lower = lwr, upper = upp, object = object, mult = mult)
        #ans = list(fit = etapr, lower = loweta, upper = uppeta, object = object, mult = mult)
      }
    }
  }
  class(ans) = "cgamp"
  return (ans)
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



##############
#predict.wps#
# >= pairs + additive + z
##############
predict.wps = function(object, newData, interval = c("none", "confidence", "prediction"), type = c("response", "link"), level = 0.95, ...) {
    #new:
    family = object$family
    cicfamily = CicFamily(family)
    muhat.fun = cicfamily$muhat.fun
    if (!inherits(object, "wps")) {
        warning("calling predict.wps(<fake-wps-object>) ...")
    }
    if (missing(newData) || is.null(newData)) {
        #if (missing(newData)) {
        #etahat = object$etahat
        #muhat = muhat.fun(etahat, fml = family$family)
        #ans = list(fit = muhat, etahat = etahat, newbigmat = object$bigmat)
        ans = list(fit = object$muhat)
        return (ans)
    }
    if (!is.data.frame(newData)) {
        #newData = as.data.frame(newData)
        stop ("newData must be a data frame!")
    }
    #shapes = object$shapes
    #new: used for ci
    ahat = object$coefs
    coef_wp = object$coef_wp
    prior.w = object$prior.w
    pen = object$pen
    dmat = object$dmat
    amat = object$amat
    y = object$y
    ##zmat includes one vector and also conc, conv, shape=17, additive
    #no more
    zmat = object$zmat
    #vmat includes one vector and also conc, conv, shape=17
    #vmat = object$vmat
    #delta includes everything
    delta = object$delta
    tt = object$tms
    Terms = delete.response(tt)
    #model.frame will re-organize newData in the original order in formula
    m = model.frame(Terms, newData)

    newdata = m
    #print (head(newdata))
    #new:
    newx0 = NULL; newxv = NULL; newx = NULL; newx_s = NULL; newu = NULL; newt = NULL; newz = NULL; newv = NULL
    #newz = list();
    ztb = object$ztb; zid1 = object$zid1; zid2 = object$zid2; iz = 1
    rn = nrow(newdata)
    #at least one pair in decrs and ks_wps
    decrs = object$decrs
    ks_wps = object$ks_wps
    sps_wps = list()
    dcss = list()
    x1_wps = x2_wps = list()
    iwps = 0
    varlist = object$varlist_wps
    #varlist = varlist[-1]
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
                    klvls = length(ztbi)
                    if (klvls > 1) {
                        newimat = matrix(0, nrow = rn, ncol = klvls-1)
                        for (i1 in 1:rn) {
                            if (newdata[i1,i] != ztbi[1]) {
                                id_col = which(ztbi %in% newdata[i1,i]) - 1
                                newimat[i1,id_col] = 1
                            }
                        }
                        newdatai = newimat
                    }
                }
            } else {
                newdatai = newdata[,i]
            }
            newz = cbind(newz, newdatai)
            iz = iz + 1
        }
        if (length(attributes(newdata[,i])$decreasing)>1) {
            iwps = iwps + 1
            #sps_wp = attributes(newdata[, i])$space
            #sps_wps[[iwps]] = sps_wp
            dcs = attributes(newdata[, i])$decreasing
            dcss[[iwps]] = dcs
            x1_wp = (newdata[, i])[, 1]
            x1_wps[[iwps]] = x1_wp
            x2_wp = (newdata[, i])[, 2]
            x2_wps[[iwps]] = x2_wp
        }
        if (is.numeric(attributes(newdata[,i])$shape)) {
            newx0 = cbind(newx0, newdata[,i])
            if ((attributes(newdata[,i])$shape > 2 & attributes(newdata[,i])$shape < 5) | (attributes(newdata[,i])$shape > 10 & attributes(newdata[,i])$shape < 13)) {
                newxv = cbind(newxv, newdata[,i])
            }
        }
        #if (is.character(attributes(newdata[,i])$shape)) {
        #    if (attributes(newdata[,i])$shape == "tree") {
        #        newt = cbind(newt, newdata[,i])
        #    }
        #    if (attributes(newdata[,i])$shape == "umbrella") {
        #        newu = cbind(newu, newdata[,i])
        #    }
        #}
    }
    #conc, or conv + shape=17 + other shapes + zmat (including one) + wps
    m_acc_add = 0
    spl_add = splpr_add = NULL
    if (!is.null(newx0)) {
        shapes = object$shapes
        #put s=17 first as in wps.fit
        if (any(shapes == 17)) {
            kshapes = length(shapes)
            obs = 1:kshapes
            idx_s = obs[which(shapes == 17)]; idx = obs[which(shapes != 17)]
            newx1 = newx0
            shapes0 = 1:kshapes*0
            newx1[,1:length(idx_s)] = newx0[,idx_s]
            shapes0[1:length(idx_s)] = shapes[idx_s]
            if (length(idx) > 0) {
                newx1[,(1 + length(idx_s)):kshapes] = newx0[,idx]
                shapes0[(1 + length(idx_s)):kshapes] = shapes[idx]
            }
            newx0 = newx1; shapes = shapes0
        }
        nc = ncol(newx0)
        ms = object$ms
        knotsuse_add = object$knotsuse_add
        #17 is first
        xmat0_add = object$xmat0_add

        for (i in 1:nc) {
            msi = ms[[i]]
            shpi = shapes[i]
            ki = knotsuse_add[[i]]
            xi = xmat0_add[,i]
            xipr = newx0[,i]
            if (any(xipr > max(ki)) | any(xipr < min(ki))) {
                stop ("No extrapolation is allowed in cgam prediction!")
            }
            #if (shpi >= 9) {
            dpri = makedelta(xipr, shpi, knots = ki, suppre = TRUE, interp = TRUE)
            splpri = dpri$amat
            if (shpi > 10 & shpi < 13) {
                xs = sort(xi)
                ord = order(xi)
                #nx = length(xi)
                #nx is n
                obs = 1:n
                nr = nrow(splpri)
                #nc = length(xipr)
                #rn is the length of a new vector
                ms2 = matrix(0, nrow = nr, ncol = rn)
                for (i1 in 1:rn) {
                    for (i2 in 1:nr) {
                        ms2[i2, i1] = my_line(xp = xipr[i1], y = msi[i2, ][ord], x = xs, end = n, start = 1)$yp
                    }
                }
                splpri = splpri - ms2
            } else {
                splpri = splpri - msi
            }
            #} else {splpri = pred_del(xi, shpi, xipr, msi)}
            mi = dim(splpri)[1]
            m_acc_add = m_acc_add + mi
            #splpr_add  = rbind(splpr_add, splpri)
            splpr_add  = cbind(splpr_add, t(splpri))
        }
    }
    #temp
    #newv = cbind(1:rn*0 + 1, newz, newxv)
    #newv = cbind(newxv, splpr_add, 1:rn*0 + 1, newz)
    newv = cbind(newxv, splpr_add, newz)
    np = ncol(newv)
    n = length(y)
    splpr = NULL
    splpr_lst = mus = list()
    m_acc = 0
    #loop through total number of wps pairs
    for (i in 1:iwps) {
        k1 = ks_wps[[i]][[1]]
        k2 = ks_wps[[i]][[2]]
        nxi1 = x1_wps[[i]]
        nxi2 = x2_wps[[i]]
        if (any(nxi1 > max(k1)) | any(nxi1 < min(k1)) | any(nxi2 > max(k2)) | any(nxi2 < min(k2))) {
            stop ("No extrapolation is allowed in cgam prediction!")
        }
        #sps_wp = sps_wps[[i]]
        #space, m1_0, m2_0 don't matter
        deli = makedelta_wps(nxi1, nxi2, m1_0 = length(k1), m2_0 = length(k2), k1, k2, space = c("E",
        "E"), decreasing = decrs[[i]], interp = TRUE)
        splpri = deli$delta
        #remove the constant vector
        #if (i > 1) {
        #splpri = splpri[,-1]
        #}
        mi = dim(splpri)[2]
        m_acc = m_acc + mi
        #spl = rbind(spl, spli)
        splpr = cbind(splpr, splpri)
        splpr_lst[[i]] = splpri
        #ignore z for now
        #mui = cbind(1, splpri) %*% c(ahat[1], ahat[which(varlist == i)+np])
        #print (dim(splpri))
        #print (length(ahat))
        mui = splpri %*% ahat[which(varlist == i)]+np
        mus[[i]] = muhat.fun(mui, fml = family$family)
    }
    m_acc = m_acc + m_acc_add
    xmatpr = cbind(newv, splpr)
    muhatpr = xmatpr %*% ahat
    muhatpr = muhat.fun(muhatpr, fml = family$family)
    if ("none" %in% interval) {
        ans = list(fit = muhatpr, xmatpr = xmatpr, spls = splpr_lst, mus = mus)
        #new: add prediction interval
    } else if (interval == "confidence" | interval == "prediction") {
        if (family$family == "gaussian") {
            #bmat includes constant + zmat, not weighted
            bmat = delta
            np = object$d0
            nloop = 100
            #nsec will be too big for ordinal, ignore for now
            #nsec = 2^m_acc
            #print (dim(zmat))
            if (is.null(prior.w)) {
                prior.w = rep(1, n)
            }
            if (!all(prior.w == 1)) {
                for (i in 1:n) {
                    bmat[,i] = bmat[,i] * sqrt(prior.w[i])
                    #spl[,i] = spl[,i] * sqrt(prior.w[i])
                    zmat[i,] = zmat[i,] * sqrt(prior.w[i])
                    vmat[i,] = vmat[i,] * sqrt(prior.w[i])
                }
            }
            #spl is constrained splines: additive + wps
            #spl = bmat[,-c(1:np),drop=FALSE]
            #spl = bmat
            #if (np >= 1) {
            #  spl = bmat[,-c(1:np),drop=FALSE]
            #}
            #muhat0 = (object$muhat)/sqrt(prior.w)
            #muhat0: unweighted
            muhat0 = object$muhat
            sighat = (object$sig2hat)^(1/2)
            nv = np
            ## estimate sector probabilties
            #new: not use sector = 1:nsec*0; simulate for 1,000 times and record the faces used more than once
            # if times / nloop < 1e-3, then delete
            sector = NULL
            times = NULL
            df = NULL
            faces = list()
            #set.seed(1)
            #first column shows sector; second column shows times
            for (iloop in 1:nloop) {
                ysim = muhat0 + rnorm(n)*sighat
                ysim = ysim * sqrt(prior.w)
                #ans = coneB(ysim, t(spl), zmat, face = face)
                #cf = round(ans$coefs[(nv+1):(m_acc+nv)], 10)
                qv0 = crossprod(bmat)
                dv0 = crossprod(dmat)
                #pen2 = find_pen(aims = edf, Q = qv0, B = bmat, D = dv0, PNT = TRUE, Y = ysim, D0 = dmat)

                #nxi1 = x1_wps[[1]]
                #nxi2 = x2_wps[[1]]
                #mat = cbind(1, nxi1, nxi2, nxi1*nxi2)
                #mu_para = mat %*% solve(t(mat) %*% mat) %*% t(mat) %*% ysim
                #ssr = sum((ysim - mu_para)^2)
                #sc = ssr / (n - ncol(mat))
                #pen = max(1e-6, 1 * sc)

                qv = qv0 + pen * dv0
                cv = crossprod(bmat, ysim)
                #amat = diag(m)
                #amat[1, ] = 0
                #ans = qprog(qv, cv, amat, 1:nrow(amat)*0, msg = FALSE)
                #cf = round(ans$thetahat[(np+1):(m_acc+np)],10)
                #sec = 1:m_acc*0
                #sec[cf > 0] = 1
                umat = chol(qv)
                uinv = solve(umat)
                atil = amat %*% uinv
                cvec = t(uinv) %*% t(bmat) %*% ysim
                ans = coneA(cvec, atil, msg = FALSE)
                #phihat = ans$thetahat
                #ahat = uinv %*% phihat

                #new: use polar cone instead
                face = ans$face
                faces[[iloop]] = face
                sec = 1:nrow(amat)*0
                sec[face] = 1

                #sec = 1:nrow(amat)*0+1
                #cf = round(ans$thetahat,10)
                #sec[cf > 0] = 0

                r = makebin(sec) + 1
                if (iloop == 1) {
                    df = rbind(df, c(r, 1))
                    sector = rbind(sector, sec)
                } else {
                    if (r %in% df[,1]) {
                        ps = which(df[,1] %in% r)
                        df[ps,2] = df[ps,2] + 1
                    } else {
                        df = rbind(df, c(r, 1))
                        sector = rbind(sector, sec)
                    }
                }
            }
            #remove times/sum(times) < 1e-3??
            sm_id = which((df[,2]/nloop) < 1e-3)
            if (any(sm_id)) {
                df = df[-sm_id, ,drop=FALSE]
                sector = sector[-sm_id, ,drop=FALSE]
            }
            #new:
            ns = nrow(df)
            bsec = df
            bsec[,2] = bsec[,2] / sum(bsec[,2])
            ord = order(bsec[,1])
            bsec = bsec[ord, ,drop=FALSE]
            sector = sector[ord, ,drop=FALSE]

            ### calculate the mixture cov(alpha) matrix:
            #spl = t(spl)
            #obs = 1:m_acc;oobs = 1:(m_acc+nv)
            #acov = matrix(0, nrow = m_acc+nv, ncol = m_acc+nv)
            #for (is in 1:ns) {
            #    if (bsec[is,2] > 0) {
            #        jvec = getbin(bsec[is,1], m_acc)
            #        if (sum(jvec) == 1) {
            #            smat = cbind(vmat, spl[,which(jvec==1),drop=FALSE])
            #        } else if (sum(jvec) == 0) {
            #            smat = vmat
            #        } else {
            #            smat = cbind(vmat, spl[,which(jvec==1),drop=FALSE])
            #        }
            #        acov1 = bsec[is,2]*solve(t(smat)%*%smat)
            #        acov2 = matrix(0,nrow=m_acc+nv,ncol=m_acc+nv)
            #        jobs = 1:(m_acc+nv)>0
            #        jm = 1:m_acc>0
            #        jm[obs[jvec==0]] = FALSE
            #        jobs[(nv+1):(m_acc+nv)] = jm
            #        nobs = oobs[jobs==TRUE]
            #        for (i in 1:sum(jobs)) {
            #            acov2[nobs[i],jobs] = acov1[i,]
            #        }
            #        acov = acov+acov2
            #    }
            #}
            #get covariance of phihat first
            #umat = chol(qv)
            #uinv = solve(umat)
            #atil = amat %*% uinv
            dp = -atil
            for (i in 1:nrow(dp)) {
                dpi = dp[i,,drop=F]
                dpi = dpi/sqrt(tcrossprod(dpi)[1])
                dp[i,]=dpi
            }
            dp = t(dp)
            imat = diag(nrow(qv))
            acov0 = matrix(0, nrow=nrow(qv), ncol=ncol(qv))
            for (is in 1:ns) {
                if (bsec[is,2] > 0) {
                    #jvec = getbin(bsec[is,1], ncol(dp))
                    jvec = sector[is, ]
                    if (sum(jvec) == 0) {
                        pmat_is = imat
                    } else {
                        smat = dp[,which(jvec==1),drop=FALSE]
                        #if (qr(smat)$rank < ncol(smat)) {
                        #save(smat,file='smat.Rda')
                        #fc = bsec[is,1]
                        #pj = bsec[is,2]
                        #save(jvec, file='jvec.Rda')
                        #save(fc, file='fc.Rda')
                        #save(pj, file='pj.Rda')
                        #save(faces, file='faces.Rda')
                        #save(bsec, file='bsec.Rda')
                        #smat = qr.Q(qr(smat), complete = TRUE)[, 1:(qr(smat)$rank), drop = FALSE]
                        #}
                        pmat_is_p = smat %*% solve(crossprod(smat), t(smat))
                        pmat_is = (imat-pmat_is_p)
                    }
                    acov0 = acov0 + bsec[is,2]*pmat_is%*%t(uinv)%*%qv0%*%uinv%*%pmat_is
                }
            }
            acov = uinv%*%acov0%*%t(uinv)*sighat^2
            #only xmatpr and splpr have interpolation points
            #xmatpr = cbind(newv, splpr)
            #muhatpr = xmatpr %*% object$coefs
            #temp:
            #muhatpr = xmatpr %*% object$coefs[1:ncol(xmatpr), ,drop=FALSE]
            #new: C.I. level
            mult = qnorm((1 - level)/2, lower.tail=FALSE)
            #hl = 2*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
            #hl = mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
            #ans = list(fit = muhatpr, lower = muhatpr - hl, upper = muhatpr + hl)
        }
        if (interval == "confidence") {
            hl = mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
        }
        if (interval == "prediction") {
            hl = mult*sqrt(sighat^2+diag(xmatpr%*%acov%*%t(xmatpr)))
        }
        ans = list(fit = muhatpr, lower = muhatpr - hl, upper = muhatpr + hl, object = object, acov = acov, mult = mult, newz = newz)
    }
    class(ans) = "wpsp"
    return (ans)
}

##############
#predict.trispl
#not finished: 1 pair + z only
##############
predict.trispl = function(object, newData, interval = c("none", "confidence", "prediction"), type = c("response", "link"), level = 0.95, ...) {
    #new:
    family = object$family
    cicfamily = CicFamily(family)
    muhat.fun = cicfamily$muhat.fun
    if (!inherits(object, "trispl")) {
        warning("calling predict.trispl(<fake-trispl-object>) ...")
    }
    if (missing(newData) || is.null(newData)) {
        #if (missing(newData)) {
        #etahat = object$etahat
        #muhat = muhat.fun(etahat, fml = family$family)
        #ans = list(fit = muhat, etahat = etahat, newbigmat = object$bigmat)
        ans = list(fit = object$muhat)
        return (ans)
    }
    if (!is.data.frame(newData)) {
        #newData = as.data.frame(newData)
        stop ("newData must be a data frame!")
    }
    #shapes = object$shapes
    #new: used for ci
    ahat = object$coefs
    coef_tri = object$coef_tri

    prior.w = object$prior.w
    pen = object$pen

    dmat = object$pmatc
    amat = object$amatc
    #delta includes everything
    delta = object$dmatc
    capk_lst = object$capk_lst
    trimat_lst = object$trimat_lst

    y = object$y
    n = length(y)
    #zmat does't include one vector in trispl
    zmat = object$zmat
    #vmat includes one vector and also conc, conv, shape=17
    #vmat = object$vmat

    tt = object$tms
    Terms = delete.response(tt)
    #model.frame will re-organize newData in the original order in formula
    m = model.frame(Terms, newData)

    newdata = m
    #print (head(newdata))
    #new:
    newx0 = NULL; newxv = NULL; newx = NULL; newx_s = NULL; newu = NULL; newt = NULL; newz = NULL; newv = NULL
    #newz = list();
    ztb = object$ztb; zid1 = object$zid1; zid2 = object$zid2; iz = 1
    rn = nrow(newdata)
    #at least one pair in decrs and ks_wps
    #decrs = object$decrs
    #ks_wps = object$ks_wps
    #sps_wps = list()
    #dcss = list()
    #x1_wps = x2_wps = list()
    iwps = 0
    #varlist = object$varlist_wps
    #varlist = varlist[-1]

    itri = 0
    convexity = object$cvss
    cvss = list()
    x1_tris = x2_tris = list()
    icvs = 0
    varlist_tri = object$varlist_tri
    ks_tri = object$knots_lst

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
                    klvls = length(ztbi)
                    if (klvls > 1) {
                        newimat = matrix(0, nrow = rn, ncol = klvls-1)
                        for (i1 in 1:rn) {
                            if (newdata[i1,i] != ztbi[1]) {
                                id_col = which(ztbi %in% newdata[i1,i]) - 1
                                newimat[i1,id_col] = 1
                            }
                        }
                        newdatai = newimat
                    }
                }
            } else {
                newdatai = newdata[,i]
            }
            newz = cbind(newz, newdatai)
            iz = iz + 1
        } else if (length(attributes(newdata[,i])$decreasing)>1) {
            iwps = iwps + 1
            #sps_wp = attributes(newdata[, i])$space
            #sps_wps[[iwps]] = sps_wp
            dcs = attributes(newdata[, i])$decreasing
            dcss[[iwps]] = dcs
            x1_wp = (newdata[, i])[, 1]
            x1_wps[[iwps]] = x1_wp
            x2_wp = (newdata[, i])[, 2]
            x2_wps[[iwps]] = x2_wp
        } else if (is.numeric(attributes(newdata[,i])$shape)) {
            newx0 = cbind(newx0, newdata[,i])
            if ((attributes(newdata[,i])$shape > 2 & attributes(newdata[,i])$shape < 5) | (attributes(newdata[,i])$shape > 10 & attributes(newdata[,i])$shape < 13)) {
                newxv = cbind(newxv, newdata[,i])
            }
        } else if (is.character(attributes(newdata[,i])$shape) && (attributes(newdata[,i])$categ == "tri")) {
            itri = itri + 1
            cvs = attributes(newdata[,i])$cvs
            cvss[[itri]] = cvs
            x1_tri = (newdata[,i])[, 1]
            x1_tris[[itri]] = x1_tri
            x2_tri = (newdata[,i])[, 2]
            x2_tris[[itri]] = x2_tri
        }
    }
    #conc, or conv + shape=17 + other shapes + zmat (including one) + wps
    m_acc_add = 0
    spl_add = splpr_add = NULL
    if (!is.null(newx0)) {
        shapes = object$shapes
        #put s=17 first as in wps.fit
        if (any(shapes == 17)) {
            kshapes = length(shapes)
            obs = 1:kshapes
            idx_s = obs[which(shapes == 17)]; idx = obs[which(shapes != 17)]
            newx1 = newx0
            shapes0 = 1:kshapes*0
            newx1[,1:length(idx_s)] = newx0[,idx_s]
            shapes0[1:length(idx_s)] = shapes[idx_s]
            if (length(idx) > 0) {
                newx1[,(1 + length(idx_s)):kshapes] = newx0[,idx]
                shapes0[(1 + length(idx_s)):kshapes] = shapes[idx]
            }
            newx0 = newx1; shapes = shapes0
        }
        nc = ncol(newx0)
        ms = object$ms
        knotsuse_add = object$knotsuse_add
        #17 is first
        xmat0_add = object$xmat0_add

        for (i in 1:nc) {
            msi = ms[[i]]
            shpi = shapes[i]
            ki = knotsuse_add[[i]]
            xi = xmat0_add[,i]
            xipr = newx0[,i]
            if (any(xipr > max(ki)) | any(xipr < min(ki))) {
                stop ("No extrapolation is allowed in cgam prediction!")
            }
            #if (shpi >= 9) {
            dpri = makedelta(xipr, shpi, knots = ki, suppre = TRUE, interp = TRUE)
            splpri = dpri$amat
            if (shpi > 10 & shpi < 13) {
                xs = sort(xi)
                ord = order(xi)
                #nx = length(xi)
                #nx is n
                obs = 1:n
                nr = nrow(splpri)
                #nc = length(xipr)
                #rn is the length of a new vector
                ms2 = matrix(0, nrow = nr, ncol = rn)
                for (i1 in 1:rn) {
                    for (i2 in 1:nr) {
                        ms2[i2, i1] = my_line(xp = xipr[i1], y = msi[i2, ][ord], x = xs, end = n, start = 1)$yp
                    }
                }
                splpri = splpri - ms2
            } else {
                splpri = splpri - msi
            }
            #} else {splpri = pred_del(xi, shpi, xipr, msi)}
            mi = dim(splpri)[1]
            m_acc_add = m_acc_add + mi
            #splpr_add  = rbind(splpr_add, splpri)
            splpr_add  = cbind(splpr_add, t(splpri))
        }
    }
    #temp
    #newv = cbind(1:rn*0 + 1, newz, newxv)
    #for tri.fit, don't include the constant vector
    #newv = cbind(newxv, splpr_add, 1:rn*0 + 1, newz)
    #newv = cbind(newxv, splpr_add, newz)
    #np = ncol(newv)
    #m_acc = 0

    splpr = NULL
    splpr_lst = mus = list()
    #loop through total number of wps pairs
    for (i in 1:itri) {
        k1 = ks_tri[[i]][,1]
        m1 = length(k1)
        k2 = ks_tri[[i]][,2]
        m2 = length(k2)
        nxi1 = x1_tris[[i]]
        nxi2 = x2_tris[[i]]
        trimat = trimat_lst[[i]]
        capk = capk_lst[[i]]
        if (any(nxi1 > max(k1)) | any(nxi1 < min(k1)) | any(nxi2 > max(k2)) | any(nxi2 < min(k2))) {
            stop ("No extrapolation is allowed in cgam prediction!")
        }
        #sps_wp = sps_wps[[i]]
        #space, m1_0, m2_0 don't matter
        deli = makedelta_tri(nxi1, nxi2, m1, m2, k1, k2, trimat, capk, space = c("E",
        "E"), cvs = convexity[[i]], interp = TRUE)
        splpri = deli
        #mi = dim(splpri)[2]
        #m_acc = m_acc + mi
        #spl = rbind(spl, spli)
        splpr = cbind(splpr, splpri)
        splpr_lst[[i]] = splpri
        #ignore z for now
        #mui = cbind(1, splpri) %*% c(ahat[1], ahat[which(varlist == i)+np])
        #mus[[i]] = muhat.fun(mui, fml = family$family)
    }
    # m_acc = m_acc + m_acc_add

    #in tri.fit: additive + trispl + z + wps (ignore)
    xmatpr = cbind(splpr_add, splpr, newz)
    muhatpr = xmatpr %*% ahat
    muhatpr = muhat.fun(muhatpr, fml = family$family)
    if ("none" %in% interval) {
        ans = list(fit = muhatpr, xmatpr = xmatpr, spls = splpr_lst, mus = mus, object = object)
        #new: add prediction interval
    } else if (interval == "confidence" | interval == "prediction") {
        if (family$family == "gaussian") {
            #bmat includes constant + zmat, not weighted
            bmat = delta
            np = object$d0
            nloop = 500
            #nsec will be too big for ordinal, ignore for now
            #nsec = 2^m_acc
            #print (dim(zmat))
            if (is.null(prior.w)) {
                prior.w = rep(1, n)
            }
            if (!all(prior.w == 1)) {
                for (i in 1:n) {
                    bmat[,i] = bmat[,i] * sqrt(prior.w[i])
                    #spl[,i] = spl[,i] * sqrt(prior.w[i])
                    #zmat[i,] = zmat[i,] * sqrt(prior.w[i])
                    #vmat[i,] = vmat[i,] * sqrt(prior.w[i])
                }
            }
            #spl is constrained splines: additive + wps
            #spl = bmat[,-c(1:np),drop=FALSE]
            #muhat0 = (object$muhat)/sqrt(prior.w)
            #muhat0: unweighted
            muhat0 = object$muhat
            sighat = (object$sig2hat)^(1/2)
            #no use:
            nv = np
            ## estimate sector probabilties
            #new: not use sector = 1:nsec*0; simulate for 1,000 times and record the faces used more than once
            # if times / nloop < 1e-3, then delete
            sector = NULL
            times = NULL
            df = NULL
            #faces = list()
            #set.seed(1)
            #first column shows sector; second column shows times
            for (iloop in 1:nloop) {
                #print (iloop)
                ysim = muhat0 + rnorm(n)*sighat
                ysim = ysim * sqrt(prior.w)
                #ans = coneB(ysim, t(spl), zmat, face = face)
                #cf = round(ans$coefs[(nv+1):(m_acc+nv)], 10)
                qv0 = crossprod(bmat)
                dv0 = crossprod(dmat)
                #pen2 = find_pen(aims = edf, Q = qv0, B = bmat, D = dv0, PNT = TRUE, Y = ysim, D0 = dmat)

                qv = qv0 + pen * dv0
                cv = crossprod(bmat, ysim)
                umat = chol(qv)
                uinv = solve(umat)
                atil = amat %*% uinv
                cvec = t(uinv) %*% t(bmat) %*% ysim
                ans = coneA(cvec, atil, msg = FALSE)
                #phihat = ans$thetahat
                #ahat = uinv %*% phihat

                #new: use polar cone instead
                face = ans$face
                #faces[[iloop]] = face
                sec = 1:nrow(amat)*0
                sec[face] = 1

                #sec = 1:nrow(amat)*0+1
                #cf = round(ans$thetahat,10)
                #sec[cf > 0] = 0

                r = makebin(sec) + 1
                if (iloop == 1) {
                    df = rbind(df, c(r, 1))
                    sector = rbind(sector, sec)
                } else {
                    if (r %in% df[,1]) {
                        ps = which(df[,1] %in% r)
                        df[ps,2] = df[ps,2] + 1
                    } else {
                        df = rbind(df, c(r, 1))
                        sector = rbind(sector, sec)
                    }
                }
            }
            #remove times/sum(times) < 1e-3??
            sm_id = which((df[,2]/nloop) < 1e-3)
            if (any(sm_id)) {
                df = df[-sm_id, ,drop=FALSE]
                sector = sector[-sm_id, ,drop=FALSE]
            }
            #new:
            ns = nrow(df)
            bsec = df
            bsec[,2] = bsec[,2] / sum(bsec[,2])
            ord = order(bsec[,1])
            bsec = bsec[ord, ,drop=FALSE]
            sector = sector[ord, ,drop=FALSE]

            ### calculate the mixture cov(alpha) matrix:
            dp = -atil
            for (i in 1:nrow(dp)) {
                dpi = dp[i,,drop=F]
                dpi = dpi/sqrt(tcrossprod(dpi)[1])
                dp[i,]=dpi
            }
            dp = t(dp)
            imat = diag(nrow(qv))
            acov0 = matrix(0, nrow=nrow(qv), ncol=ncol(qv))
            for (is in 1:ns) {
                if (bsec[is,2] > 0) {
                    #jvec = getbin(bsec[is,1], ncol(dp))
                    jvec = sector[is, ]
                    if (sum(jvec) == 0) {
                        pmat_is = imat
                    } else {
                        smat = dp[,which(jvec==1),drop=FALSE]
                        pmat_is_p = smat %*% solve(crossprod(smat), t(smat))
                        pmat_is = (imat-pmat_is_p)
                    }
                    acov0 = acov0 + bsec[is,2]*pmat_is%*%t(uinv)%*%qv0%*%uinv%*%pmat_is
                }
            }
            acov = uinv%*%acov0%*%t(uinv)*sighat^2
            #only xmatpr and splpr have interpolation points
            #xmatpr = cbind(newv, splpr)
            #muhatpr = xmatpr %*% object$coefs
            #temp:
            #muhatpr = xmatpr %*% object$coefs[1:ncol(xmatpr), ,drop=FALSE]
            #new: C.I. level
            mult = qnorm((1 - level)/2, lower.tail=FALSE)
            #hl = 2*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
            #hl = mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
            #ans = list(fit = muhatpr, lower = muhatpr - hl, upper = muhatpr + hl)
        }
        if (interval == "confidence") {
            hl = mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
        }
        if (interval == "prediction") {
            hl = mult*sqrt(sighat^2+diag(xmatpr%*%acov%*%t(xmatpr)))
        }
        ans = list(fit = muhatpr, lower = muhatpr - hl, upper = muhatpr + hl, object = object, acov = acov, mult = mult, newz = newz)
    }
    class(ans) = "trisplp"
    return (ans)
}


###
#trispl makedelta
###
makedelta_tri = function(x1, x2, m1 = 0, m2 = 0, k1 = NULL, k2 = NULL, trimat = NULL, capk = NULL, space = c("E",
"E"), cvs = c(TRUE, TRUE), interp = TRUE) {
    n = length(x1)
    xmat = cbind(x1, x2)
    delta = NULL
    #### now determine which triangle contains each point
    xtri = 1:n*0;still = 1:n>0
    ntri = (m2 - 1) * (2 * m1 - 1)
    knots = cbind(k1, k2)
    for (j in 1:ntri) {
        #print (j)
        for (i in 1:n) {
            if (still[i]) {
                if (intri(xmat[i,],knots[trimat[j,1],],knots[trimat[j,2],],knots[trimat[j,3],])) {
                    still[i] = FALSE
                    xtri[i] = j
                }
            }
        }
    }
    ### make the design matrix
    dmat = matrix(0, nrow = n, ncol = capk)
    bmat = matrix(1, nrow = 3, ncol = 3)
    a = 1:3
    for (i in 1:n) {
        for (j in 1:3) {
            #print (c(i,j))
            a[j] = trimat[xtri[i],j]
            bmat[j,2] = knots[a[j],1]
            bmat[j,3] = knots[a[j],2]
        }
        binv = solve(bmat)
        for (j in 1:3) {
            dmat[i,a[j]] = binv[1,j] + binv[2,j]*x1[i] + binv[3,j]*x2[i]
        }
    }
    delta = cbind(delta, dmat)
    return (delta)
}

###########################################
#create a 3D plot for a cgam or wps object#
###########################################
plotpersp <- function(object,...) {
  UseMethod("plotpersp", object)
}
# plotpersp <- function(object, x1 = NULL, x2 = NULL,...) {
#     x1nm <- deparse(substitute(x1))
#     x2nm <- deparse(substitute(x2))
#     if (inherits(object, "cgamp")) {
#         xnms_add <- object$object$xnms_add
#     } else {
#         xnms_add <- object$xnms_add
#     }
#     if (inherits(object, "wpsp")) {
#         xnms_wp <- object$object$xnms_wp
#     } else {
#         xnms_wp <- object$xnms_wp
#     }
#     if (inherits(object, "trisplp")) {
#         xnms_tri <- object$object$xnms_tri
#     } else {
#         xnms_tri <- object$xnms_tri
#     }
#     is_add <- all(c(any(grepl(x1nm, xnms_add, fixed = TRUE)), any(grepl(x2nm, xnms_add, fixed = TRUE))))
#     #print (is_add)
#     is_wps <- all(c(any(grepl(x1nm, xnms_wp, fixed = TRUE)), any(grepl(x2nm, xnms_wp, fixed = TRUE))))
#     #print (is_wps)
#     is_tri <- all(c(any(grepl(x1nm, xnms_tri, fixed = TRUE)), any(grepl(x2nm, xnms_tri, fixed = TRUE))))
#     #print (is_tri)
#     if (missing(x1) | missing(x2)) {
#         UseMethod("plotpersp")
#     } else {
#         cs = class(object)
#         if (length(cs) == 1 & is.null(x1nm) & is.null(x2nm)) {
#             UseMethod("plotpersp")
#         } else {
#             if (is_wps) {
#                 #print (x1)
#                 #print (x2nm)
#                 if (inherits(object, "wpsp")) {
#                     #print (x1nm)
#                     #print (x2nm)
#                     #print (head(x1))
#                     plotpersp.wpsp(object, x1, x2, x1nm, x2nm,...)
#                 } else {
#                     plotpersp.wps(object, x1, x2, x1nm, x2nm,...)
#                 }
#             } else if (is_add) {
#                 if (inherits(object, "cgamp")){
#                     plotpersp.cgamp(object, x1, x2, x1nm, x2nm,...)
#                 } else {
#                     plotpersp.cgam(object, x1, x2, x1nm, x2nm,...)
#                 }
#             } else if (is_tri) {
#                 if (inherits(object, "trisplp")) {
#                     plotpersp.trisplp(object, x1, x2, x1nm, x2nm,...)
#                 } else {
#                     plotpersp.trispl(object, x1, x2, x1nm, x2nm,...)
#                 }
#             } else {
#                 stop ("Nonparametric components must be from the same class!")
#             }
#         }
#     }
# }


################
#plotpersp.cgam#
################
plotpersp.cgam <- function(object, x1 = NULL, x2 = NULL, x1nm = NULL, x2nm = NULL,
                           data = NULL, surface = "mu", categ = NULL,
                           col = NULL, random = FALSE, ngrid = 12,
                           xlim = range(x1), ylim = range(x2), zlim = NULL,
                           xlab = NULL, ylab = NULL, zlab = NULL, th = NULL,
                           ltheta = NULL, main = NULL, ticktype = "simple",...) {
	#if (class(object) == "list") {
	#	object <- object[[1]]
	#}
    #print (class(object))
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
	#x1nm <- deparse(substitute(x1))
	#x2nm <- deparse(substitute(x2))
	#print (x1nm)
	#print (x2nm)
#stop (print (x1nm))
	xnms <- object$xnms_add
	xmat <- object$xmat_add
	#if (x1nm == "NULL" | x2nm == "NULL") {
	if (is.null(x1nm) | is.null(x2nm)) {
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
	if (any(class(object) == "wps")) {
		d0 <- object$d0
		np_add <- object$np_add
		p <- d0 - np_add
		pb <- object$pb
		zmat <- zmat[, (pb+1):(pb+p), drop = FALSE]
#remove the one
		zmat <- zmat[, -1, drop = FALSE]
		#print (head(zmat))
	}
	#zcoefs = object$zcoefs
#new: exclude the coef for the one vector
#temp:trispl not include one
if (fml != "ordered" & all(class(object) != "trispl")) {
	zcoefs <- object$zcoefs[-1]
} else {
	zcoefs <- object$zcoefs
}
#print (zcoefs)
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
	x_grid <- ngrid
	y_grid <- ngrid
	x1g <- 0:x_grid / x_grid * .95 * (max(xm[,1]) - min(xm[,1])) + min(xm[,1]) + .025 * (max(xm[,1]) - min(xm[,1]))
 	n1 <- length(x1g)
	x2g <- 0:y_grid / y_grid * .95 * (max(xm[,2]) - min(xm[,2])) + min(xm[,2]) + .025 * (max(xm[,2]) - min(xm[,2]))
	n2 <- length(x2g)
	xgmat <- matrix(nrow = n1, ncol = n2)
	eta0 <- object$coefs[1]
	thvecs <- object$etacomps
	#print ('thvecs')
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
		#xgmat <- muhat.fun(xgmat, fml = fml)
		if (fml != "gaussian" & fml != "ordered") {
			for (i2 in 1:n2) {
				for (i1 in 1:n1) {
					xgmat[i1, i2] <- muhat.fun(xgmat[i1, i2], fml = fml)
				}
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
		#zi <- zmat[, pos1:pos2, drop = FALSE]
		#z_add <- 1:nrow(zi)*0
		#zcoefsi <- zcoefs[pos1:pos2]
#print (zcoefsi)
		zcoefsi = zcoefs[pos1:pos2]
#include the base level
		zcoefsi = c(0, zcoefsi)
		z_add = sort(zcoefsi)
		kz_add <- length(z_add)
		#for (j in 1:ncol(zi)) {
		#	zij <- zi[,j]
		#	zijhat <- zij * zcoefsi[j]
		#	z_add <- z_add + zijhat
		#}
		#z_add <- unique(z_add)
		#kz_add <- length(z_add)
#new: plot the smallest one first:
		#z_add <- z_add[order(z_add)]
#print (z_add)
		for (iz in 1:kz_add) {
			xgmats[[iz]] <- xgmat + as.numeric(x3_add) + z_add[iz]
			#mins <- c(mins, min(xgmats[[iz]]))
			#maxs <- c(maxs, max(xgmats[[iz]]))
			if (fml != "gaussian" & fml != "ordered") {
				for (i2 in 1:n2) {
					for (i1 in 1:n1) {
						xgmats[[iz]][i1, i2] <- muhat.fun(xgmats[[iz]][i1, i2], fml = fml)
					}
				}
			}
			#xgmat[[iz]] <- muhat.fun(xgmat[[iz]], fml = fml)
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
			} else if (fml == "poisson" | fml == "gaussian" | fml == "Gamma") {
				zlab <- paste("Est mean of", ynm)
			}
		}
		if (surface == "eta") {
			if (fml == "binomial") {
				zlab <- paste("Est log odds ratio of", ynm)
			}  else if (fml == "poisson" | fml == "Gamma") {
				zlab <- paste("Est log mean of", ynm)
			} else if (fml == "gaussian") {
				zlab <- paste("Est mean of", ynm)
			}
		}
	}
    #if (is.null(zlim)) {
    #	zlim <- range(xgmat, na.rm = TRUE)
    #}
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
			col0 <- col
			if (col0 == "heat" | col0 == "topo" | col0 == "terrain" | col0 == "cm") {
				#col0 <- col
				ncol <- 100
				facets <- facetcols <- list()
				col <- list()
				for (i in 1:kxgm) {
					nr <- nrow(xgmats[[i]])
					nc <- ncol(xgmats[[i]])
					facets[[i]] <- (xgmats[[i]])[-1,-1] + (xgmats[[i]])[-1,-nc] + (xgmats[[i]])[-nr,-1] + (xgmats[[i]])[-nr,-nc]
					facetcols[[i]] <- cut(facets[[i]], ncol)
					#print (head(facetcols[[i]]))
					if (col0 == "heat") {
						col[[i]] <- (heat.colors(ncol))[facetcols[[i]]]
						#print (head(col[[i]]))
					} else if (col0 == "topo") {
						col[[i]] <- (topo.colors(ncol))[facetcols[[i]]]
					} else if (col0 == "terrain") {
						col[[i]] <- (terrain.colors(ncol))[facetcols[[i]]]
					} else {
						col[[i]] <- (cm.colors(ncol))[facetcols[[i]]]
					}
				}
			} else if (length(col0) < kxgm) {
				rem <- kxgm - length(col0)
				nrem <- length(rem)
				rem_col <- palette[1:nrem]
				col <- c(col0, rem_col)
#new:
				#nr <- nrow(xgmat)
				#nc <- ncol(xgmat)
				#ncol <- 100
				#facet <- xgmat[-1,-1] + xgmat[-1,-nc] + xgmat[-nr,-1] + xgmat[-nr,-nc]
				#facetcol <- cut[facet, ncol]
				#col <- topo.colors[facetcol]

			} else if (length(col0) > kxgm) {
				col <- col0[1:kxgm]
				#print (paste("The first", kxgm, "colors are used!"))
			}
			if (random) {
				print ("User defined colors are used!")
			}
		}
		#print (col[[1]][1:10])
#print (kxgm)
#new: set th for decr or incr
		decrs = shapes[c(x1id, x2id)] %in% c(2, 6, 8, 10, 15, 16)
		incrs = shapes[c(x1id, x2id)] %in% c(1, 5, 7, 9, 13, 14)
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
		#print (col[[1]][1:10])
		for (i in 1:kxgm) {
			#print (col[i])
#i = 1
#print (i)
#print (length(col))
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
				xlab = ylab = zlab = " "
				box = FALSE
				axes = FALSE
			}
#print (length(col))
			if (is.list(col)) {
				coli = unlist(col[[i]])
				#print (head(coli))
			} else {coli = col[i]}
			#print (head(coli))
			#print (head(xgmat[[i]]))
			#print (head(col[[1]]))
            if (is.null(zlim)) {
                lwr = min(mins)
                upp = max(maxs)
                zlim0 = c(lwr - (upp-lwr)/3, upp + (upp-lwr)/3)
            } else {
                zlim0 = zlim
            }
			persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab,
			      col = coli, main = main, theta = ang, ltheta = ltheta, ticktype = ticktype, box = box, axes = axes,...)
			par(new = TRUE)
		}
	par(new = FALSE)
	} else {
		if (is.null(col)) {
			if (random) {
				col <- sample(palette, size = 1, replace = FALSE)
			} else {
				#col <- "white"
				#col <- color[facetcol]
				nr <- nrow(xgmat)
				nc <- ncol(xgmat)
				ncol <- 100
				facet <- xgmat[-1,-1] + xgmat[-1,-nc] + xgmat[-nr,-1] + xgmat[-nr,-nc]
				#print (facet)
				facetcol <- cut(facet, ncol)
				col <- heat.colors(ncol)[facetcol]
			}
		} else {
			#if (length(col) > 1) {
			#	col <- col[1]
			#	print ("The first color is used!")
			#	col <- heat.colors(x_grid*y_grid)
			#}
			if (col == "heat" | col == "topo" | col == "terrain" | col == "cm") {
				nr <- nrow(xgmat)
				nc <- ncol(xgmat)
				ncol <- 100
				facet <- xgmat[-1,-1] + xgmat[-1,-nc] + xgmat[-nr,-1] + xgmat[-nr,-nc]
				facetcol <- cut(facet, ncol)
				if (col == "heat") {
					col <- heat.colors(ncol)[facetcol]
				} else if (col == "topo") {
					col <- topo.colors(ncol)[facetcol]
				} else if (col == "terrain") {
					col <- terrain.colors(ncol)[facetcol]
				} else {
					col <- cm.colors(ncol)[facetcol]
				}
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
		incrs = shapes[c(x1id, x2id)] %in% c(1, 5, 7, 9, 13, 14)
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
        if (is.null(zlim)) {
            lwr = min(xgmat)
            upp = max(xgmat)
            zlim0 = c(lwr - (upp-lwr)/3, upp + (upp-lwr)/3)
        } else {
            zlim0 = zlim
        }
		persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab, col = col, main = main, theta = ang, ltheta = ltheta, ticktype = ticktype,...)
        rslt = list(zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab, theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype, z_add = z_add, x3_add = x3_add)
        invisible(rslt)
	}
#print (col)
}


#############################################################
#apply plotpersp on a predict.cgam or predict.cgamm object
#############################################################
plotpersp.cgamp = function(object, x1=NULL, x2=NULL, x1nm=NULL, x2nm=NULL,
                           data=NULL, up = TRUE, main=NULL, cex.main=.8,
                           xlab = NULL, ylab = NULL, zlab = NULL, zlim = NULL,
                           th = NULL, ltheta = NULL, ticktype = "detailed",...) {
    #obj is prediction for cgam or cgamm
    if (!inherits(object, "cgamp")) {
        warning("calling plotpersp(<fake-cgam-prediction-object>) ...")
    }
    t_col = function(color, percent = 50, name = NULL) {
        rgb.val <- col2rgb(color)
        ## Make new color using input color as base and alpha set by transparency
        t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
        maxColorValue = 255,
        alpha = (100-percent)*255/100,
        names = name)
        ## Save the color
        invisible(t.col)
    }
    if (up) {
        mycol = t_col("green", perc = 90, name = "lt.green")
    } else {
        mycol = t_col("pink", perc = 80, name = "lt.pink")
    }

    acov = object$acov
    mult = object$mult
    #obj is the cgam or cgamm fit
    obj = object$object
    #ah = obj$coef_wp

    xnms = obj$xnms_add
    xmat = obj$xmat_add
    bigmat = obj$bigmat
    #decrs = obj$decrs

    #if (x1nm == "NULL" | x2nm == "NULL") {
    if (is.null(x1nm) | is.null(x2nm)) {
        if (length(xnms) >= 2) {
            x1nm = xnms[1]
            x2nm = xnms[2]
            x1id = 1
            x2id = 2
            x1 = xmat[, 1]
            x2 = xmat[, 2]
        } else {stop ("Number of non-parametric predictors must >= 2!")}
    }
    #labels = obj$labels
    #labels = labels[which(grepl("warp", labels, fixed = TRUE))]

    #is_fac = obj$is_fac
    ynm = obj$ynm
    #varlist = obj$varlist_wps
    #varlist = varlist[-1]
    kts = obj$knots
    np = obj$d0

    knms = length(xnms)
    obs = 1:knms

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
        x1id <- obs[xnms == x1nm]

        if (length(x2) != nrow(xmat)) {
            warning ("Number of observations in the data set is not the same as the number of elements in x2!")
        }
        x2id <- obs[xnms == x2nm]
    } else {
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
    #xm = xmat[, c(x1id, x2id)]
    xm0 = object$newx0
    xm = xm0[,c(x1id,x2id)]
    #print (x1id)
    #print (x2id)
    #print (head(x1))
    #print (head(x2))
    thvs_upp = object$thvs_upp
    thvs_lwr = object$thvs_lwr
    mins = min(thvs_lwr)+obj$coefs[1]
    maxs = max(thvs_upp)+obj$coefs[1]
    #print (mins)
    #print (maxs)
    #print (obj$coefs[1])
    #if (up) {
    if (is.null(zlim)) {
        zlim = c(mins-(maxs-mins)/2.2, maxs+(maxs-mins)/2.2)
    }
    #} else {
    #    zlim = c(mins-(maxs-mins)/3, maxs)
    #}
    res = plotpersp.cgam(obj, x1=x1, x2=x2, x1nm=x1nm, x2nm=x2nm, zlim=zlim, col='white', xlab=xlab, ylab=ylab, zlab=zlab, th=th, ltheta=ltheta, ticktype=ticktype)
    #print ('check')
    ngrid = res$ngrid
    #print (res$zlim)

    x_grid = ngrid
    y_grid = ngrid
    x1g = 0:x_grid / x_grid * .95 * (max(xm[,1]) - min(xm[,1])) + min(xm[,1]) + .025 * (max(xm[,1]) - min(xm[,1]))
    n1 = length(x1g)
    x2g = 0:y_grid / y_grid * .95 * (max(xm[,2]) - min(xm[,2])) + min(xm[,2]) + .025 * (max(xm[,2]) - min(xm[,2]))
    n2 = length(x2g)
    xgmat = matrix(nrow = n1, ncol = n2)
    eta0 = obj$coefs[1]
    if (up) {
        thvecs = thvs_upp
    } else {
        thvecs = thvs_lwr
    }
    for (i2 in 1:n2) {
        for (i1 in 1:n1) {
            x1a = max(xm[xm[,1] <= x1g[i1], 1])
            x1b = min(xm[xm[,1] > x1g[i1], 1])
            v1a = min(thvecs[x1id, xm[,1] == x1a])
            v1b = min(thvecs[x1id, xm[,1] == x1b])
            alp = (x1g[i1] - x1a) / (x1b - x1a)
            th1add = (1 - alp) * v1a + alp * v1b
            x2a = max(xm[xm[,2] <= x2g[i2],2])
            x2b =min(xm[xm[,2] > x2g[i2],2])
            v2a = min(thvecs[x2id, xm[,2] == x2a])
            v2b = min(thvecs[x2id, xm[,2] == x2b])
            alp = (x2g[i2] - x2a) / (x2b - x2a)
            th2add = (1 - alp) * v2a + alp * v2b
            xgmat[i1,i2] = eta0 + th1add + th2add
        }
    }

    z_add = res$z_add
    x3_add = res$x3_add
    xgmat = xgmat + z_add + x3_add
    fml = obj$family$family
    #if (fml != "gaussian" & fml != "ordered") {
    #    for (i2 in 1:n2) {
    #        for (i1 in 1:n1) {
    #            xgmat[i1, i2] <- muhat.fun(xgmat[i1, i2], fml = fml)
    #        }
    #    }
    #}
    if (up) {
        if (is.null(main)) {
            main = "Cgam Surface with Upper 95% Confidence Surface"
        }
    }
    if (!up) {
        if (is.null(main)) {
            main = "Cgam Surface with Lower 95% Confidence Surface"
        }
    }
    par(new = TRUE)
    persp(x1g, x2g, xgmat, zlim = res$zlim, xlab = "", ylab = "", zlab = "", theta = res$theta,
          ltheta = res$ltheta, cex.axis = res$cex.axis, main = main, cex.main = cex.main,
          ticktype = res$ticktype, col=mycol, box=FALSE, axes=FALSE,...)
    par(new=FALSE)
}



#################
#plotpersp.wps#
################
plotpersp.wps = function(object, x1 = NULL, x2 = NULL, x1nm = NULL,
                         x2nm = NULL, data = NULL, surface = "C",
                         categ = NULL, col = NULL, random = FALSE,
                         xlim = range(x1), ylim = range(x2), zlim = NULL,
                         xlab = NULL, ylab = NULL, zlab = NULL, th = NULL,
                         ltheta = NULL, main = NULL, ticktype = "simple",...) {
  #print ('wps')
  if (!inherits(object, "wps")) {
    warning("calling plotpersp(<fake-wps-object>) ...")
  }
  #x1nm = deparse(substitute(x1))
  #x2nm = deparse(substitute(x2))
  #print (x1nm)
  #print (x2nm)
  xnms = object$xnms_wp
  xmat = object$xmat_wp
  #if (x1nm == "NULL" | x2nm == "NULL") {
  if (is.null(x1nm) | is.null(x2nm)) {
    if (length(xnms) >= 2) {
      x1nm = xnms[1]
      x2nm = xnms[2]
      x1id = 1
      x2id = 2
      x1 = xmat[, 1]
      x2 = xmat[, 2]
    } else {stop ("Number of non-parametric predictors must >= 2!")}
  }
  #xnms = object$xnms
  #xmat = object$xmat
  labels = object$labels
  labels = labels[which(grepl("warp", labels, fixed = TRUE))]
  #new:
  is_fac = object$is_fac
  ynm = object$ynm
  #xmat is delta
  #delta = object$delta
  znms = object$znms
  decrs0 = object$decrs
  kznms = length(znms)
  #zmat include 1 vector if only wps + add
  zmat = object$zmat
  np = d0 = object$d0
  pb = object$pb
  np_add = object$np_add
  p = d0 - np_add
  #print ('call wps')
  #zmat = zmat[, (pb+1):(pb+d0), drop = FALSE]
  #check!
  #zmat = zmat[, (pb+1):(pb+p), drop = FALSE]
  #print ('zmat')
  #print (head(zmat))
  if (any(class(object) == "trispl")) {
    zmat0 = zmat
  } else {
    #test more...
    if (!is.null(zmat)) {
      zmat = zmat[, (pb+1):(pb+p), drop = FALSE]
      #new constant in zmat now
      #if (d0 > 1) {
      #  zmat0 = zmat[, -1, drop = FALSE]
      #} else {zmat0 = NULL}
      zmat0 = zmat
    } else {
      zmat0 = zmat
    }
  }
  #print (head(zmat0))
  if (all(class(object) != "trispl")) {
    #zcoefs = object$zcoefs[-1]
    #no more constant vector in zmat
    zcoefs = object$zcoefs
  } else {zcoefs = object$zcoefs}
  #print (zcoefs)
  zid1 = object$zid1
  zid2 = object$zid2
  ah = object$coef_wp
  #ahu = object$coefsu
  varlist = object$varlist_wps
  #varlist = varlist[-1]
  kts = object$ks_wps
  #k1 = object$k1
  #k2 = object$k2
  #p = dim(zmat)[2]
  #p = d0
  family = object$family
  fml = family$family
  cicfamily = CicFamily(family)
  muhat.fun = cicfamily$muhat.fun
  #additive
  thvecs = object$etacomps
  xnms_add = object$xnms_add
  xmat_add = object$xmat_add
  knms = length(xnms_add)
  x3_add = 0
  if (knms >= 1) {
    #x3id <- obs[-c(x1id, x2id)]
    #kx3 <- length(x3id)
    for (i in 1:knms) {
      x3i = xmat_add[, i]
      x3i_use = max(x3i[x3i <= median(x3i)])
      x3i_add = min(thvecs[i, x3i == x3i_use])
      x3_add = x3_add + x3i_add
    }
  }
  #if (!is.null(categ)) {
  #	if (!is.character(categ)) {
  #		warning("categ must be a character argument!")
  #	} else if (!any(znms == categ)) {
  #print ('TRUE')
  #		warning(paste(categ, "is not an exact character name defined in the cgam fit!"))
  #		categ = NULL
  #	} else {
  #		obsz = 1:kznms
  #		zid = obsz[znms == categ]
  #		if (!(is_fac[zid])) {
  #			categ = NULL
  #		}
  #	}
  #}
  if (!is.null(categ)) {
    if (!is.character(categ)) {
      warning("categ must be a character argument!")
    } else if (!any(znms == categ)) {
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
      #if (!(is_fac[zid])) {
      #  categ = NULL
      #}
    }
  }
  #new: switch xnms
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop ("User need to make the data argument a data frame with names for each variable!")
    }
    datnms = names(data)
    if (!any(datnms == x1nm) | !any(datnms == x2nm)) {
      stop ("Check the accuracy of the names of x1 and x2!")
    }
    x1 = data[ ,which(datnms == x1nm)]
    x2 = data[ ,which(datnms == x2nm)]
  } else {
    if (all(xnms != x1nm)) {
      #stop (paste(paste("'", x1nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
      #new: in case of wrong data fame
      if (length(x1) != nrow(xmat)) {
        stop ("Number of observations in the data set is not the same as the number of elements in x1!")
      }
      bool = apply(xmat, 2, function(x) all(x1 == x))
      if (any(bool)) {
        #x1id = obs[bool]
        x1nm = xnms[bool]
      } else {
        stop (paste(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the wps fit!"))
      }
    }
    if (all(xnms != x2nm)) {
      #stop (paste(paste("'", x2nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
      if (length(x2) != nrow(xmat)) {
        stop ("Number of observations in the data set is not the same as the number of elements in x2!")
      }
      bool = apply(xmat, 2, function(x) all(x2 == x))
      if (any(bool)) {
        #x2id = obs[bool]
        x2nm = xnms[bool]
      } else {
        stop (paste(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the wps fit!"))
      }
    }
  }
  xnm12 = c(x1nm, x2nm)
  id_lab = which(xnms %in% xnm12)
  xnm12_lab = labels[id_lab]
  xnm_other = xnms[-id_lab]
  id1 = id2 = ipr = NULL
  #print ('id_lab')
  #print (xnm12_lab)
  if (length(unique(xnm12_lab)) > 1 | length(id_lab) != 2) {
    stop ("Two non-parametric predictors do not form a warped-plane surface!")
  } else {
    id1 = sort(id_lab)[1]
    id2 = sort(id_lab)[2]
    ipr = id2 / 2
  }
  decrs = decrs0[[ipr]]
  ktsi = kts[[ipr]]
  k1 = ktsi[[1]]
  k2 = ktsi[[2]]
  m1 = length(k1)
  m2 = length(k2)
  #if (x1nm0 != xnms[1] & x2nm0 != xnms[2]) {
  #	x1nm = x2nm0
  #	x2nm = x1nm0
  #	tmp = x1
  #	x1 = x2
  #	x2 = tmp
  #} else {x1nm = x1nm0; x2nm = x2nm0}
  #print (paste('id1: ', id1))
  #print (paste('id2: ', id2))
  if (x1nm != xnms[id1] & x2nm != xnms[id2]) {
    nm = x1nm
    x1nm = x2nm
    x2nm = nm
    tmp = x1
    x1 = x2
    x2 = tmp
  }
  # apl = 1:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
  # #aplu = 1:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
  # apl[1] = ah[1]
  # #aplu[1] = ahu[1]
  # #apl[2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))] = ah[(p + 1):(m1 + m2 - 1 + (m1 - 1) * (m2 - 1) + p - 1)]
  # #aplu[2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))] = ahu[(p + 1):(m1 + m2 - 1 + (m1 - 1) * (m2 - 1) + p - 1)]
  # apl[2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))] = (ah[-c(1:p)])[which(varlist == ipr)]
  # #aplu[2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))] = (ahu[-1])[which(varlist == ipr)]
  # mupl = matrix(apl[1], nrow = m1, ncol = m2)
  # #muplu = matrix(aplu[1], nrow = m1, ncol = m2)
  # for (i1 in 2:m1) {
  #     mupl[i1, ] = mupl[i1, ] + apl[i1]
  #     #muplu[i1, ] = muplu[i1, ] + aplu[i1]
  # }
  # for (i2 in 2:m2) {
  #     mupl[, i2] = mupl[, i2] + apl[m1 - 1 + i2]
  #     #muplu[, i2] = muplu[, i2] + aplu[m1 - 1 + i2]
  # }
  # for (i1 in 2:m1) {
  #     for (i2 in 2:m2) {
  #         mupl[i1, i2] = mupl[i1, i2] + apl[m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2]
  #         #muplu[i1, i2] = muplu[i1, i2] + aplu[m1 + m2 - 2 + (i1 - 2) * (m2 - 1) + i2]
  #     }
  # }
  #plot the estimates at knots
  newData = expand.grid(k1,k2)
  colnames(newData) = xnm12
  npr = round(length(xnms) / 2, 0L)
  new_other = NULL
  if (npr > 1) {
    new_other = matrix(0, nrow=nrow(newData), ncol=(2*(npr-1)))
    kts_other = kts[-ipr]
    for (i in 1:(npr-1)) {
      for (j in 1:2) {
        newi = mean(kts_other[[i]][[j]])
        new_other[,(i-1)*2+j] = newi
      }
    }
    newd = cbind(newData, new_other)
    colnames(newd) = c(xnm12, xnm_other)
    newData = as.data.frame(newd)
  }

  if (length(xnms_add) > 0) {
    #new_add = matrix(0, nrow=nrow(newData), ncol=ncol(xmat_add))
    means = apply(xmat_add, 2, mean)
    new_add = matrix(rep(means, nrow(newData)), ncol=ncol(xmat_add), byrow=T)
    nms = colnames(newData)
    newd = cbind(newData, new_add)
    colnames(newd) = c(nms, xnms_add)
    newData = as.data.frame(newd)
  }
  #test!
  Mode = function(x) {
    ux = unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  newz = NULL
  if (!is.null(zmat)) {
    newz = matrix(0, nrow=nrow(newData), ncol=ncol(zmat))
    for(j in 1:(ncol(zmat))) {
      newz[,j] = Mode(zmat[,j])
    }
  }
  if (!is.null(newz)) {
    znms = object$znms
    nz = matrix(0, nrow=nrow(newData), ncol=ncol(newz))
    nznm = gsub("[\\(\\)]", "", regmatches(znms, gregexpr("\\(.*?\\)", znms))[[1]])
    #new: zmat has > 1 columns -> not work, will make predict.wps give an error
    #nznm = paste0(nznm, 1:ncol(newz), sep="")
    nms = colnames(newData)
    newData = cbind(newData, nz)
    #check more:
    if (length(nznm) > 0) {
      colnames(newData) = c(nms, rep(nznm, ncol(newz)))
    }
    if (length(nznm) == 0) {
      colnames(newData) = c(nms, znms)
    }
  }

  pfit.knots = predict.wps(object, newData, interval='none')

  xmatpr = pfit.knots$xmatpr
  spls = pfit.knots$spls
  #ignore z
  spl_use = spls[[ipr]]

  #get the fit for each pair
  mus = pfit.knots$mus
  muhat_use = spl_use%*%ah[which(varlist == ipr)+np]

  mupl = matrix(0, m1, m2)
  for(i2 in 1:m2) {
    for(i1 in 1:m1) {
      mupl[i1,i2] = muhat_use[(i2-1)*m1 + i1]
    }
  }
  #new: more families
  if (fml != "gaussian") {
    for (i1 in 1:m1) {
      for (i2 in 1:m2) {
        mupl[i1, i2] = muhat.fun(mupl[i1, i2], fml = fml)
        #muplu[i1, i2] = muhat.fun(muplu[i1, i2], fml = fml)
      }
    }
    #mupl = muhat.fun(mupl, fml = fml)
    #muplu = muhat.fun(muplu, fml = fml)
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
      #mupl = mupl[m1:1, 1:m2];
      #muplu = muplu[m1:1, 1:m2]
      if (is.null(ang)) {
        ang = 40
      }
    } else if (!decrs[1] & decrs[2]) {
      #k2 = -k2[m2:1];
      #mupl = mupl[1:m1, m2:1];
      #muplu = muplu[1:m1, m2:1]
      if (is.null(ang)) {
        ang = 240
      }
    } else if (decrs[1] & decrs[2]) {
      #k1 = -k1[m1:1]; k2 = -k2[m2:1];
      #mupl = mupl[m1:1, m2:1];
      #muplu = muplu[m1:1, m2:1]
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
    mupl = mupl + as.numeric(z_add) + as.numeric(x3_add)
    #muplu = muplu + as.numeric(z_add)
    #new:
    if (fml != "gaussian") {
      for (i1 in 1:m1) {
        for (i2 in 1:m2) {
          mupl[i1, i2] = muhat.fun(mupl[i1, i2], fml = fml)
          #muplu[i1, i2] = muhat.fun(muplu[i1, i2], fml = fml)
        }
      }
      #mupl = muhat.fun(mupl, fml = fml)
      #muplu = muhat.fun(muplu, fml = fml)
    }
    mins = min(mupl); maxs = max(mupl)
    #minsu = min(muplu); maxsu = max(muplu)
  } else {
    mupls = muplus = list()
    mins = maxs = NULL
    minsu = maxsu = NULL
    obsz = 1:kznms
    zid = obsz[znms == categ]
    pos1 = zid1[zid]; pos2 = zid2[zid]
    zcoefsi = zcoefs[pos1:pos2]
    #?include the base level
    zcoefsi = c(0, zcoefsi)
    z_add = sort(zcoefsi)
    kz_add = length(z_add)
    #print (kz_add)
    for (iz in 1:kz_add) {
      mupls[[iz]] = mupl + z_add[iz] + as.numeric(x3_add)
      mins = c(mins, min(mupls[[iz]]))
      maxs = c(maxs, max(mupls[[iz]]))
      #muplus[[iz]] = muplu + z_add[iz]
      #minsu = c(minsu, min(muplus[[iz]]))
      #maxsu = c(maxsu, max(muplus[[iz]]))
    }
    if (fml != "gaussian") {
      for (iz in 1:kz_add) {
        mupli = mupls[[iz]]
        #muplui = muplus[[iz]]
        for (i1 in 1:m1) {
          for (i2 in 1:m2) {
            mupli[i1, i2] = muhat.fun(mupli[i1, i2], fml = fml)
            #muplui[i1, i2] = muhat.fun(muplui[i1, i2], fml = fml)
          }
        }
        #mupli = muhat.fun(mupli, fml = fml)
        #muplui = muhat.fun(muplui, fml = fml)
        mupls[[iz]] = mupli
        #muplus[[iz]] = muplui
      }
    }
  }
  palette = c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
  if (is.null(xlab)) {
    xlab = x1nm
  }
  if (is.null(ylab)) {
    ylab = x2nm
  }
  if (is.null(zlab)) {
    if (fml == "binomial") {
      zlab = paste("Pr(", ynm, ")")
    } else if (fml == "poisson" | fml == "gaussian" | fml == "Gamma") {
      zlab = paste("Est mean of", ynm)
    }
  }
  if (is.null(categ)) {
    if (is.null(col)) {
      if (random) {
        col = sample(palette, size = 1, replace = FALSE)
      }  else {
        #col = "white"
        if (surface == 'C') {
          musurf = mupl
          #} else if (surface == 'U') {
          #	musurf = muplu
        }
        nr = nrow(musurf)
        nc = ncol(musurf)
        ncol = 100
        facet = musurf[-1,-1] + musurf[-1,-nc] + musurf[-nr,-1] + musurf[-nr,-nc]
        facetcol = cut(facet, ncol)
        col = heat.colors(ncol)[facetcol]
      }
    } else if (col == "heat" | col == "topo" | col == "terrain" | col == "cm") {
      if (surface == 'C') {
        musurf = mupl
        #} else if (surface == 'U') {
        #	musurf = muplu
      }
      nr = nrow(musurf)
      nc = ncol(musurf)
      ncol = 100
      facet = musurf[-1,-1] + musurf[-1,-nc] + musurf[-nr,-1] + musurf[-nr,-nc]
      facetcol = cut(facet, ncol)
      if (col == "heat") {
        col = heat.colors(ncol)[facetcol]
      } else if (col == "topo") {
        col = topo.colors(ncol)[facetcol]
      } else if (col == "terrain") {
        col = terrain.colors(ncol)[facetcol]
      } else {
        col = cm.colors(ncol)[facetcol]
      }
    }

    if (surface == 'C') {
      musurf = mupl
      #if (is.null(main)) {
      #	main = 'Constrained Warped-Plane Spline Surface'
      #}
      if (is.null(zlim)) {
        lwr = min(mins)
        upp = max(maxs)
        #print (lwr)
        #print (upp)
        zlim0 = c(lwr - (upp-lwr)/5, upp + (upp-lwr)/5)
      } else {
        zlim0 = zlim
      }
      #} else if (surface == 'U') {
      #	musurf = muplu
      #if (is.null(main)) {
      #	main =  'Unconstrained Warped-Plane Spline Surface'
      #}
      #zlim0 = c(min(minsu), max(maxsu))
    }
    #print (head(musurf))
    xlim = range(x1)
    ylim = range(x2)
    #persp(k1, k2, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype,...)
    persp(k1, k2, musurf, zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab,
          theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype,...)
    rslt = list(zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab, theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype, z_add = z_add, x3_add = x3_add)
    invisible(rslt)
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
      col0 = col
      if (col0 == "heat" | col0 == "topo" | col0 == "terrain" | col0 == "cm") {
        ncol = 100
        facets = facetcols = list()
        if (surface == "C") {
          xgmats = mupls
        }
        #if (surface == "U") {
        #	xgmats = muplus
        #}
        col = list()
        for (i in 1:kxgm) {
          nr = nrow(xgmats[[i]])
          nc = ncol(xgmats[[i]])
          facets[[i]] = (xgmats[[i]])[-1,-1] + (xgmats[[i]])[-1,-nc] + (xgmats[[i]])[-nr,-1] + (xgmats[[i]])[-nr,-nc]
          facetcols[[i]] = cut(facets[[i]], ncol)
          if (col0 == "heat") {
            col[[i]] = (heat.colors(ncol))[facetcols[[i]]]
          } else if (col0 == "topo") {
            col[[i]] = (topo.colors(ncol))[facetcols[[i]]]
          } else if (col0 == "terrain") {
            col[[i]] = (terrain.colors(ncol))[facetcols[[i]]]
          } else {
            col[[i]] = (cm.colors(ncol))[facetcols[[i]]]
          }
        }
      } else if (length(col0) < kxgm) {
        #rem = kxgm - length(col)
        #nrem = length(rem)
        #rem_col = palette[1:nrem]
        #col = c(col, rem_col)
        #new:
        col = topo.colors(kxgm)
      } else if (length(col0) > kxgm) {
        col = col[1:kxgm]
        #print (paste("The first", kxgm, "colors are used!"))
      }
      #if (random) {
      #print ("User defined colors are used!")
      #}
    }
    for (i in 1:kxgm) {
      mupli = mupls[[i]]
      #muplui = muplus[[i]]
      if (surface == 'C') {
        musurf = mupli
        #if (is.null(main)) {
        #	main = 'Constrained Warped-Plane Spline Surface'
        #}
        if (is.null(zlim)) {
          lwr = min(mins)
          upp = max(maxs)
          zlim0 = c(lwr - (upp-lwr)/5, upp + (upp-lwr)/5)
        } else {
          zlim0 = zlim
        }
        #} else if (surface == 'U') {
        #	musurf = muplui
        #if (is.null(main)) {
        #	main = 'Unconstrained Warped-Plane Spline Surface'
        #}
        #zlim0 = c(min(minsu), max(maxsu))
      }
      #par(mar = c(4, 2, 2, 2))
      #print (sub)
      #persp(k1, k2, musurf,  sub = sub)
      if (is.list(col)) {
        coli = unlist(col[[i]])
      } else {coli = col[i]}
      #xlim = range(x1)
      #ylim = range(x2)
      #persp(k1, k2, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = coli, cex.axis = .75, main = main, ticktype = ticktype,...)
      persp(k1, k2, musurf, zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab, theta = ang,
            ltheta = ltheta, col = coli, cex.axis = .75, main = main, ticktype = ticktype,...)
      par(new = TRUE)
    }
    par(new = FALSE)
  }
}

##########################################
#apply plotpersp on a wps.predict object
#>=1 wps pair + z + additive #done
##########################################
###############################################################
#plotpersp for a predict.wps object
#check the zlim when there's a smooth additive component more
###############################################################
plotpersp.wpsp = function(object, x1=NULL, x2=NULL, x1nm=NULL, x2nm=NULL,
                          data=NULL, up = TRUE, main=NULL, cex.main=.8, xlab = NULL,
                          ylab = NULL, zlab = NULL, th = NULL, ltheta = NULL,
                          ticktype = "simple",...) {
  #obj is prediction for wps
  if (!inherits(object, "wpsp")) {
    warning("calling plotpersp(<fake-wpsp-object>) ...")
  }
  t_col = function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 maxColorValue = 255,
                 alpha = (100-percent)*255/100,
                 names = name)
    ## Save the color
    invisible(t.col)
  }
  if (up) {
    mycol = t_col("green", perc = 90, name = "lt.green")
  } else {
    mycol = t_col("pink", perc = 80, name = "lt.pink")
  }

  acov = object$acov
  mult = object$mult
  #obj is the wps fit
  obj = object$object
  ah = obj$coef_wp

  xnms = obj$xnms_wp
  xmat = obj$xmat_wp
  delta = obj$delta
  decrs = obj$decrs

  xnms_add = obj$xnms_add
  xmat_add = obj$xmat_add
  #if (x1nm == "NULL" | x2nm == "NULL") {
  if (is.null(x1nm) | is.null(x2nm)) {
    if (length(xnms) >= 2) {
      x1nm = xnms[1]
      x2nm = xnms[2]
      x1id = 1
      x2id = 2
      x1 = xmat[, 1]
      x2 = xmat[, 2]
    } else {stop ("Number of non-parametric predictors must >= 2!")}
  }

  labels = obj$labels
  labels = labels[which(grepl("warp", labels, fixed = TRUE))]

  is_fac = obj$is_fac
  ynm = obj$ynm

  varlist = obj$varlist_wps
  #varlist = varlist[-1]

  kts = obj$ks_wps
  np = obj$d0
  #switch xnms
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop ("User need to make the data argument a data frame with names for each variable!")
    }
    datnms = names(data)
    if (!any(datnms == x1nm) | !any(datnms == x2nm)) {
      stop ("Check the accuracy of the names of x1 and x2!")
    }
    x1 = data[ ,which(datnms == x1nm)]
    x2 = data[ ,which(datnms == x2nm)]
  } else {
    if (all(xnms != x1nm)) {
      #stop (paste(paste("'", x1nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
      #new: in case of wrong data fame
      if (length(x1) != nrow(xmat)) {
        stop ("Number of observations in the data set is not the same as the number of elements in x1!")
      }
      bool = apply(xmat, 2, function(x) all(x1 == x))
      if (any(bool)) {
        #x1id = obs[bool]
        x1nm = xnms[bool]
      } else {
        stop (paste(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the wps fit!"))
      }
    }
    if (all(xnms != x2nm)) {
      #stop (paste(paste("'", x2nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
      if (length(x2) != nrow(xmat)) {
        stop ("Number of observations in the data set is not the same as the number of elements in x2!")
      }
      bool = apply(xmat, 2, function(x) all(x2 == x))
      if (any(bool)) {
        #x2id = obs[bool]
        x2nm = xnms[bool]
      } else {
        stop (paste(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the wps fit!"))
      }
    }
  }
  xnm12 = c(x1nm, x2nm)
  id_lab = which(xnms %in% xnm12)
  xnm12_lab = labels[id_lab]
  #new:
  xnm_other = xnms[-id_lab]
  id1 = id2 = ipr = NULL

  if (length(unique(xnm12_lab)) > 1 | length(id_lab) != 2) {
    stop ("Two non-parametric predictors do not form a warped-plane surface!")
  } else {
    id1 = sort(id_lab)[1]
    id2 = sort(id_lab)[2]
    ipr = id2 / 2
  }

  #find the pairs to be plotted
  decrs_use = decrs[[ipr]]
  kts_use = kts[[ipr]]
  x1p = kts_use[[1]]
  x2p = kts_use[[2]]

  k1 = length(x1p)
  k2 = length(x2p)
  if (x1nm != xnms[id1] & x2nm != xnms[id2]) {
    nm = x1nm
    x1nm = x2nm
    x2nm = nm
    tmp = x1
    x1 = x2
    x2 = tmp
    #new:
    xnm12 = c(x1nm, x2nm)
  }
  newData = expand.grid(x1p,x2p)

  #find the names of the pair
  colnames(newData) = xnm12

  #ignore z if there's any
  #temp:
  newz = object$newz
  npr = round(length(xnms) / 2, 0L)

  new_other = NULL
  if (npr > 1) {
    new_other = matrix(0, nrow=nrow(newData), ncol=(2*(npr-1)))
    kts_other = kts[-ipr]
    for (i in 1:(npr-1)) {
      for (j in 1:2) {
        newi = mean(kts_other[[i]][[j]])
        new_other[,(i-1)*2+j] = newi
      }
    }
    newd = cbind(newData, new_other)
    colnames(newd) = c(xnm12, xnm_other)
    newData = as.data.frame(newd)
  }

  if (length(xnms_add) > 0) {
    #new_add = matrix(0, nrow=nrow(newData), ncol=ncol(xmat_add))
    means = apply(xmat_add, 2, mean)
    new_add = matrix(rep(means, nrow(newData)), ncol=ncol(xmat_add), byrow=T)
    nms = colnames(newData)
    newd = cbind(newData, new_add)
    colnames(newd) = c(nms, xnms_add)
    newData = as.data.frame(newd)
  }

  if (!is.null(newz)) {
    znms = obj$znms
    nz = matrix(0, nrow=nrow(newData), ncol=ncol(newz))
    nznm = gsub("[\\(\\)]", "", regmatches(znms, gregexpr("\\(.*?\\)", znms))[[1]])
    #new: zmat has > 1 columns -> not work, will make predict.wps give an error
    #nznm = paste0(nznm, 1:ncol(newz), sep="")
    nms = colnames(newData)
    newData = cbind(newData, nz)
    colnames(newData) = c(nms, nznm)
  }

  #print (head(x1))
  #print (x1nm)
  #print (x2nm)
  #print (colnames(newData))
  #determine zlim beforehand
  pfit.knots = predict.wps(obj, newData, interval='none')

  #get the spline for each pair
  xmatpr = pfit.knots$xmatpr
  spls = pfit.knots$spls
  #ignore z
  #spl_use = cbind(1, spls[[ipr]])
  spl_use = spls[[ipr]]

  #get the fit for each pair
  mus = pfit.knots$mus
  #muhat_use = mus[[ipr]]
  #test more:
  #muhat_use = spl_use%*%ah[c(1, which(varlist == ipr)+np)]
  muhat_use = spl_use%*%ah[which(varlist == ipr)+np]

  #get the acov for each pair
  if (length(xnms_add) > 0) {
    k_add = length(obj$coef_add)
    k_conv = sum(obj$shapes == 11 || obj$shapes == 12)
    #acov_use = acov[c(1, which(varlist == ipr)+np+k_conv+k_add), c(1, which(varlist == ipr)+np+k_conv+k_add)]
    keep_id = which(varlist == ipr)+np+k_conv+k_add
    acov_use = acov[keep_id, keep_id]
  } else {
    #acov_use = acov[c(1, which(varlist == ipr)+np), c(1, which(varlist == ipr)+np)]
    keep_id = which(varlist == ipr)+np
    acov_use = acov[keep_id, keep_id]
  }

  lower = muhat_use - mult*sqrt(diag(spl_use%*%acov_use%*%t(spl_use)))
  upper = muhat_use + mult*sqrt(diag(spl_use%*%acov_use%*%t(spl_use)))
  #lower = pfit.knots$lower
  #upper = pfit.knots$upper

  #zlwr = min(object$lower) - (max(object$fit) - min(object$fit)) / 2 + x3_add
  #zupp = max(object$upper) + (max(object$fit) - min(object$fit)) / 2 + x3_add
  #check the limits more
  zlwr = min(lower) - (max(object$fit) - min(object$fit)) / 10
  zupp = max(upper) + (max(object$fit) - min(object$fit)) / 10

  res = plotpersp.wps(obj, x1=x1, x2=x2, x1nm, x2nm, col='white', zlim=c(zlwr, zupp), xlab=xlab, ylab=ylab, zlab=zlab, th=th, ltheta=ltheta, ticktype=ticktype)
  #res = plotpersp.wps(obj, x1=x1, x2=x2, x1nm, x2nm, col='white', xlab=xlab, ylab=ylab, zlab=zlab, th=th, ltheta=ltheta, ticktype=ticktype)

  surf = matrix(0, k1, k2)
  for(i2 in 1:k2) {
    for(i1 in 1:k1) {
      if (up) {
        #surf[i1,i2] = muhat_use[(i2-1)*k1 + i1]
        surf[i1,i2] = upper[(i2-1)*k1 + i1]
      } else {
        surf[i1,i2] = lower[(i2-1)*k1 + i1]
      }
    }
  }
  z_add = res$z_add
  x3_add = res$x3_add
  surf = surf + z_add + x3_add
  if (up) {
    if (is.null(main)) {
      main = "Warped-Plane Surface with Upper 95% Confidence Surface"
    }
  }
  if (!up) {
    if (is.null(main)) {
      main = "Warped-Plane Surface with Lower 95% Confidence Surface"
    }
  }
  par(new = TRUE)
  #persp(x1p, x2p, surf, zlim = res$zlim, xlab = res$xlab, ylab = res$ylab, zlab = res$zlab, theta = res$theta, ltheta = res$ltheta, cex.axis = res$cex.axis, main = main, cex.main = cex.main, ticktype = res$ticktype, col=mycol)
  persp(x1p, x2p, surf, zlim = res$zlim, xlab = "", ylab = "", zlab = "",
        theta = res$theta, ltheta = res$ltheta, cex.axis = res$cex.axis,
        main = main, cex.main = cex.main, ticktype = res$ticktype, col=mycol, box=FALSE, axes=FALSE,...)
  par(new=FALSE)
}

#####
#wps#
#####
wps_getedf = function(ahati, sm, amat, amat0, xw0, xmat0, qmat0, p) {
	nz = 1:dim(amat)[1] < sm
	nz[amat %*% ahati > sm] = TRUE
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
			pa = gamj %*% solve(t(gamj) %*% gamj) %*% t(gamj)
			dj1 = dim(gamj)[1]; dj2 = dim(gamj)[2]
			#imat = matrix(0, nrow = dj1, ncol = dj1)
			#for (i in 1:dj1) {imat[i, i] = 1}
			imat = diag(dj1)
			uvecs = matrix(rnorm(dj1 * (dj1 - dj2)), nrow = dj1)
			wmatt = (imat - pa) %*% uvecs
			aqt = qr(wmatt); wmat = qr.Q(aqt)
			b0mat = xw0 %*% wmat
		}
	} else {
		b0mat = xw0
		inmat = diag(dim(xmat0)[2])
		wmat = inmat
	}
	p0mat = xmat0 %*% wmat %*% solve(t(wmat) %*% qmat0 %*% wmat) %*% t(wmat) %*% t(xmat0)
	edfi = sum(diag(p0mat)) + p - 1
	return (edfi)
}

########################################
#new wps.fit, with new amat and dmat
########################################
wps.fit = function(x1t, x2t, y, zmat = NULL, xmat_add = NULL, delta_add = NULL, delta_ut = NULL, varlist_add = NULL, shapes_add = NULL, np_add = 0, shapes = NULL, w = NULL, pen = 0, pnt = TRUE, cpar = 1.5, decrs = c(FALSE, FALSE), delta = NULL, kts = NULL, wt.iter = FALSE, family = gaussian(), cic = FALSE, nsim = 100, nprs = 1, idx_s = NULL, idx = NULL, gcv = FALSE, pvf = FALSE) {
	cicfamily = CicFamily(family)
	linkfun = cicfamily$linkfun
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
	#new: scale penalty term later
	#scy = var(y)
	one = 1:n*0 + 1
# 	if (is.null(zmat)) {
# 		zmat = matrix(one, ncol = 1)
# 	} else {
# 		if (dim(zmat)[1] != n) {
# 			stop ("Error: # rows of zmat must equal length of y")
# 		}
# 		zproj = zmat %*% solve(crossprod(zmat), t(zmat))
# 		onep = one - zproj %*% one
# # if the one vector is not in the space of zmat, then include the one vector in zmat
# 		if (sum(onep^2) > 1e-12) {
# 			zmat = cbind(one, zmat)
# 		}
# 	}
#p doesn't include the additive x
	p = 0
	if (!is.null(zmat)){
	  p = dim(zmat)[2]
	  dimnames(zmat)[[2]] = NULL
	}
#new: we have additive x's
#we already get warp.delta (delta); the 1st col of xmat is one
#bmat include conv, conc, shape = 17 and additive
#include bmat in zmat because we don't have penalty for bmat
	imat = diag(n)
  #delta0 = delta
	#dproj = delta%*%solve(crossprod(delta))%*%t(delta)
	#zproj = zmat %*% solve(crossprod(zmat), t(zmat))
	#zmat = (imat - dproj) %*% zmat
	#delta = (imat - zproj) %*% delta
	#save(delta,file='delta.Rda')
	xmat0 = delta
	bmat = NULL
	pb = 0
	if (!is.null(delta_add)) {
		bmat = t(delta_add)
		pb = ncol(bmat)
		zmat = cbind(bmat, zmat)
		xmat0 = cbind(xmat0, bmat)
	}
	#stop (print (xmat0[1,]))
#print (all(xmat0[, 1] == 1))
#xmat is additive, z, surface
	if (p >= 1) {
		#xmat = cbind(zmat, delta[, -1])
	  #imat = diag(n)
	  #zproj = zmat %*% solve(crossprod(zmat), t(zmat))
	  #delta = (imat - zproj) %*% delta
	  #zproj = zmat %*% solve(crossprod(zmat), t(zmat))
	  #onep = one - zproj %*% one
	  #onep = one %*% solve(crossprod(one), t(one))
	  #delta = (imat - onep) %*% delta
	  #zmat = (imat - onep)%*%zmat
	  xmat = cbind(zmat, delta)
	  #xmat = cbind(delta,zmat)
	} else {
		#xmat = xmat0
		xmat = delta
	}
# constraint matrix
if (nprs >= 1) {
	amat_lst = list()
	dmat_lst = list()
	varlist = NULL
	for (ipr in 1:nprs) {
		ktsi = kts[[ipr]]
		k1 = ktsi[[1]]
		k2 = ktsi[[2]]
		m1 = length(k1)
		m2 = length(k2)
		amat = matrix(0, nrow = 2*m1*m2 - m1 - m2, ncol = m1*m2)
		irow = 0
		for(j in 1:(m2-1)){
		  for(i in 1:m1){
		    irow = irow+1
		    amat[irow,m2*(i-1)+j]=-1
		    amat[irow,m2*(i-1)+j+1]=1
		  }
		}
		for(j in 1:m2){
		  for(i in 1:(m1-1)){
		    irow = irow+1
		    amat[irow,m2*(i)+j]=1
		    amat[irow,m2*(i-1)+j]=-1
		  }
		}
		amat_lst[[ipr]] = amat
# penalty matrix
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
		dmat_lst[[ipr]] = dmat
		#if (ipr > 1) {
#no constant for the 2nd pair
		# 	vari = 1:(m1*m2-1)*0 + ipr
		# } else {
		# 	vari = 1:m1*m2*0 + ipr
		# }
		vari = 1:m1*m2*0 + ipr
#this varlist is for surface
		varlist = c(varlist, vari)
	}
	amat = as.matrix(bdiag(amat_lst))
	dmat = as.matrix(bdiag(dmat_lst))
}
#amat0 is for the wps surface now
	amat0 = amat
#print (amat0[,1])
	nr0 = nrow(amat0); nc0 = ncol(amat0)
#amat include the original zmat
	if (p >= 1) {
		#amatz =  matrix(0, nrow = 2 * m1 * m2 - m1 - m2, ncol = p)
		#amat = cbind(amatz, amat[, 2:(m1 * m2)])
#p doesn't include additive components but include the one vector
		amatz = matrix(0, nrow = nrow(amat), ncol = p)
		#amat = cbind(amatz, amat[, -1])
		amat = cbind(amatz, amat)
		#amat = cbind(amat, amatz)
	}
#now make amat and amat0 to include additive components
	nr = nrow(amat); nc = ncol(amat)
	if (pb >= 1) {
		tmp = matrix(0, nrow = (nr+pb), ncol = (nc+pb))
		tmp0 = matrix(0, nrow = (nr0+pb), ncol = (nc0+pb))
		amatb = diag(pb)
#np_add is shape == 17 and conc (4,12); conv (3,11)
		if (np_add > 0) {
			amatb[1:np_add,1:np_add] = 0
		}
		tmp[1:pb,1:pb] = amatb
		tmp0[1:pb,1:pb] = amatb
		tmp[(pb+1):(nr+pb), (pb+1):(nc+pb)] = amat
		tmp0[(pb+1):(nr0+pb), (pb+1):(nc0+pb)] = amat0
		amat = tmp
		amat0 = tmp0
	}
	dmat0 = dmat
	#print (paste('dmat', dim(dmat)))
#not penalize the addtive part and z, zmat include both
	if (p >= 1) {
		dmatz = matrix(0, nrow = dim(dmat)[1], ncol = p)
		#dmat = cbind(dmatz, dmat[, -1])
		dmat = cbind(dmatz, dmat)
		#dmat = cbind(dmat, dmatz)
	}
	if (pb >= 1) {
		dmatb = matrix(0, nrow = dim(dmat)[1], ncol = pb)
		dmat = cbind(dmatb, dmat)
		dmat0 = cbind(dmatb, dmat0)
	}
#}
# weight
# xmat0 not include zmat;only used in edf
	#print (dim(xmat))
	#print (dim(amat))
	#print (dim(dmat))
	if (is.null(w)) {
		xw0 = xmat0
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
				zw[i, ] = zmat[i, ] * sqrt(w[i])
			}
		} else {
			xw0 = xmat0; xw = xmat; yw = y; zw = zmat
		}
	}
# transform to cone projection
	sc = 1
#print (pen)
	sm = 0
	if (round(pen, 6) > sm) {
		ps = pen
		if (pen > 1) {
		  warning('Large penalty term will make the inference wrong!')
		}
	} else if (pnt & (round(pen, 6) == 0)) {
#new: penalty interplolation
	  if (!gcv) {
	    ps = make_pen(n)
	    #print (ps)
	  } else {
	    #ps = make_pen(n, xw=xmat, xmat=xmat, dmat=dmat, y=y, amat=amat, gcv=gcv)
	    ng = 9
	    lams = 2^(0:8)
	    #lams = lams/2^8/n^(1/3)
	    #lams = 10*lams/2^8/n^(3/4)
	    lams = lams/2^8/n^(3/4)
	    #print (lams)
	    gcvs = 1:ng*0
	    for(i in 1:ng) {
	      # pen = lams[ipen]
	      # qv = t(xmat)%*%xmat+pen*t(dmat)%*%dmat
	      # umat = chol(qv)
	      # uinv = solve(umat)
	      # atil = amat%*%uinv
	      # zvec = t(uinv)%*%t(xmat)%*%y
	      # ans = coneA(zvec,atil)
	      # ahatc = uinv%*%ans$thetahat
	      # dimp = nrow(qv)
	      # if (length(ans$face)==0) {
	      #   pmat = diag(dimp)
	      # } else if (length(ans$face)==1){
	      #   aj = matrix(atil[ans$face,],nrow=1)
	      #   pmat = -t(aj)%*%solve(aj%*%t(aj))%*%aj
	      #   for(i in 1:dimp){pmat[i,i] = 1+pmat[i,i]}
	      # } else {
	      #   aj = atil[ans$face,]
	      #   pmat = -t(aj)%*%solve(aj%*%t(aj))%*%aj
	      #   for(i in 1:dimp){pmat[i,i] = 1+pmat[i,i]}
	      # }
	      # ## get df
	      # bigp = xmat%*%uinv%*%pmat%*%t(uinv)%*%t(xmat)
	      # degf = sum(diag(bigp))
	      # muhat = xmat%*%ahatc
	      # sse = sum((y-muhat)^2)
	      # gcvs[ipen] = sse/(1-degf/n)^2

	      pen = lams[i]
	      qv0 = crossprod(xw)
	      dv0 = crossprod(dmat)
	      qv = qv0 + pen * dv0
	      cv = crossprod(xw, y)
	      umat = chol(qv)
	      uinv = solve(umat)
	      atil = amat %*% uinv
	      cvec = t(uinv) %*% t(xw) %*% y
	      ansi = coneA(cvec, atil, msg = FALSE)
	      face = ansi$face
	      phihat = ansi$thetahat
	      ahat = uinv %*% phihat
	      muhat = xmat %*% ahat

	      sse = sum((y-muhat)^2)
	      dp = -atil
	      dp = t(dp)
	      imat = diag(nrow(qv))
	      if (length(face) == 0) {
	        pm = imat
	      } else {
	        smat = dp[,face,drop=FALSE]
	        pmat_polar = smat %*% solve(crossprod(smat), t(smat))
	        pm = (imat-pmat_polar)
	      }
	      bigp = xmat%*%uinv%*%pm%*%t(uinv)%*%t(xmat)
	      edfi = sum(diag(bigp))
	      gcvi = sse/(1-edfi/n)^2
	      gcvs[i] = gcvi
	    }
	    ps = min(lams[gcvs == min(gcvs)])
	  }
	} else if (!pnt) {
	  ps = 1e-6
	}
if (!wt.iter) {
	qmat = t(xw) %*% xw + ps * t(dmat) %*% dmat
	qmat0 = t(xw0) %*% xw0 + ps * t(dmat0) %*% dmat0
	umat = chol(qmat)
	uinv = solve(umat)
	atil = amat %*% uinv
	cvec = t(uinv) %*% t(xw) %*% yw
	ans = coneA(cvec, atil, msg = FALSE)
	face = ans$face
	phihat = ans$thetahat
	ahat = uinv %*% phihat
	#muhat = xw %*% ahat
	muhat = xmat %*% ahat
	muhatkeep = muhat
	etahatkeep = muhat
	coefkeep = ahat
	llh = llh.fun(y, muhatkeep, etahatkeep, phihat=NULL, n, w, fml = family$family)
} else {
	etahat = etahat.fun(n, y, fml = family$family)
	gr = gr.fun(y, etahat, weights = w, fml = family$family)
	wt = wt.fun(y, etahat, n, weights = w, fml = family$family)
	#cvec = crossprod(xw, (wt * etahat - gr))
	#qmat = t(xw) %*% diag(wt) %*% xw + ps * t(dmat) %*% dmat
	#qmat0 = t(xw0) %*% diag(wt) %*% xw0 + ps * t(dmat0) %*% dmat0
	cvec = crossprod(xmat, (wt * etahat - gr))
	qmat = t(xmat) %*% diag(wt) %*% xmat + ps * t(dmat) %*% dmat
	qmat0 = t(xmat0) %*% diag(wt) %*% xmat0 + ps * t(dmat0) %*% dmat0
	ans = qprog(qmat, cvec, amat, 1:nrow(amat)*0, msg = FALSE)
	face = ans$face
	ahat = ans$thetahat
	#etahat = xw %*% ahat
	etahat = xmat %*% ahat
	muhat = muhat.fun(etahat, fml = family$family)
	diff = 1
	nrep = 0
	sm = 1e-8
	while (diff > sm & nrep < 100) {
		oldmu = muhat
		nrep = nrep + 1
		gr = gr.fun(y, etahat, weights = w, fml = family$family)
		wt = wt.fun(y, etahat, n, weights = w, fml = family$family)
		#cvec = crossprod(xw, (wt * etahat - gr))
		#qmat = t(xw) %*% diag(wt) %*% xw + ps * t(dmat) %*% dmat
		#qmat0 = t(xw0) %*% diag(wt) %*% xw0 + ps * t(dmat0) %*% dmat0
		cvec = crossprod(xmat, (wt * etahat - gr))
		qmat = t(xmat) %*% diag(wt) %*% xmat + ps * t(dmat) %*% dmat
		qmat0 = t(xmat0) %*% diag(wt) %*% xmat0 + ps * t(dmat0) %*% dmat0
		ans = qprog(qmat, cvec, amat, 1:nrow(amat)*0, msg = FALSE)
		ahat = ans$thetahat
		#etahat = xw %*% ahat
		etahat = xmat %*% ahat
		muhat = muhat.fun(etahat, fml = family$family)
		diff = mean((muhat - oldmu)^2)
    }
    muhatkeep = muhat
	  etahatkeep = etahat
	  coefkeep = ahat
    llh = llh.fun(y, muhatkeep, etahatkeep, phihat=NULL, n, w, fml = family$family)
}
# get trace of "proj" matrix
# include the additive part
	sm = 1e-8
	nz = 1:dim(amat)[1] < sm
#only additive and wps will give nz = TRUE
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
			#imat = matrix(0, nrow = dj1, ncol = dj1)
			#for (i in 1:dj1) {imat[i, i] = 1}
			imat = diag(dj1)
			uvecs = matrix(rnorm(dj1 * (dj1 - dj2)), nrow = dj1)
			wmatt = (imat - pa) %*% uvecs
			aqt = qr(wmatt); wmat = qr.Q(aqt)
			b0mat = xw0 %*% wmat
		}
	} else {
		b0mat = xw0
		inmat = diag(dim(xmat0)[2])
		wmat = inmat
	}
	p0mat = xmat0 %*% wmat %*% solve(t(wmat) %*% qmat0 %*% wmat) %*% t(wmat) %*% t(xmat0)
	#edf = sum(diag(p0mat)) + p - 1
#new: to find edf differently
	#print (ps)
	dp = -atil
	#for (i in 1:nrow(dp)) {
	#  dpi = dp[i,,drop=F]
	#  dpi = dpi/sqrt(tcrossprod(dpi)[1])
	#  dp[i,]=dpi
	#}
	dp = t(dp)
	imat = diag(nrow(qmat))
	#print (dim(imat))
	#print (dim(dp))
	if (length(face) == 0) {
	  pmat = imat
	} else {
	  smat = dp[,face,drop=FALSE]
	  pmat_polar = smat %*% solve(crossprod(smat), t(smat))
	  pmat = (imat-pmat_polar)
	}
	#print (dim(smat))
	bigp = xmat%*%uinv%*%pmat%*%t(uinv)%*%t(xmat)
  #print (pmat)
	edf = sum(diag(bigp))

	#inmat = matrix(0, nrow = n, ncol = n)
	#diag(inmat) = 1
	#sse1 = sum((yw - xw %*% ahat)^2)
	#zcoefs = ahat[1:p]
	#sse0 = sum((yw - zw %*% zcoefs)^2)
#new: more than gaussian
#zcoefs don't include shape = 17
	thvecs = NULL
	thvecs_ut = NULL
	coef_add = 0
	coef_ut = 0
	coef_wp = ahat
	if (pb > 0) {
		#print (pb)
		coef_add = ahat[1:pb]
		coef_wp = ahat[-(1:pb)]
		capl = ncol(xmat_add)
		if (!is.null(xmat_add)){
			#print (paste('capl: ', capl))
			thvecs = matrix(0, nrow = capl, ncol = n)
	    	ncon = 0
	   	 	vcoef_add = coef_add[1:np_add]
	   	 	#print (vcoef_add)
	    	lconv = sum(shapes_add > 2 & shapes_add < 5 | shapes_add > 10 & shapes_add < 13)
	    	if (lconv > 0) {
	    		dcoefs = coef_add[-c(1:lconv)]
	    		delta_add2 = delta_add[-c(1:lconv), ,drop = FALSE]
	    	} else {
	   	 		dcoefs = coef_add
	    		delta_add2 = delta_add
	    	}
	    	for (i in 1:capl) {
	    		#thvecs[i,] = t(delta_add[varlist_add == i,]) %*% coef_add[varlist_add == i]
	    		thvecs[i,] = t(delta_add2[varlist_add == i,]) %*% dcoefs[varlist_add == i]
				if (shapes_add[i] > 2 & shapes_add[i] < 5 | shapes_add[i] > 10 & shapes_add[i] < 13) {
            		ncon = ncon + 1
					thvecs[i,] = thvecs[i,] + vcoef_add[ncon] * xmat_add[,i]
					#print (vcoef_add[ncon])
            	}
	    	}
	    	#if (!is.null(idx_s)) {
	    	if (length(idx_s) > 0) {
				thvecs0 = thvecs
				thvecs0[idx_s,] = thvecs[1:length(idx_s), ]
				#if (!is.null(idx)) {
				if (length(idx) > 0) {
					thvecs0[idx,] = thvecs[(1+length(idx_s)):capl, ]
				}
				thvecs = thvecs0
			}
		}
		if (!is.null(delta_ut)) {
			#print (length(coef_add))
			#print (nrow(delta_ut))
			#print ((length(coef_add) - nrow(delta_ut) + 1):length(coef_add))
				coef_ut = coef_add[(length(coef_add) - nrow(delta_ut) + 1):length(coef_add)]
				thvecs_ut = t(delta_ut) %*% coef_ut
		}
		if (!is.null(thvecs_ut)) {
				thvecs = rbind(thvecs, t(thvecs_ut))
		}
	}# else {coef_add = 0; coef_wp = ahat; thvecs = NULL}
	zcoefs = NULL
	if (!is.null(zmat)) {
	  zcoefs = ahat[(pb+1):(pb+p)]
	}
#vcoefs include zcoefs and shape = 17, conv and conc
	# if (np_add > 0) {
	# 	vcoefs = ahat[c(1:np_add, (pb+1):(pb+p))]
	# 	vmat = zw[,c(1:np_add, (pb+1):(pb+p)), drop=F]
	# } else {
	# 	vcoefs = zcoefs
	# 	vmat = zw[,(pb+1):(pb+p), drop=F]
	# }
#muvhat includes the weight part
	# muvhat = muhat.fun(vmat %*% vcoefs, fml = family$family)
	#sse1 = sum((yw - muhat)^2)
	sse1 = sum((yw - xw %*% ahat)^2)
	# sse0 = sum((yw - muvhat)^2)
#new use edf instead if (n - cpar * edf) < 0
	# np = ncol(vmat)
	np = 0
	if (!is.null(zmat)) {
	  np = ncol(zmat)
	}
	#if ((n - np - cpar * edf) <= 0) {
	if ((n - cpar * edf) <= 0) {
		sig2hat = sse1 / edf
	} else {
		sig2hat = sse1 / (n - cpar * edf)
		#print (sse1)
	}
	#print (paste('sig2hat', sig2hat))
	if (p >= 1) {
		#w2 = as.vector(w / deriv.fun(muhatkeep, fml = family$family))
	  inmat = diag(n)
		#pm = one %*% solve(crossprod(one)) %*% t(one)
		#covmat = solve(t(zmat[,(pb+1):(pb+p), drop=F]) %*% (inmat - p0mat + pm) %*% zmat[,(pb+1):(pb+p), drop=F])

	  #new:
	  # qmat0 = t(xw0) %*% xw0 + ps * t(dmat0) %*% dmat0
	  # umat0 = chol(qmat0)
	  # uinv0 = solve(umat0)
	  # atil0 = amat0 %*% uinv0
	  #
	  # dp0 = -t(atil0)
	  # imat0 = diag(nrow(qmat0))
	  # if (length(face) == 0) {
	  #   pmat0 = imat0
	  # } else {
	  #   smat0 = dp0[,face,drop=FALSE]
	  #   pmat_polar0 = smat0 %*% solve(crossprod(smat0), t(smat0))
	  #   pmat0 = (imat0-pmat_polar0)
	  # }
	  # bigp0 = xmat0%*%uinv0%*%pmat0%*%t(uinv0)%*%t(xmat0)
	  # zm = zmat[,(pb+1):(pb+p), drop=F]
	  # covmat = solve(t(zm) %*% (inmat - bigp0) %*% zm)

	  #test:
	  #one = 1:nrow(pmat)*0 + 1
	  #pm = one %*% solve(crossprod(one)) %*% t(one)
	  #print (dim(pmat))
	  #print (dim(pm))
	  mat1 = uinv %*% pmat %*% t(uinv) %*% t(xmat)
	  covmat0 = sig2hat * mat1 %*% t(mat1)
	  covmat = covmat0[(pb+1):(pb+p), (pb+1):(pb+p), drop=FALSE]
	  #print (covmat)
		#if (wt.iter) {
		#	sez = sqrt(diag(covmat))
		#} else {
		#sez = sqrt(diag(covmat) * sig2hat)
	  sez = sqrt(diag(covmat))
		#print (paste('sig:',sig2hat))
		#print (sez)
		#}
		tz = zcoefs / sez
		#print (zcoefs)
		#print (tz)
#new use edf instead if (n - cpar * edf) < 0
		if ((n - cpar * edf) <= 0) {
			pz = 2 * (1 - pt(abs(tz), edf))
			if (np > 1) {
				warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
			}
#print ('Check pz!')
		} else {
			pz = 2 * (1 - pt(abs(tz), n - cpar * edf))
			#print (edf)
			#print (cpar)
			#print (pz)
		}
	}
# get unconstrained penalized estimator
#if (!wt.iter) {
#	prmatu = xw %*% solve(qmat) %*% t(xw)
#new:
#	etahatu = prmatu %*% yw
#	muhatu = muhat.fun(etahatu, fml = family$family)
#	sseu = sum((yw - muhatu)^2)
#	edfu = sum(diag(prmatu))
#} else {
#	if (is.null(w)) {
#		w = 1:n*0 + 1
#	}
#	prior.w = w
#	prmatu = xw %*% solve(qmat) %*% t(xw)
#	nrep = 0
#	muhatu = mean(y) + 1:n*0
#	etahatu = linkfun(muhatu)
#	diff = 1
#	if (family$family == "binomial") {
#		mdiff = abs(max(muhatu) - 1) > sm
#	} else {mdiff = TRUE}
#	while (diff > sm & mdiff & nrep < n^2) {
#		nrep = nrep + 1
#		oldmu = muhatu
#		zhat = etahatu + (y - muhatu) * deriv.fun(muhatu, fml = family$family)
#		w = as.vector(prior.w * (deriv.fun(muhatu, fml = family$family))^(-1))
#		#b = solve(tvmat %*% vmat) %*% tvmat %*% zhat
#		#etahat = xmat %*% b
#		etahatu = prmatu %*% zhat
#		muhatu = muhat.fun(etahatu, fml = family$family)
#		diff = mean((muhatu - oldmu)^2)
#		mdiff = abs(max(muhatu) - 1)
#		if (family$family == "binomial") {
#			mdiff = abs(max(muhatu) - 1) > sm
#		} else {mdiff = TRUE}
#	}
#	sseu = sum((y - muhatu)^2)
#	edfu = sum(diag(prmatu))
#}

#new: get cic
	cic_val = NULL
	if (is.null(w)) {
		w = rep(1, n)
	}
	if (cic & nsim > 0) {
		edfs = 1:nsim*0
		if (!wt.iter) {
			for (isim in 1:nsim) {
				ysim = rnorm(n)
				ysimw = ysim * sqrt(w)
				cveci = t(uinv) %*% t(xw) %*% ysimw
				ansi = coneA(cveci, atil, msg = FALSE)
				phihati = ansi$thetahat
				ahati = uinv %*% phihati
				edfi = wps_getedf(ahati, sm, amat, amat0, xw0, xmat0, qmat0, np)
				edfs[isim] = edfi
			}
			#cic = llh + log(2 * (mean(edfs) + np) / (n - np - 1.5 * mean(edfs)) + 1)
			#cic = log(sse1) + log(2 * (mean(edfs) + np) / (n - np - 1.5 * mean(edfs)) + 1)
			cic_val = llh + log(2 * (mean(edfs)) / (n - np - 1.5 * (mean(edfs) - np)) + 1)
		} else {
			if (family$family == "poisson") {
				mu0 = mean(y)
			} else {mu0 = NULL}
			for (isim in 1:nsim) {
				ysim = ysim.fun(n, mu0, fml = family$family)
				etahat = etahat.fun(n, ysim, fml = family$family)
				gr = gr.fun(ysim, etahat, weights = w, fml = family$family)
				wt = wt.fun(ysim, etahat, n, weights = w, fml = family$family)
				cvec = crossprod(xmat, (wt * etahat - gr))
				qmat = t(xmat) %*% diag(wt) %*% xmat + ps * t(dmat) %*% dmat
				qmat0 = t(xmat0) %*% diag(wt) %*% xmat0 + ps * t(dmat0) %*% dmat0
				ans = qprog(qmat, cvec, amat, 1:nrow(amat)*0, msg = FALSE)
				ahat = ans$thetahat
				etahat = xmat %*% ahat
				muhat = muhat.fun(etahat, fml = family$family)
				diff = 1
				if (family$family == "binomial") {
					mdiff = abs(max(muhat) - 1) > 1e-8
				} else {mdiff = TRUE}
				nrep = 0
				while (diff > sm & nrep < 100 & mdiff) {
					oldmu = muhat
					nrep = nrep + 1
					gr = gr.fun(ysim, etahat, weights = w, fml = family$family)
					wt = wt.fun(ysim, etahat, n, weights = w, fml = family$family)
					cvec = crossprod(xmat, (wt * etahat - gr))
					qmat = t(xmat) %*% diag(wt) %*% xmat + ps * t(dmat) %*% dmat
					qmat0 = t(xmat0) %*% diag(wt) %*% xmat0 + ps * t(dmat0) %*% dmat0
					ansi = qprog(qmat, cvec, amat, 1:nrow(amat)*0, msg = FALSE)
					ahati = ansi$thetahat
					etahat = xmat %*% ahati
					muhat = muhat.fun(etahat, fml = family$family)
					diff = mean((muhat - oldmu)^2)
					if (family$family == "binomial") {
						mdiff = abs(max(muhat) - 1) > 1e-8
					} else {mdiff = TRUE}
    			 }
			 	edfi = wps_getedf(ahati, sm, amat, amat0, xw0, xmat0, qmat0, np)
			 	edfs[isim] = edfi
			}
			#cic = llh + log(2 * (mean(edfs) + np) / (n - np - 1.5 * mean(edfs)) + 1)
			cic_val = llh + log(2 * (mean(edfs)) / (n - np - 1.5 * (mean(edfs) - np)) + 1)
		}
	}
	#print (edf)
#new: test for constant vs wps
	#vmat = qr.Q(qr(t(atil)), complete = TRUE)[, -(1:(qr(t(atil))$rank)), drop = FALSE]
	vmat = cbind(1:n*0+1)
	if (!is.null(zmat)) {
	  vmat = cbind(vmat, zmat)
	}
	np = ncol(vmat)
	if (!is.null(w)) {
	  pvmat = vmat %*% solve(crossprod(vmat), t(vmat))
	} else {
	  vw = vmat
	  for (i in 1:nrow(vmat)) {
	    vw[i, ] = vmat[i, ] * sqrt(w[i])
	  }
	  pvmat = vmat %*% solve(crossprod(vw), t(vw))
	}
	#phihat0 = pvmat %*% cvec
	#ahat0 = uinv %*% phihat0
	#muhat0 = xmat %*% ahat0
	muhat0 = pvmat %*% y

	if (!is.null(w)) {
	  sse0 = sum(w * (y - muhat0)^2)
	  sse1 = sum(w * (y - muhat)^2)
	} else {
	  sse0 = sum((y - muhat0)^2)
	  sse1 = sum((y - muhat)^2)
	}
	bval = (sse0 - sse1) / sse0
	#print (bval)

	#g = qr(amat)
	#m = ncol(amat)
	#dim0 = m - g$rank

	#test = TRUE
	pval = NULL
	if (pvf) {
	  nloop = 200
	  sm = 1e-5
	  if (bval > sm) {
	    #bdist = 0:nloop*0
	    bdist = NULL
	    for (iloop in 1:nloop) {
	      #ysim = muhat0 + rnorm(n)
	      #new: add (sig2hat)^(1/2) when there's z
	      ysim = muhat0 + rnorm(n)*(sig2hat)^(1/2)
	      ysimw = ysim
	      if (!is.null(w)) {
	        ysimw = ysim * sqrt(w)
	      }
	      cveci = t(uinv) %*% t(xw) %*% ysimw
	      ansi = try(coneA(cveci, atil))
	      if (any(class(ansi) %in% 'try-error')) {
	        next
	      }
	      phi = ansi$thetahat
	      ah = uinv %*% phi
	      muh = xmat %*% ah

	      #phi0 = pvmat %*% cveci
	      #ah0 = uinv %*% phi0
	      #muh0 = xmat %*% ah0
	      #muh0 = pvmat %*% ysimw
	      muh0 = pvmat %*% ysim

	      if (!is.null(w)) {
	        sse0 = sum(w * (ysim - muh0)^2)
	        sse1 = sum(w * (ysim - muh)^2)
	      } else {
	        sse0 = sum((ysim - muh0)^2)
	        sse1 = sum((ysim - muh)^2)
	      }
	      bstat = (sse0 - sse1) / sse0
	      bdist = c(bdist, bstat)
	      #bdist[iloop] = bstat
	    }
	    pval = sum(bdist > bval)/nloop
	  } else {pval = 1}
	}

	ans = new.env()
	ans$pval = pval
	ans$bval = bval
	#ans$k1 = k1
	#ans$k2 = k2
	for (ipr in 1:nprs) {
		decri = decrs[[ipr]]
		ktsi = kts[[ipr]]
		k1 = ktsi[[1]]
		k2 = ktsi[[2]]
		if (decri[1]) {
			k1 = -rev(k1)
		}
		if (decri[2]) {
			k2 = -rev(k2)
		}
		ktsi[[1]] = k1
		ktsi[[2]] = k2
		kts[[ipr]] = ktsi
	}
	ans$kts = kts
	ans$muhat = muhatkeep
	ans$etahat = etahatkeep
	#ans$muplot = mupl
	#ans$muhatu = muhatu
	#ans$muplotu = muplu
	ans$gcv = sse1 / (1 - edf / n)^2
	#ans$gcvu = sseu / (1 - edfu/n)^2
	#ans$ssr = sse1
	ans$sse1 = sse1
	#ans$sse0 = sse0
	ans$edf = edf
	if (cic & nsim > 0) {
		ans$edf0 = mean(edfs) #+ np
	} else {ans$edf0 = np}
	#ans$nz = amat %*% ahat
	#ans$edfu = edfu
	ans$coef_add = coef_add
	ans$coef_ut = coef_ut
	ans$coef_wp = coef_wp
	ans$coefs = coefkeep
	#ans$coefsu = solve(qmat) %*% t(xw) %*% yw
#include coef for one vector
	#ans$zcoefs = ahat[1:p]
	ans$zcoefs = zcoefs
	ans$zmat = zmat
	#ans$zmat_0 = zmat_0
	ans$sig2hat = sig2hat
	ans$delta = xmat
	ans$pen = ps
	ans$d0 = p + np_add
	#print (paste('p: ', p))
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
	ans$cic = cic_val
	ans$varlist = varlist
	ans$etacomps = thvecs
  #new:
  ans$amat = amat
  ans$dmat = dmat
  #ans$vmat = vmat
	return (ans)
}

####################################################################
#four monotonicity functions for warped-plane fit
####################################################################
s.incr.incr <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
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
    attr(xm, "categ") <- "warp"
    #class(xm) <- "warp"
    #warp <<- TRUE
    return (xm)
}

s.incr.decr <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "wps_id"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "decreasing") <- c(FALSE, TRUE)
    attr(xm, "categ") <- "warp"
    #warp <<- TRUE
    #class(xm) <- "warp"
    return (xm)
}

s.decr.incr <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
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
    attr(xm, "categ") <- "warp"
    #warp <<- TRUE
    #class(xm) <- "warp"
    return (xm)
}

s.decr.decr <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
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
    attr(xm, "categ") <- "warp"
    #warp <<- TRUE
    #class(xm) <- "warp"
    return (xm)
}

###############################################################
#makedelta_wps: make delta to a pair of warped-plane variables#
###############################################################
makedelta_wps <- function(x1t, x2t, m1_0 = 0, m2_0 = 0, k1 = 0, k2 = 0, space = c("E", "E"), decreasing = c(FALSE, FALSE), interp = FALSE)
{
    # x1 and x2 no need to sort
    # if decreasing  (all calculations done for doubly-increasing case)
    n = length(x1t)
    if (decreasing[1]) {
        #print (k1 != 0)
        x1 = -x1t
        #if (!is.null(k1)) {
        if (!all(k1 == 0)) {
            m1 = length(k1); k1 = -k1[m1:1]
        }
    } else {x1 = x1t}
    if (decreasing[2]) {
        x2 = -x2t
        #if (!is.null(k2)) {
        if (!all(k2 == 0)) {
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
            #m1 = 2 * round(n^(1/6)) + 4
            #new:
            #m1 = round(6*n^(1/6))
            #m1 = round(4*n^(1/6))
            m1 = round(5*n^(1/6))
            #print (m1)
        } else {m1 = m1_0}
        if (space[1] == "Q") {
            k1 = quantile(unique(x1), probs = seq(0, 1, length = m1), names = FALSE)
        }
        if (space[1] == "E") {
            k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) + min(x1)
            #print (0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) +  min(x1) )
            #k1 = 0:(m1 - 1) / (m1 - 1)
        }
        #k1 = 0:(m1 - 1) / (m1 - 1) * (max(x1) - min(x1)) + min(x1)
    } else { m1 = length(k1) }
    if (make2) {
        #new:
        if (m2_0 < 4 | round(m2_0, 0) != m2_0) {
            if (m2_0 != 0) {
                warning ('At least four knots should be used! Number of knots is re-defined!')
            }
            #m2 = 2 * round(n^(1/6)) + 4
            #m2 = round(6*n^(1/6))
            #m2 = round(4*n^(1/6))
            m2 = round(5*n^(1/6))
        } else {m2 = m2_0}
        if (space[2] == "Q") {
            k2 = quantile(unique(x2), probs = seq(0, 1, length = m2), names = FALSE)
        }
        if (space[2] == "E") {
            k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2) - min(x2)) + min(x2)
            #k2 = 0:(m2 - 1) / (m2 - 1)
        }
        #k2 = 0:(m2 - 1) / (m2 - 1) * (max(x2) - min(x2)) + min(x2)
    } else { m2 = length(k2) }
    ## check to see if empty knot intervals
    #new: if it's for prediction, then skip
    if (!interp) {
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
    }
    # make the basis functions
    #b1 = matrix(0, nrow = n, ncol = m1)
    #b2 = matrix(0, nrow = n, ncol = m2)
    #for (i in 2:(m1 - 1)) {
    #    i1 = x1 >= k1[i-1] & x1 <= k1[i]
    #    b1[i1, i] = (x1[i1] - k1[i-1]) / (k1[i] - k1[i-1])
    #    i2 = x1 > k1[i] & x1 <= k1[i+1]
    #    b1[i2, i] = (k1[i+1] - x1[i2]) / (k1[i+1] - k1[i])
    #}
    # i1 = x1 >= k1[1] & x1 <= k1[2]
    # b1[i1, 1] = (k1[2] - x1[i1]) / (k1[2] - k1[1])
    # i2 = x1 > k1[m1-1] & x1 <= k1[m1]
    # b1[i2, m1] = (x1[i2] - k1[m1-1]) / (k1[m1] - k1[m1-1])
    #
    # for (i in 2:(m2 - 1)) {
    #     i1 = x2 >= k2[i-1] & x2 <= k2[i]
    #     b2[i1, i] = (x2[i1] - k2[i-1]) / (k2[i] - k2[i-1])
    #     i2 = x2 > k2[i] & x2 <= k2[i+1]
    #     b2[i2, i] = (k2[i+1] - x2[i2]) / (k2[i+1] - k2[i])
    # }
    # i1 = x2 >= k2[1] & x2 <= k2[2]
    # b2[i1, 1] = (k2[2] - x2[i1]) / (k2[2] - k2[1])
    # i2 = x2 > k2[m2-1] & x2 <= k2[m2]
    # b2[i2, m2] = (x2[i2] - k2[m2-1]) / (k2[m2] - k2[m2-1])
    # ## design matrix
    # xmat0 = matrix(nrow = n, ncol = m1 + m2 - 1 + (m1 - 1) * (m2 - 1))
    # xmat0[ ,1] = 1:n*0 + 1
    # xmat0[ ,2:m1] = b1[ ,2:m1]
    # xmat0[ ,(m1 + 1):(m1 + m2 - 1)] = b2[ ,2:m2]
    # for (i in 1:(m1 - 1)) {
    #     xmat0[ ,(m1 + m2 + (i - 1) * (m2 - 1)):(m1 + m2 - 1 + i * (m2 - 1))] = b1[ ,i + 1] * b2[ ,2:m2]
    # }
    #if (dim(zmat)[2] > 1) {
    #	xmat = cbind(zmat, xmat0[ ,2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))])
    #} else {xmat = xmat0}
    #bmat = xmat0[ ,2:(m1 + m2 - 1 + (m1 - 1) * (m2 - 1))]
    # ignore the zmat for now
    # xmat = xmat0
    # columns of delta are edges, different from other make_delta in cgam
    bmat = matrix(0, nrow=n, ncol=m1*m2)
    for(i in 1:(m1-1)){
      for(j in 1:(m2-1)){
        ii=x1>=k1[i]&x1<=k1[i+1]&x2>=k2[j]&x2<=k2[j+1]
        #rg_ij = (k1[i+1]-k1[i])*(k2[j+1]-k2[j])
        bmat[ii,m2*(i-1)+j]=(k1[i+1]-x1[ii])*(k2[j+1]-x2[ii])/((k1[i+1]-k1[i])*(k2[j+1]-k2[j]))
        bmat[ii,m2*(i-1)+j+1]=(k1[i+1]-x1[ii])*(x2[ii]-k2[j])/((k1[i+1]-k1[i])*(k2[j+1]-k2[j]))
        bmat[ii,m2*i+j]=(x1[ii]-k1[i])*(k2[j+1]-x2[ii])/((k1[i+1]-k1[i])*(k2[j+1]-k2[j]))
        bmat[ii,m2*i+j+1]=(x1[ii]-k1[i])*(x2[ii]-k2[j])/((k1[i+1]-k1[i])*(k2[j+1]-k2[j]))
      }
    }

    ans = new.env()
    ans$delta = bmat
    ans$k1 = k1
    ans$k2 = k2
    return (ans)
    #attr(delta, "shape") = "warp"
    #delta
}

##############################
#penalty interpolation
##############################
make_pen = function(n, xw=NULL, xmat=NULL, dmat=NULL, y=NULL, amat=NULL, gcv=FALSE) {
  if (!gcv) {
    if (n <= 100) {
      lambda = 0.06
    }
    if (n >= 5000) {
      lambda = 0.026
    }
    if (n > 100 & n < 5000) {
      ns = c(1, 2, 4, 8, 10, 20, 50)*100
      lams = c(6, 5, 4, 3.4, 3.2, 2.8, 2.6)/100
      diff_vec = sign(n - ns)
      st = rev(which(diff_vec == 1L))[1]
      ed = which(diff_vec == -1L)[1]

      my_line = function(xp = NULL, y, x, end=2, start=1) {
        slope = NULL
        intercept = NULL
        yp = NULL
        slope = (y[end] - y[start]) / (x[end] - x[start])
        intercept = y[end] - slope * x[end]
        yp = intercept + slope * xp
        return (yp)
      }
      yvec = c(lams[st], lams[ed])
      xvec = c(ns[st], ns[ed])
      lambda = my_line(xp = n, y = yvec, x = xvec)
    }
  } else {
    #use gcv to find lambda
    #ng = 20
    #lams = seq(1e-4, 1, length=20)
    #lams = 2^(1:ng)
    #lams = 2*lams/max(lams)
    #gcvs = 1:ng*0

    ng = 9
    lams = 2^(0:8)
    lams = lams/2^8/n^(1/3)
    #print (lams)
    gcvs = 1:ng*0

    for(i in 1:ng) {
      pen = lams[i]
      qv0 = crossprod(xw)
      dv0 = crossprod(dmat)
      qv = qv0 + pen * dv0
      cv = crossprod(xw, y)
      umat = chol(qv)
      uinv = solve(umat)
      atil = amat %*% uinv
      cvec = t(uinv) %*% t(xw) %*% y
      ansi = coneA(cvec, atil, msg = FALSE)
      face = ansi$face
      phihat = ansi$thetahat
      ahat = uinv %*% phihat
      muhat = xmat %*% ahat

      sse = sum((y-muhat)^2)
      dp = -atil
      dp = t(dp)
      imat = diag(nrow(qv))
      if (length(face) == 0) {
        pm = imat
      } else {
        smat = dp[,face,drop=FALSE]
        pmat_polar = smat %*% solve(crossprod(smat), t(smat))
        pm = (imat-pmat_polar)
      }
      bigp = xmat%*%uinv%*%pm%*%t(uinv)%*%t(xmat)
      #tst = bigp%*%y
      #print (all.equal(tst, muhat))
      edfi = sum(diag(bigp))
      #edfi = sum(diag(xw %*% solve(qv) %*% t(xw)))
      gcvi = sse/(1-edfi/n)^2
      gcvs[i] = gcvi

      #print (sse)
    }
    #plot(lams, gcvs, type='o')
    #lambda = (lams[which.min(gcvs)])[1]
    lambda = min(lams[gcvs == min(gcvs)])
  }
  return (lambda)
}

#######################
#amat and penalty bmat#
#######################
makeamat_wps <- function(kts, nprs) {
#if (nprs >= 1) {
	amat_lst = list()
	dmat_lst = list()
	varlist = NULL
	for (ipr in 1:nprs) {
		ktsi = kts[[ipr]]
		k1 = ktsi[[1]]
		k2 = ktsi[[2]]
		m1 = length(k1)
		m2 = length(k2)
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
		if (ipr > 1) {
			amat = amat[,-1]
		}
		amat_lst[[ipr]] = amat
# penalty matrix
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
		if (ipr > 1) {
			dmat = dmat[,-1]
		}
		#Dmat = cbind(Dmat, dmat)
		dmat_lst[[ipr]] = dmat
		if (ipr > 1) {
#no constant for the 2nd pair
			vari = 1:(m1*m2-1)*0 + ipr
		} else {
			vari = 1:m1*m2*0 + ipr
		}
		varlist = c(varlist, vari)
	}
	amat = as.matrix(bdiag(amat_lst))
	#dmat = Dmat
	dmat = as.matrix(bdiag(dmat_lst))
	ans = list(amat = amat, dmat = dmat, varlist = varlist)
	return (ans)
#}
}
#################
#new coef method#
#################
coef.cgam <- function(object,...) {
  ans <- object$coefs
  ans
}

coef.wps <- function(object,...) {
  ans <- object$coefs
  ans
}

coef.trispl <- function(object,...) {
  ans <- object$coefs
  ans
}

coef.cgam.polr <- function(object,...) {
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

fitted.trispl <- function(object,...) {
  ans <- object$muhat
  ans
}

#fitted.cgam.polr <- function(object,...) {
#  ans <- object$muhat
#  ans
#}
fitted.cgam.polr <- function(object,...) {
	a <- object$zeta
	eta <- object$muhat
	lev <- object$lev
	nc <- length(a) + 1
	n <- length(eta)
	ps <- matrix(0, nrow = nc, ncol = n)
	for (i in 1:nc) {
		if (i == 1) {
			ps[i,] <- pfun(a[1] - eta)
		} else if (i == nc) {
			ps[i,] <- 1 - pfun(a[nc-1] - eta)
		} else {
			ps[i,] <- pfun(a[i] - eta) - pfun(a[i-1] - eta)
		}
	}
	ps <- t(ps)
	dimnames(ps) <- list(1:n, lev)
	return (ps)
}

#############################
#shape selection part
#############################
#get(x = "s", pos = "package:cgam")
#new: add npop; per.mutate
ShapeSelect <- function(formula, family = gaussian, cpar = 2, data = NULL, weights = NULL, npop = 200, per.mutate = 0.05, genetic = FALSE) {
	#if (exists("s", parent.frame()) & class(get("s", envir = parent.frame())) == "function") {
  if (exists("s", parent.frame()) & inherits(get("s", envir = parent.frame()), "function")) {
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
		#if (class(y) == "factor") {
  	if(inherits(y, "factor")){
			y <- ifelse(y == levels(y)[1], 0, 1)
		}
#new: test
		#if (class(y) == "character") {
		if(inherits(y, "character")){
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
		ans <- ConstrGA(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, nmod = nmod, zfacs = zfacs, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat, npop = npop, per.mutate = per.mutate)
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
			invisible(ConstrGA(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, nmod = nmod, zfacs = zfacs, time.est = TRUE, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat, npop = npop, per.mutate = per.mutate))
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
				ans <- ConstrGA(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, nmod = nmod, zfacs = zfacs, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat, npop = npop, per.mutate = per.mutate)
				tm <- proc.time() - ptm
			}
#too many models or memory problem
		} else {
			print ("Go genetic! Too many models!")
			ptm <- proc.time()
			ans <- ConstrGA(y, xmat, zmat, trmat, family = family, shpsx = shpsx, shpsvx = shpsvx, shpsz = shpsz, shpst = shpst, cpar = cpar, nmod = nmod, zfacs = zfacs, weights = weights, vzmat = vzmat, vzfacs = vzfacs, vxmat = vxmat, vtrmat = vtrmat, npop = npop, per.mutate = per.mutate)
			tm <- tm + proc.time() - ptm
		}
	}
#print (vznms)
	colnames(ans$pop2)[1:(capl+capk+capt)] = c(xnms, znms, trnms)
	rslt <- list(pop = ans$pop2, top = ans$pop2[1,], fitness = ans$fitness, tm = tm, xnms = xnms, znms = znms, trnms = trnms, zfacs = zfacs, mxf = ans$mxf, mnf = ans$mnf, GA = ans$GA, vzcat = ans$vzcat, vzmat = vzmat, shpsx = shpsx, vxnms = vxnms, vznms = vznms, vtrnms = vtrnms)
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
ConstrGA = function(y, xmat, zmat, trmat, family = gaussian, shpsx = NULL, shpsvx = NULL, shpsz = NULL, shpst = NULL, cpar = 1.2, nmod = 2e+4, zfacs = NULL, time.est = FALSE, weights = NULL, vzmat = NULL, vzfacs = NULL, vxmat = NULL, vtrmat = NULL, npop = NULL, per.mutate = NULL) {
	#linkfun = family$linkfun
	cicfamily = CicFamily(family)
	linkfun = cicfamily$linkfun
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
        #if (nmod < 2e+6) {
        #    npop = 200
        #} else {
        #    npop = 200 + round((nmod - 2e+6) / nmod * 20)
        #}
        #new: now we allow user to speficy npop
        if (npop <= 0 | round(npop) != npop) {
            warning('npop must be a positive integer! npop=200 will be used!')
            if (nmod < 2e+6) {npop = 200} else {npop = 200 + round((nmod - 2e+6) / nmod * 20)}
        }
        #if (npop > nmod) {
        #    warning('npop must be smaller than the total number of possible models! npop=200 will be used!')
        #    if (nmod < 2e+6) {npop = 200} else {npop = 200 + round((nmod - 2e+6) / nmod * 20)}
        #}
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
				llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), phihat=NULL, n = n, weights = weights, fml = family$family)
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
        #imut = trunc(runif(nm) * (npop - 1) + 2)
        #new for imut
        if (per.mutate <= 0 | per.mutate >= 0.5) {
            warning('per.mutate should be a percentage > 0 and < 0.5! per.mute = 0.05 will be used!')
            per.mutate = 0.05

        }
        nm = round(npop * per.mutate)
        imut = sample(1:npop, size = nm, replace = FALSE)
        #print (nm)
        #print (paste('imute:',imut))
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
					llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), phihat=NULL, n = n, weights = weights, fml = family$family)
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
					llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), phihat=NULL, n = n, weights = weights, fml = family$family)
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
					llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), phihat=NULL, n = n, weights = weights, fml = family$family)
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
								llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), phihat=NULL, n = n, weights = weights, fml = family$family)
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
ConstrALL = function(y, xmat, zmat, trmat, family = gaussian, shpsx = NULL, shpsvx = NULL, shpsz = NULL, shpst = NULL, cpar = 1.2, zfacs = NULL, weights = NULL, vzmat = NULL, vzfacs = NULL, vxmat = NULL, vtrmat = NULL) {
	#linkfun = family$linkfun
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
		llh = llh.fun(y = y, muhat = rep(mean(y), n), etahat = rep(mean(y), n), phihat=NULL, n = n, weights = weights, fml = family$family)
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
irls = function(y, bigmat, np, nsim = 100, family = gaussian, weights = NULL, cpar = 1.2) {
	#linkfun = family$linkfun
	cicfamily = CicFamily(family)
	linkfun = cicfamily$linkfun
	llh.fun = cicfamily$llh.fun
	etahat.fun = cicfamily$etahat.fun
	gr.fun = cicfamily$gr.fun
	wt.fun = cicfamily$wt.fun
	zvec.fun = cicfamily$zvec.fun
	muhat.fun = cicfamily$muhat.fun
	ysim.fun = cicfamily$ysim.fun
	deriv.fun = cicfamily$deriv.fun
	dev.fun = cicfamily$dev.fun
	if (family$family == "binomial" | family$family == "poisson" | family$family == "Gamma") {
		wt.iter = TRUE
	} else {wt.iter = FALSE}
	n = length(y)
	if (is.null(weights)) {
		weights = 1:n*0 + 1
	}
	sm = 1e-8
	m = dim(bigmat)[1] - np
# new: initialize cvec
	cvec = NULL
	if (wt.iter) {
		etahat = etahat.fun(n, y, fml = family$family)
		gr = gr.fun(y, etahat, weights, fml = family$family)
		wt = wt.fun(y, etahat, n, weights, fml = family$family)
		cvec = wt * etahat - gr
	} else {wt = wt.fun(y, etahat, n, weights, fml = family$family)}
	zvec = zvec.fun(cvec, wt, y, fml = family$family)
    gmat = t(bigmat %*% sqrt(diag(wt)))
	if (m > 0) {
#np always >= 1
		dsend = gmat[, (np + 1):(np + m), drop = FALSE]
        zsend = gmat[, 1:np, drop = FALSE]
        #ans = coneB(zvec, t(dsend), zsend, msg = FALSE)
        ans = coneB(zvec, dsend, zsend, msg = FALSE)
        face = ans$face
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
				wt = wt.fun(y, etahat, n, weights, fml = family$family)
				cvec = wt * etahat - gr
				zvec = zvec.fun(cvec, wt, y, fml = family$family)
				gmat = t(bigmat %*% sqrt(diag(wt)))
				dsend = gmat[, (np + 1):(np + m), drop = FALSE]
        		zsend = gmat[, 1:np, drop = FALSE]
                #ans = coneB(zvec, t(dsend), zsend, msg = FALSE)
                ans = coneB(zvec, dsend, zsend, msg = FALSE, face = face)
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
	llh = llh.fun(y, muhat, etahat, phihat=NULL, n, weights, fml = family$family)
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
				wt = wt.fun(ysim, etahat, n, weights, fml = family$family)
				cvec = wt * etahat - gr
			} else {wt = wt.fun(ysim, etahat, n, weights, fml = family$family)}
			zvec = zvec.fun(cvec, wt, ysim, fml = family$family)
            gmat = t(bigmat %*% sqrt(diag(wt)))
           	dsend = gmat[, (np + 1):(np + m), drop = FALSE]
            zsend = gmat[, 1:np, drop = FALSE]
            #ans = try(coneB(zvec, t(dsend), zsend, msg = FALSE))
            ans = try(coneB(zvec, dsend, zsend, msg = FALSE))
            face = ans$face
			#if (class(ans) == "try-error") next
      if(inherits(ans, "try-error")) next
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
					wt = wt.fun(ysim, etahat, n, weights, fml = family$family)
					cvec = wt * etahat - gr
					#zvec <- cvec / sqrt(wt)
					zvec = zvec.fun(cvec, wt, y, fml = family$family)
					gmat = t(bigmat %*% sqrt(diag(wt)))
					dsend = gmat[, (np + 1):(np + m), drop = FALSE]
					zsend = gmat[, 1:np, drop = FALSE]
                    #ans = try(coneB(zvec, t(dsend), zsend, msg = FALSE))
                    ans = try(coneB(zvec, dsend, zsend, msg = FALSE, face = face))
					#if (class(ans) == "try-error") next
          if(inherits(ans, "try-error")) next
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

#### TRIANGLE SPLINE ######
##  y is response
##  x1, x2 are continuous predictors
##  zmat is matrix of parametrically modeled covariates
##  for fixed covariate values, E(y) is convex in (x1,x2) -- nonadditive!
####
trispl.fit = function(x1t, x2t, y, zmat = NULL, xmat_add = NULL, delta_add = NULL, varlist_add = NULL, shapes_add = NULL, np_add = 0, xmat_wp = NULL, delta_wp = NULL, varlist_wps = NULL, amat_wp = NULL, dmat_wp = NULL, w = NULL, lambda = 0, pnt = TRUE, cpar = 1.2, cvss = c(TRUE, TRUE), delta = NULL, kts = NULL, nkts = NULL, wt.iter = FALSE, family = gaussian(), nsim = 0, nprs = 1) {
#additive + tri surface + z(exclude one) + warp surface
	cicfamily = CicFamily(family)
	linkfun = cicfamily$linkfun
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
#zmat not includes one
	if (!is.null(zmat)){
		if (length(zmat) == n) {
			zmat = matrix(zmat, ncol = 1)
		}
		p = dim(zmat)[2]
	} else{p = 0}
	#print (p)
	#one = 1:n*0 + 1
	#if (is.null(zmat)) {
	#	zmat = matrix(one, ncol = 1)
	#} else {
	#	if (dim(zmat)[1] != n) {
	#		stop ("Error: # rows of zmat must equal length of y")
	#	}
	#	zproj = zmat %*% solve(crossprod(zmat), t(zmat))
	#	onep = one - zproj %*% one
# if the one vector is not in the space of zmat, then include the one vector in zmat
	#	if (sum(onep^2) > 1e-12) {
	#		zmat = cbind(one, zmat)
	#	}
	#}
	#p = dim(zmat)[2]
	dimnames(zmat)[[2]] = NULL
	if (nprs >= 1) {
		amat_lst = list()
#penalty matrix
		#Pmat = NULL
		pmat_lst = list()
		knots_lst = list()
		m12_lst = list()
		trimat_lst = list()
		bmat_lst = list()
		capk_lst = list()
#design matrix
		Dmat = NULL
		varlist = NULL
		#print (nkts)
		for (ipr in 1:nprs) {
			cvs = cvss[[ipr]]
			x1 = x1t[,ipr]
			x2 = x2t[,ipr]
			nktsi = nkts[[ipr]]
			if (length(x1) != n | length(x2) != n) {
				stop ("ERROR -- x1, x2, and y must have same lengths")
			}
			xmat = cbind(x1, x2)
#temp:
			#if (nktsi[1] > 0) {
			#	m1 = nktsi[i]
			#} else {
				#m1 = 2 + round(n^(1/3))
			  m1 = round(5*n^(1/6))
			#}
			s = (max(x1) - min(x1)) / (m1-1)
			h = s*sqrt(3) / 2
			m2 = trunc((max(x2) - min(x2)) / h) + 2
			m0 = trunc(m2/2)
			capk = 2*m0*m1+m0
			if (m0 < m2/2) {
				capk = capk+m1+1
			}
			knots = matrix(0, nrow = capk, ncol = 2)
			l1 = min(x1)-s/2
			u1 = max(x1)+s/2
			l2 = min(x2)
			u2 = l2 + (m2-1)*h
			rw = 0
			for (i in 1:m0) {
				for (j in 1:(m1 + 1)) {
					rw = rw + 1
					knots[rw,1] = l1 + (j - 1) * s
					knots[rw,2] = l2 + (i - 1) * 2 * h
				}
				for (j in 1:m1) {
					rw = rw + 1
					knots[rw,1] = l1 + s/2 + (j - 1) * s
					knots[rw,2] = l2 + (i - 1) * 2 * h + h
				}
			}
			if (m0 < m2/2) {
				i = m0 + 1
				for (j in 1:(m1 + 1)) {
					rw = rw + 1
					knots[rw,1] = l1 + (j - 1) * s
					knots[rw,2] = l2 + (i - 1) * 2 * h
				}
			}
			knots[,2] = knots[,2] - (u2 - max(x2))/2
			u2 = max(knots[,2])
			l2 = min(knots[,2])
#	plot(c(l1,u1),c(l2,u2),pch = "")
#	points(x1,x2,col = "slategray")
#	points(knots[,1],knots[,2],pch = 20)
### now number the triangles and give knots for each
			ntri = (m2 - 1) * (2 * m1 - 1)
			trimat = matrix(0, nrow = ntri, ncol = 3)
			if (m0 == m2/2) {
				mlim = m0 - 1
			} else {mlim = m0}
			rt = 0
			for (i in 1:mlim) {
				for (j in 1:m1) {
					rt = rt + 1
					trimat[rt,1] = (2*m1 + 1)*(i - 1) + j
					trimat[rt,2] = (2*m1 + 1)*(i - 1) + j + 1
					trimat[rt,3] = (2*m1 + 1)*(i - 1) + m1 + 1 + j
				}
				for (j in 1:(m1 - 1)) {
					rt = rt + 1
					trimat[rt,1] = (2*m1 + 1)*(i - 1) + j + 1
					trimat[rt,2] = (2*m1 + 1)*(i - 1) + m1 + 1 + j
					trimat[rt,3] = (2*m1 + 1)*(i - 1) + m1 + 2 + j
				}
				for (j in 1:(m1 - 1)) {
					rt = rt + 1
					trimat[rt,1] = (2*m1 + 1)*i + j + 1
					trimat[rt,2] = (2*m1 + 1)*(i - 1) + m1 + 1 + j
					trimat[rt,3] = (2*m1 + 1)*(i - 1) + m1 + 2 + j
				}
				for (j in 1:m1) {
					rt = rt + 1
					trimat[rt,1] = (2*m1 + 1)*(i - 1) + m1 + 1 + j
					trimat[rt,2] = (2*m1 + 1)*i + j
					trimat[rt,3] = (2*m1 + 1)*i + j + 1
				}
			}
			if (m0 == m2/2) {
				i = m0
				for (j in 1:m1) {
					rt = rt + 1
					trimat[rt,1] = (2*m1 + 1)*(i - 1) + j
					trimat[rt,2] = (2*m1 + 1)*(i - 1) + j + 1
					trimat[rt,3] = (2*m1 + 1)*(i - 1) + m1 + 1 + j
				}
				for (j in 1:(m1 - 1)) {
					rt = rt + 1
					trimat[rt,1] = (2*m1 + 1)*(i - 1) + j + 1
					trimat[rt,2] = (2*m1 + 1)*(i - 1) + m1 + 1 + j
					trimat[rt,3] = (2*m1 + 1)*(i - 1) + m1 + 2 + j
				}
			}
#### for each triangle, list adjacent knots (up to three)
#print (trimat)
			adjk = matrix(0,nrow = ntri,ncol = 3)
			kij = 1:6
			s1 = 1:3
			s2 = 1:3
			found = 1:ntri*0
			for (i in 1:(ntri - 1)) {
				for (j in (i + 1):ntri) {
					kij[1:3] = trimat[i,]
					kij[4:6] = trimat[j,]
					if (length(unique(kij)) == 4) {
						found[i] = found[i] + 1
						found[j] = found[j] + 1
						s1[1] = sum(kij == trimat[i,1])
						s1[2] = sum(kij == trimat[i,2])
						s1[3] = sum(kij == trimat[i,3])
						s2[1] = sum(kij == trimat[j,1])
						s2[2] = sum(kij == trimat[j,2])
						s2[3] = sum(kij == trimat[j,3])
						adjk[i,found[i]] = (kij[4:6])[s2 == 1]
						adjk[j,found[j]] = (kij[1:3])[s1 == 1]
					}
				}
			}
#### now determine which triangle contains each point
			xtri = 1:n*0;still = 1:n>0
			for (j in 1:ntri) {
				for (i in 1:n) {
					if (still[i]) {
						if (intri(xmat[i,],knots[trimat[j,1],],knots[trimat[j,2],],knots[trimat[j,3],])) {
							still[i] = FALSE
							xtri[i] = j
						}
					}
				}
			}
### make the design matrix
			dmat = matrix(0, nrow = n, ncol = capk)
			bmat = matrix(1, nrow = 3, ncol = 3)
			a = 1:3
			for (i in 1:n) {
				for (j in 1:3) {
					a[j] = trimat[xtri[i],j]
					bmat[j,2] = knots[a[j],1]
					bmat[j,3] = knots[a[j],2]
				}
				binv = solve(bmat)
				for (j in 1:3) {
					dmat[i,a[j]] = binv[1,j] + binv[2,j]*x1[i] + binv[3,j]*x2[i]
				}
			}
			#dmatc = cbind(dmat,zmat)
			Dmat = cbind(Dmat, dmat)
			vari = 1:capk*0 + ipr
			varlist = c(varlist, vari)
#print (head(dmatc))
### make the constraint matrix
			nconstr = sum(adjk > 0)
			amat = matrix(0, nrow = nconstr, ncol = capk)
			nr = 0
			for (itri in 1:ntri) {
				a = trimat[itri,]
				for (j in 1:3) {
					bmat[j,2] = knots[a[j],1]
					bmat[j,3] = knots[a[j],2]
				}
				binv = solve(bmat)
				for (i in 1:3) {
					k = adjk[itri,i]
					if(k > 0) {
						nr = nr + 1
						amat[nr,k] = 1
						for (j in 1:3) {
							amat[nr,a[j]] = -(binv[1,j] + binv[2,j]*knots[k,1] + binv[3,j]*knots[k,2])
						}
					}
				}
			}
			#m = dim(amat)[1]
			#a0 = matrix(0, nrow = m, ncol = p)
#new: s.conc.conc
			if (all(!cvs)) {
				amat = -amat
			} else if (any(cvs) & any(!cvs)) {
				stop ('Convex-concave or concave-convex fit is not implemented!')
			}
			amat_lst[[ipr]] = amat
			#amatc = cbind(amat, a0)
## make the penalty matrix
# find common edges & penalize change in slope
			pmat = matrix(0, nrow = 4 * m1 * m2,ncol = capk)
			six = 1:6;nu = 1:4;obs = 1:4
			nr = 0
			for (i in 1:(ntri-1)) {
				for (j in (i + 1):ntri) {
					six[1:3] = trimat[i,]
					six[4:6] = trimat[j,]
					sixu = unique(six)
					if (length(sixu) == 4) {
						nr = nr + 1
						for (ii in 1:4) {nu[ii] = sum(six == sixu[ii])}
						a = sixu[min(obs[nu == 2])]
						b = sixu[max(obs[nu == 2])]
						c = sixu[min(obs[nu == 1])]
						d = sixu[max(obs[nu == 1])]
						bmat[1,2] = knots[a,1]
						bmat[1,3] = knots[a,2]
						bmat[2,2] = knots[b,1]
						bmat[2,3] = knots[b,2]
						bmat[3,2] = knots[c,1]
						bmat[3,3] = knots[c,2]
						binv = solve(bmat)
						pmat[nr,a] = binv[1,1] + binv[2,1] * knots[d,1] + binv[3,1] * knots[d,2]
						pmat[nr,b] = binv[1,2] + binv[2,2] * knots[d,1] + binv[3,2] * knots[d,2]
						pmat[nr,c] = binv[1,3] + binv[2,3] * knots[d,1] + binv[3,3] * knots[d,2]
						pmat[nr,d] = -1
					}
				}
			}
			pmat = pmat[1:nr,]
			#p0 = matrix(0, nrow = nr, ncol = p)
			#Pmat = cbind(Pmat, pmat)
			#pmatc = cbind(pmat, p0)
			pmat_lst[[ipr]] = pmat
			#print (paste('nr: ', nr))
			knots_lst[[ipr]] = knots
			m12_lst[[ipr]] = c(m1,m2)
			trimat_lst[[ipr]] = trimat
			bmat_lst[[ipr]] = bmat
			capk_lst[[ipr]] = capk
		}
		amat = as.matrix(bdiag(amat_lst))
		m = dim(amat)[1]
		if (p > 0) {
			a0 = matrix(0, nrow = m, ncol = p)
		} else {a0 = NULL}
		amatc = cbind(amat, a0)
#one = 1:n*0 + 1
#pm = one %*% solve(crossprod(one), t(one))
#Dmat = Dmat - pm %*% Dmat
		dmatc = cbind(Dmat, zmat)
		pmat = as.matrix(bdiag(pmat_lst))
		if (p > 0) {
			p0 = matrix(0, nrow = nrow(pmat), ncol = p)
		} else {p0 = NULL}
		pmatc = cbind(pmat, p0)
#capk is for all columns of tri-surfaces pairs
		capk = sum(unlist(capk_lst))
	}
	#print (dim(amatc))
	nr = nrow(amatc); nc = ncol(amatc)
	pb = 0
	if (!is.null(delta_add)) {
		pb = nrow(delta_add)
		dmatc = cbind(t(delta_add), dmatc)
	}
	if (pb >= 1) {
		tmp = matrix(0, nrow = (nr+pb), ncol = (nc+pb))
		amatb = diag(pb)
#np_add is shape == 17
		if (np_add > 0) {
			amatb[1:np_add,1:np_add] = 0
		}
		tmp[1:pb,1:pb] = amatb
		tmp[(pb+1):(nr+pb), (pb+1):(nc+pb)] = amatc
		amatc = tmp
	}
	if (pb >= 1) {
		pmatb = matrix(0, nrow = dim(pmatc)[1], ncol = pb)
		pmatc = cbind(pmatb, pmatc)
	}
	if (!is.null(amat_wp)) {
		pwp = ncol(delta_wp)
		dmatc = cbind(dmatc, delta_wp)
		amatc = as.matrix(bdiag(amatc, amat_wp))
		pmatc = as.matrix(bdiag(pmatc, dmat_wp))
	} else {pwp = 0}
#new: always let lambda > 0
	if (round(lambda, 6) > 1e-6) {
		ps = lambda
	#} else if (pnt & (round(lambda, 6) == 0)) {
		#mat = cbind(1, x1, x2, x1*x2)
		#mat = cbind(1, x1t[,1], x2t[,1], x1t[,1]*x2t[,1])
		#if (nprs >= 2) {
		#	for (ipr in 2:nprs) {
		#		mat = cbind(mat, x1t[,ipr], x2t[,ipr], x1t[,ipr]*x2t[,ipr])
		#	}
		#}
		#mu_para = mat %*% solve(t(mat) %*% mat) %*% t(mat) %*% y
		#ssr = sum((y - mu_para)^2)
		#sc = ssr / (n - ncol(mat))
		#ps = max(1e-6, 1 * sc)
	} else {
	  ps = 10*1/n^(1/3)
	  #ps = 1e-6
	}
	lambda = ps
	#print (paste('pen:',lambda))
if (!wt.iter) {
	# weight
	yw = y
	dw = dmatc
	if (!is.null(w)) {
		if (min(w) > 1e-8) {
			yw = y * sqrt(w)
			for (i in 1:n) {
				dw[i, ] = dmatc[i, ] * sqrt(w[i])
			}
		} else {print ("check the user-defined weights!")}
	}
	#print (dim(dmatc))
	umatc = chol(t(dw) %*% dw + lambda * t(pmatc) %*% pmatc)
	uinv = solve(umatc)
	atil = amatc %*% uinv
	zvec = t(uinv) %*% t(dw) %*% yw
	ans = coneA(zvec, atil)
	thhat = uinv %*% ans$thetahat
	muhat = dmatc %*% thhat
	etahat = muhat
} else {
	etahat = etahat.fun(n, y, fml = family$family)
	gr = gr.fun(y, etahat, weights = w, fml = family$family)
	wt = wt.fun(y, etahat, n, weights = w, fml = family$family)
	cvec = crossprod(dmatc, (wt * etahat - gr))
	umatc = chol(t(dmatc) %*% diag(wt) %*% dmatc + lambda * t(pmatc) %*% pmatc)
	uinv = solve(umatc)
	atil = amatc %*% uinv
	zvec = t(uinv) %*% t(dmatc) %*% y
	ans = coneA(zvec, atil)
	thhat = uinv %*% ans$thetahat
	etahat = dmatc %*% thhat
	muhat = muhat.fun(etahat, fml = family$family)
	diff = 1
	nrep = 0
	sm = 1e-6
	while (diff > sm & nrep < 100) {
		oldmu = muhat
		nrep = nrep + 1
		gr = gr.fun(y, etahat, weights = w, fml = family$family)
		wt = wt.fun(y, etahat, n, weights = w, fml = family$family)
		cvec = crossprod(dmatc, (wt * etahat - gr))
		umatc = chol(t(dmatc) %*% diag(wt) %*% dmatc + lambda * t(pmatc) %*% pmatc)
		uinv = solve(umatc)
		atil = amatc %*% uinv
		zvec = t(uinv) %*% t(dmatc) %*% y
		ans = coneA(zvec, atil)
		thhat = uinv %*% ans$thetahat
		etahat = dmatc %*% thhat
		muhat = muhat.fun(etahat, fml = family$family)
		diff = mean((muhat - oldmu)^2)
    }
    #llh = llh.fun(y, muhat, etahat, n, w, fml = family$family)
}
	llh = llh.fun(y, muhat, etahat, phihat=NULL, n, w, fml = family$family)
	coefkeep = ans$thetahat
	muhatkeep = muhat
	etahatkeep = etahat
	thhatkeep = thhat
	if (pb > 0) {
		#print (dim(delta_add))
		#print (varlist_add)
		coef_add = thhat[1:pb]
		capl = ncol(xmat_add)
		thvecs = matrix(0, nrow = capl, ncol = n)
	    ncon = 0
	    #vcoef_add = coef_add[1:capl]
	    vcoef_add = coef_add[1:np_add]
	    lconv = sum(shapes_add > 2 & shapes_add < 5 | shapes_add > 10 & shapes_add < 13)
	    if (lconv > 0) {
	    	dcoefs = coef_add[-c(1:lconv)]
	    	delta_add2 = delta_add[-c(1:lconv), ,drop = FALSE]
	    } else {
	    	dcoefs = coef_add
	    	delta_add2 = delta_add
	    }
	    for (i in 1:capl) {
	    	thvecs[i,] = t(delta_add2[varlist_add == i,]) %*% dcoefs[varlist_add == i]
			if (shapes_add[i] > 2 & shapes_add[i] < 5 | shapes_add[i] > 10 & shapes_add[i] < 13) {
            	ncon = ncon + 1
				thvecs[i,] = thvecs[i,] + vcoef_add[ncon] * xmat_add[,i]
            }
	    }
	} else {coef_add = 0; thvecs = NULL}
	if (p > 0) {
		zcoefs = thhat[(pb+capk+1):(pb+capk+p)]
	} else {zcoefs = 0}
	if (pwp > 0) {
		nth = length(thhat)
		coef_wp = thhat[(nth-pwp+1):nth]
	} else {coef_wp = 0}
#vcoefs include zcoefs and shape = 17
	vcoefs = zcoefs
	#vmat = dmatc[,(pb+capk+1):(pb+capk+p), drop=F]
	#vmat = cbind(1:n*0+1, zmat)
	vmat = zmat
	if (np_add > 0) {
		vcoefs = c(thhat[1:np_add], vcoefs)
		vmat = cbind(dmatc[, 1:np_add, drop=F], vmat)
	}
	#pv = vmat %*% solve(crossprod(vmat), t(vmat))
	if (is.null(vmat)) {
		muvhat = 0
	} else {
		muvhat = muhat.fun(vmat %*% vcoefs, fml = family$family)
	}
	if (is.null(w)) {
		w = rep(1, n)
	}
	sse1 = sum(w*(y - muhat)^2)
	sse0 = sum(w*(y - muvhat)^2)
#new use edf instead if (n - cpar * edf) < 0
if (!is.null(vmat)) {
	np = ncol(vmat)
} else {np = 0}
### find the GCV
	atil0 = atil[round(atil %*% coefkeep, 8) == 0,]
	ans1 = qr(t(atil0))
	ared = qr.Q(ans1)[,1:ans1$rank]
	prmat = ared %*% solve(t(ared) %*% ared) %*% t(ared)
	#imat = matrix(0, nrow = capk+p, ncol = capk+p)
	#for (i in 1:(capk+p)) {imat[i,i] = 1}
	#diag(imat) = 1
	imat = diag(ncol(dmatc))
	bigpr = dmatc %*% uinv %*% (imat - prmat) %*% t(uinv) %*% t(dmatc)
	edf = sum(diag(bigpr))
	gcv = sse1 / (1 - edf/n)^2
	#if ((n - np - cpar * edf) <= 0) {
	if ((n - cpar * edf) <= 0) {
		sig2hat = sse1 / edf
	} else {
		sig2hat = sse1 / (n - cpar * edf)
	}
	#print (sig2hat)
	#if (p >= 1) {
	#	#inmat = matrix(0, nrow = n, ncol = n)
	#	#diag(inmat) = 1
	#	inmat = diag(n)
	#	one = 1:n*0+1
	#	pm = one %*% solve(crossprod(one)) %*% t(one)
		#print (dim(zmat))
		#print (dim(bigpr))
	#	covmat = solve(t(zmat) %*% (inmat - bigpr) %*% zmat)
	#	print (covmat)
	#	print (sig2hat)
	#	sez = sqrt(diag(covmat) * sig2hat)
	#	tz = zcoef / sez
#new use edf instead if (n - cpar * edf) < 0
	#	if ((n - p - cpar * edf) <= 0) {
	#		pz = 2 * (1 - pt(abs(tz), edf))
	#		if (p > 1) {
	#			warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
	#		}
#print ('Check pz!')
	#	} else {
	#		pz = 2 * (1 - pt(abs(tz), n - p - cpar * edf))
	#	}
	#} else {pz = NULL; sez = NULL}
	pz = NULL; sez = NULL
	#new: get cic
	cic = NULL
	if (is.null(w)) {
		w = rep(1, n)
	}
	dw = dmatc
	for (i in 1:n) {
		dw[i, ] = dmatc[i, ] * sqrt(w[i])
	}
	if (nsim > 0) {
		edfs = 1:nsim*0
		if (!wt.iter) {
			for (isim in 1:nsim) {
				ysim = rnorm(n)
				ysim = ysim * sqrt(w)
				zvec = t(uinv) %*% t(dw) %*% ysim
				ansi = coneA(zvec, atil, msg=FALSE)
				thhati = uinv %*% ansi$thetahat
				muhati = dmatc %*% thhati
				etahati = muhati
				edfi = tri_getedf(cf = ansi$thetahat, atil, dmatc, uinv, ncol(dmatc))
				edfs[isim] = edfi
			}
			#cic = log(sse1) + log(2 * (mean(edfs) + p) / (n - p - 1.5 * mean(edfs)) + 1)
			#cic = llh + log(2 * (mean(edfs) + p) / (n - p - 1.5 * mean(edfs)) + 1)
			cic = llh + log(2 * (mean(edfs)) / (n - np - 1.5 * (mean(edfs) - np)) + 1)
		} else {
			if (family$family == "poisson") {
				mu0 = mean(y)
			} else {mu0 = NULL}
			for (isim in 1:nsim) {
				ysim = ysim.fun(n, mu0, fml = family$family)
				etahat = etahat.fun(n, ysim, fml = family$family)
				gr = gr.fun(ysim, etahat, weights = w, fml = family$family)
				wt = wt.fun(ysim, etahat, n, weights = w, fml = family$family)
				cvec = crossprod(dmatc, (wt * etahat - gr))
				umatc = chol(t(dmatc) %*% diag(wt) %*% dmatc + lambda * t(pmatc) %*% pmatc)
				uinv = solve(umatc)
				atil = amatc %*% uinv
				zvec = t(uinv) %*% t(dmatc) %*% ysim
				ansi = coneA(zvec, atil, msg=FALSE)
				thhat = uinv %*% ansi$thetahat
				etahat = dmatc %*% thhat
				muhat = muhat.fun(etahat, fml = family$family)
				diff = 1
				nrep = 0
				while (diff > sm & nrep < 100) {
					oldmu = muhat
					nrep = nrep + 1
					gr = gr.fun(ysim, etahat, weights = w, fml = family$family)
					wt = wt.fun(ysim, etahat, n, weights = w, fml = family$family)
					cvec = crossprod(dmatc, (wt * etahat - gr))
					umatc = chol(t(dmatc) %*% diag(wt) %*% dmatc + lambda * t(pmatc) %*% pmatc)
					uinv = solve(umatc)
					atil = amatc %*% uinv
					zvec = t(uinv) %*% t(dmatc) %*% ysim
					ansi = coneA(zvec, atil)
					thhat = uinv %*% ansi$thetahat
					etahat = dmatc %*% thhat
					muhat = muhat.fun(etahat, fml = family$family)
					diff = mean((muhat - oldmu)^2)
    			}
				edfi = tri_getedf(cf = ansi$thetahat, atil, dmatc, uinv, ncol(dmatc))
			 	edfs[isim] = edfi
			}
			#cic = llh + log(2 * (mean(edfs) + p) / (n - p - 1.5 * mean(edfs)) + 1)
			cic = llh + log(2 * (mean(edfs)) / (n - np - 1.5 * (mean(edfs) - np)) + 1)
		}
	}
	ans = new.env()
	ans$family = family
	ans$edf = edf
	if (nsim > 0) {
		ans$edf0 = mean(edfs) #+ p
	} else {ans$edf0 = p}
	ans$zcoefs = zcoefs
	ans$pvals.beta = pz
	ans$se.beta = sez
	ans$muhat = muhatkeep
	ans$etahat = etahatkeep
	ans$gcv = gcv
	#ans$knots = knots
	ans$thhat = thhatkeep
	ans$coef_tri = thhatkeep[(pb+1):(pb+capk)]
	ans$trimat = trimat
	ans$x1 = x1
	ans$x2 = x2
	ans$d0 = p
	ans$zmat = zmat
	ans$capk = capk
	ans$cic = cic
	ans$sse1 = sse1
	#ans$sse0 = sse0
	ans$varlist = varlist
	ans$coef_add = coef_add
	ans$coef_wp = coef_wp
	ans$etacomps = thvecs
	ans$pen = lambda
	ans$knots_lst = knots_lst
	#new: in trispl, knots' length is not m1 or m2, it's much larger
	ans$m12_lst = m12_lst
	ans$trimat_lst = trimat_lst
	ans$bmat_lst = bmat_lst
	ans$capk_lst = capk_lst

    #new
    ans$dmatc = dmatc
    ans$pmatc = pmatc
    ans$amatc = amatc
    ans$sig2hat = sig2hat
	ans
}

tri_getedf = function(cf = NULL, atil, dmatc, uinv, nc) {
	atil0 = atil[round(atil %*% cf, 8) == 0,]
	ans1 = qr(t(atil0))
	ared = qr.Q(ans1)[,1:ans1$rank]
	prmat = ared %*% solve(t(ared) %*% ared) %*% t(ared)
	imat = diag(nc)
	bigpr = dmatc %*% uinv %*% (imat - prmat) %*% t(uinv) %*% t(dmatc)
	edf = sum(diag(bigpr))
	return (edf)
}



###
#trispl makedelta: used in predict.trispl
###
makedelta_tri = function(x1, x2, m1 = 0, m2 = 0, k1 = NULL, k2 = NULL, trimat = NULL, capk = NULL, space = c("E",
"E"), cvs = c(TRUE, TRUE), interp = TRUE) {
    n = length(x1)
    xmat = cbind(x1, x2)
    delta = NULL
    #### now determine which triangle contains each point
    xtri = 1:n*0;still = 1:n>0
    ntri = (m2 - 1) * (2 * m1 - 1)
    knots = cbind(k1, k2)
    for (j in 1:ntri) {
        for (i in 1:n) {
            if (still[i]) {
                if (intri(xmat[i,],knots[trimat[j,1],],knots[trimat[j,2],],knots[trimat[j,3],])) {
                    still[i] = FALSE
                    xtri[i] = j
                }
            }
        }
    }
    ### make the design matrix
    dmat = matrix(0, nrow = n, ncol = capk)
    bmat = matrix(1, nrow = 3, ncol = 3)
    a = 1:3
    for (i in 1:n) {
        for (j in 1:3) {
            a[j] = trimat[xtri[i],j]
            bmat[j,2] = knots[a[j],1]
            bmat[j,3] = knots[a[j],2]
        }
        binv = solve(bmat)
        for (j in 1:3) {
            dmat[i,a[j]] = binv[1,j] + binv[2,j]*x1[i] + binv[3,j]*x2[i]
        }
    }
    delta = cbind(delta, dmat)
    return (delta)
}


##
s.conv.conv <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "tri_cvs"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "cvs") <- c(TRUE, TRUE)
    attr(xm, "categ") <- "tri"
    #warp <<- TRUE
    #class(xm) <- "tri"
    return (xm)
}

s.conc.conc <- function(x1, x2, numknots = c(0, 0), knots = list(k1 = 0, k2 = 0), space = c("E", "E"))
{
    cl <- match.call()
    pars1 <- match.call()[2]
    pars2 <- match.call()[3]
    xm <- cbind(x1, x2)
    attr(xm, "name") <- c(deparse(pars1$x1), deparse(pars2$x2))
    attr(xm, "shape") <- "tri_ccs"
    attr(xm, "numknots") <- numknots
    attr(xm, "knots") <- knots
    attr(xm, "space") <- space
    attr(xm, "cvs") <- c(FALSE, FALSE)
    attr(xm, "categ") <- "tri"
    #warp <<- TRUE
    #class(xm) <- "tri"
    return (xm)
}
################
#summary.trispl#
################
#summary.trispl <- function(object,...) {
#	if (!is.null(object$zcoefs)) {
#		coefs <- object$zcoefs
#		se <- object$se.beta
#		#tval <- object$tz
#		pvalbeta <- object$pvals.beta
#		tval <- coefs / se
#		n <- length(coefs)
#		#sse0 <- object$SSE0
#		#sse1 <- object$SSE1
#		zid <- object$zid
#new: zid1, zid2 just index zmat not bigmat
#		zid1 <- object$zid1
#		zid2 <- object$zid2
#		#tms <- object$tms
#		#zmat <- object$zmat
#		#is_mat <- object$is_mat
#		is_param <- object$is_param
#		is_fac <- object$is_fac
#		vals <- object$vals
#		tms <- object$tms
#new:
#		cic <- object$cic
#		rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "t.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
#		rownames(rslt1)[1] <- "(Intercept)"
#		if (n > 1) {
#			lzid <- length(zid1)
#			for (i in 1:lzid) {
#				pos1 <- zid1[i]; pos2 <- zid2[i]
#				for (j in pos1:pos2) {
#					if (!is_param[i]) {
#						rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")
#					} else {
#						rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")
#					}
#				}
#			}
#		}
#		rslt1 <- as.matrix(rslt1)
#		#if (!is.null(sse0) & !is.null(sse1)) {
#		#	rslt2 <- data.frame("SSE.Linear" = sse0, "SSE.Full" = sse1)
#new:
#		#	rownames(rslt2)[1] <- ""
#			#rslt2 <- as.matrix(rslt2)
#		#	ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs)
#		#	class(ans) <- "summary.wps"
#		#	ans
#		#} else {
#			ans <- list(call = object$call, coefficients = rslt1, zcoefs = coefs, cic = cic)
#			class(ans) <- "summary.trispl"
#			ans
#		#}
#	} else {
#		ans <- list(zcoefs = object$zcoefs)
#		class(ans) <- "summary.trispl"
#		ans
#	}
#}


#######################
#print.summary.trispl #
#######################
#print.summary.trispl <- function(x,...) {
#	if (!is.null(x$zcoefs)) {
#	#if (!is.null(x$se.beta)) {
#		cat("Call:\n")
#		print(x$call)
#		cat("\n")
#		cat("Coefficients:")
#		cat("\n")
#		printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
#		cat("\n")
#		if (!is.null(x$cic)) {
#			cat("CIC: ", round(x$cic,4), "\n", sep = "")
#		}
#	} else {
#		print ("No linear predictor is defined")
#	}
#}

#################
#plotpersp.tri#
################
plotpersp.trispl = function(object, x1 = NULL, x2 = NULL, x1nm = NULL, x2nm = NULL,
                            data = NULL, surface = "C", categ = NULL, col = NULL,
                            random = FALSE, ngrid = 12, xlim = range(x1),
                            ylim = range(x2), zlim = NULL, xlab = NULL,
                            ylab = NULL, zlab = NULL, th = NULL, ltheta = NULL,
                            main = NULL, ticktype = "simple",...) {
    if (!inherits(object, "trispl")) {
        warning("calling plotpersp(<fake-trisp-object>) ...")
    }
    #x1nm = deparse(substitute(x1))
    #x2nm = deparse(substitute(x2))
    xnms = object$xnms_tri
    xmat = object$xmat_tri
    #if (x1nm == "NULL" | x2nm == "NULL") {
    if (is.null(x1nm) | is.null(x2nm)) {
        if (length(xnms) >= 2) {
            x1nm = xnms[1]
            x2nm = xnms[2]
            x1id = 1
            x2id = 2
            x1 = xmat[, 1]
            x2 = xmat[, 2]
        } else {stop ("Number of non-parametric predictors must >= 2!")}
    }
    #x1nm0 = xnms[1]
    #x2nm0 = xnms[2]
    #x1 = xmat[, 1]
    #x2 = xmat[, 2]
    is_fac = object$is_fac
    ynm = object$ynm
    #xmat is delta
    #delta = object$delta
    znms = object$znms
    kznms = length(znms)
    #zmat not include 1 vector
    zmat = object$zmat
    #zmat0 = zmat[, -1, drop = FALSE]
    zmat0 = zmat
    #zcoefs = object$zcoefs[-1]
    zcoefs = object$zcoefs
    zid1 = object$zid1
    zid2 = object$zid2
    p = object$d0
    cvss = object$cvss
    labels = object$labels
    labels = labels[which(grepl("tri", labels, fixed = TRUE))]
    varlist = object$varlist_tri
    #varlist = varlist[-1]
    #knots = object$kts
    #trimat = object$trimat
    family = object$family
    fml = family$family
    #print (family)
    cicfamily = CicFamily(family)
    muhat.fun = cicfamily$muhat.fun
    #if (!is.null(categ)) {
    #	if (!is.character(categ)) {
    #		warning("categ must be a character argument!")
    #	} else if (!any(znms == categ)) {
    #print ('TRUE')
    #		warning(paste(categ, "is not an exact character name defined in the cgam fit!"))
    #		categ = NULL
    #	} else {
    #		obsz = 1:kznms
    #		zid = obsz[znms == categ]
    #		if (!(is_fac[zid])) {
    #			categ = NULL
    #		}
    #	}
    #}
    if (!is.null(categ)) {
        if (!is.character(categ)) {
            warning("categ must be a character argument!")
        } else if (!any(znms == categ)) {
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
    #new: switch xnms
    if (!is.null(data)) {
        if (!is.data.frame(data)) {
            stop ("User need to make the data argument a data frame with names for each variable!")
        }
        datnms = names(data)
        if (!any(datnms == x1nm) | !any(datnms == x2nm)) {
            stop ("Check the accuracy of the names of x1 and x2!")
        }
        x1 = data[ ,which(datnms == x1nm)]
        x2 = data[ ,which(datnms == x2nm)]
    } else {
        if (all(xnms != x1nm)) {
            #stop (paste(paste("'", x1nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
            #new: in case of wrong data fame
            if (length(x1) != nrow(xmat)) {
                stop ("Number of observations in the data set is not the same as the number of elements in x1!")
            }
            bool = apply(xmat, 2, function(x) all(x1 == x))
            if (any(bool)) {
                #x1id = obs[bool]
                x1nm = xnms[bool]
            } else {
                stop (paste(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the wps fit!"))
            }
        }
        if (all(xnms != x2nm)) {
            #stop (paste(paste("'", x2nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
            if (length(x2) != nrow(xmat)) {
                stop ("Number of observations in the data set is not the same as the number of elements in x2!")
            }
            bool = apply(xmat, 2, function(x) all(x2 == x))
            if (any(bool)) {
                #x2id = obs[bool]
                x2nm = xnms[bool]
            } else {
                stop (paste(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the wps fit!"))
            }
        }
    }
    xnm12 = c(x1nm, x2nm)
    #print (xnms)
    #print(xnm12)
    id_lab = which(xnms %in% xnm12)
    #print (labels)
    #print (xnm12)
    xnm12_lab = labels[id_lab]
    id1 = id2 = ipr = NULL
    if (length(unique(xnm12_lab)) > 1 | length(id_lab) != 2) {
        stop ("Two non-parametric predictors do not form a triangle-spline surface!")
    } else {
        id1 = sort(id_lab)[1]
        id2 = sort(id_lab)[2]
        ipr = id2 / 2
    }
    #decrs = decrs0[[ipr]]
    #ktsi = kts[[ipr]]
    #k1 = ktsi[[1]]
    #k2 = ktsi[[2]]
    knots_lst = object$knots_lst
    trimat_lst = object$trimat_lst
    bmat_lst = object$bmat_lst
    thhat_all = object$coef_tri
    capk_lst = object$capk_lst

    knots = knots_lst[[ipr]]
    trimat = trimat_lst[[ipr]]
    bmat = bmat_lst[[ipr]]
    capk = capk_lst[[ipr]]
    cvs = cvss[[ipr]]
    thhat = thhat_all[which(varlist == ipr)]
    if (p > 0) {
        thhat = c(thhat, zcoefs)
    }
    #additive
    thvecs = object$etacomps
    xnms_add = object$xnms_add
    xmat_add = object$xmat_add
    knms = length(xnms_add)
    x3_add = 0
    if (knms >= 1) {
        #x3id <- obs[-c(x1id, x2id)]
        #kx3 <- length(x3id)
        for (i in 1:knms) {
            x3i = xmat_add[, i]
            x3i_use = max(x3i[x3i <= median(x3i)])
            x3i_add = min(thvecs[i, x3i == x3i_use])
            x3_add = x3_add + x3i_add
        }
    }
    ntri = nrow(trimat)
    if (x1nm != xnms[1] & x2nm != xnms[2]) {
        nm = x1nm
        x1nm = x2nm
        x2nm = nm
        tmp = x1
        x1 = x2
        x2 = tmp
    } else {x1nm = x1nm; x2nm = x2nm}

    if (is.null(th) | !is.numeric(th)) {
        ang = NULL
        if (cvs[1] & cvs[2]) {
            if (is.null(ang)) {
                ang = 140
            }
        } else if (!cvs[1] & !cvs[2]) {
            if (is.null(ang)) {
                ang = 25
            }
        }
    } else {ang = th}
    if (is.null(ltheta) | !is.numeric(ltheta)) {
        ltheta <- -135
    }
    #make the surface
    g = ngrid
    ng = g^2
    x1g = 0:(g-1)/(g-1) * (max(x1)-min(x1)) + min(x1)
    x2g = 0:(g-1)/(g-1) * (max(x2)-min(x2)) + min(x2)
    gg = matrix(0, nrow = ng, ncol = 2)
    for (i in 1:g) {
        gg[((i-1) * g + 1):(i * g),1] = x1g
        gg[((i-1) * g + 1):(i * g),2] = 1:g * 0 + x2g[i]
    }
    gtri = 1:ng * 0
    still = 1:ng>0
    for (j in 1:ntri) {
        for (i in 1:ng) {
            if (still[i]) {
                if (intri(gg[i,],knots[trimat[j,1],],knots[trimat[j,2],],knots[trimat[j,3],])) {
                    still[i] = FALSE
                    gtri[i] = j
                }
            }
        }
    }
    gmat = matrix(0,nrow = ng,ncol = capk)
    #print (gmat)
    bmat = matrix(1,nrow = 3,ncol = 3)
    a = 1:3
    for (i in 1:ng) {
        for (j in 1:3) {
            #if (gtri[i] > 0) {
            a[j] = trimat[gtri[i],j]
            bmat[j,2] = knots[a[j],1]
            bmat[j,3] = knots[a[j],2]
            #}
        }
        binv = solve(bmat)
        for (j in 1:3) {
            gmat[i,a[j]] = binv[1,j] + binv[2,j] * gg[i,1] + binv[3,j] * gg[i,2]
        }
    }
    #?
    #if (p > 0) {
    g0 = matrix(0,nrow = dim(gmat)[1],ncol = p)
    #} else {g0 = NULL}
    #print (g0)
    #print (dim(gmat))
    #print (dim(thhat))
    gvals = cbind(gmat, g0) %*% thhat
    #print (gvals)
    mupl = matrix(gvals, nrow = g)
    m1 = g
    m2 = ncol(mupl)
    if (fml != "gaussian") {
        for (i1 in 1:m1) {
            for (i2 in 1:m2) {
                mupl[i1, i2] = muhat.fun(mupl[i1, i2], fml = fml)
            }
        }
        #mupl = muhat.fun(mupl, fml = fml)
    }
    if (is.null(categ)) {
        z_add = 0
        if (!is.null(znms)) {
            #if (p > 0) {
            #print ('skip')
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
        mupl = mupl + as.numeric(z_add) + as.numeric(x3_add)
        #print (mupl)
        m1 = nrow(mupl)
        m2 = ncol(mupl)
        if (fml != "gaussian") {
            for (i1 in 1:m1) {
                for (i2 in 1:m2) {
                    mupl[i1, i2] = muhat.fun(mupl[i1, i2], fml = fml)
                }
            }
            #mupl = muhat.fun(mupl, fml = fml)
        }
        mins = min(mupl); maxs = max(mupl)
    } else {
        mupls = list()
        mins = maxs = NULL
        obsz = 1:kznms
        zid = obsz[znms == categ]
        pos1 = zid1[zid]; pos2 = zid2[zid]
        zcoefsi = zcoefs[pos1:pos2]
        #include the base level
        zcoefsi = c(0, zcoefsi)
        z_add = sort(zcoefsi)
        kz_add = length(z_add)
        for (iz in 1:kz_add) {
            mupls[[iz]] = mupl + z_add[iz] + as.numeric(x3_add)
            mins = c(mins, min(mupls[[iz]]))
            maxs = c(maxs, max(mupls[[iz]]))
        }
        m1 = nrow(mupl)
        m2 = ncol(mupl)
        if (fml != "gaussian") {
            for (iz in 1:kz_add) {
                mupli = mupls[[iz]]
                for (i1 in 1:m1) {
                    for (i2 in 1:m2) {
                        mupli[i1, i2] = muhat.fun(mupli[i1, i2], fml = fml)
                    }
                }
                #mupli = muhat.fun(mupli, fml = fml)
                mupls[[iz]] = mupli
            }
        }
    }
    palette = c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
    if (is.null(xlab)) {
        xlab = x1nm
    }
    if (is.null(ylab)) {
        ylab = x2nm
    }
    if (is.null(zlab)) {
        if (fml == "binomial") {
            zlab = paste("Pr(", ynm, ")")
        } else if (fml == "poisson" | fml == "gaussian" | fml == "Gamma") {
            zlab = paste("Est mean of", ynm)
        }
    }
    if (is.null(categ)) {
        if (is.null(col)) {
            if (random) {
                col = sample(palette, size = 1, replace = FALSE)
            } else {
                #col = "white"
                musurf = mupl
                nr = nrow(musurf)
                nc = ncol(musurf)
                ncol = 100
                facet = musurf[-1,-1] + musurf[-1,-nc] + musurf[-nr,-1] + musurf[-nr,-nc]
                #print (facet)
                facetcol = cut(facet, ncol)
                col = heat.colors(ncol)[facetcol]
            }
        } else {
            if (col == "heat" | col == "topo" | col == "terrain" | col == "cm") {
                nr = nrow(mupl)
                nc = ncol(mupl)
                ncol = 100
                facet = mupl[-1,-1] + mupl[-1,-nc] + mupl[-nr,-1] + mupl[-nr,-nc]
                facetcol = cut(facet, ncol)
                if (col == "heat") {
                    col = heat.colors(ncol)[facetcol]
                } else if (col == "topo") {
                    col = topo.colors(ncol)[facetcol]
                } else if (col == "terrain") {
                    col = terrain.colors(ncol)[facetcol]
                } else {
                    col = cm.colors(ncol)[facetcol]
                }
            }
            if (random) {
                print ("User defined color is used!")
            }
        }
        #if (surface == 'C') {
        musurf = mupl
        #if (is.null(main)) {
        #	main = 'Constrained Warped-Plane Spline Surface'
        #}
        if (is.null(zlim)) {
            lwr = min(mins)
            upp = max(maxs)
            zlim0 = c(lwr - (upp-lwr)/5, upp + (upp-lwr)/5)
        } else {
            zlim0 = zlim
        }
        #}
        #persp(k1, k2, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype,...)
        persp(x1g, x2g, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype,...)
        res = list(musurf = mupl, gg = gg, x1g=x1g, x2g=x2g, z_add = z_add, x3_add = x3_add, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype)
        invisible(res)
    } else {
        kxgm = length(mupls)
        if (is.null(col)) {
            if (random) {
                #new:
                col = topo.colors(kxgm)
            } else {
                if (kxgm > 1 & kxgm < 11) {
                    col = palette[1:kxgm]
                } else {
                    #new:
                    col = topo.colors(kxgm)
                }
            }
        } else {
            col0 = col
            if (col0 == "heat" | col0 == "topo" | col0 == "terrain" | col0 == "cm") {
                #col0 <- col
                ncol = 100
                facets = facetcols = list()
                col = list()
                for (i in 1:kxgm) {
                    nr = nrow(mupls[[i]])
                    nc = ncol(mupls[[i]])
                    facets[[i]] = (mupls[[i]])[-1,-1] + (mupls[[i]])[-1,-nc] + (mupls[[i]])[-nr,-1] + (mupls[[i]])[-nr,-nc]
                    facetcols[[i]] = cut(facets[[i]], ncol)
                    #print (head(facetcols[[i]]))
                    if (col0 == "heat") {
                        col[[i]] = (heat.colors(ncol))[facetcols[[i]]]
                        #print (head(col[[i]]))
                    } else if (col0 == "topo") {
                        col[[i]] = (topo.colors(ncol))[facetcols[[i]]]
                    } else if (col0 == "terrain") {
                        col[[i]] = (terrain.colors(ncol))[facetcols[[i]]]
                    } else {
                        col[[i]] = (cm.colors(ncol))[facetcols[[i]]]
                    }
                }
            } else if (length(col0) < kxgm) {
                #new:
                col = topo.colors(kxgm)
            } else if (length(col0) > kxgm) {
                col = col0[1:kxgm]
            }
        }
        for (i in 1:kxgm) {
            mupli = mupls[[i]]
            #if (surface == 'C') {
            musurf = mupli
            #if (is.null(main)) {
            #	main = 'Constrained Warped-Plane Spline Surface'
            #}
            if (is.null(zlim)) {
                lwr = min(mins)
                upp = max(maxs)
                zlim0 = c(lwr - (upp-lwr)/5, upp + (upp-lwr)/5)
            } else {
                zlim0 = zlim
            }
            #}
            if (is.list(col)) {
                coli = unlist(col[[i]])
                #print (head(coli))
            } else {coli = col[i]}
            #persp(k1, k2, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = col[i], cex.axis = .75, main = main, ticktype = ticktype,...)
            persp(x1g, x2g, musurf, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = x1nm, ylab = x2nm, zlab = ynm, theta = ang, ltheta = ltheta, col = coli, cex.axis = .75, main = main, ticktype = ticktype,...)
            par(new = TRUE)
        }
        par(new = FALSE)
    }
}


##########################################
#apply plotpersp on a trispl.predict object
#not done: >= 1 tri pair + additive + z, ignore wps
#one pair + z only
##########################################
plotpersp.trisplp = function(object, x1=NULL, x2=NULL, x1nm=NULL, x2nm=NULL, data=NULL, up = TRUE, main=NULL, cex.main=.8, xlab = NULL, ylab = NULL, zlab = NULL, th = NULL, ltheta = NULL, ticktype = "simple",...) {
    #obj is prediction for trispl
    if (!inherits(object, "trisplp")) {
        warning("calling plotpersp(<fake-trisplp-object>) ...")
    }
    t_col = function(color, percent = 50, name = NULL) {
        rgb.val <- col2rgb(color)
        ## Make new color using input color as base and alpha set by transparency
        t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
        maxColorValue = 255,
        alpha = (100-percent)*255/100,
        names = name)
        ## Save the color
        invisible(t.col)
    }
    if (up) {
        mycol = t_col("green", perc = 90, name = "lt.green")
    } else {
        mycol = t_col("pink", perc = 80, name = "lt.pink")
    }

    acov = object$acov
    mult = object$mult
    #obj is the wps fit
    obj = object$object

    xnms = obj$xnms_tri
    xmat = obj$xmat_tri
    xnms_add = obj$xnms_add
    xmat_add = obj$xmat_add
    delta = obj$dmatc
    #if (x1nm == "NULL" | x2nm == "NULL") {
    if (is.null(x1nm) | is.null(x2nm)) {
        if (length(xnms) >= 2) {
            x1nm = xnms[1]
            x2nm = xnms[2]
            x1id = 1
            x2id = 2
            x1 = xmat[, 1]
            x2 = xmat[, 2]
        } else {stop ("Number of non-parametric predictors must >= 2!")}
    }

    labels = obj$labels
    labels = labels[which(grepl("tri", labels, fixed = TRUE))]

    is_fac = obj$is_fac
    ynm = obj$ynm

    varlist = obj$varlist_tri

    kts = obj$knots_lst
    np = obj$d0
    #switch xnms
    if (!is.null(data)) {
        if (!is.data.frame(data)) {
            stop ("User need to make the data argument a data frame with names for each variable!")
        }
        datnms = names(data)
        if (!any(datnms == x1nm) | !any(datnms == x2nm)) {
            stop ("Check the accuracy of the names of x1 and x2!")
        }
        x1 = data[ ,which(datnms == x1nm)]
        x2 = data[ ,which(datnms == x2nm)]
    } else {
        if (all(xnms != x1nm)) {
            #stop (paste(paste("'", x1nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
            #new: in case of wrong data fame
            if (length(x1) != nrow(xmat)) {
                stop ("Number of observations in the data set is not the same as the number of elements in x1!")
            }
            bool = apply(xmat, 2, function(x) all(x1 == x))
            if (any(bool)) {
                #x1id = obs[bool]
                x1nm = xnms[bool]
            } else {
                stop (paste(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the wps fit!"))
            }
        }
        if (all(xnms != x2nm)) {
            #stop (paste(paste("'", x2nm0, "'", sep = ''), "is not an exact predictor name defined in the cgam fit!"))
            if (length(x2) != nrow(xmat)) {
                stop ("Number of observations in the data set is not the same as the number of elements in x2!")
            }
            bool = apply(xmat, 2, function(x) all(x2 == x))
            if (any(bool)) {
                #x2id = obs[bool]
                x2nm = xnms[bool]
            } else {
                stop (paste(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the wps fit!"))
            }
        }
    }
    xnm12 = c(x1nm, x2nm)
    id_lab = which(xnms %in% xnm12)
    xnm12_lab = labels[id_lab]
    #new:
    xnm_other = xnms[-id_lab]
    id1 = id2 = ipr = NULL

    if (length(unique(xnm12_lab)) > 1 | length(id_lab) != 2) {
        stop ("Two non-parametric predictors do not form a triangle-spline surface!")
    } else {
        id1 = sort(id_lab)[1]
        id2 = sort(id_lab)[2]
        ipr = id2 / 2
    }

    #find the pairs to be plotted
    if (x1nm != xnms[1] & x2nm != xnms[2]) {
        nm = x1nm
        x1nm = x2nm
        x2nm = nm
        tmp = x1
        x1 = x2
        x2 = tmp
        #new:
        xnm12 = c(x1nm, x2nm)
    } else {x1nm = x1nm; x2nm = x2nm}

    newz = object$newz
    npr = round(length(xnms) / 2, 0L)

    #zlwr = min(object$lower) - (max(object$fit) - min(object$fit)) / 4
    #zupp = max(object$upper) + (max(object$fit) - min(object$fit)) / 4

    res = plotpersp.trispl(obj, x1, x2, ngrid = 7, col='white', xlab=xlab, ylab=ylab, zlab=zlab, th=th, ltheta=ltheta, ticktype=ticktype)

    newData = as.data.frame(res$gg)

    #if (!is.null(newz)) {
    #    znms = obj$znms
    #    nz = matrix(0, nrow=nrow(gg), ncol=ncol(newz))
    #    nznm = gsub("[\\(\\)]", "", regmatches(znms, gregexpr("\\(.*?\\)", znms))[[1]])
    #    gg = cbind(gg, nz)
    #    colnames(gg) = c(x1nm, x2nm, nznm)
    #} else {
    #    colnames(gg) = c(x1nm, x2nm)
    #}

    if (npr > 1) {
        new_other = matrix(0, nrow=nrow(newData), ncol=(2*(npr-1)))
        kts_other = kts[,-c(ipr*2-1,ipr*2)]
        for (i in 1:(npr-1)) {
            for (j in 1:2) {
                newi = mean(kts_other[,(i-1)*2+j])
                new_other[,(i-1)*2+j] = newi
            }
        }
        newd = cbind(newData, new_other)
        colnames(newd) = c(xnm12, xnm_other)
        newData = as.data.frame(newd)
    }

    if (length(xnms_add) > 0) {
        #new_add = matrix(0, nrow=nrow(newData), ncol=ncol(xmat_add))
        means = apply(xmat_add, 2, mean)
        new_add = matrix(rep(means, nrow(newData)), ncol=ncol(xmat_add), byrow=T)
        nms = colnames(newData)
        newd = cbind(newData, new_add)
        colnames(newd) = c(nms, xnms_add)
        newData = as.data.frame(newd)
    }

    if (!is.null(newz)) {
        znms = obj$znms
        nz = matrix(0, nrow=nrow(newData), ncol=ncol(newz))
        nznm = gsub("[\\(\\)]", "", regmatches(znms, gregexpr("\\(.*?\\)", znms))[[1]])
        nms = colnames(newData)
        newData = cbind(newData, nz)
        colnames(newData) = c(nms, nznm)
    }

    gg = newData
    pfit.gg = predict.trispl(obj, gg, interval='none')
    #upper = pfit.gg$upper
    #lower = pfit.gg$lower

    xmatpr = pfit.gg$xmatpr
    muhat = pfit.gg$fit
    lower = muhat - mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
    upper = muhat + mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))

    #spls = pfit.bb$spls
    #ignore z
    #spl_use = cbind(1, spls[[ipr]])
    #get the fit for each pair
    #mus = pfit.knots$mus
    #muhat_use = mus[[ipr]]
    #get the acov for each pair
    #acov_use = acov[c(1, which(varlist == ipr)+np), c(1, which(varlist == ipr)+np)]

    #lower = muhat_use - mult*sqrt(diag(spl_use%*%acov_use%*%t(spl_use)))
    #upper = muhat_use + mult*sqrt(diag(spl_use%*%acov_use%*%t(spl_use)))

    x1g = res$x1g
    x2g = res$x2g
    k1 = length(x1g)
    k2 = length(x2g)

    surf = matrix(0, k1, k2)
    for(i2 in 1:k2) {
        for(i1 in 1:k1) {
            if (up) {
                surf[i1,i2] = upper[(i2-1)*k1 + i1]
            } else {
                surf[i1,i2] = lower[(i2-1)*k1 + i1]
            }
        }
    }
    #gaussian only
    z_add = res$z_add
    x3_add = res$x3_add
    surf = surf + z_add + x3_add
    if (up) {
        if (is.null(main)) {
            main = "Triangle-Spline Surface with Upper 95% Confidence Surface"
        }
    }
    if (!up) {
        if (is.null(main)) {
            main = "Triangle-Spline Surface with Lower 95% Confidence Surface"
        }
    }
    par(new = TRUE)
    persp(x1g, x2g, surf, zlim = res$zlim, xlab = "", ylab = "", zlab = "", theta = res$theta,
          ltheta = res$ltheta, cex.axis = res$cex.axis, main = main, cex.main = cex.main,
          ticktype = res$ticktype, col=mycol, box=FALSE, axes=FALSE,...)
    par(new=FALSE)
}


#####################################################################################
#### SUBROUTINES
#####################################################################################
### subroutine to determine if point p is inside the triangle formed by a,b,c
intri = function(p, a, b, c) {
	if (sameside(p, a, b, c) & sameside(p, b, a, c) & sameside(p, c, a, b)) {
		ins = TRUE
	} else {ins = FALSE}
	ins
}

### subroutine to determine if p1 and p2 are on the same side of the line formed by a,b
sameside = function(p1, p2, a, b) {
	cp1 = cpfun(b-a, p1-a)
	cp2 = cpfun(b-a, p2-a)
	if (cp1 * cp2 >= 0) {ss = TRUE} else {ss = FALSE}
	ss
}

## subroutine to find 3rd value of cross product of a and b, where a=(a1,a2,0) and b=(b1,b2,0)
cpfun = function(a, b) {
	c = a[1] * b[2] - a[2] * b[1]
	c
}

###########
#prop odds#
###########
cgam.polr <- function(formula, data = NULL, weights = NULL, family = NULL, nsim = 0, cpar = 1.2)
{
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  #print (mf)
  mf[[1L]] <- as.name("model.frame")
  #print (class(eval.parent(mf$data)))
  #if(is.matrix(eval.parent(mf$data))) {
  #	mf$data <- as.data.frame(data)
  #} else {mf$data <- data}
  mf <- eval(mf, parent.frame())
  ynm <- names(mf)[1]
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  shapes1 <- NULL; shapes2 <- NULL
  xmat <- NULL; xnms <- NULL
  nums <- NULL; ks <- list(); sps <- NULL; xid <- 1
  zmat <- NULL; zid <- NULL; zid0 <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_param <- NULL; is_fac <- NULL; vals <- NULL; st <- 1; ed <- 1
  ztb <- list(); iztb <- 1
  pl <- NULL
  tree.delta <- NULL
  tid1 <- NULL; tid2 <- NULL; tpos2 <- 0
  tr <- NULL
  umbrella.delta <- NULL
  uid1 <- NULL; uid2 <- NULL; upos2 <- 0
  umb <- NULL
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
  ans <- cgam.polr.fit(y, xmat0, zmat, umbrella.delta = umbrella.delta, tree.delta = tree.delta, shapes = shapes1, nums0, ks0, sps0, family, idx_s, idx, nsim, cpar, weights)
  rslt <- list(muhat = ans$muhat, zeta = ans$ac, wta = ans$wta, lev = ans$lev, vcoefs = ans$vcoefs, xcoefs = ans$xcoefs, zcoefs = ans$zcoefs, coefs = ans$coefs, ucoefs = ans$ucoefs, tcoefs = ans$tcoefs, cic = ans$cic, d0 = ans$d0, edf0 = NULL, etacomps = ans$etacomps, xmat = xmat, zmat = zmat, ztb = ztb, tr = tr, umb = umb, tree.delta = tree.delta, umbrella.delta = umbrella.delta, bigmat = ans$bigmat, shapes = shapes1, shapesx = shapes1, wt = ans$wt, wt.iter = TRUE, family = family, SSE0 = NULL, SSE1 = NULL, pvals.beta = NULL, se.beta = NULL, df.null = ans$df.null, df = ans$df, df.residual = ans$df.residual, null.deviance = ans$dev.null, deviance = ans$dev, tms = mt, capm = ans$capm, capms = ans$capms, capk = ans$capk, capt = ans$capt, capu = ans$capu, xid1 = ans$xid1, xid2 = ans$xid2, tid1 = tid1, tid2 = tid2, uid1 = uid1, uid2 = uid2, zid = zid, vals = vals, zid1 = zid1, zid2 = zid2, nsim = nsim, xnms = xnms,  ynm = ynm, znms = znms, is_param = is_param, is_fac = is_fac, knots = ans$knots, numknots = ans$numknots, sps = sps, ms = NULL, cpar = cpar, pl = pl, idx_s = idx_s, idx = idx, xmat0 = ans$xmat2, knots0 = ans$knots2, numknots0 = ans$numknots2, sps0 = ans$sps2, ms0 = ans$ms2)
  rslt$call <- cl
  class(rslt) <- "cgam.polr"
  return (rslt)
}

#
cgam.polr.fit = function(y = NULL, xmat = NULL, zmat = NULL, umbrella.delta = NULL, tree.delta = NULL, shapes = NULL, numknots = NULL, knots = NULL, space = NULL, family = NULL, idx_s = NULL, idx = NULL, nsim = 0, cpar = 1.2, weights = NULL)  {
	fmin = function(a) {
		nc = length(a) + 1
		y_id = t(col(matrix(0, length(y), nc)) == y)
		llh = 0
		for (i in 1:nc) {
			if (i == 1) {
				llhi = sum(w * y_id[i,] * log(pfun(a[1] - eta)))
			} else if (i == nc) {
				llhi = sum(w * y_id[i,] * log(1 - pfun(a[nc-1] - eta)))
			} else {
				llhi = sum(w * y_id[i,] * log(pfun(a[i] - eta) - pfun(a[i-1] - eta)))
			}
			llh = llh + llhi
		}
		-2*llh
	}
	gmin = function(a) {
		nc = length(a) + 1
		y_id = t(col(matrix(0, n, nc)) == y)
		gr = NULL
		for (i in 1:(nc-1)) {
			if (i == 1) {
				gri = -sum(y_id[1,] * w * (1 - pfun(a[1] - eta))) + sum(y_id[2,] * w * dfun(a[1] - eta) / (pfun(a[2] - eta) - pfun(a[1] - eta)))
			} else if (i == (nc-1)) {
				gri = -sum(y_id[nc-1,] * w * dfun(a[nc-1] - eta) / (pfun(a[nc-1] - eta) - pfun(a[nc-2] - eta))) + sum(I(y == nc) * w * dfun(a[nc-1] - eta) / (1 - pfun(a[nc-1] - eta)))
			#gri = -sum(y_id[nc-1,id] * w[id] * dfun(a[nc-1] - eta[id]) / (pfun(a[nc-1] - eta[id]) - pfun(a[nc-2] - eta[id]))) + sum(y_id[nc,id]  * w[id] * dfun(a[nc-1] - eta[id]) / (1 - pfun(a[nc-1] - eta[id])))
			} else {
				gri = -sum(y_id[i,] * w * dfun(a[i] - eta) / (pfun(a[i] - eta) - pfun(a[i-1] - eta))) + sum(y_id[i+1,] * w * dfun(a[i] - eta) / (pfun(a[i+1] - eta) - pfun(a[i] - eta)))
			}
			gr = c(gr, gri)
		}
		return (gr)
	}
	if (is.factor(y)) {
		lev = levels(y)
    	y = unclass(y)
	} else {
		lev = attributes(factor(y))$levels
	}
	n = length(y)
	w = weights
	if(is.null(w)) {
		w = rep(1, n)
	}
	sc_x = FALSE
	capl = length(xmat) / n
	#print (dim(xmat))
	if (capl < 1) {capl = 0}
	if (round(capl, 8) != round(capl, 1)) {stop ("Incompatible dimensions for xmat!")}
	if (capl > 0 & sc_x) {
		for (i in 1:capl) {xmat[,i] = (xmat[,i] - min(xmat[,i])) / (max(xmat[,i]) - min(xmat[,i]))}
	}
	capk = length(zmat) / n
	#print (capk)
	if (capk < 1) {capk = 0}
	if (round(capk, 8) != round(capk, 1)) {stop ("Incompatible dimensions for zmat!")}
#smooth only
#if (!is.null(shapes)) {
	capls = sum(shapes == 17)
#} else {capls = 0}
	bigmat = NULL; delta = NULL
	varlist = NULL
	xid1 = NULL; xid2 = NULL; xpos2 = 0
	knotsuse = list(); numknotsuse = NULL
	mslst = list()
	capm = 0
	capms = 0
	if (capl - capls > 0) {
		#print (class(xmat[,1]))
		#print (head(xmat))
		del1_ans = makedelta(xmat[, 1], shapes[1], numknots[1], knots[[1]], space = space[1], interp=T)
		#del1_ans = makedelta(xmat[, 1], shapes[1])
		del1 = del1_ans$amat
		knotsuse[[1]] = del1_ans$knots
		mslst[[1]] = del1_ans$ms
		numknotsuse = c(numknotsuse, length(del1_ans$knots))
        m1 = length(del1) / n
#new code: record the number of columns of del1 if shapes0[1] == 17:
		if (shapes[1] == 17) {capms = capms + m1}
        var1 = 1:m1*0 + 1
		xpos1 = xpos2 + 1
		xpos2 = xpos2 + m1
		xid1 = c(xid1, xpos1)
		xid2 = c(xid2, xpos2)
		if (capl == 1) {
        	delta = del1
         	varlist = var1
        } else {
	    	for (i in 2:capl) {
#new code:
	        	del2_ans = makedelta(xmat[,i], shapes[i], numknots[i], knots[[i]], space = space[i], interp=T)
				#del2_ans = makedelta(xmat[,i], shapes[i])
				del2 = del2_ans$amat
				knotsuse[[i]] = del2_ans$knots
				mslst[[i]] = del2_ans$ms
				numknotsuse = c(numknotsuse, length(del2_ans$knots))
				m2 = length(del2) / n
#new code: record the number of columns of del2 if shapes0[i] == 17:
				if (shapes[i] == 17) {capms = capms + m2}
				xpos1 = xpos2 + 1
				xpos2 = xpos2 + m2
				xid1 = c(xid1, xpos1)
				xid2 = c(xid2, xpos2)
				delta = rbind(del1, del2)
				varlist = 1:(m1 + m2)*0
				varlist[1:m1] = var1
				varlist[(m1 + 1):(m1 + m2)] = (1:m2)*0 + i
				var1 = varlist
				m1 = m1 + m2
				del1 = delta
	      	}
	    }
	    np = 0
		if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk > 0) {
			bigmat = rbind(t(zmat), t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
			np = capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)  + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk == 0) {
			bigmat = rbind(t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
			np = sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk > 0) {
			bigmat = rbind(t(zmat), delta)
			np = capk + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk == 0) {
			bigmat = delta
			np = capms
		} else {
			print ("error in capk, shapes!")
		}
#new:
	capm = length(delta) / n - capms
	} else {
	  	if (capk + capls > 0) {
#new:
			if (capls  < 1 & capk > 0) {
          			#bigmat = rbind(1:n*0 + 1, t(zmat))
          			bigmat = t(zmat)
          			#np = 1 + capk
          			np = capk
			} else if (capls > 0) {
				delta = NULL; varlist = NULL
				del1_ans = makedelta(xmat[,1], 17, numknots[1], knots[[1]], space = space[1], interp=T)
				#del1_ans = makedelta(xmat[,i], 9)
				del1 = del1_ans$amat
				knotsuse[[1]] = del1_ans$knots
				mslst[[1]] = del1_ans$ms
				numknotsuse = c(numknotsuse, length(del1_ans$knots))
				m1 = length(del1) / n
				var1 = 1:m1*0 + 1
				xpos1 = xpos2 + 1
				xpos2 = xpos2 + m1
				xid1 = c(xid1, xpos1)
				xid2 = c(xid2, xpos2)
				if (capls == 1) {
        			delta = del1
         			varlist = var1
          		} else {
					for (i in 2:capls) {
	        			del2_ans = makedelta(xmat[,i], 17, numknots[i], knots[[i]], space = space[i], interp=T)
						#del2_ans = makedelta(xmat[,i], 9)
						del2 = del2_ans$amat
						knotsuse[[i]] = del2_ans$knots
						mslst[[i]] = del2_ans$ms
						numknotsuse = c(numknotsuse, length(del2_ans$knots))
						m2 = length(del2) / n
						xpos1 = xpos2 + 1
						xpos2 = xpos2 + m2
						xid1 = c(xid1, xpos1)
						xid2 = c(xid2, xpos2)
						delta = rbind(del1, del2)
						varlist = 1:(m1 + m2)*0
						varlist[1:m1] = var1
						varlist[(m1 + 1):(m1 + m2)] = (1:m2)*0 + i
						var1 = varlist
						m1 = m1 + m2
						del1 = delta
	      			}
				}
				if (capk < 1){
					#bigmat = rbind(1:n*0 + 1, delta)
					bigmat = delta
					capms = length(delta) / n
					#np = 1 + capms
					np = capms
				} else {
					#bigmat = rbind(1:n*0 + 1, t(zmat), delta)
					bigmat = rbind(t(zmat), delta)
					capms = length(delta) / n
					#np = 1 + capk + capms
					np =  capk + capms
				}
			}
        } else {bigmat = NULL; capm = 0; capms = 0; np = 0}
	}
	#print (head(zmat))
	#print (head(t(bigmat)))
#new:
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
#coneB:
	#n = length(y)
	nc = length(unique(y)) - 1
	#lev = attributes(factor(y))$levels
#initialize a
if (all(w == 1)) {
	a = sapply(1:nc, function(ct) log(sum(y <= ct) / (n - sum(y <= ct))))
} else {
	#a = sapply(1:nc, function(ct) log(w * sum(y <= ct) / (n - sum(y <= ct))))
	a = 1:nc*0
	for (i in 1:nc) {
		#idi = which(y <= i)
		a[i] = log(sum((y <= i)* w) / (sum(w) - sum((y <= i)* w)))
	}
}
	gr = 1:n*0
	wt = 1:n*0
	eta = 1:n*0
	diff = 100
	nrep = 0
	olddev = 100
	oldmu = 1:n*0+1
    face = NULL
	while(diff > 1e-8 & nrep < 400){
		nrep = nrep+1
#eta step:
		gr = fgr(a, y, eta, w)
		wt = fwt(a, y, eta, w)
		q0 = diag(as.vector(wt))
		cvec = wt * eta - gr
		zvec = 1:n*0
     	zvec[wt == 0] = 1 / 1e-4
     	zvec[wt > 0] = cvec[wt > 0] / sqrt(wt[wt > 0])
		gmat = t(bigmat)
		for (j in 1:n) {gmat[j,] = bigmat[,j] * sqrt(wt[j])}
		if (np > 0) {
			zsend = gmat[, 1:np]
			if ((capm + capt + capu) > 0) {
				dsend = t(gmat[, (np+1):(nrow(bigmat))])
			} else {dsend = NULL}
		} else {
			dsend = t(gmat)
			zsend = NULL
		}
		if ((capm + capt + capu) > 0) {
            #ans = coneB(zvec, dsend, vmat = zsend)
            if (nrep > 1) {
                ans = coneB(zvec, t(dsend), vmat = zsend, face = face)
            } else {
                ans = coneB(zvec, t(dsend), vmat = zsend)
                face = ans$face
            }
			#ans = coneB(zvec, t(dsend))
			bh = coef(ans)
		} else {
			bh = solve(t(zsend) %*% zsend, t(zsend) %*% zvec)
		}
		eta = t(bigmat) %*% bh
#a step:
#res = optim(a, fmin, gr = gmin, hessian = TRUE)
#a = res$par
        gr_a = fgr_a(a, y, eta, w)
		wt_a = fwt_a(a, y, eta, w)
		qmat_a = wt_a
		a = a - solve(qmat_a, gr_a)
#a = atil + bh[1]
		#if (a[1] >= a[2]) {
		if (any(a != sort(a))) {
			print ('a wrong!')
			break
		}
		#dev = res$value
		#diff = (dev-olddev)^2
		#olddev = dev
		mu = predict_polr(a, eta)$mu
		diff = mean((mu-oldmu)^2)
		oldmu = mu
	}
	names(a) = paste(lev[-length(lev)], lev[-1L], sep="|")
	#print (diff)
	yhat = eta
	coefskeep = bh
	#df_obs = sum(abs(coefskeep) > 0) + length(a)
	#print (coefskeep)
########################
#if capk > 0, we have:#
########################
	zcoefs = NULL
	if (capk > 0) {
		zcoefs = coefskeep[1:capk]
	}
	#print (zcoefs)
	vcoefs = NULL
	if (np > 0) {
		vcoefs = coefskeep[1:np]
	}
#######################
#if capm > 0, we have:#
#######################
	xcoefs = NULL
#new:
	if (capl > 0) {
		xcoefs = coefskeep[(np - capms + 1):(np + capm)]
	}
#######################
#if capu > 0, we have:#
#######################
	ucoefs = NULL
	if (capu > 0) {
		ucoefs = coefskeep[(np + 1 + capm):(np + capm + capu)]
	}
#######################
#if capt > 0, we have:#
#######################
	tcoefs = NULL
	if (capt > 0) {
		tcoefs =
		 coefskeep[(np + 1 + capm + capu):(np + capm + capu + capt)]
	}
#########################################################
#if we have at least one constrained predictor, we have:#
#########################################################
	thvecs = NULL
	if (capl > 0) {
#new code:
		dcoefs = coefskeep[(np - capms + 1):(np + capm)]
#####################################################
#thvecs is f(x), where x has one of the eight shapes#
#####################################################
		thvecs = matrix(0, nrow = capl, ncol = n)
	    ncon = 1
	    for (i in 1:capl) {
	    	thvecs[i,] = t(delta[varlist == i,]) %*% dcoefs[varlist == i]
			if (shapes[i] > 2 & shapes[i] < 5 | shapes[i] > 10 & shapes[i] < 13) {
            	ncon = ncon + 1
				thvecs[i,] = thvecs[i,] + vcoefs[capk + ncon] * xmat[,i]
            }
	    }
	}
#new:order thvecs back
	#if (!is.null(idx_s)) {
	if (length(idx_s) > 0) {
		thvecs0 = thvecs
		thvecs0[idx_s,] = thvecs[1:length(idx_s), ]
		#if (!is.null(idx)) {
		if (length(idx) > 0) {
			thvecs0[idx,] = thvecs[(1+length(idx_s)):capl, ]
		}
		thvecs = thvecs0
	}
	thvecs_ut = NULL
	if (capu + capt > 0) {
		thvecs_ut = t(delta_ut) %*% coefskeep[(np + 1 + capm):(np + capm + capu + capt)]
	}
	if (!is.null(thvecs_ut)) {
		thvecs = rbind(thvecs, t(thvecs_ut))
	}
	etakeep = eta
	muhatkeep = eta
	df_obs = sum(abs(coefskeep) > 0) + length(a)
#deviance
	dev = dev_fun(y, a = a, eta = eta, flat = FALSE, w = w)
	#print (dev)
	llh = dev/n
#null deviance
	dev.null = dev_fun(y, a = a, eta = eta, flat = TRUE, w = w)
	#print (dev.null)
	cic = NULL
	dfmean = 0
	if (nsim > 0) {
		if (capm + capms > 0) {
			dfs = 1:nsim*0
			for (isim in 1:nsim) {
				#print (isim)
				ysti = rlogis(n, scale = 1)
				#print (nc)
				cts = quantile(ysti, probs = seq(0, 1, length = (nc+2)))
				#print (cts)
				yordi = cut(ysti, breaks=cts, include.lowest = TRUE, labels = c(1:(nc+1)), ordered = TRUE)
				ysim = as.numeric(levels(yordi))[yordi]
				edfi = polr_getedf(ysim, bigmat, np, capm, shapes, w)
		  		#dfs[isim] = sum(abs(bhi) > 0)
		  		dfs[isim] = edfi
			}
			if (any(shapes == 11) | any(shapes == 12)) {
				dfmean = mean(dfs) - sum(shapes > 10 & shapes < 13)
			} else {dfmean = mean(dfs)}
		} else {dfmean = 0}
		#cic = llh + log(1 + 2 * (dfmean + np) / (n - np - cpar * dfmean))
#check! length(a)?
		if ((n - np - 1.5 * (dfmean - np)) <= 0) {
			cic = llh + log(1 + 2 * dfmean / (dfmean - np))
		} else {
        	cic = llh + log(1 + 2 * dfmean / (n - np - 1.5 * (dfmean - np)))
		}
	}
	#print (cic)
	env = new.env()
	#tm=proc.time()-tm0
	#env$tm=tm
	env$dev=dev
	env$dev.null=dev.null
	#print (lev)
	env$lev=lev
	env$etacomps=thvecs
	env$zcoefs=zcoefs
	env$ac=a
	env$etahat=eta
	env$muhat=eta
	env$vcoefs=vcoefs
	env$xcoefs=xcoefs
	env$zcoefs=zcoefs
	env$ucoefs=ucoefs
	env$tcoefs=tcoefs
	env$coefs=coefskeep
	env$cic=cic
	env$d0=np
	env$capm = capm
	env$capms = capms
	env$capk = capk
	env$capu = capu
	env$capt = capt
	env$edf=df_obs
	env$edf0=dfmean
	env$bigmat=bigmat
	#env$family=family
	env$sse0=NULL
	env$sse1=NULL
	env$pvals.beta=NULL
	env$se.beta=NULL
	env$wt=wt
	env$gr=gr
	env$mu=mu
	env$oldmu=oldmu
	env$diff=diff
	env$bh=bh
	env$thvecs=thvecs
	#gr_a=NULL
	env$gra=gr_a
	env$wta=wt_a
	knotsuse2=knotsuse
	numknotsuse2=numknotsuse
	mslst2=mslst
	xmat2=xmat
	#if (!is.null(idx_s)) {
	if (length(idx_s) > 0) {
		knotsuse0 = knotsuse
		numknotsuse0 = numknotsuse
		mslst0 = mslst
		knotsuse0[idx_s] = knotsuse[1:length(idx_s)]
		numknotsuse0[idx_s] = numknotsuse[1:length(idx_s)]
		mslst0[idx_s] = mslst[1:length(idx_s)]
		#if (!is.null(idx)) {
		if (length(idx) > 0) {
			knotsuse0[idx] = knotsuse[(1+length(idx_s)):capl]
			numknotsuse0[idx] = numknotsuse[(1+length(idx_s)):capl]
			mslst0[idx] = mslst[(1+length(idx_s)):capl]
		}
		knotsuse = knotsuse0
		numknotsuse = numknotsuse0
		mslst = mslst0
	}
	env$knots = knotsuse
	env$numknots = numknotsuse
	env$xid1 = xid1
	env$xid2 = xid2
	env$xmat2 = xmat2
	env$knotsuse2 = knotsuse2
	env$numknotsuse2 = numknotsuse2
	env$df = n - length(a)
	env$df.null = n - length(a)
	env$df.residual = n - cpar * df_obs
	#env$df.residual = n - np - length(a) - cpar * df_obs
	#env$sps2 = sps2
	#env$ms2 = ms2
	return(env)
}

###################
#summary.cgam.polr#
###################
summary.cgam.polr = function(object,...) {
	coefs = object$zcoefs
	n = length(coefs)
	zid = object$zid
	zid1 = object$zid1
	zid2 = object$zid2
	tms = object$tms
	is_param = object$is_param
	is_fac = object$is_fac
	vals = object$vals
	zeta = object$zeta
	hess = object$wta
    pc = length(zeta)
    np = object$d0
    df_obs = object$df_obs
    cpar = object$cpar
    #pvals = object$pvals
    #if ((n - np - cpar * df_obs) <= 0) {
	#			pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  df_obs))
	#			warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
	#} else {
	#	pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  n - np - cpar * df_obs))
	#}
	if (length(coefs) >= 1) {
		rslt1 = data.frame("Estimate" = round(coefs, 4))
		#rownames(rslt1)[1] <- "(Intercept)"
		if (n >= 1) {
			lzid = length(zid1)
			for (i in 1:lzid) {
				pos1 = zid1[i]
				pos2 = zid2[i]
				for (j in pos1:pos2) {
					if (!is_param[i]) {
						rownames(rslt1)[j] = paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j], sep = "")
					} else {
						rownames(rslt1)[j] = paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")
					}
				}
			}
		}
		rslt1 = as.matrix(rslt1)
	} else {rslt1 = NULL}
	rslt = matrix(0, pc, 3L, dimnames = list(names(zeta), c("Estimate", "StdErr", "t.value")))#"p.value"
    rslt[, 1L] = round(zeta, 4)
    hess = object$wta
    covmat = solve(hess)
    rslt[, 2L] = round(sqrt(diag(covmat)), 4)
    rslt[, 3L] = round(rslt[, 1L]/rslt[, 2L], 4)
    #tstat = rslt[, 3L]
	#pvals = 2 * (1 - pt(abs(tstat),  n - np - pc - cpar * (nzeta+df_obs)))
    #rslt[, 4L] = pvals
    #object$coefficients = coefs
    #object$pc = pc
    #if(correlation)
    #    object$correlation <- (vc/sd)/rep(sd, rep(pc+q, pc+q))
    ans = list(call = object$call, coefficients = rslt, cic = object$cic, pc=pc, deviance=object$deviance, df.residual=object$df.residual, rslt1 = rslt1)
    class(ans) = "summary.cgam.polr"
    return (ans)
}


#########################
#print.summary.cgam.polr#
#########################
print.summary.cgam.polr = function(x,...) {
    #rslt = x$coefficients
    #pc = x$pc
    #if(pc > 0) {
    #    cat("\nCoefficients:\n")
    #    print(x$coefficients[seq_len(pc), , drop=FALSE], quote = FALSE,
    #          digits = digits, ...)
    #} else {
    #    cat("\nNo coefficients\n")
    #}
    cat("Call:\n")
	print(x$call)
	if (!is.null(x$rslt1)) {
		cat("\n")
		cat("Coefficients:")
		cat("\n")
		print (x$rslt1)
	}
    cat("\nIntercepts:\n")
    #print(rslt[(pc+1L):nrow(rslt), , drop=FALSE], quote = FALSE)
    printCoefmat(x$coefficients, P.values = FALSE, has.Pvalue = FALSE)
    #cat("\nResidual Deviance:", format(x$deviance, nsmall=2L), "\n")
    cat("\nResidual deviance: ", round(x$deviance, 4), " ", "on ", x$df.residual, " ", "observed degrees of freedom", sep="", "\n")
    if (!is.null(x$cic)) {
    	cat("CIC:", format(x$cic), "\n")
    }
    #if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    #if(!is.null(correl <- x$correlation)) {
    #   cat("\nCorrelation of Coefficients:\n")
    #    ll <- lower.tri(correl)
    #    correl[ll] <- format(round(correl[ll], digits))
    #    correl[!ll] <- ""
    #    print(correl[-1L, -ncol(correl)], quote = FALSE, ...)
    #}
    invisible(x)
}

dev_fun = function(y, a = NULL, eta = NULL, flat = TRUE, w = NULL) {
	n = length(y)
	if (is.null(w)) {
		w = rep(1, n)
	}
	if (flat) {
		eta = rep(0, n)
		luy = length(unique(y))-1
		if (all(w == 1)) {
			a = sapply(1:luy, function(ct) log(sum(y <= ct)/(n-sum(y <= ct))))
		} else {
			a = 1:luy*0
			for (i in 1:luy) {
				a[i] = log(sum((y <= i)* w) / (sum(w) - sum((y <= i)* w)))
			}
		}
	}
	nc = length(a) + 1
	y_id = t(col(matrix(0, n, nc)) == y)
	llh = 0
	for (i in 1:nc) {
		if (i == 1) {
			llhi = sum(w * y_id[i,] * log(pfun(a[1] - eta)))
		} else if (i == nc) {
			llhi = sum(w * y_id[i,] * log(1 - pfun(a[nc-1] - eta)))
		} else {
			llhi = sum(w * y_id[i,] * log(pfun(a[i] - eta) - pfun(a[i-1] - eta)))
		}
		llh = llh + llhi
	}
	#llh = -2/n*llh
	dev = -2*llh
	return (dev)
}

##
predict_polr = function(a, eta) {
	nc = length(a) + 1
	n = length(eta)
	ps = matrix(0, nrow = nc, ncol = n)
	cums = matrix(0, nrow = nc, ncol = n)
	for (i in 1:nc) {
		if (i == 1) {
			ps[i,] = pfun(a[1] - eta)
			#cums[, i] = ps[,  i]
		} else if (i == nc) {
			ps[i,] = 1 - pfun(a[nc-1] - eta)
		} else {
			ps[i,] = pfun(a[i] - eta) - pfun(a[i-1] - eta)
			#cums[, i] = ps[,  i] + ps[, i-1]
		}
	}
	cums = apply(ps, 2, cumsum)
	rhs = a[nc-1] - eta
	mu = exp(rhs) / (1 + exp(rhs))
	rslt = new.env()
	rslt$ps = ps
	rslt$mu = mu
	rslt$cums = cums
	return (rslt)
}

polr_getedf = function(ysim, bigmat, np, capm, shapes, w = NULL) {
	fmin = function(a) {
		nc = length(a) + 1
		y_id = t(col(matrix(0, length(ysim), nc)) == ysim)
		llh = 0
		for (i in 1:nc) {
			if (i == 1) {
				llhi = sum(w * y_id[i,] * log(pfun(a[1] - eta)))
			} else if (i == nc) {
				llhi = sum(w * y_id[i,] * log(1 - pfun(a[nc-1] - eta)))
			} else {
				llhi = sum(w * y_id[i,] * log(pfun(a[i] - eta) - pfun(a[i-1] - eta)))
			}
			llh = llh + llhi
		}
		-2*llh
	}
	gmin = function(a) {
		nc = length(a) + 1
		y_id = t(col(matrix(0, n, nc)) == ysim)
		gr = NULL
		for (i in 1:(nc-1)) {
			if (i == 1) {
				gri = -sum(y_id[1,] * w * (1 - pfun(a[1] - eta))) + sum(y_id[2,] * w * dfun(a[1] - eta) / (pfun(a[2] - eta) - pfun(a[1] - eta)))
			} else if (i == (nc-1)) {
				gri = -sum(y_id[nc-1,] * w * dfun(a[nc-1] - eta) / (pfun(a[nc-1] - eta) - pfun(a[nc-2] - eta))) + sum(I(ysim == nc) * w * dfun(a[nc-1] - eta) / (1 - pfun(a[nc-1] - eta)))
			#gri = -sum(y_id[nc-1,id] * w[id] * dfun(a[nc-1] - eta[id]) / (pfun(a[nc-1] - eta[id]) - pfun(a[nc-2] - eta[id]))) + sum(y_id[nc,id]  * w[id] * dfun(a[nc-1] - eta[id]) / (1 - pfun(a[nc-1] - eta[id])))
			} else {
				gri = -sum(y_id[i,] * w * dfun(a[i] - eta) / (pfun(a[i] - eta) - pfun(a[i-1] - eta))) + sum(y_id[i+1,] * w * dfun(a[i] - eta) / (pfun(a[i+1] - eta) - pfun(a[i] - eta)))
			}
			gr = c(gr, gri)
		}
		return (gr)
	}
	nc = length(unique(ysim)) - 1
	n = length(ysim)
	if (is.null(w)) {
		w = rep(1, n)
	}
	if (all(w == 1)) {
		a = sapply(1:nc, function(ct) log(sum(ysim <= ct) / (n - sum(ysim <= ct))))
	} else {
		a = 1:nc*0
		for (i in 1:nc) {
			a[i] = log(sum((ysim <= i)* w) / (sum(w) - sum((ysim <= i)* w)))
		}
	}
	gr = 1:n*0
	wt = 1:n*0
	eta = 1:n*0
	diff = 100
	nrep = 0
	oldmu = 1:n*0+1
    face = NULL
	while(diff > 1e-8 & nrep < 400){
		nrep = nrep+1
#eta step:
		gr = fgr(a, ysim, eta, w)
        wt = fwt(a, ysim, eta, w)
		q0 = diag(as.vector(wt))
		cvec = wt * eta - gr
		zvec = 1:n*0
     	zvec[wt == 0] = 1 / 1e-4
     	zvec[wt > 0] = cvec[wt > 0] / sqrt(wt[wt > 0])
		gmat = t(bigmat)
		for (j in 1:n) {gmat[j,] = bigmat[,j] * sqrt(wt[j])}
		if (np > 0) {
			zsend = gmat[, 1:np]
			if (capm > 0) {
				dsend = t(gmat[, (np+1):(nrow(bigmat))])
			} else {dsend = NULL}
		} else {
			dsend = t(gmat)
			zsend = NULL
		}
		#print (dsend)
		if (capm > 0) {
            #ans = coneB(zvec, dsend, vmat = zsend)
            if (nrep > 1) {
                ans = coneB(zvec, t(dsend), vmat = zsend, face = face)
            } else {
                ans = coneB(zvec, t(dsend), vmat = zsend)
                face = ans$face
            }
			bh = coef(ans)
		} else {
			bh = solve(t(zsend) %*% zsend, t(zsend) %*% zvec)
		}
		#print (bh)
		eta = t(bigmat) %*% bh
#a step:
#res = optim(a, fmin, gr = gmin)
#a = res$par
#print (a)
		gr_a = fgr_a(a, ysim, eta,  w)
		wt_a = fwt_a(a, ysim, eta, w)
		qmat_a = wt_a
		a = a - solve(qmat_a, gr_a)
#a = atil + bh[1]
		#if (a[1] >= a[2]) {
		if (any(a != sort(a))) {
			print ('a wrong!')
			break
		}
		mu = predict_polr(a, eta)$mu
		diff = mean((mu-oldmu)^2)
		oldmu = mu
	}
	edf = sum(abs(bh) > 0)
	return (edf)
}

#big phi
pfun = function(x, sc = 1) {
    if (sc == 1) {
		#return (exp(x) / (1 + exp(x)))
		return (1 - 1 / (1 + exp(x)))
	} else {
		xtil = x / sc
		return (exp(xtil) / (1 + exp(xtil)))
	}
}

#small phi
dfun = function(x, sc = 1) {
     #d1 = exp(x) / (1 + exp(x))
     d1 = 1 - 1 / (1 + exp(x))
     d2 = 1 / (1 + exp(x))
     if (sc == 1) {
     	#print (range(d1*d2))
	     return (d1 * d2)
     } else {
		xtil = x / sc
		d1 = exp(xtil) / (1 + exp(xtil))
     	d2 = 1 / (1 + exp(xtil))
		return (1 / sc * d1 * d2)
     }
}

#derivative of small phi
ddfun = function(x, sc = 1) {
    if (sc == 1) {
      	#dd = dfun(x) * (1 - exp(x)) / (1 + exp(x))
      	dd = dfun(x) * (2 / (1 + exp(x)) - 1)
      	return (dd)
    } else {
		xtil = x / sc
		dd = 1 / sc * dfun(xtil, sc = sc) * (1 - exp(xtil)) / (1 + exp(xtil))
		return (dd)
    }
}

#k = 4: gr w.r.t eta
fgr = function(a, y, eta, w = NULL) {
	n = length(y)
	if (is.null(w)) {
		w = rep(1, n)
	}
	nc = length(a) + 1
	y_id = t(col(matrix(0, n, nc)) == y)
	gr = 1:n*0
	for (i in 1:nc) {
		if (i == 1) {
			gri = y_id[i,] * (1 - pfun(a[1] - eta)) * w
		} else if (i == nc) {
			gri = -y_id[i,] * (pfun(a[nc-1] - eta)) * w
		} else {
			gri = y_id[i,] * (1 - pfun(a[i-1] - eta) - pfun(a[i] - eta)) * w
		}
		gr = gr + gri
	}
	gr = round(gr, 6)
	return (gr)
}

#k = 4: Hess w.r.t eta
fwt = function(a, y, eta, w = NULL) {
	n = length(y)
	if (is.null(w)) {
		w = rep(1, n)
	}
	nc = length(a) + 1
	y_id = t(col(matrix(0, n, nc)) == y)
	wt = 1:n*0
	for (i in 1:nc) {
		if (i == 1) {
			wti = y_id[i,] * dfun(a[1] - eta) * w
		} else if (i == nc) {
			wti = y_id[i,] * dfun(a[nc-1] - eta) * w
		} else {
			wti = y_id[i,] * (dfun(a[i-1] - eta) + dfun(a[i] - eta)) * w
		}
		wt = wt + wti
	}
	wt = round(wt, 6)
	return (wt)
}

#k = 4: gr w.r.t a1, a2, a3
fgr_a = function(a, y, eta, w = NULL) {
	n = length(y)
	if (is.null(w)) {
		w = rep(1, n)
	}
	gr = NULL
	nc = length(a) + 1
	y_id = t(col(matrix(0, n, nc)) == y)
	for (i in 1:(nc-1)) {
		if (i == 1) {
			id = (pfun(a[2] - eta) - pfun(a[1] - eta)) != 0
			gri = -sum(y_id[1,id] * w[id] * (1 - pfun(a[1] - eta[id]))) + sum(y_id[2,id] * w[id] * dfun(a[1] - eta[id]) / (pfun(a[2] - eta[id]) - pfun(a[1] - eta[id])))
		} else if (i == (nc-1)) {
			id = (pfun(a[nc-1] - eta) - pfun(a[nc-2] - eta)) != 0
			gri = -sum(y_id[nc-1,id] * w[id] * dfun(a[nc-1] - eta[id]) / (pfun(a[nc-1] - eta[id]) - pfun(a[nc-2] - eta[id]))) + sum(I(y[id] == nc) * w[id] * dfun(a[nc-1] - eta[id]) / (1 - pfun(a[nc-1] - eta[id])))
			#gri = -sum(y_id[nc-1,id] * w[id] * dfun(a[nc-1] - eta[id]) / (pfun(a[nc-1] - eta[id]) - pfun(a[nc-2] - eta[id]))) + sum(y_id[nc,id]  * w[id] * dfun(a[nc-1] - eta[id]) / (1 - pfun(a[nc-1] - eta[id])))
		} else {
			id1 = (pfun(a[i] - eta) - pfun(a[i-1] - eta)) != 0
			id2 = (pfun(a[i+1] - eta) - pfun(a[i] - eta)) != 0
			id12 = id1 & id2
			gri = -sum(y_id[i,id12] * w[id12] * dfun(a[i] - eta[id12]) / (pfun(a[i] - eta[id12]) - pfun(a[i-1] - eta[id12]))) + sum(y_id[i+1,id12] * w[id12] * dfun(a[i] - eta[id12]) / (pfun(a[i+1] - eta[id12]) - pfun(a[i] - eta[id12])))
		}
		gr = c(gr, gri)
	}
	return (gr)
}

#k = 4: Hess diagonal w.r.t a1, a2, a3 and off-diagonal
fwt_a = function(a, y, eta, w = NULL) {
	n = length(y)
	if (is.null(w)) {
		w = rep(1, n)
	}
	nc = length(a) + 1
	y_id = t(col(matrix(0, n, nc)) == y)
	wt = matrix(0, (nc-1), (nc-1))
	for (i in 1:(nc-1)) {
		if (i == 1) {
			id = (pfun(a[2] - eta) - pfun(a[1] - eta)) != 0
			wt[1,1] = sum(y_id[1,id] * w[id] * dfun(a[1] - eta[id])) + sum(y_id[2,id] * w[id] * (ddfun(a[1] - eta[id]) / (pfun(a[2] - eta[id]) - pfun(a[1] - eta[id])) + (dfun(a[1] - eta[id]) / (pfun(a[2] - eta[id]) - pfun(a[1] - eta[id])))^2))
			wt[1,2] = wt[2,1] = -sum(y_id[2,id] * w[id] * (dfun(a[2] - eta[id]) * dfun(a[1] - eta[id])) / (pfun(a[2] - eta[id]) - pfun(a[1] - eta[id]))^2)
		} else if (i == (nc-1)) {
			id = (pfun(a[nc-1] - eta) - pfun(a[nc-2] - eta)) != 0
			wt[nc-1, nc-1] = -sum(y_id[nc-1,id] * w[id] * (ddfun(a[nc-1] - eta[id]) / (pfun(a[nc-1] - eta[id]) - pfun(a[nc-2] - eta[id])) - (dfun(a[nc-1] - eta[id]) / (pfun(a[nc-1] - eta[id]) - pfun(a[nc-2] - eta[id])))^2)) + sum(y_id[nc,id] * w[id] * dfun(a[nc-1] - eta[id]))
		} else {
			id1 = (pfun(a[i] - eta) - pfun(a[i-1] - eta)) != 0
			id2 = (pfun(a[i+1] - eta) - pfun(a[i] - eta)) != 0
			id12 = id1 & id2
			wt[i,i] = -sum(y_id[i,id12] * w[id12] * (ddfun(a[i] - eta[id12]) / (pfun(a[i] - eta[id12]) - pfun(a[i-1] - eta[id12])) - (dfun(a[i] - eta[id12]) / (pfun(a[i] - eta[id12]) - pfun(a[i-1] - eta[id12])))^2)) + sum(y_id[i+1,id12] * w[id12] * (ddfun(a[i] - eta[id12]) / (pfun(a[i+1] - eta[id12]) - pfun(a[i] - eta[id12])) + (dfun(a[i] - eta[id12]) / (pfun(a[i+1] - eta[id12]) - pfun(a[i] - eta[id12])))^2))
			wt[i,i+1] = wt[i+1,i] = -sum(y_id[i+1,id12] * w[id12] * (dfun(a[i+1] - eta[id12]) * dfun(a[i] - eta[id12])) / (pfun(a[i+1] - eta[id12]) - pfun(a[i] - eta[id12]))^2)
		}
	}
	return (wt)
}

####
#estimate the probability of the latent variable lying
#between cut-points
####
#predict.cgam.polr = function(object,...) {
#fitted.cgam.polr = function(object) {	#
#	a = object$zeta
#	eta = object$muhat
#	lev = object$lev
#	nc = length(a) + 1
#	n = length(eta)
#	ps = matrix(0, nrow = nc, ncol = n)
#	for (i in 1:nc) {
#		if (i == 1) {
#			ps[i,] = pfun(a[1] - eta)
#		} else if (i == nc) {
#			ps[i,] = 1 - pfun(a[nc-1] - eta)
#		} else {
#			ps[i,] = pfun(a[i] - eta) - pfun(a[i-1] - eta)
#		}
#	}
#	ps = t(ps)
#	dimnames(ps) = list(1:n, lev)
#	return (ps)
#}

#####
Ord = function(link = "identity") {
	linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("identity"))
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(linktemp, " link not available for ordered categorical family; available links are \"identity\"")
    }
	structure(list(family = "ordered", link = linktemp), class = "family")
}

#########################################
#subroutines for confidence interval
#########################################
makebin = function(x) {
	k = length(x)
	r = 0
	for (i in 1:k) {
		r = r + x[k-i+1]*2^(i-1)
	}
	r
}

getvars = function(num) {
	i = num
	digits = 0
	power = 0
	while (digits == 0) {
		if (i < 2^power) {
			digits = power
		}
		power = power+1
	}
	binry = 1:digits*0
	if (num > 0) {
		binry[1] = 1
	}
	i = i - 2^(digits - 1)
	power = digits - 2
	for (p in power:0) {
		if (i >= 2^p) {
			i = i - 2^p
			binry[digits-p] = 1
		}
	}
	binry
}

getbin = function(num, capl) {
	br = getvars(num-1)
	digits = length(br)
	binrep = 1:capl*0
	binrep[(capl-digits+1):capl] = br
	binrep
}


########################################################
#cgam.pv:get p-values for smooth-constrained components#
#when there is only one x and one z, and test for x
#then it's not a sub-cone test and we don't use amat0 and coneA
#include mixed-effect for gaussian: weights is uinv_mat
###############################################################
cgam.pv = function(y, xmat, zmat, shapes, delta=NULL, np=NULL, capms=NULL,
                   numknotsuse=NULL, varlist=NULL, family=gaussian(),
                   weights=NULL, test_id=1, nsims=1000, skip=TRUE) {
  n = length(y)
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

  capl = length(xmat) / n
  if (capl < 1) {capl = 0}
  capk = length(zmat) / n
  if (capk < 1) {capk = 0}
  capls = sum(shapes == 17)
  one = 1:n*0 + 1
  if (!skip) {
    #numknots = c(0,0,0)
    #knots = list(); for (i in 1:3) knots[[i]] = 0
    #space = rep('E',3)
    numknots <- rep(0, capl)
    knots <- list(); for (i in 1:capl) knots[[i]] <- 0
    space <- rep('E',capl)

    delta <- NULL
    varlist <- NULL
    xid1 <- NULL; xid2 <- NULL; xpos2 <- 0
    knotsuse <- list(); numknotsuse <- NULL
    mslst <- list()
    #new:
    capm <- 0
    capms <- 0

    del1_ans <- makedelta(xmat[, 1], shapes[1], numknots[1], knots[[1]], space = space[1])
    del1 <- del1_ans$amat
    knotsuse[[1]] <- del1_ans$knots
    mslst[[1]] <- del1_ans$ms
    if(shapes[1] >= 9 & shapes[1] <= 17) {
      numknotsuse <- c(numknotsuse, length(del1_ans$knots))
    } else{
      numknotsuse <- c(numknotsuse, nrow(del1))
    }
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
        if(shapes[i] >= 9 & shapes[i] <= 17) {
          numknotsuse <- c(numknotsuse, length(del2_ans$knots))
        } else{
          numknotsuse <- c(numknotsuse, nrow(del2))
        }
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
  }
  umbrella.delta = tree.delta = NULL
  #ignore ut for now
  delta_ut = NULL; caput = 0
  if (!is.null(umbrella.delta) | !is.null(tree.delta)) {
    delta_ut = rbind(umbrella.delta, tree.delta)
    caput = nrow(delta_ut)
  }
  #del1 is edges for the component being tested; del0 is the cone with del1 removed
  del1 = t(delta[varlist == test_id, ])
  if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk > 0) {
    del0 = cbind(t(delta[varlist != test_id, ]), one, zmat, xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13])
    del = cbind(t(delta), one, zmat, xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13])
    np = 1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)  + capms
  } else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk == 0) {
    del0 = cbind(t(delta[varlist != test_id, ]), one, xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13])
    del = cbind(t(delta), one, xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13])
    np = 1 + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) + capms
  } else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk > 0) {
    del0 = cbind(t(delta[varlist != test_id, ]), one, zmat)
    del = cbind(t(delta), one, zmat)
    np = 1 + capk + capms
  } else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk == 0) {
    del0 = cbind(t(delta[varlist != test_id, ]), one)
    del = cbind(t(delta), one)
    np = 1 + capms
  } else {
    print ("error in capk, shapes!")
  }

  m0 = dim(del0)[2]
  m = dim(del)[2]

  #nkts is the col number of edges for each component
  nkts = numknotsuse
  id_add = which(shapes >= 13 & shapes <= 17)
  if (length(id_add) >= 1) {
    nkts[id_add] = nkts[id_add] + 1
  }
  nr = sum(nkts)
  nc = nr + np
  amat = matrix(0, nrow=nr, ncol=nc)
  for(i in 1:nr){amat[i,i]=1}

  #if (capl>1) {
  nr0 = sum(nkts) - nkts[test_id]
  #} else {nr0 = sum(nkts)}
  #nc0 = nr0 + np
  #debugged
  if (any(shapes == 17)) {
    nc0 = nr0 + np - capms
  } else {
    nc0 = nr0 + np
  }
  amat0 = matrix(0, nrow=nr0, ncol=nc0)

  #when there is only one component, flat vs s.incr
  if (nr0 > 0) {
    for(i in 1:nr0){amat0[i,i]=1}
  } #else {'Use shapereg in coneproj'}
  wt.iter = FALSE
  if (family$family %in% c("binomial", "poisson")) {
    wt.iter = TRUE
  }
  if (!wt.iter) {
    if (is.null(weights)) {
      weights = 1:n*0 + 1
    }
    w = weights
    #new: to check if w is a matrix or not to include the mixed-effect case
    ### null hyp fit
    if(!is.matrix(w)){
      ytil = y*sqrt(w)
      deltil = del0
      for(i in 1:m0){deltil[,i] = deltil[,i]*sqrt(w)}
    } else {
      #w is uinv_mat
      ytil = w %*% y
      deltil = w %*% del0
    }
    #print (m0)
    #print (np)
    #print (head(del0))
    if (m0 > np) {
      umat = chol(crossprod(deltil))
      uinv = solve(umat)
      ytiltil = t(uinv) %*% crossprod(deltil, ytil)
      atil = amat0 %*% uinv

      ans = coneA(ytiltil, atil)
      chat = uinv %*% ans$thetahat
      mu1 = del0 %*% chat
      #d0 = cbind(del2, del3)
    } else {
      cvec = t(deltil) %*% ytil
      chat = solve(crossprod(deltil), cvec)
      mu1 = del0 %*% chat
    }

    if (m0 > np) {
      d0 = t(delta[varlist != test_id, ])
      pr0 = -d0 %*% solve(crossprod(d0), t(d0))
      for(i in 1:n){pr0[i,i] = 1+pr0[i,i]}
      d1tr = pr0 %*% del1
      #} else {
      #d1tr = del1
      #}
      ##  do weighted projection of xi=z, weights w, onto d1tr cone
      d1trw = d1tr
      nk1 = nkts[test_id]
      #new
      if(!is.matrix(w)){
        for(i in 1:nk1){d1trw[,i] = d1tr[,i]*sqrt(w)}
      } else {
        d1trw = w %*% d1trw
      }
      ans1 = coneB(ytil, d1trw)
      num = sum((d1trw%*%ans1$coef)^2)

      ## do weighted projection of xi=z, weights w, onto big cone
      delw = del
      #new
      if(!is.matrix(w)){
        for(i in 1:m){delw[,i] = del[,i]*sqrt(w)}
      } else {
        delw = w %*% delw
      }

      ans2 = coneB(ytil,delw[,1:(m-np)],delw[,(m-np+1):m])
      coef0 = ans2$coef
      coef = coef0
      coef[1:(m-np)] = coef0[(np+1):m]
      coef[(m-np+1):m] = coef0[1:np]
      ss2 = sum((ytil-delw%*%coef)^2)
      #ss2 = sum((ytil-delw[,1:(m-np)]%*%ans2$coef[(np+1):m]-delw[,(m-np+1):m]*ans2$coef[1:np])^2)
      df2 = n-sum(abs(ans2$coef)>1e-8)
      #df2 = n - m
      bstat = num/(num+ss2)
      ## get mixing distn
      mdist = 0:nk1*0
      for(isims in 1:nsims){
        ysim = rnorm(n)
        anss = coneB(ysim, d1trw)
        d = sum(anss$coef>1e-8)
        mdist[d+1] = mdist[d+1]+1
      }
      mdist = mdist/nsims

      pval = mdist[1]
      for(i in 1:nk1){
        pval = pval+pbeta(bstat,i/2,df2/2)*mdist[i+1]
      }
      pv = 1-pval
    } else {
      #print (dim(del))
      del2 = del
      #new
      if(!is.matrix(w)){
        for(i in 1:m){del2[,i] = del2[,i] * sqrt(w)}
      } else {
        del2 = w %*% del2
      }
      nk = ncol(del1)
      ans2 = coneB(ytil, del2[,1:nk], deltil)
      chat2 = ans2$coefs
      #nd = np?
      nd = ncol(deltil)
      that1 = deltil%*%chat2[1:nd]+del2[,1:nk]%*%chat2[(nd+1):m]
      ss1 = sum((ytil-that1)^2)
      ss0 = sum((ytil-deltil%*%chat)^2)

      bstat = (ss0-ss1)/ss0
      ## get mixing distribution
      mdist = 0:nk*0
      keepd = 1:nsims
      for(isims in 1:nsims){
        ysim = rnorm(n)
        anss = coneB(ysim, del2[,1:nk], deltil)
        d = sum(abs(round(anss$coef,8))>0)
        #mdist[d-2] = mdist[d-2]+1
        mdist[d+1-np] = mdist[d+1-np]+1
        #mdist[d] = mdist[d]+1
        keepd[isims] = d
        if(d < 3){zkeep=ysim}
      }
      mdist = mdist/nsims

      pval = mdist[1]
      for(i in 1:nk) {
        pval = pval + pbeta(bstat,i/2,(n-np-i)/2)*mdist[i+1]
      }
      pv = 1-pval
    }
  } else {
    #print (head(del0))
    ### null hyp fit
    eta0 = 1:n*0
    z = 1:n
    #mu0 = muhat.fun(eta0, fml=family$family)
    mu0 = 1:n*0 + 1/2
    mu1 = mu0
    if (is.null(weights)) {
      weights = 1:n*0 + 1
    }
    prior.w = weights
    end = FALSE; nrep = 0; face = NULL
    while (!end & nrep < 100) {
      nrep = nrep+1
      if (family$famil == "binomial") {
        ch = mu0>1e-5 & mu0<1-1e-5
      }
      if (family$famil == "poisson") {
        ch = mu0>1e-5
      }
      z[ch] = eta0[ch] + (y[ch]-mu0[ch]) * deriv.fun(mu0[ch], fml = family$family)
      z[!ch] = eta0[!ch]
      w = as.vector(prior.w * (deriv.fun(mu0, fml = family$family))^(-1))
      #w = mu0*(1-mu0)
      #to avoid deltil to be singular; test!
      w[which(round(w, 7) == 0)] = 1e-7
      ztil = z*sqrt(w)
      deltil = del0
      #print (m0)
      #print (any(round(w, 7) == 0))
      #print (w[round(w, 7) == 0])
      for(i in 1:m0){deltil[,i] = deltil[,i]*sqrt(w)}
      if (m0 > np) {
        #print ('True')
        #print (qr(crossprod(deltil))$rank)
        umat = chol(crossprod(deltil))
        uinv = solve(umat)
        ztiltil = t(uinv) %*% crossprod(deltil, ztil)
        atil = amat0 %*% uinv
        #if (nrep > 1) {
        ans = coneA(ztiltil, atil, face=face)
        #} else {
        #    ans = coneA(ztiltil, atil)
        #    face = ans$face
        #if (nrep == 107 | nrep == 108) {
        #    print (face)
        #}
        #}
        chat = uinv %*% ans$thetahat
      } else {
        #cvec = t(deltil) %*% ztil
        #chat = solve(crossprod(deltil), cvec)
        chat = solve(t(deltil)%*%deltil)%*%t(deltil)%*%ztil
      }
      eta1 = del0 %*% chat
      mu1 = muhat.fun(eta1, fml=family$family)
      #tb=eta1>50
      #mu1=eta1
      #mu1[!tb]=exp(eta1[!tb])/(1+exp(eta1[!tb]))
      #mu1[tb]=1
      ch2 = mu1 > 2*max(y)
      if (length(ch2) >= 1) {
        mu1[ch2] = 2*max(y)
      }
      dist = sqrt(sum(mu0-mu1)^2/sum(mu0^2))
      #print (dist)
      if (dist < 1e-5){end = TRUE} else {eta0=eta1;mu0=mu1}
    }
    #print (head(mu1,10))
    #if m0 > np, del0 is a cone containing a linear space
    if (m0 > np) {
      d0 = t(delta[varlist != test_id, ])
      #pr0 = -d0%*%solve(crossprod(d0), t(d0))
      #for(i in 1:n){pr0[i,i] = 1+pr0[i,i]}
      #d1tr = pr0%*%del1
      ##  do weighted projection of xi=z, weights w, onto d1tr cone
      #d1trw = d1tr
      #nk1 = nkts[test_id]
      #for(i in 1:nk1){d1trw[,i] = d1tr[,i]*sqrt(w)}
      #to make sure edges are orthogonal:
      m0a = dim(d0)[2]
      d0w = d0
      for(i in 1:m0a){d0w[,i] = d0[,i]*sqrt(w)}
      pr0 = -d0w %*% solve(t(d0w) %*% d0w) %*% t(d0w)
      for(i in 1:n){pr0[i,i]=1+pr0[i,i]}
      #}
      del1w = del1
      nk1 = nkts[test_id]
      for(i in 1:nk1){del1w[,i] = del1[,i]*sqrt(w)}
      #if (capl > 1) {
      d1trw = pr0%*%del1w
      #} else {
      #d1trw = del1w
      ans1 = coneB(ztil, d1trw)
      num = sum((d1trw %*% ans1$coef)^2)

      ## do weighted projection of xi=z, weights w, onto big cone
      delw = del
      for(i in 1:m){delw[,i] = del[,i]*sqrt(w)}

      ans2 = coneB(ztil,delw[,1:(m-np)],delw[,(m-np+1):m])
      coef0 = ans2$coef
      coef = coef0
      coef[1:(m-np)] = coef0[(np+1):m]
      coef[(m-np+1):m] = coef0[1:np]
      #ss2 = sum((ztil-delw[,1:(m-np)]%*%ans2$coef[(np+1):m]-delw[,(m-np+1):m]*ans2$coef[1:np])^2)
      ss2 = sum((ztil-delw%*%coef)^2)
      #df2=n-sum(abs(ans2$coef)>1e-8)
      df2 = n - m
      bstat = num/(num+ss2)
      #print (bstat)
      #print (dim(amat))
      #print (dim(amat0))
      mdist = 0:nk1*0
      for(isims in 1:nsims){
        zsim = rnorm(n)
        anss = coneB(zsim, d1trw)
        d = sum(anss$coef>1e-8)
        mdist[d+1] = mdist[d+1]+1
      }
      mdist = mdist/nsims

      pval = mdist[1]
      for(i in 1:nk1){
        pval = pval+pbeta(bstat,i/2,df2/2)*mdist[i+1]
      }
      pv = 1-pval
    } else {
      #print (head(del))
      #only one x and >= one z; del2 is a copy of del
      del2 = del
      #print (head(del))
      #print (head(w))
      for(i in 1:m){del2[,i] = del2[,i]*sqrt(w)}
      nk = ncol(del1)
      #print (nk)
      #print (head(del2))
      #print (nk)
      #print (head(deltil))
      #print (nrep)
      #print (head(ztil))
      #print (head(del2))
      ans2 = coneB(ztil, del2[,1:nk], deltil)
      chat2 = ans2$coefs
      #print (chat2)
      #nd = np?
      nd = ncol(deltil)
      that1 = deltil%*%chat2[1:nd]+del2[,1:nk]%*%chat2[(nd+1):m]

      ss1 = sum((ztil - that1)^2)
      ss0 = sum((ztil - deltil %*% chat)^2)
      #print (ss0)
      #print (ss1)
      bstat = (ss0 - ss1)/ss0
      #print (bstat)
      ## get mixing distribution
      #set.seed(123)
      mdist = 0:nk*0
      keepd = 1:nsims
      for(isims in 1:nsims){
        zsim = rnorm(n)
        ans = coneB(zsim, del2[,1:nk], deltil)
        d = sum(abs(round(ans$coef,8))>0)
        #mdist[d-2] = mdist[d-2]+1
        mdist[d+1-np] = mdist[d+1-np]+1
        #mdist[d] = mdist[d]+1
        keepd[isims] = d
        if(d < 3){zkeep=zsim}
      }
      mdist = mdist/nsims

      pval = mdist[1]
      for(i in 1:nk) {
        pval = pval + pbeta(bstat,i/2,(n-np-i)/2)*mdist[i+1]
      }
      pv = 1-pval
      #print (pv)
    }
  }
  ans = list(pv = pv, edf = sum(abs(ans2$coef)>1e-8), coef=ans2$coef, bstat=bstat)
  return (ans)
}

############################################################
#cgam.pvz is to test categorical predictors, not every level
#for mixed-effect model, use uinv_mat as w
############################################################
cgam.pvz = function(y, bigmat, df_obs, sse1 = NULL, np = 1, zid = 1, zid1 = 1, zid2 = 1, muhat = NULL,
                    etahat = NULL, coefskeep = NULL, wt.iter=FALSE, family=gaussian(), weights=NULL,
                    uinv_mat=NULL) {
  n = length(y)
  cicfamily = CicFamily(family)
  llh.fun = cicfamily$llh.fun
  #new: use log link in gamma
  linkfun = cicfamily$linkfun
  etahat.fun = cicfamily$etahat.fun
  gr.fun = cicfamily$gr.fun
  wt.fun = cicfamily$wt.fun
  zvec.fun = cicfamily$zvec.fun
  muhat.fun = cicfamily$muhat.fun
  ysim.fun = cicfamily$ysim.fun
  deriv.fun = cicfamily$deriv.fun
  dev.fun = cicfamily$dev.fun

  m = nrow(bigmat)
  df_full = min(m, 1.2*df_obs)
  if (is.null(weights)) {
    weights = 1:n*0 + 1
  }
  prior.w = weights
  #n = ncol(bigmat)
  lz = length(zid)
  #we only need to get ssers in the following code
  pvs = ssers = ssefs = fstats = edfs = rep(0, lz)
  #sse_f = NULL
  for(iz in 1:lz){
    pos1 = zid1[iz]
    pos2 = zid2[iz]
    bigmat0 = bigmat[-c((pos1+1):(pos2+1)), ,drop = FALSE]
    #np0 is the # of columns remaining in the vmat; it's for H_0
    np0 = np - (pos2-pos1+1)
    #lvs is the # of levels in a z, including the reference level
    lvs = (pos2-pos1+1)
    edfs[iz] = lvs
    #initialize
    cvec = NULL
    etahat = NULL
    if (wt.iter) {
      etahat = etahat.fun(n, y, fml = family$family)
      gr = gr.fun(y, etahat, weights, fml = family$family)
      wt = wt.fun(y, etahat, n, weights, fml = family$family)
      cvec = wt * etahat - gr
    } else {wt = wt.fun(y, etahat, n, weights, fml = family$family)}

    #new:
    if(is.null(uinv_mat)){
      zvec = zvec.fun(cvec, wt, y, fml = family$family)
      gmat = t(bigmat0)
      for (i in 1:n) {gmat[i,] = bigmat0[,i] * sqrt(wt[i])}
    } else {
      #only used for gaussian
      zvec = uinv_mat %*% y
      gmat = t(bigmat0)
      gmat = uinv_mat %*% gmat
    }
    dsend = gmat[, -c(1:np0), drop = FALSE]
    zsend = gmat[,1:np0 , drop = FALSE]
    ans = coneB(zvec, dsend, zsend)
    face = ans$face
    etahat = t(bigmat0) %*% ans$coefs
    #coefs = ans$coefs
    #muhat = t(bigmat0) %*% coefs
    muhat = etahat
    if(is.null(uinv_mat)){
      sse_r = sum(prior.w * (y - muhat)^2)
    }else {
      sse_r = sum(uinv_mat %*% (y - muhat)^2)
    }
    sse_f = sse1
    fstati = (sse_r - sse_f) / lvs / sse_f*(n-df_full)
    #print (sse_r)
    #print (sse_f)
    #print (lvs)
    #print (n-df_full)
    #print (dim(bigmat0))
    #print (dim(bigmat))
    pvi = 1 - pf(fstati, lvs, n-df_full)
    if (wt.iter) {
      sm = 1e-7
      muhat = muhat.fun(etahat, fml = family$family)
      diff = 1
      if (family$family == "binomial") {
        mdiff = abs(max(muhat) - 1) > sm
      } else {mdiff = TRUE}
      nrep = 0
      ##########
      #iterate!#
      ##########
      while (diff > sm & mdiff & nrep < n^2){
        oldmu = muhat
        nrep = nrep + 1
        gr = gr.fun(y, etahat, weights, fml = family$family)
        wt = wt.fun(y, etahat, n, weights, fml = family$family)
        cvec = wt * etahat - gr
        #zvec <- cvec / sqrt(wt)
        zvec = zvec.fun(cvec, wt, y, fml = family$family)
        #gmat <- t(bigmat %*% sqrt(diag(wt)))
        gmat = t(bigmat0)
        for (i in 1:n) {gmat[i,] = bigmat0[,i] * sqrt(wt[i])}
        dsend = gmat[, -c(1:np0), drop = FALSE]
        zsend = gmat[,1:np0 , drop = FALSE]
        #ans <- coneB(zvec, t(dsend), zsend)
        ans = coneB(zvec, dsend, zsend, face = face)
        etahat = t(bigmat0) %*% ans$coefs
        muhat = muhat.fun(etahat, fml = family$family)
        diff = mean((muhat - oldmu)^2)
        mdiff = abs(max(muhat) - 1)
        if (family$family == "binomial") {
          mdiff = abs(max(muhat) - 1) > sm
        } else {mdiff = TRUE}
      }
      z = 1:n*0
      if (family$famil == "binomial") {
        ch = muhat>1e-7 & muhat<1-1e-7
      }
      if (family$famil == "poisson") {
        ch = muhat>1e-7
      }
      z[ch] = etahat[ch] + (y[ch]-muhat[ch]) * deriv.fun(muhat[ch], fml = family$family)
      z[!ch] = etahat[!ch]
      w = as.vector(prior.w * (deriv.fun(muhat, fml = family$family))^(-1))
      #to avoid deltil to be singular
      w[which(round(w, 7) == 0)] = 1e-7
      ztil = z*sqrt(w)
      deltil = t(bigmat0)
      m0 = nrow(bigmat0)
      for(i in 1:m0){deltil[,i] = deltil[,i]*sqrt(w)}
      sse_r = sum((ztil - deltil %*% ans$coefs)^2)

      deltil = t(bigmat)
      m = nrow(bigmat)
      for(i in 1:m){deltil[,i] = deltil[,i]*sqrt(w)}
      umat = chol(t(deltil)%*%deltil)
      #umat = chol(crossprod(deltil))
      uinv = solve(umat)
      ztiltil = t(uinv) %*% t(deltil) %*% ztil
      #amat = matrix(0, nrow = m, ncol = m)
      #nk = m - np
      #for(ik in 1:nk){amat[ik,ik] = 1}
      amat = diag(m-np)
      zerom = matrix(0, nrow=nrow(amat), ncol=np)
      amat = cbind(zerom, amat)
      atil = amat %*% uinv
      ans = coneA(ztiltil, atil)
      sse_f = sum((ztil - deltil %*% uinv %*% ans$thetahat)^2)
      df_full = min(m, 1.2*ans$df)
      fstati = (sse_r - sse_f) / lvs / sse_f*(n-df_full)
      pvi = 1 - pf(fstati, lvs, n-df_full)
    }
    ssers[iz] = sse_r
    ssefs[iz] = sse_f
    fstats[iz] = fstati
    #pvi = 1 - pf(fstati, lvs, n-df_full)
    pvs[iz] = pvi
  }
  rslt = list(pvs = pvs, ssers = ssers, ssefs = ssefs, fstats = fstats, edfs = edfs, mu1 = muhat)
  return (rslt)
}

#############
#anova.cgam#
############
anova.cgam <- function(object,...){
    family <- object$family
    call <- object$call
    pvs <- object$pvs
    pvsz <- object$pvsz
    capl <- object$capl
    capk <- object$capk
    s.edf <- object$s.edf
    z.edf <- object$z.edf
    bstats <- object$bstats
    fstats <- object$fstats
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

    #new:
    rslt1 <- pTerms.pv <- NULL
    if (!is.null(pvsz)) {
        rslt1 <- data.frame("df" = z.edf, "F-stat" = round(fstats, 4), "p.value" = round(pvsz, 4))
        lzid <- length(zid1)
        for (i in 1:lzid) {
            rownames(rslt1)[i] <- attributes(tms)$term.labels[zid[i] - 1]
        }
        rslt1 <- as.matrix(rslt1)
        pTerms.pv <- pvsz
        #rownames(pTerms.pv) <- rownames(rslt1)
    }
    rslt2 <- s.pv <- NULL
    if (!is.null(pvs)) {
        rslt2 <- data.frame("edf" = round(s.edf, 4), "mixture of Beta" = round(bstats, 4), "p.value" = round(pvs, 4))
        #rownames(rslt2) <- attributes(tms)$term.labels
        #debugged: check more
        if (!is.null(zid)) {
          if(inherits(object, "cgamm")){
            rownames(rslt2) <- rev(rev((attributes(tms)$term.labels)[-(zid-1)])[-1])
          } else {
            rownames(rslt2) <- (attributes(tms)$term.labels)[-(zid-1)]
          }
        } else {
            if(inherits(object, "cgamm")){
              rownames(rslt2) <- rev(rev((attributes(tms)$term.labels))[-1])
            } else{
              rownames(rslt2) <- (attributes(tms)$term.labels)
            }
        }
        s.pv <- pvs
        #rownames(s.pv) <- rownames(rslt2)
    }
    ans <- list(call = call, coefficients1 = rslt1, coefficients2 = rslt2, family = family, s.pv = s.pv, pTerms.pv = pTerms.pv)
    class(ans) <- "anova.cgam"
    return(ans)
}


print.anova.cgam <- function(x,...){
    print(x$family)
    cat("Formula:\n")
    print(x$call)
    if (!is.null(x$coefficients1)) {
        cat("\n")
        cat("Parametric terms:\n")
        printCoefmat(x$coefficients1, P.values = TRUE, has.Pvalue = TRUE)
    }
    #cat("\n")
    if (!is.null(x$coefficients2)) {
        cat("\n")
        cat("Approximate significance of smooth terms: \n")
        printCoefmat(x$coefficients2, P.values = TRUE, has.Pvalue = TRUE)
    }
}

###############################################
#subroutines for monotonic variance estimation
###############################################
varest = function(y, x, muhat=NULL, shape=9, var.knots=0, db.exp=FALSE){
    n = length(y)
    order.id = order(x)
    #new:rk to order var back
    #rk = rank(x)
    #x = x[order.id]
    xs = sort(x)
    y = y[order.id]

    ndegree = 2
    #if (length(muhat) == 1){
    #  muhat = rep(0,n)
    #}
    #if (length(muhat) == n){
    #  muhat = muhat
    #}
    if (n < 20) {
        print("ERROR: must have at least 20 observations")
    }
    if (length(x) != length(y)) {
        print("ERROR: length of x must be length of y")
    }
    if (length(var.knots) > 1) {
        var.kint = var.knots[-c(1, length(var.knots))]
    } else {
        br = c(30, 100, 200, 400, 700, 1000, 1e+10)
        obs = 1:7
        var.nk = min(obs[n < br]) + 2
        var.knots = 0:(var.nk - 1)/(var.nk - 1) * (max(x) - min(x)) + min(x)
        var.kint = var.knots[-c(1, length(var.knots))]
    }
    nknots = length(var.kint)
    nk = nknots + 2 + 1

    #initial coefficients
    theta0 = seq(1, 0.5, length = nk)
    dif.theta = 1
    amat = matrix(0, ncol = nk, nrow = nk)

    #if (increasing == 1){
    if (shape == 9) {
        for (i in 1:(nk - 1)) {
            amat[i, i] = 1
            amat[i, i + 1] = -1
        }
        amat[nk, nk] = 1
    }

    #if (increasing == 2){
    if (shape == 10) {
        for (i in 1:(nk - 1)) {
            amat[i, i] = -1
            amat[i, i + 1] = 1
        }
        amat[nk, 1] = 1
    }

    bvec = c(rep(0,nk-1), 10^-10)
    Bint = bs(xs, knots = var.kint, degree = 2, intercept = T)

    if (!is.null(muhat)) {
        r = y - muhat
    } else {
        r = y
    }
    sm = 1e-6
    nrep = 0
    while (dif.theta > sm) {
        nrep = nrep + 1
        dl = db.penal(theta0, Bint, r, db.exp)
        d2l = ddb.penal(theta0, Bint)
        res = qprog(d2l, c(t(theta0) %*% d2l) - dl, amat, bvec)
        #if (class(res) == "try-error") {
        #    break
        #}
        theta1.new = res$thetahat
        var1 = 1/(Bint %*% theta1.new)
        var0 = 1/(Bint %*% theta0)
        dif.theta = sum(abs(var1 - var0))/sum(abs(var0))
        theta0 = theta1.new + 1e-10
        #print (dif.theta)
    }
    res0 = theta1.new

    if (!db.exp) {
        vhat = 1/Bint %*% theta1.new
    } else if (db.exp) {
        vhat = (1/Bint %*% theta1.new)^2
    }
    #new
    #print (vhat)
    #vhat = vhat[rk]
    #print (vhat)
    vhat=vhat[rank(x)]
    res = list(vhat = vhat, muhat = muhat, x = x, y = y, var.knots = var.knots)
    return(res)
}

#functions used in the varest function for normal errors
db.penal <- function(betas,B,r, db.exp=FALSE){
    nk= ncol(B)
    dl=vector()
    if (!db.exp) {
        for (i in 1:nk){
            dl[i]=sum(-(1/(B%*%betas))*B[,i] + (r^2)*B[,i])
        }
    } else {
        for (i in 1:nk){
            dl[i]=sum(-(1/(B%*%betas))*B[,i] + sqrt(2)*abs(r)*B[,i])
        }
    }
    return(dl)
}

ddb.penal <- function(betas,B){
    nk=ncol(B)
    d2l=diag(0,nk)
    for (i in 1:nk){
        for (j in 1:nk){
            d2l[i,j]= sum((1/((B%*%betas)^2))*(B[,i]*B[,j]))
        }
    }
    return(d2l)
}


ddb.penal_dexp <- function(betas,B){
    nk=ncol(B)
    d2l=diag(0,nk)
    for (i in 1:nk){
        for (j in 1:nk){
            d2l[i,j] = sum((1/((B%*%betas)^2))*(B[,i]*B[,j]))
        }
    }
    return(d2l)
}


