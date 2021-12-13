cgamm = function(formula, nsim = 0, family = gaussian(), cpar = 1.2, data = NULL, weights = NULL, sc_x = FALSE, sc_y = FALSE, bisect = TRUE, reml = TRUE) {
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0L)
  	mf <- mf[c(1L, m)]
    mf[[1]] <- quote(lme4::lFormula)
    #print (mf)
    mf <- eval(mf, parent.frame(1L))$fr
    #print (mf)
#the following is the same as cgam
#print (names(mf))
    ynm <- names(mf)[1]
    mt <- attr(mf, "terms")
    #print (names(mt))
    y <- model.response(mf, "any")
    #print (mt)
    #print (ynm)
    #print (head(y))
    #stop (print (head(mf)))
    #save(mf, file='mf.Rda')
    #stop ('check')
    shapes1 <- NULL; shapes2 <- NULL
    xmat <- NULL; xnms <- NULL
    tr <- NULL; pl <- NULL; umb <- NULL
    tree.delta <- NULL; umbrella.delta <- NULL
    tid1 <- NULL; tid2 <- NULL; tpos2 <- 0
    uid1 <- NULL; uid2 <- NULL; upos2 <- 0
    nums <- NULL; ks <- list(); sps <- NULL; xid <- 1
    zmat <- NULL; zid <- NULL; zid0 <- NULL; zid1 <- NULL; zid2 <- NULL; znms <- NULL; is_param <- NULL; is_fac <- NULL; vals <- NULL; st <- 1; ed <- 1
    ztb <- list(); iztb <- 1
#new: exclude the column for id
    #nc <- ncol(mf) - 1
    nc <- ncol(mf)
    #print (head(mf))
#test:
    id <- mf[,nc]
    #print (nc)
    #print (id)
    szs <- unname(table(id))
    for (i in 2:(nc-1)) {
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
    #if (family$family == "binomial" | family$family == "poisson") {
    #    wt.iter = TRUE
    #} else {wt.iter = FALSE}
    wt.iter = FALSE
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
    #print (bisect)
    ans <- cgamm.fit(y = y, xmat = xmat0, zmat = zmat, id = id, shapes = shapes0, numknots = nums0, knots = ks0, space = sps0, nsim = nsim, family = family, cpar = cpar, wt.iter = wt.iter, umbrella.delta = umbrella.delta, tree.delta = tree.delta, weights = weights, sc_x = sc_x, sc_y = sc_y, idx_s = idx_s, idx = idx, bisect = bisect, reml = reml)
    rslt <- list(muhat = ans$muhat, coefs = ans$coefs, bh = ans$bh, zcoefs = ans$zcoefs, pvals.beta = ans$pvals.beta, se.beta = ans$se.beta, vcoefs = ans$vcoefs, ahat = ans$ahat, sig2hat = ans$sig2hat, siga2hat = ans$siga2hat, thhat = ans$thhat, bigmat = ans$bigmat, gtil=ans$gtil, dd2=ans$dd2, id = id, szs = szs, shapes = shapes0, numknots = ans$numknots, knots = ans$knots, space = sps0, d0 = ans$np, xmat_add = xmat, xmat0 = ans$xmat2, knots0 = ans$knots2, numknots0 = ans$numknots2, sps0 = ans$sps2, ms0 = ans$ms2, etacomps = ans$etacomps, xnms_add = xnms, xid1 = ans$xid1, xid2 = ans$xid2, ynm = ynm, y = y, znms = znms, zmat = zmat, ztb = ztb, zid = zid, zid1 = zid1, zid2 = zid2, vals = vals, family = family, is_fac = is_fac, is_param = is_param, tms = mt, capm = ans$capm, capms = ans$capms, capk = ans$capk, capt = ans$capt, capu = ans$capu, mod.lmer = ans$mod.lmer, pv.siga2 = ans$pv.siga2, ci.siga2 = ans$ci.siga2, ci.siga2.bi = ans$ci.siga2.bi, ci.th = ans$ci.th, ci.rho = ans$ci.rho, ci.sig2 = ans$ci.sig2, ones = ans$ones, resid_df_obs = ans$resid_df_obs, edf = ans$edf)
    rslt$call <- cl
    class(rslt) <- c("cgamm", "cgam")
    return (rslt)
}

############
#cgamm.fit
############
cgamm.fit = function(y, xmat, zmat, id, shapes, numknots, knots, space, nsim, family = gaussian(), cpar = 1.2, wt.iter = FALSE, umbrella.delta = NULL, tree.delta = NULL, weights = NULL, sc_x = FALSE, sc_y = FALSE, idx_s = NULL, idx = NULL, bisect = FALSE, modlmer = FALSE, reml = TRUE) {
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
	
	n = length(y)
#new for random-effect
	szs = unname(table(id))
    #print (id)
    ncl = length(szs)
	balanced = FALSE
	if (length(unique(szs)) == 1) {balanced = TRUE}
	ycl = f_ecl(y, ncl, szs)
	sm = 1e-7 
	capl = length(xmat) / n
	if (capl < 1) {capl = 0}
	if (round(capl, 8) != round(capl, 1)) {stop ("Incompatible dimensions for xmat!")}
#new:
	if (capl > 0 & sc_x) {
		for (i in 1:capl) {xmat[,i] = (xmat[,i] - min(xmat[,i])) / (max(xmat[,i]) - min(xmat[,i]))}
	}
#new:
	if (sc_y) {
		sc = sd(y)
		y = y / sc	
	}
	capk = length(zmat) / n
    #print (head(zmat,20))
	if (capk < 1) {capk = 0}
	if (round(capk, 8) != round(capk, 1)) {stop ("Incompatible dimensions for zmat!")}
#new:
	capls = sum(shapes == 17)
####################################################
#get basis functions for the constrained components#
####################################################	
	delta = NULL
	varlist = NULL
	xid1 = NULL; xid2 = NULL; xpos2 = 0
	knotsuse = list(); numknotsuse = NULL
	mslst = list()
#new:
	capm = 0
	capms = 0
    #if (capl - capls > 0) {
#test: no need for capk + capls > 0
    if (capl > 0) {
        del1_ans = makedelta(xmat[, 1], shapes[1], numknots[1], knots[[1]], space = space[1])
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
	        	del2_ans = makedelta(xmat[,i], shapes[i], numknots[i], knots[[i]], space = space[i])
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
        xvec = NULL
		if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk > 0) {
            xvec = t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13])
			bigmat = rbind(1:n*0 + 1, t(zmat), xvec, delta)
			np = 1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)  + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk == 0) {
            xvec = t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13])
			bigmat = rbind(1:n*0 + 1, xvec, delta)
			np = 1 + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk > 0) {
			bigmat = rbind(1:n*0 + 1, t(zmat), delta)
			np = 1 + capk + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) == 0 & capk == 0) {
			bigmat = rbind(1:n*0 + 1, delta)
			np = 1 + capms
		} else {
			print ("error in capk, shapes!")
		} 
#new: capm is the number of columns of edges for constrained x's
		capm = length(delta) / n - capms
	}
	if (!is.null(umbrella.delta)) {
		bigmat = rbind(bigmat, umbrella.delta)
		capu = length(umbrella.delta) / n
	} else {capu = 0}
	if (!is.null(tree.delta)) {
		bigmat = rbind(bigmat, tree.delta)
		capt = length(tree.delta) / n
	} else {capt = 0}
	if (!is.null(umbrella.delta) | !is.null(tree.delta)) 
		delta_ut = rbind(umbrella.delta, tree.delta)
########	
#step 0
########
	if (!wt.iter) {
        if (is.null(weights)) {
            weights = 1:n*0 + 1
        }
        wt = weights
#new:
        zvec = wt^(1/2) * y
        gmat = t(bigmat)
        for (i in 1:n) {gmat[i,] = bigmat[,i] * sqrt(wt[i])}
        if (any(shapes != 17)) {
            dsend = gmat[, (np + 1):(np + capm + capu + capt), drop = FALSE]
            zsend = gmat[, 1:np, drop = FALSE]
            #ans = coneB(y, dsend, zsend)
            ans = coneB(zvec, dsend, zsend)
            edf = ans$df
            face = ans$face
            bh = coef(ans)
            #test:
            if (any(round(bh[1:np],6) < 0)) {
                pos = (1:np)[which(round(bh[1:np],6) < 0)]
                face = unique(c(pos, face))
            }
        } else {
            bh = solve(crossprod(gmat), t(gmat)) %*% zvec
            edf = nrow(bigmat)
            face = 1:edf
        }
#new: reml; find the edges with positive coefficients
#print (bh)
#print (edf)
#print (face)
        xtx = xtx2 = NULL
        dd = t(bigmat[face, ,drop = FALSE])
        #if (balanced) {
        #    xtx = t(dd)%*%dd
        #}
        #print (dim(dd))
        xms = ones = list()
        st = 1
        ed = 0
        for (icl in 1:ncl) {
            sz = szs[icl]
            ed = ed + sz
            xms[[icl]] = dd[st:ed, ,drop=F]
            onevec = 1:sz*0+1
            onemat = onevec%*%t(onevec)
            ones[[icl]] = onemat
            st = ed + 1
        }
        #print (dim(ones[[96]]))
        #if (balanced) {
        #     oneMat = as.matrix(bdiag(ones))
        #    xtx2 = crossprod(dd, oneMat) %*% dd
        #}
		muhat = t(bigmat) %*% bh
		oldmu = muhat
        #evec = y - muhat
        #ecl = f_ecl(evec, ncl, szs)
        #ebars = sapply(ecl, mean)
        #ansi = try(ansi0<-uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf, type='ub', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), silent=TRUE)
        #thhat = ansi$root
        #sig2hat = fsig(thhat, szs, ecl, ncl, N=n, edf=edf, D=nrow(bigmat), type=type)
        #siga2hat = sig2hat * thhat
        #ahat = ebars*szs*thhat/(1+szs*thhat)

		diff = 10
		nrep = 0
		while (diff > 1e-7 & nrep < 10) {
##########
#step 1: update covmat
##########
            nrep = nrep + 1
			evec = y - muhat
			ecl = f_ecl(evec, ncl, szs)
            mod.lmer = NULL
            if (!balanced) {
                #print (n)
                #ansi = try(uniroot(fth2, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n))
                if (modlmer) {
                    mod.lmer = lmer(evec~-1+(1|id), REML=reml)
                    thhat = summary(mod.lmer)$optinfo$val^2
                } else {
                    #ansi = try(uniroot(fth2, c(-10, 1e+3), szs=szs, ycl=ecl, N=n))
                    if (reml) {
                      ansi = try(ansi0<-uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf, type='ub', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), silent=TRUE)
                      #ansia = try(ansi0a<-uniroot(fsigarm, c(1e-10, 1e+3), sige = sig2hat^.5, szs=szs, ycl=ecl, N=n, xcl=xms, p=edf, type='ub', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), silent=TRUE)
                    } else {
                      ansi = try(ansi0<-uniroot(fth2, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n), silent=TRUE)
                    }
                    if (class(ansi) == "try-error") {
                        thhat = 0
                    } else {
                        thhat = ansi$root
                    }
                }
                type = "ub"
            } else {
                if (modlmer) {
                    mod.lmer = lmer(evec~-1+(1|id), REML=reml)
                    thhat = summary(mod.lmer)$optinfo$val^2
                } else {
                    #sometimes return a negative number; use uniroot instead.
                    #thhat = fth(ecl, ncl, N=n)
                    #ansi = try(uniroot(fth2, c(-10, 1e+3), szs=szs, ycl=ecl, N=n))
                    if (reml) {
                        #print (edf)
                        #print (dim(xms[[1]]))
                        #ansi = try(uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf,  type='b', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones))
                        #ansi = tryCatch(uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf,  type='b', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), error=return(0))
                        #print (class(ansi))
                        #print (ansi)
                        ansi = try(ansi0<-uniroot(fth2rm, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n, xcl=xms, p=edf,  type='b', xtx=xtx, xtx2=xtx2, xmat_face=dd, ones=ones), silent=TRUE)
                    } else {
                        ansi = try(ansi0<-uniroot(fth2, c(1e-10, 1e+3), szs=szs, ycl=ecl, N=n), silent=TRUE)
                    }
                    if (class(ansi) == "try-error") {
                        thhat = 0
                    } else {
                        thhat = ansi$root
                    }
                }
                type = "b"
            }
            #mod.lmer = lmer(evec~-1+(1|id), REML=F)
            #modsm = summary(mod.lmer)
            #vs = modsm$var
            #nvs = names(vs)
            #siga2hat = (attributes(vs$nvs)$stddev)^2
            #sig2hat = modsm$sigma^2
            #thhat = siga2hat/sig2hat
##################
#step 2: update muhat
##################
			ytil = NULL 
#gtil is edges
			gtil = NULL
			st = 1
			ed = 0
            if (balanced) {
                oneMat = ones[[1]]
                sz = szs[1]
            } else {
                sz = max(szs)
                pos = which(szs == sz)[1]
                oneMat = ones[[pos]]
            }
            #print (sz)
            #print (dim(onemat))
            vi = diag(sz) + oneMat*thhat
            covi = vi
            umat = t(chol(covi))
            uinv = solve(umat)
            #uinv0 is used for unbalanced
            uinv0 = uinv
            #print (szs)
			for (icl in 1:ncl) {
				sz = szs[icl]
                #one = matrix(rep(1, sz), ncol=1)
                #onemat = tcrossprod(one)
                if (!balanced) {
                    #onemat = ones[[icl]]
                    #vi = diag(sz) + onemat*thhat
                    #covi = vi
                    #umat = t(chol(covi))
                    #uinv = solve(umat)
                    uinv = uinv0[1:sz, 1:sz, drop=FALSE]
                }
                yi = ycl[[icl]]
				ytil = c(ytil, uinv %*% as.matrix(yi, ncol=1))
				ed = ed + sz
				gtil = rbind(gtil, uinv %*% gmat[st:ed, ,drop=F])
				st = ed + 1
			}
            if (any(shapes != 17)) {
                dsend = gtil[, (np + 1):(np + capm + capu + capt), drop = FALSE]
                zsend = gtil[, 1:np, drop = FALSE]
                ans = coneB(ytil, dsend, vmat = zsend, face=face)
                edf = ans$df
                face = ans$face
                bh = coef(ans)
                #test:
                if (any(round(bh[1:np],6) < 0)) {
                    pos = (1:np)[which(round(bh[1:np],6) < 0)]
                    face = unique(c(pos, face))
                    #print (bh[1:np])
                    #print (face)
                }
            } else {
                bh = solve(crossprod(gtil), t(gtil)) %*% ytil
                edf = nrow(bigmat)
                face = 1:edf
            }
            muhat = t(bigmat) %*% bh
			diff = mean((oldmu - muhat)^2)
			oldmu = muhat
            dd = t(bigmat[face, ,drop = FALSE])
            dd2 = gtil[,face,drop=FALSE]
            if (reml) {
                #dd = t(bigmat[face, ,drop = FALSE])
                #if (balanced) {
                #    xtx = t(dd)%*%dd
                #}
                xms = list()
                st = 1
                ed = 0
                for (icl in 1:ncl) {
                    sz = szs[icl]
                    ed = ed + sz
                    xms[[icl]] = dd[st:ed, ,drop=F]
                    #onevec = 1:sz*0+1
                    #onemat = onevec%*%t(onevec)
                    #ones[[icl]] = onemat
                    st = ed + 1
                }
                #if (balanced) {
                #    #oneMat = as.matrix(bdiag(ones))
                #    xtx2 = crossprod(dd, oneMat) %*% dd
                #}
            }
		}
		ebars = sapply(ecl, mean)
		sig2hat = fsig(thhat, szs, ecl, ncl, N=n, edf=edf, D=nrow(bigmat), type=type)
		siga2hat = sig2hat * thhat 
		ahat = ebars*szs*thhat/(1+szs*thhat)
	}
    #vectors for inference about param coefs:
    onevw = NULL; zmatw = NULL; xvecw = NULL; dusew = NULL
    st = 1
    ed = 0
    dd_nv = dd[,-c(1:np),drop=FALSE]
    #uinv0 is the one in the last iteration
    for (icl in 1:ncl) {
        sz = szs[icl]
        onevi = 1:sz*0+1
        if (!balanced) {
            uinv = uinv0[1:sz, 1:sz, drop=FALSE]
        }
        ed = ed + sz
        onevw = rbind(onevw, uinv %*% onevi)
        if (!is.null(xvec)) {
            xvecw = rbind(xvecw, uinv %*% xvec[st:ed, ,drop=F])
        }
        if (capk > 0) {
            zmatw = rbind(zmatw, uinv %*% zmat[st:ed, ,drop=F])
        }
        if (!is.null(dd_nv)) {
            dusew = rbind(dusew, uinv %*% dd_nv[st:ed, ,drop=F])
        }
        #print (c(st,ed))
        st = ed + 1
    }
    #dusew = dusew[,-c(1:np),drop=FALSE]
    #print (dim(dusew))
    df_obs = sum(abs(bh) > 0)
#new: tests and cis
#F-test for siga2
    pv.siga2 = ranef.test(ecl, szs, n, ncl)
#c.i. for siga2
    ci1 = NULL
    #if (type == 'b') {
        ci1 = ranef.ci(ecl, szs, n, ncl, level = 0.95)
    #}
    ci.siga2 = ci1
#c.i. for th and rho
    ci.th = ci.rho = ci.siga2.bi = NULL
    if (bisect) {
        #print (edf)
        ci2 = ranef.cith(thhat, sig2hat, siga2hat, ahat, ecl, szs, n, ncl, level = 0.95, xms=xms, p=edf, reml=reml)
        ci.th = ci2$ci.th
        ci.rho = ci2$ci.rho
        
        cia = ranef.cisiga(sig2hat, siga2hat, ahat, ecl, szs, n, ncl, level = 0.95, xms=xms, p=edf, reml=reml)
        ci.siga2.bi = cia$ci.siga2
    }
#c.i. for sig2
    ci.sig2 = ranef.cisig2(ecl, n, ncl, level = 0.95)
    coefskeep = bh
    thvecs = NULL
    if (capl > 0) {
        #new code:
        dcoefs = coefskeep[(np - capms + 1):(np + capm)]
        #dcoefs <- coefskeep[(np + 1):(np + capm)]
        #####################################################
        #thvecs is f(x), where x has one of the eight shapes#
        #####################################################
        thvecs = matrix(nrow = capl, ncol = n)
        ncon = 1
        for (i in 1:capl) {
            thvecs[i,] = t(delta[varlist == i,]) %*% dcoefs[varlist == i]
            #new:
            #if (shapes[i] > 2 & shapes[i] < 5) {
            if (shapes[i] > 2 & shapes[i] < 5 | shapes[i] > 10 & shapes[i] < 13) {
                ncon = ncon + 1
                #thvecs[i,] <- thvecs[i,] + zcoefs[ncon] * xmat[,i]
                thvecs[i,] = thvecs[i,] + vcoefs[capk + ncon] * xmat[,i]
            }
        }
    }
    #new:order thvecs back
    if (length(idx_s) > 0) {
        thvecs0 = thvecs
        thvecs0[idx_s,] = thvecs[1:length(idx_s), ]
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
#temp
    ncl = length(szs)
    bhmt = matrix(rep(bh[-1], each = ncl), nrow = ncl)
    coefs = cbind(ahat - bh[1], bhmt)
    colnames(coefs) = c("(Intercept)", paste("edge", 1:ncol(bhmt)))
    #zcoefs include one
    se.beta = 1:(capk + 1)*0
    tstat = 1:(capk + 1)*0
    pvals.beta = 1:(capk + 1)*0
    #zcoefs include onevec
####################################
#inference abour parametric coefs; need the weighted version
####################################
    zcoefs = bh[1:(1+capk)]
    imat = diag(n)
    pj = 0
    if (ncol(dusew) >= 1) {
        if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0) {
            #xvec = xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]
            #pm = cbind(1:n*0 + 1, xvec-mean(xvec) , t(bigmat_nv[duse, , drop = FALSE]))
            pm = cbind(xvecw , dusew)
        } else {
            #pm = cbind(1:n*0 + 1, t(bigmat_nv[duse, , drop = FALSE]))
            pm = dusew
        }
        pj = pm %*% solve(crossprod(pm), t(pm))
    }
    #pone = onevec %*% solve(crossprod(onevec), t(onevec))
    se2 = solve(t(cbind(onevw, zmatw)) %*% (imat - pj) %*% cbind(onevw, zmatw))*sig2hat
    se.beta = sqrt(as.vector(diag(se2)))
    tstat = zcoefs / se.beta
    if (n<=200){cpar=1.5}
    if ((n - cpar * df_obs) <= 0) {
        pvals.beta=2 * (1 - pt(abs(tstat),  df_obs))
        warning ('Effective degrees of freedom is close to the number of observations! Inference about parametric covariates is not reliable!')
    } else {
        #pvals.beta[i] <- 2 * (1 - pt(abs(tstat[i]),  n - np - cpar * df_obs))
        pvals.beta=2 * (1 - pt(abs(tstat),  n - cpar * df_obs))
    }
#############################
#end of inference abour parametric coefs
#############################
    if (capl > 0) {
        xcoefs = bh[(capk + 2):np]
    } else {xcoefs = NULL}
    if (np > 0) {
        vcoefs = bh[1:np]
    } else {vcoefs = NULL}
#new: order back
    knotsuse2 = knotsuse
    numknotsuse2 = numknotsuse
    mslst2 = mslst
    xmat2 = xmat
    if (length(idx_s) > 0) {
        knotsuse0 = knotsuse
        numknotsuse0 = numknotsuse
        mslst0 = mslst
        knotsuse0[idx_s] = knotsuse[1:length(idx_s)]
        numknotsuse0[idx_s] = numknotsuse[1:length(idx_s)]
        mslst0[idx_s] = mslst[1:length(idx_s)]
        if (length(idx) > 0) {
            knotsuse0[idx] = knotsuse[(1+length(idx_s)):capl]
            numknotsuse0[idx] = numknotsuse[(1+length(idx_s)):capl]
            mslst0[idx] = mslst[(1+length(idx_s)):capl]
        }
        knotsuse = knotsuse0
        numknotsuse = numknotsuse0
        mslst = mslst0
    }
    rslt = list(muhat = muhat, coefs = coefs, bh = bh, zcoefs = zcoefs, pvals.beta = pvals.beta, se.beta = se.beta, vcoefs = vcoefs, ahat = ahat, sig2hat = sig2hat, siga2hat = siga2hat, thhat = thhat, bigmat = bigmat, gtil=gtil, dd2=dd2, np = np, knots = knotsuse, knots2 = knotsuse2, numknots = numknotsuse, numknots2 = numknotsuse2, ms = mslst, ms2 = mslst2, xmat2 = xmat2, xid1 = xid1, xid2 = xid2, capm = capm, capms = capms, capk = capk, capt = capt, capu = capu, etacomps = thvecs, mod.lmer = mod.lmer, pv.siga2 = pv.siga2, ci.siga2 = ci.siga2, ci.siga2.bi = ci.siga2.bi, ci.th = ci.th, ci.rho = ci.rho, ci.sig2 = ci.sig2, ones = ones, resid_df_obs = n - cpar * df_obs, edf = df_obs)
    return (rslt)
}

#####
#mle#
#####
#balanced:
fth = function(ycl, ncl, N) {
	ni = N/ncl
	ybar = sapply(ycl, mean)
	y = unlist(ycl)
	num = ni^2*sum(ybar^2) - sum(y^2)
	den = ni*sum(y^2) - ni^2*sum(ybar^2)
	return (num/den)
}

#unbalanced (uniroot):
fth2 = function(th, szs, ycl, N) {
	ybar = sapply(ycl, mean)
	y = unlist(ycl)
	num = sum(szs^2*ybar^2/(1+szs*th)^2)
	den = sum(y^2) - sum(th*szs^2*ybar^2/(1+szs*th))
	obj = N/2*num/den - 1/2*sum(szs/(1+szs*th))
	return (obj)
}

######
#reml#
######
fth2rm = function(th, szs, ycl, N, xcl, p=2, type='b', xtx=NULL, xtx2=NULL, xmat_face=NULL, ones=NULL) {
    ybar = sapply(ycl, mean)
    y = unlist(ycl)
    num = sum(szs^2*ybar^2/(1+szs*th)^2)
    den = sum(y^2) - sum(th*szs^2*ybar^2/(1+szs*th))
    ncl = length(ycl)
    hmat = matrix(0, p, p)
    xtils = list()
    #if (type == 'b') {
    #    n = szs[1]
    #    hmat = xtx - th/(1+n*th)*xtx2
    #} else {
    ones2 = list()
    for(icl in 1:ncl){
        ni = szs[icl]
        xi = xcl[[icl]]
        xm = xi
        #onevec = 1:ni*0+1
        #onemat = onevec%*%t(onevec)
        onemat = ones[[icl]]
        #ones[[icl]] = onemat*th/(1+ni*th)
        ones2[[icl]] = onemat/(1+ni*th)^2
        rinv = diag(ni) - th/(1+ni*th)*onemat
        hmat = hmat + t(xm) %*% rinv %*% xm
        #xtil = t(onevec)%*%xm
        #xtils[[icl]] = xtil
    }
    #oneMat = as.matrix(bdiag(ones))
    #oneMat2 = as.matrix(bdiag(ones2))
    #hmat = xtx - crossprod(xmat_face, oneMat) %*% xmat_face
    
    #ones=ones2=list()
    #    st = 1
    #    ed = 0
    #    for (icl in 1:ncl) {
    #        sz=szs[icl]
    #        ed = ed + sz
    #        onevec = 1:sz*0+1
    #        onemat = onevec%*%t(onevec)
    #        ones[[icl]] = onemat*th/(1+sz*th)
    #        ones2[[icl]] = onemat/(1+sz*th)^2
    #        st = ed + 1
    #    }
    #    oneMat = as.matrix(bdiag(ones))
    #    oneMat2 = as.matrix(bdiag(ones2))
    #    hmat = xtx - crossprod(xmat_face, oneMat) %*% xmat_face
    #}
    #hinv = solve(hmat)
    lmat = chol(hmat)
    hinv = chol2inv(lmat)
    tr = 0
    #if (type == 'b') {
    #    n = szs[1]
        #onevec = 1:n*0+1
        for (icl in 1:ncl) {
            ni = szs[icl]
            xi = xcl[[icl]]
            onevec = 1:ni*0+1
            xtil = t(onevec)%*%xi
            tr = tr + sum(diag(hinv %*% crossprod(xtil)/(1+ni*th)^2))
            #tr = tr + 1/(1+n*th)^2*sum(diag(hinv %*% crossprod(xtil)))
        }
    #    tr = 1/(1+n*th)^2*sum(diag(hinv %*% xtx2))
    #} else {
        # for(icl in 1:ncl) {
        #    ni = szs[icl]
        #    xtil = xtils[[icl]]
        #    tr = tr + sum(diag(hinv %*% crossprod(xtil)/(1+ni*th)^2))
        #}
        #    xtx2ub = crossprod(xmat_face, oneMat2) %*% xmat_face
        #    tr = sum(diag(hinv %*% xtx2ub))
    #}
    rml = 1/2*tr
    obj = (N-p)/2*num/den - 1/2*sum(szs/(1+szs*th)) + rml
    return (obj)
}


########################################
#sig^2 (use sig2hat in with N-d-1 now.)
########################################
fsig = function(thhat, szs, ycl, ncl, N, edf, D, type='b') {
	ybars = sapply(ycl, mean)
    d = min(1.5*edf, D)
	if (type == 'b') {
		sz = N/ncl
		sig2hat = (sum(unlist(ycl)^2) - sz^2*thhat/(1+sz*thhat) * sum(ybars^2))/(N-d-1)
	} else {
		sig2hat = (sum(unlist(ycl)^2) - sum(thhat*szs^2*ybars^2/(1 + szs*thhat)))/(N-d-1)
	}
	return (sig2hat)
}

########################################
#function to get the residual vector for each cluster
#return a list
########################################
f_ecl = function(evec, ncl, sz) {
	ecl = list()
	st = 1
	ed = 0
	for (icl in 1:ncl) {
		if (length(sz) > 1) {
			szi = sz[icl]
		} else {szi = sz}
		ed = ed + szi
		ecl[[icl]] = evec[st:ed]
		st = ed + 1
	}
	return (ecl)
}

#########################
#F-test siga2           #
#balanced and unbalanced#
#########################
ranef.test = function(ecl, szs, N, ncl) {
    evec = unlist(ecl)
    ebars = sapply(ecl, mean)
    sse = 0
    ssb = 0
    ebar = mean(evec)
    for(icl in 1:ncl) {
        sz = szs[icl]
        eibar = ebars[[icl]]
        ei = ecl[[icl]]
        ssei = sum((ei - eibar)^2)
        sse = sse + ssei
        ssbi = sz*(eibar - ebar)^2
        ssb = ssb + ssbi
    }
    mse = sse/(N-ncl)
    msb = ssb/(ncl-1)
    fstat = msb/mse
    pv = 1-pf(fstat, df1=ncl-1, df2=N-ncl)
    return (pv)
}

##########################
#C.I. for siga2          #
#balanced only           #
##Satterthwaite procedure#
##########################
ranef.ci = function(ecl, szs, N, ncl, level = 0.95) {
    evec = unlist(ecl)
    ebars = sapply(ecl, mean)
    sse = 0
    ssb = 0
    ebar = mean(evec)
    for(icl in 1:ncl) {
        sz = szs[icl]
        eibar = ebars[[icl]]
        ei = ecl[[icl]]
        ssei = sum((ei - eibar)^2)
        sse = sse + ssei
        ssbi = sz*(eibar - ebar)^2
        ssb = ssb + ssbi
    }
    mse = sse/(N-ncl)
    msb = ssb/(ncl-1)
    sz2 = 1/(ncl-1)*(N - sum(szs^2)/N)
    siga2hat = (msb-mse)/sz2
    dfsiga = (sz2*siga2hat)**2/(msb**2/(ncl-1)+mse**2/(N-ncl))
    alpha = 1-level
    lwr = dfsiga*siga2hat/qchisq(1-alpha/2, df=dfsiga)
    upp = dfsiga*siga2hat/qchisq(alpha/2, df=dfsiga)
    ci = c(lwr, upp)
    return (ci)
}


##################################################
#C.I. for thhat  and siga2/(siga2+sig2)          #
#balanced and unbalanced                         #
#use profile-log-like and bisection              #
##################################################
ranef.cith = function(thhat, sig2hat, siga2hat, ahat, ecl, szs, N, ncl, level = 0.95, xms, p, reml=TRUE) {
    N = sum(szs)
    evec = unlist(ecl)
    thval = fmin(thhat, ncl, ecl, N, xms, p, reml=reml)
    ebars = sapply(ecl, mean)
    #sigaval = fmin2(siga2hat^.5, sig2hat, ncl, ecl, N)
    ans = try(ans0<-uniroot(fn2, ncl=ncl, ycl=ecl, N=N, thval=thval, xms=xms, p=p, reml=reml, interval=c(1e-10, thhat^1), tol=.Machine$double.eps),silent=TRUE)
    if (class(ans) == 'try-error') {
        lwr = 0#1e-7
    } else {lwr = ans$root}
    #upp = uniroot(fn2, ncl=ncl, ycl=ecl, N=N, thval=thval, xms=xms, p=p, reml=reml, interval=c(thhat^1, 1e+4), tol=.Machine$double.eps)$root
    ans2 = try(ans20<-uniroot(fn2, ncl=ncl, ycl=ecl, N=N, thval=thval, xms=xms, p=p, reml=reml, interval=c(thhat^1, 1e+4), tol=.Machine$double.eps),silent=TRUE)
    if (class(ans2) == 'try-error') {
        upp = 1e+4
    } else {upp = ans2$root}
    ci = c(lwr, upp)
    ci2 = c(lwr/(1+lwr), upp/(1+upp))
    ans = list(ci.th = ci, ci.rho = ci2)
    return (ans)
}

##################################################
#C.I. for siga2; sig2hat being mle               #
#balanced and unbalanced                         #
#use profile-log-like and bisection              #
##################################################
ranef.cisiga = function(sig2hat, siga2hat, ahat, ecl, szs, N, ncl, level = 0.95, xms, p, reml=TRUE) {
    N = sum(szs)
    evec = unlist(ecl)
    thval = fmin2(siga2=siga2hat, sig2hat, ncl, ecl, N, xms, p, reml=reml)
    ebars = sapply(ecl, mean)
    #sigaval = fmin2(siga2hat^.5, sig2hat, ncl, ecl, N)
    ans = try(ans0<-uniroot(fn2a, sig2hat=sig2hat, ncl=ncl, ycl=ecl, N=N, thval=thval, xms=xms, p=p, reml=reml, interval=c(1e-10, siga2hat^1), tol=.Machine$double.eps),silent=TRUE)
    if (class(ans) == 'try-error') {
        lwr = 0#1e-7
    } else {lwr = ans$root}
    #upp = uniroot(fn2, ncl=ncl, ycl=ecl, N=N, thval=thval, xms=xms, p=p, reml=reml, interval=c(thhat^1, 1e+4), tol=.Machine$double.eps)$root
    ans2 = try(ans20<-uniroot(fn2a, sig2hat=sig2hat, ncl=ncl, ycl=ecl, N=N, thval=thval, xms=xms, p=p, reml=reml, interval=c(siga2hat^1, 1e+4), tol=.Machine$double.eps),silent=TRUE)
    if (class(ans2) == 'try-error') {
        upp = 1e+4
    } else {upp = ans2$root}
    ci = c(lwr, upp)
    #ci2 = c(lwr/(1+lwr), upp/(1+upp))
    #ans = list(ci.th = ci, ci.rho = ci2)
    ans = list(ci.siga2 = ci)
    return (ans)
}

##################################################
#C.I. for sig2                                   #
#balanced and unbalanced                         #
##################################################
ranef.cisig2 = function(ecl, N, ncl, level = 0.95) {
    evec = unlist(ecl)
    ebars = sapply(ecl, mean)
    sse = 0
    ebar = mean(evec)
    for(icl in 1:ncl) {
        eibar = ebars[[icl]]
        ei = ecl[[icl]]
        ssei = sum((ei - eibar)^2)
        sse = sse + ssei
    }
    alpha = 1-level
    lwr = sse/qchisq(1-alpha/2, df=N-ncl)
    upp = sse/qchisq(alpha/2, df=N-ncl)
    ci = c(lwr, upp)
    return (ci)
}

###################################
#subroutines for tests and c.i.'s#
###################################
#neg p: negative profile-loglikelihood of theta
fmin = function(theta, ncl, ycl, N, xms=NULL, p=2, reml=TRUE) {
    if (reml) {
        acc1 = acc2 = acc3 = 0
        hmat = matrix(0, p, p)
        for (i in 1:ncl) {
            yi = ycl[[i]]
            ni = length(yi)
            one = matrix(rep(1, ni), ncol=1)
            onemat = tcrossprod(one)
            viinv = diag(ni) - theta / (1+ni*theta) * onemat
            detvi = (1+ni*theta)
            acc1 = acc1 + t(yi) %*% viinv %*% yi
            acc2 = acc2 + log(detvi)
        
    
            xm = xms[[i]]
            rinv = diag(ni) - theta/(1+ni*theta)*onemat
            hmat = hmat + t(xm) %*% rinv %*% xm
        }
        obj = (N-p)/2 * log(acc1) + 1/2 * acc2 + 1/2 * log(det(hmat))
    } else {
        acc1 = acc2 = 0
        for (i in 1:ncl) {
            yi = ycl[[i]]
            ni = length(yi)
            one = matrix(rep(1, ni), ncol=1)
            onemat = tcrossprod(one)
            viinv = diag(ni) - theta / (1+ni*theta) * onemat
            #detvi = 1 / (1+ni*theta)
            detvi = (1+ni*theta)
            acc1 = acc1 + t(yi) %*% viinv %*% yi
            acc2 = acc2 + log(detvi)
            #vi = diag(ni) + theta * tcrossprod(one)
            #acc1 = acc1 + t(yi) %*% solve(vi, yi)
            #acc2 = acc2 + log(det(vi))
        }
        obj = N/2 * log(acc1) + 1/2 * acc2
    }
    return (obj)
}


#the theta profile-loglike curve for bisection
fn2 = function(x, ncl, ycl, N, thval, level=0.95, xms=NULL, p=2, reml=TRUE) {
    thval2 = -thval
    cts = thval2 - 1/2*qchisq(level, df=1)
    obj = -fmin(x, ncl, ycl, N, xms, p, reml) - cts[1]
    return (obj)
}

#neg p: negative profile-loglikelihood of siga2; sig2hat being fixed as mle
fmin2 = function(siga2, sig2hat, ncl, ycl, N, xms=NULL, p=2, reml=TRUE) {
    if (reml) {
        acc1 = acc2 = acc3 = 0
        hmat = matrix(0, p, p)
        for (i in 1:ncl) {
            yi = ycl[[i]]
            ni = length(yi)
            one = matrix(rep(1, ni), ncol=1)
            onemat = tcrossprod(one)
            viinv = diag(ni) - siga2 / (sig2hat+ni*siga2) * onemat
            detvi = (1+ni*siga2/sig2hat)
            acc1 = acc1 + t(yi) %*% viinv %*% yi
            acc2 = acc2 + log(detvi)
            
            
            xm = xms[[i]]
            rinv = diag(ni) - siga2/(sig2hat+ni*siga2)*onemat
            hmat = hmat + t(xm) %*% rinv %*% xm
        }
        obj = (N-p)/2 * log(acc1) + 1/2 * acc2 + 1/2 * log(det(hmat))
    } else {
        acc1 = acc2 = 0
        for (i in 1:ncl) {
            yi = ycl[[i]]
            ni = length(yi)
            one = matrix(rep(1, ni), ncol=1)
            onemat = tcrossprod(one)
            viinv = diag(ni) - siga2 / (sig2hat+ni*siga2) * onemat
            #detvi = 1 / (1+ni*theta)
            detvi = (1+ni*siga2/sig2hat)
            acc1 = acc1 + t(yi) %*% viinv %*% yi
            acc2 = acc2 + log(detvi)
            #vi = diag(ni) + theta * tcrossprod(one)
            #acc1 = acc1 + t(yi) %*% solve(vi, yi)
            #acc2 = acc2 + log(det(vi))
        }
        obj = N/2 * log(acc1) + 1/2 * acc2
    }
    return (obj)
}

#the siga2 profile-loglike curve for bisection
fn2a = function(siga2, sig2hat, ncl, ycl, N, thval, level=0.95, xms=NULL, p=2, reml=TRUE) {
    thval2 = -thval
    cts = thval2 - 1/2*qchisq(level, df=1)
    obj = -fmin2(siga2, sig2hat, ncl, ycl, N, xms, p, reml) - cts[1]
    return (obj)
}

#################
#new coef method#
#################
coef.cgamm <- function(object,...) {
    ans <- object$coefs
    ans
}

ranef.cgamm <- function(object,...) {
    ans <- object$ahat
    ans
}

fixef.cgamm <- function(object,...) {
    ans <- object$bh
    ans
}


###############
#predict.cgamm#
#smooth only  #
#guasisan only#
###############
predict.cgamm = function(object, newData, interval = c("none", "confidence", "prediction"), type = c("response", "link"), level = 0.95, n.mix = 500, var.f = NULL,...) {
    #print (is.data.frame(newData))
    #print (newData)
    #new:
    family = object$family
    cicfamily = CicFamily(family)
    muhat.fun = cicfamily$muhat.fun
    if (!inherits(object, "cgamm")) {
        warning("calling predict.cgamm(<fake-cgam-object>) ...")
    }
    if (missing(newData) || is.null(newData)) {
        #if (missing(newData)) {
        #etahat = object$etahat
        #muhat = muhat.fun(etahat, fml = family$family)
        #ans = list(fit = muhat, etahat = etahat, newbigmat = object$bigmat)
        #ans = list(fit = muhat)
        muhat = object$muhat
        ahat = object$ahat
        ans = list(fix_effect = muhat, random_effect = ahat)
        return (ans)
    }
    if (!is.data.frame(newData)) {
        #newData = as.data.frame(newData)
        stop ("newData must be a data frame!")
    }
    #shapes = object$shapes
    #new: used for ci
    prior.w = object$prior.w
    y = object$y
    muhat = object$muhat
    #shapes = object$shapesx
    shapes = object$shapes
    np = object$d0; capm = object$capm; capk = object$capk; capt = object$capt; capu = object$capu
    #new:
    xid10 = object$xid1; xid20 = object$xid2;
    #uid1 = object$uid1; uid2 = object$uid2; tid1 = object$tid1; tid2 = object$tid2
    #new:
    xmat0 = object$xmat0; knots0 = object$knots0; numknots0 = object$numknots0; sps0 = object$sps0; ms0 = object$ms0
    #xmat0 = object$xmat2; knots0 = object$knots2; numknots0 = object$numknots2; sps0 = object$sps2; ms0 = object$ms2
    zmat = object$zmat; umb = object$umb; tr = object$tr
    #new:
    ztb = object$ztb; zid1 = object$zid1; zid2 = object$zid2; iz = 1
    bigmat = object$bigmat; umbrella.delta = object$umbrella.delta; tree.delta = object$tree.delta
    #coefs = object$coefs; zcoefs = object$zcoefs; vcoefs = object$vcoefs; xcoefs0 = object$xcoefs; ucoefs = object$ucoefs; tcoefs = object$tcoefs
    #temp:
    coefs = object$bh; zcoefs = object$zcoefs; vcoefs = object$vcoefs; xcoefs0 = object$xcoefs; ucoefs = object$ucoefs; tcoefs = object$tcoefs
    tt = object$tms
    Terms = delete.response(tt)
    #model.frame will re-organize newData in the original order in formula
    #temp: need a better method, create fake id
    lbs = attributes(object$tms)$term.labels
    idlb = rev(lbs)[1]
    group = object$id
    if (!any(names(newData) %in% lbs)){
    #used for confidence interval
        newData[[idlb]] = rep(1, nrow(newData))
        m = model.frame(Terms, newData)
        #use real group for prediction interval
        #group = object$id
        group_new = NULL
        if (interval == "prediction") {
            stop ("Group information is missing for the newData!")
        }
    } else {
     #used for prediction interval
        m = model.frame(Terms, newData)
        nmsm = names(m)
        #group = object$id
        group_new = m[,which(nmsm%in%idlb)]
    }
    #print (m)
    #new: don't include id in newdata
    #rm_id = NULL
    #for (i in 1:ncol(m)) {
    #    mi = m[,i]
    #    if (all(mi == object$id)) {
    #        rm_id = i
    #        idnew = mi
    #    }
        #if (is.null(attributes(mi)$shape)) {
        #   rm_id = i
        #}
    #}
    nmsm = names(m)
    rm_id = which(nmsm%in%idlb)
    newdata = m[, -rm_id, drop=F]
    #new: need id or group for prediction interval
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
#different from cgam: simpler
    n = ncol(bigmat)
    nc = ncol(xmat0)
    spl = splpr = NULL
    m_acc = 0
    object$ms0 -> ms0
    #sh = 17 first
    object$shapes -> shapes
    for (i in 1:nc) {
        msi = ms0[[i]]
        shpi = shapes[i]
        ki = knots0[[i]]
        xi = xmat0[,i]
        xipr = newx0[,i]
        if (any(xipr > max(xi)) | any(xipr < min(xi))) {
            stop ("No extrapolation is allowed in cgamm prediction!")
        }
        deli = makedelta(xi, shpi, knots = ki)
        #print (i)
        #if (i == 2) {
        #	stop (print (msi))
        #}
        spli = deli$amat
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
    #only xmatpr and splpr have interpolation points
    xmatpr = cbind(newv, t(splpr))
    #temp:
    muhatpr = xmatpr %*% coefs[1:ncol(xmatpr), ,drop=FALSE]
    if ("none" %in% interval) {
        ans = list(fit = muhatpr)
    } else if (interval == "confidence" | interval == "prediction") {
        #new:
        ones = object$ones
        capk = object$capk
        p = 1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)
        #zmat = t(bigmat[1:p, ,drop = FALSE])
        if (family$family == "gaussian") {
            sig2hat = object$sig2hat
            siga2hat = object$siga2hat
            thhat = object$thhat
            szs = object$szs
            ahat = object$ahat
            muhat = object$muhat
            balanced = FALSE
            if (length(unique(szs)) == 1) {
                balanced = TRUE
            }
            ncl = length(szs)
            edges = t(bigmat)
#new:
            nd = ncol(edges)
            #face = 1:nd
#var.f is variance of muhat
if (is.null(var.f)) {
                
            face = NULL
            #nloop = 1000
            nloop = n.mix
            nsec = 2^m_acc
       
            if (balanced) {
                oneMat = ones[[1]]
                sz = szs[1]
            } else {
                sz = max(szs)
                pos = which(szs == sz)[1]
                oneMat = ones[[pos]]
            }
            wi = 1*(diag(sz) + thhat*oneMat)
            covi = wi
            umat = t(chol(covi))
            #uinv0 is used for unbalanced
            uinv0 = uinv = solve(umat)
            #vms = vector('list', ncl)
            #for(icl in 1:ncl) {vms[[icl]] = uinv}
            #uinvmat = as.matrix(bdiag(vms))
            if (round(thhat,6) != 0) {
            ## estimate sector probabilties:
                #sector = 1:nsec*0
                #new:
                sector = NULL
                times = NULL
                df = NULL
                
                for(iloop in 1:nloop){
                    ysim = NULL
                    etil = NULL
                    st = 1
                    ed = 0
                    for (icl in 1:ncl) {
                        sz = szs[icl]
                        ed = ed + sz
                        mui = muhat[st:ed]
                        #one = 1:sz*0+1
                        #zmi = as.matrix(one, ncol=1)
                        #wi = matrix(0, sz, sz)
                        #wi = 1*(diag(sz) + thhat*zmi%*%t(zmi))
                        if (!balanced) {
                            #oneMat = ones[[icl]]
                            #wi = 1*(diag(sz) + thhat*oneMat)
                            #covi = wi
                            #umat = t(chol(covi))
                            #uinv = solve(umat)
                            uinv = uinv0[1:sz, 1:sz, drop=FALSE]
                        }
                        #wiinv = diag(sz) - thhat/(1+sz*thhat)*oneMat
                        #upper = chol(wiinv)
                        #uinv = matrix(0, sz, sz)
                        #for(i in 1:sz){mat[i,] = rev(upper[10-i+1,])}
                        ysi = uinv %*% (rnorm(1, sd = siga2hat^.5) + as.matrix(mui, ncol=1) + rnorm(sz, sd=sig2hat^.5))
                        ysim = c(ysim, ysi)
                        emati = edges[st:ed, ,drop=F]
                        etil = rbind(etil, uinv %*% emati)
                        st = ed + 1
                    }
                    #dsend = etil[, (np+1):(capm+np)]
                    #ysim = uinvmat %*% (rep(rnorm(ncl, sd = siga2hat^.5), time=szs) + as.matrix(muhat, ncol=1) + rnorm(n, sd=sig2hat^.5))
                    #etil = uinvmat %*% edges
                    dsend = etil[, -(1:np)]
                    zsend = etil[, 1:np]
                    ans = coneB(ysim, dsend, vmat = zsend, face = NULL)
                    #face = ans$face
                    cf = round(ans$coefs[(np+1):(m_acc+np)],10)
                    #sec = 1:m_acc*0
                    #sec[cf>0] = 1
                    #r = makebin(sec)+1
                    #sector[r] = sector[r]+1
                    
                    #new:
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
            } else {
                prior.w = 1:n*0 + 1
                #sector = 1:nsec*0
                #for (iloop in 1:nloop) {
                #    ysim = muhat + rnorm(n)*sig2hat^.5
                #    ysim = ysim * sqrt(prior.w)
                #    dsend = edges[, (np+1):(m_acc+np)]
                #    zsend = edges[, 1:np]
                #    ans = coneB(ysim, dsend, zsend, face = NULL)
                #    #face = ans$face
                #    cf = round(ans$coefs[(np+1):(m_acc+np)], 10)
                #    sec = 1:m_acc*0
                #    sec[cf > 0] = 1
                #    r = makebin(sec) + 1
                #    sector[r] = sector[r] + 1
                #}
                
                #new:
                sector = NULL
                times = NULL
                df = NULL
                #first column shows sector; second column shows times
                for (iloop in 1:nloop) {
                    ysim = muhat + rnorm(n)*sig2hat^.5
                    ysim = ysim * sqrt(prior.w)
                    dsend = edges[, (np+1):(m_acc+np)]
                    zsend = edges[, 1:np]
                    ans = coneB(ysim, dsend, zsend, face = NULL)
                    cf = round(ans$coefs[(np+1):(m_acc+np)], 10)
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
            }
            ### calculate the mixture hat matrix:
            #bsec = matrix(0, nrow=nsec, ncol=2)
            #bsec[,1] = 1:nsec
            #bsec[,2] = sector / nsec
            #keep = sector > 0
            #bsec = bsec[keep,]
            #ns = dim(bsec)[1]
            #bsec[,2] = bsec[,2] / sum(bsec[,2])
            
            #new:
            ns = nrow(df)
            bsec = df
            bsec[,2] = bsec[,2] / sum(bsec[,2])
            ord = order(bsec[,1])
            bsec = bsec[ord,,drop=FALSE]
            
            # nv and np are the dim of vmat
            nv = np
            zmat = zsend
            spl = t(dsend)
            ### calculate the mixture cov(alpha) matrix:
            obs = 1:m_acc;oobs = 1:(m_acc+nv)
            acov = matrix(0, nrow = m_acc+nv, ncol = m_acc+nv)
            for (is in 1:ns) {
                if (bsec[is,2] > 0) {
                    jvec = getbin(bsec[is,1], m_acc)
                    if (sum(jvec) == 1) {
                        smat = cbind(zmat, t(spl[which(jvec==1),,drop=F]))
                    } else if (sum(jvec) == 0) {
                        smat = zmat
                    } else {
                        smat = cbind(zmat, t(spl[which(jvec==1),,drop=F]))
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
            acov = acov*sig2hat
}
            #only xmatpr and splpr have interpolation points
            #xmatpr = cbind(newv, t(splpr))
            #temp:
            #muhatpr = xmatpr %*% coefs[1:ncol(xmatpr), ,drop=FALSE]
            #new: C.I. level
            mult = qnorm((1 - level)/2, lower.tail=FALSE)
            if (interval == "confidence") {
                if (is.null(var.f)) {
                    var.f = diag(xmatpr%*%acov%*%t(xmatpr))
                    #hl = mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
                }
                hl = mult*sqrt(var.f)
            }
            if (interval == "prediction") {
                #imat = diag(n)
                bh = object$bh
                dd = t(bigmat[abs(bh)>0, ,drop=FALSE])
                #vp1 = rep(sig2hat, n)
                #vp2 = diag(xmatpr%*%acov%*%t(xmatpr))
                nnew = nrow(xmatpr)
                #imat = diag(nnew)
                vp1 = rep(sig2hat, nnew)
                if (is.null(var.f)){
                    var.f = diag(xmatpr%*%acov%*%t(xmatpr))
                }
                vp2 = var.f
    
                imat2 = diag(ncl)
                gm = siga2hat*imat2
                rinvlst = ones = list()
                if (balanced) {
                    pos = 1
                } else {
                   pos = min(which(szs == max(szs)))
                }
                n1=szs[pos]
                one1=rep(1,n1)
                onem1=tcrossprod(one1)
                im1=diag(n1)
                rinv1=im1-thhat/(1+n1*thhat)*onem1
                rinv0=rinv1
                for(i in 1:ncl){
                    if (balanced) {
                        rinvlst[[i]] = rinv1
                        ones[[i]] = one1
                    } else {
                        ni = szs[i]
                        onevec = rep(1, ni)
                        #onemat = tcrossprod(onevec)
                        onemat = onem1[1:ni,1:ni]
                        #imatn = diag(ni)
                        imatn = im1[1:ni,1:ni]
                        rinv = imatn-thhat/(1+ni*thhat)*onemat
                        #rinv = rinv0[1:ni,1:ni]
                        rinvlst[[i]] = rinv
                        ones[[i]] = onevec
                    }
                }
                rinv = as.matrix(bdiag(rinvlst))
                onem = as.matrix(bdiag(ones))
                
                mult1 = t(dd)%*%rinv
                pw = dd%*%solve(mult1%*%dd)%*%mult1
                vp30 = gm-siga2hat*thhat*t(onem)%*%rinv%*%(diag(n)-pw)%*%onem
                #vp3 = diag(onem%*%vp30%*%t(onem))
                vp3 = diag(vp30)
                
                #vp4 = -siga2hat*pw%*%onem%*%t(onem)
                #var.pred = diag(vp1 + vp2 + vp3 + 2*vp4)
                #pwpr = (xmatpr[,which(abs(bh)>0)])%*%solve(mult1%*%dd)%*%mult1
                #vp40 = -siga2hat*pwpr%*%onem
                lwr = upp = fitpr = var.pred = rep(0, nnew)
                #ugr is unique groups in the sample data
                ugr = unique(group)
                #ugr_new is unique groups in the data for prediction
                ugr_new = unique(group_new)
                ncl_new = length(ugr_new)
                for(i in 1:ncl_new){
                    gri = ugr_new[i]
                    gr_id = which(group_new%in%gri)
                    szi = length(gr_id)
                    oneveci = matrix(rep(1, szi), ncol=1)
                    evi = matrix(rep(0, ncl), nrow=1);evi[i]=1
                    
                    face = which(abs(bh)>0)
                    pwpr = (xmatpr[gr_id, face])%*%solve(mult1%*%dd)%*%mult1
                    vp40 = -siga2hat*pwpr%*%onem
                    
                    vp4i = diag(vp40%*%t(evi)%*%t(oneveci))
                    #pos_cl is the cluster position in the old sample
                    pos_cl =  which(ugr%in%gri)
                    var.predi = rep(sig2hat, szi) + vp2[gr_id] + rep(vp3[pos_cl], szi) + 2*vp4i
                    var.pred[gr_id] = var.predi
        
                    lwr[gr_id] = muhatpr[gr_id] + ahat[pos_cl] - mult*sqrt(var.predi)
                    upp[gr_id] = muhatpr[gr_id] + ahat[pos_cl] + mult*sqrt(var.predi)
                    fitpr[gr_id] = muhatpr[gr_id] + ahat[pos_cl]
                }
                #hl = mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
            }
            if (interval == "confidence") {
                ans = list(fit = muhatpr, lower = muhatpr - hl, upper = muhatpr + hl, var.f=var.f)
            }
            if (interval == "prediction") {
                ans = list(fit = fitpr, lower = lwr, upper = upp, var.f=var.f, var.pred=var.pred)
            }
        }
    }
    class(ans) = "cgamp"
    return (ans) 
}

#################################
#sub-routines for c.i.
#################################
makebin = function(x){
    k = length(x)
    r = 0
    for(i in 1:k){r = r + x[k-i+1]*2^(i-1)}
    r
}

getvars = function(num){
    i = num
    digits = 0
    power = 0
    while(digits == 0){
        if(i<2^power){digits = power}
        power = power+1
    }
    binry = 1:digits*0
    if(num>0){binry[1] = 1}
    i = i-2^(digits-1)
    power = digits-2
    for(p in power:0){
        if(i >= 2^p){
            i = i-2^p
            binry[digits-p] = 1
        }
    }
    binry
}

getbin = function(num, capl){
    br = getvars(num-1)
    digits = length(br)
    binrep = 1:capl*0
    binrep[(capl-digits+1):capl] = br
    binrep
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

###############
#summary.cgamm#
###############
summary.cgamm <- function(object,...) {
    if (!is.null(object$zcoefs)) {
        family <- object$family
        resid_df_obs <- object$resid_df_obs
        #wt.iter <- object$wt.iter
        coefs <- object$zcoefs
        se <- object$se.beta
        tval <- coefs / se
        pvalbeta <- object$pvals.beta
        n <- length(coefs)
        #sse0 <- object$SSE0
        #sse1 <- object$SSE1
        #cic <- object$cic
        #deviance <- object$deviance
        #null_deviance <- object$null_deviance
        #df <- object$df
        #null_df <- object$null_df
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
        #pvs <- object$pvs
        #s.edf <- object$s.edf
        #bstats <- object$bstats
        #if (wt.iter) {
        #    rslt1 <- data.frame("Estimate" = round(coefs, 4), "StdErr" = round(se, 4), "z.value" = round(tval, 4), "p.value" = round(pvalbeta, 4))
        #    rownames(rslt1)[1] <- "(Intercept)"
        #    if (n > 1) {
        #        lzid <- length(zid1)
        #        for (i in 1:lzid) {
        #            pos1 <- zid1[i]; pos2 <- zid2[i]
        #            for (j in pos1:pos2) {
        #                if (!is_param[i]) {
        #                    rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], rownames(rslt1)[j + 1], sep = "")
        #                } else {
        #                    rownames(rslt1)[j + 1] <- paste(attributes(tms)$term.labels[zid[i] - 1], vals[j], sep = "")
        #                }
        #            }
        #        }
        #    }
        #    rslt1 <- as.matrix(rslt1)
            #new:
            #rslt2 <- NULL
            #if (!is.null(pvs)) {
            #    rslt2 <- data.frame("edf" = round(s.edf, 4), "mixture of Beta" = round(bstats, 4), "p.value" = round(pvs, 4))
                #rownames(rslt2) <- attributes(tms)$term.labels
                #debugged: check more
            #    if (!is.null(zid)) {
            #        rownames(rslt2) <- (attributes(tms)$term.labels)[-(zid-1)]
            #    } else {
            #        rownames(rslt2) <- (attributes(tms)$term.labels)
            #    }
            #}
        #} else {
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
            #rslt2 <- NULL
            #if (!is.null(pvs)) {
            #    rslt2 <- data.frame("edf" = round(s.edf, 4), "mixture of Beta" = round(bstats, 4), "p.value" = round(pvs, 4))
            #    #debugged: check more
            #    if (!is.null(zid)) {
            #        rownames(rslt2) <- (attributes(tms)$term.labels)[-(zid-1)]
            #    } else {
            #        rownames(rslt2) <- (attributes(tms)$term.labels)
            #    }
            #}
        #}
        #if (!is.null(sse0) & !is.null(sse1)) {
        #rslt2 <- cbind(SSE.Linear = sse0, SSE.Full = sse1)
        #new:
        #	rslt2 <- data.frame("SSE.Linear" = sse0, "SSE.Full" = sse1)
        #	rownames(rslt2)[1] <- ""
        #	ans <- list(call = object$call, coefficients = rslt1, residuals = rslt2, zcoefs = coefs, cic = cic, null_deviance = null_deviance, null_df = null_df, deviance = deviance, df = df, resid_df_obs = resid_df_obs, family = family)
        #	class(ans) <- "summary.cgam"
        #	ans
        #} else {
        ans <- list(call = object$call, coefficients = rslt1, zcoefs = coefs, resid_df_obs = resid_df_obs, family = family)
        class(ans) <- "summary.cgamm"
        ans
        #}
    } else {
        ans <- list(zcoefs = object$zcoefs)
        class(ans) <- "summary.cgamm"
        ans
    }
}




######################
#print.summary.cgamm #
######################
print.summary.cgamm <- function(x,...) {
    if (!is.null(x$zcoefs)) {
        #if (!is.null(x$se.beta)) {
        cat("Call:\n")
        print(x$call)
        cat("\n")
        cat("Coefficients:")
        cat("\n")
        printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE)
        #cat("\n")
        #if (x$family$family == "binomial") {
        #    cat("(Dispersion parameter for binomial family taken to be 1)", "\n")
        #}
        #if (x$family$family == "poisson") {
        #    cat("(Dispersion parameter for poisson family taken to be 1)", "\n")
        #}
        #if (x$family$family == "gaussian") {
        #    cat("(Dispersion parameter for gaussian family taken to be ", round(x$deviance/x$df,4),")","\n", sep="")
        #}
        #cat("\n")
        #cat("Null deviance: ",round(x$null_deviance,4), "", "on", x$null_df, "", "degrees of freedom", "\n")
        #cat("Residual deviance: ",round(x$deviance,4), "", "on", x$resid_df_obs, "", "observed degrees of freedom", "\n")
        #if (is.null(x$cic)) {
        #	message("CIC value is not available when there is no shape-restricted predictor")
        #} else {message("CIC: ", round(x$cic,4))}
        #if (!is.null(x$coefficients2)) {
        #    cat("\n")
        #    cat("Approximate significance of smooth terms: \n")
        #    printCoefmat(x$coefficients2, P.values = TRUE, has.Pvalue = TRUE)
        #}
        #if (!is.null(x$cic)) {
        #    cat("CIC: ", round(x$cic,4))
        #}
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


###########################################
#create a 3D plot for a cgam or wps object#
###########################################
plotpersp <- function(object,...) {
  UseMethod("plotpersp", object)
}
# 
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
# plotpersp.cgam <- function(object, x1 = NULL, x2 = NULL, x1nm = NULL, x2nm = NULL,
#                            data = NULL, surface = "mu", categ = NULL, col = NULL,
#                            random = FALSE, ngrid = 12, xlim = range(x1), 
#                            ylim = range(x2), zlim = NULL, xlab = NULL, ylab = NULL, 
#                            zlab = NULL, th = NULL, ltheta = NULL, main = NULL, ticktype = "simple",...) {
#     #if (class(object) == "list") {
#     #	object <- object[[1]]
#     #}
#     #print (class(object))
#     if (!inherits(object, "cgam")) {
#         warning("calling plotpersp(<fake-cgam-object>) ...")
#     }
#     #xmat <- object$xmat
#     #cl <- match.call()
#     #nms <- cl[-c(1, 2)]
#     #lnms <- length(nms)
#     #x1nm <- nms[1]$x
#     #x1nm <- deparse(x1nm)
#     #x2nm <- nms[2]$x
#     #x2nm <- deparse(x2nm)
#     #new: default is plotpersp(object)
#     #x1nm <- deparse(substitute(x1))
#     #x2nm <- deparse(substitute(x2))
#     #print (x1nm)
#     #print (x2nm)
#     #stop (print (x1nm))
#     xnms <- object$xnms_add
#     xmat <- object$xmat_add
#     #if (x1nm == "NULL" | x2nm == "NULL") {
#     if (is.null(x1nm) | is.null(x2nm)) {
#         if (length(xnms) >= 2) {
#             x1nm <- xnms[1]
#             x2nm <- xnms[2]
#             x1id <- 1
#             x2id <- 2
#             x1 <- xmat[, 1]
#             x2 <- xmat[, 2]
#         } else {stop ("Number of non-parametric predictors must >= 2!")}
#     }
#     ynm <- object$ynm
#     #xmat <- object$xmat
#     #print (dim(xmat))
#     is_fac <- object$is_fac
#     is_param <- object$is_param
#     family <- object$family
#     fml <- family$family
#     cicfamily <- CicFamily(family)
#     muhat.fun <- cicfamily$muhat.fun
#     znms <- object$znms
#     kznms <- length(znms)
#     if (!is.null(categ)) {
#         if (!is.character(categ)) {
#             warning("categ must be a character argument!")
#         } else if (!any(znms == categ)) {
#             #in.or.out case
#             #if (!is.null(attr(object, "sub"))) {
#             #	if (!is.null(categ)) {
#             #		categ = paste("factor(", categ, ")", sep = "")
#             #	}
#             #} else {
#             #	warning(paste(categ, "is not an exact character name defined in the cgam fit!"))
#             #	categ = NULL
#             #}
#             if (any(grepl(categ, znms))) {
#                 id = which(grepl(categ, znms))
#                 znmi = znms[id]
#                 if (grepl("as.factor", znmi)) {
#                     categ = paste("as.factor(", categ, ")", sep = "")
#                 } else if (grepl("factor", znmi)) {
#                     categ = paste("factor(", categ, ")", sep = "")
#                 } else {print(paste(categ, "is not an exact character name defined in the cgam fit!"))}
#             } else {print(paste(categ, "is not an exact character name defined in the cgam fit!"))}
#         } else {
#             obsz = 1:kznms
#             zid = obsz[znms == categ]
#             #linear term:
#             if (!(is_fac[zid])) {
#                 categ = NULL
#             }
#         }
#     }
#     shapes <- object$shapes
#     #new:
#     #zid1 = object$zid1 - 1 - length(shapes)
#     #zid2 = object$zid2 - 1 - length(shapes)
#     zid1 <- object$zid1
#     zid2 <- object$zid2
#     kznms <- length(znms)
#     zmat <- object$zmat
#     if (any(class(object) == "wps")) {
#         d0 <- object$d0
#         np_add <- object$np_add
#         p <- d0 - np_add
#         pb <- object$pb
#         zmat <- zmat[, (pb+1):(pb+p), drop = FALSE]
#         #remove the one
#         zmat <- zmat[, -1, drop = FALSE]
#         #print (head(zmat))
#     }
#     #zcoefs = object$zcoefs
#     #new: exclude the coef for the one vector
#     #temp:trispl not include one
#     if (fml != "ordered" & all(class(object) != "trispl")) {
#         zcoefs <- object$zcoefs[-1]
#     } else {
#         zcoefs <- object$zcoefs
#     }
#     #print (zcoefs)
#     #xmatnms <- object$xmatnms
#     knms <- length(xnms)
#     obs <- 1:knms
#     #if (!any(xmatnms == x1nm)) {
#     #	warning(paste(x1nm, "is not an exact character name defined in the cgam fit!"))
#     #}
#     #if (!any(xmatnms == x2nm)) {
#     #	warning(paste(x2nm, "is not an exact character name defined in the cgam fit!"))
#     #}
#     #x1id = obs[xmatnms == x1nm]
#     #x2id = obs[xmatnms == x2nm]
#     if (!is.null(data)) {
#         if (!is.data.frame(data)) {
#             stop ("User need to make the data argument a data frame with names for each variable!")
#         }
#         datnms <- names(data)
#         if (!any(datnms == x1nm) | !any(datnms == x2nm)) {
#             stop ("Check the accuracy of the names of x1 and x2!")
#         }
#         x1 <- data[ ,which(datnms == x1nm)]
#         x2 <- data[ ,which(datnms == x2nm)]
#         if (length(x1) != nrow(xmat)) {
#             warning ("Number of observations in the data set is not the same as the number of elements in x1!")
#         }
#         #bool <- apply(xmat, 2, function(x) all(x1 == x))
#         #if (any(bool)) {
#         x1id <- obs[xnms == x1nm]
#         #}
#         if (length(x2) != nrow(xmat)) {
#             warning ("Number of observations in the data set is not the same as the number of elements in x2!")
#         }
#         #bool <- apply(xmat, 2, function(x) all(x2 == x))
#         #if (any(bool)) {
#         x2id <- obs[xnms == x2nm]
#         #}
#     } else {
#         #if (any(xmatnms == x1nm)) {
#         #	x1id <- obs[xmatnms == x1nm]
#         #} else {
#         #	bool <- apply(xmat, 2, function(x) all(x1 == x))
#         #	if (any(bool)) {
#         #		x1id <- obs[bool]
#         #	}
#         #}
#         #if (any(xmatnms == x2nm)) {
#         #	x2id <- obs[xmatnms == x2nm]
#         #} else {
#         #	bool <- apply(xmat, 2, function(x) all(x2 == x))
#         #	if (any(bool)) {
#         #		x2id <- obs[bool]
#         #	}
#         #}
#         #new: x1 and x2 are in .Globe, not in formula
#         if (all(xnms != x1nm)) {
#             if (length(x1) != nrow(xmat)) {
#                 stop ("Number of observations in the data set is not the same as the number of elements in x1!")
#             }
#             bool <- apply(xmat, 2, function(x) all(x1 == x))
#             if (any(bool)) {
#                 x1id <- obs[bool]
#                 #change x1nm to be the one in formula
#                 x1nm <- xnms[bool]
#             } else {
#                 stop (paste(paste("'", x1nm, "'", sep = ''), "is not a predictor defined in the cgam fit!"))
#             }
#         } else {
#             x1id <- obs[xnms == x1nm]
#         }
#         if (all(xnms != x2nm)) {
#             if (length(x2) != nrow(xmat)) {
#                 stop ("Number of observations in the data set is not the same as the number of elements in x2!")
#             }
#             bool <- apply(xmat, 2, function(x) all(x2 == x))
#             if (any(bool)) {
#                 x2id <- obs[bool]
#                 x2nm <- xnms[bool]
#             } else {
#                 stop (paste(paste("'", x2nm, "'", sep = ''), "is not a predictor defined in the cgam fit!"))
#             }
#         } else {
#             x2id <- obs[xnms == x2nm]
#         }
#     }
#     #xmat is not the one in fit
#     #print (length(x1))
#     #print (length(x2))
#     #xm <- cbind(x1, x2)
#     xm <- xmat[, c(x1id, x2id)]
#     #print (all(xm == cbind(x1, x2)))
#     #print (head(cbind(x1, x2)))
#     x_grid <- ngrid
#     y_grid <- ngrid
#     x1g <- 0:x_grid / x_grid * .95 * (max(xm[,1]) - min(xm[,1])) + min(xm[,1]) + .025 * (max(xm[,1]) - min(xm[,1]))
#     n1 <- length(x1g)
#     x2g <- 0:y_grid / y_grid * .95 * (max(xm[,2]) - min(xm[,2])) + min(xm[,2]) + .025 * (max(xm[,2]) - min(xm[,2]))
#     n2 <- length(x2g)
#     xgmat <- matrix(nrow = n1, ncol = n2)
#     eta0 <- object$coefs[1]
#     thvecs <- object$etacomps
#     #print ('thvecs')
#     for (i2 in 1:n2) {
#         for (i1 in 1:n1) {
#             x1a <- max(xm[xm[,1] <= x1g[i1], 1])
#             x1b <- min(xm[xm[,1] > x1g[i1], 1])
#             v1a <- min(thvecs[x1id, xm[,1] == x1a])
#             v1b <- min(thvecs[x1id, xm[,1] == x1b])
#             alp <- (x1g[i1] - x1a) / (x1b - x1a)
#             th1add <- (1 - alp) * v1a + alp * v1b
#             x2a <- max(xm[xm[,2] <= x2g[i2],2])
#             x2b <- min(xm[xm[,2] > x2g[i2],2])
#             v2a <- min(thvecs[x2id, xm[,2] == x2a])
#             v2b <- min(thvecs[x2id, xm[,2] == x2b])
#             alp <- (x2g[i2] - x2a) / (x2b - x2a)
#             th2add <- (1 - alp) * v2a + alp * v2b
#             xgmat[i1,i2] <- eta0 + th1add + th2add
#         }
#     }
#     x3_add <- 0
#     if (knms >= 3) {
#         x3id <- obs[-c(x1id, x2id)]
#         kx3 <- length(x3id)
#         for (i in 1:kx3) {
#             x3i <- xmat[, x3id[i]]
#             x3i_use <- max(x3i[x3i <= median(x3i)])
#             x3i_add <- min(thvecs[x3id[i], x3i == x3i_use])
#             x3_add <- x3_add + x3i_add
#         }
#     }
#     if (surface == "eta") {
#         xgmat <- xgmat + as.numeric(x3_add)
#     }
#     if (is.null(categ) & surface == "mu") {
#         z_add <- 0
#         if (!is.null(znms)) {
#             kzids <- length(zid1)
#             for (i in 1:kzids) {
#                 pos1 <- zid1[i]; pos2 <- zid2[i]
#                 zi <- zmat[, pos1:pos2, drop = FALSE]
#                 zcoefsi <- zcoefs[pos1:pos2]
#                 for (j in 1:ncol(zi)){
#                     uzij <- unique(zi[,j])
#                     kuzij <- length(uzij)
#                     nmodej <- sum(zi[,j] == uzij[1])
#                     zij_mode <- uzij[1]
#                     for (u in 2:kuzij) {
#                         if (sum(zi[,j] == uzij[u]) > nmodej) {
#                             zij_mode <- uzij[u]
#                             nmodej <- sum(zi[,j] == uzij[u])
#                         }
#                     }
#                     obsuzij <- 1:length(uzij)
#                     uzhatij <- uzij * zcoefsi[j]
#                     zij_add <- uzhatij[obsuzij[uzij == zij_mode]]
#                     z_add <- z_add + zij_add
#                 }
#             }
#         }
#         xgmat <- xgmat + as.numeric(x3_add) + as.numeric(z_add)
#         #xgmat <- muhat.fun(xgmat, fml = fml)
#         if (fml != "gaussian" & fml != "ordered") {
#             for (i2 in 1:n2) {
#                 for (i1 in 1:n1) {
#                     xgmat[i1, i2] <- muhat.fun(xgmat[i1, i2], fml = fml)
#                 }
#             }
#         }
#     } else if (!is.null(categ) & surface == "mu"){
#         xgmats <- list()
#         mins <- NULL; maxs <- NULL
#         obsz <- 1:kznms
#         zid <- obsz[znms == categ]
#         #print (class(znms[znms == categ]))
#         pos1 <- zid1[zid]; pos2 <- zid2[zid]
#         #print (pos1)
#         #print (pos2)
#         #zi <- zmat[, pos1:pos2, drop = FALSE]
#         #z_add <- 1:nrow(zi)*0
#         #zcoefsi <- zcoefs[pos1:pos2]
#         #print (zcoefsi)
#         zcoefsi = zcoefs[pos1:pos2]
#         #include the base level
#         zcoefsi = c(0, zcoefsi)
#         z_add = sort(zcoefsi)
#         kz_add <- length(z_add)
#         #for (j in 1:ncol(zi)) {
#         #	zij <- zi[,j]
#         #	zijhat <- zij * zcoefsi[j]
#         #	z_add <- z_add + zijhat
#         #}
#         #z_add <- unique(z_add)
#         #kz_add <- length(z_add)
#         #new: plot the smallest one first:
#         #z_add <- z_add[order(z_add)]
#         #print (z_add)
#         for (iz in 1:kz_add) {
#             xgmats[[iz]] <- xgmat + as.numeric(x3_add) + z_add[iz]
#             #mins <- c(mins, min(xgmats[[iz]]))
#             #maxs <- c(maxs, max(xgmats[[iz]]))
#             if (fml != "gaussian" & fml != "ordered") {
#                 for (i2 in 1:n2) {
#                     for (i1 in 1:n1) {
#                         xgmats[[iz]][i1, i2] <- muhat.fun(xgmats[[iz]][i1, i2], fml = fml)
#                     }
#                 }
#             }
#             #xgmat[[iz]] <- muhat.fun(xgmat[[iz]], fml = fml)
#             mins <- c(mins, min(xgmats[[iz]]))
#             maxs <- c(maxs, max(xgmats[[iz]]))
#         }
#     }
#     if (is.null(xlab)) {
#         #xlab = deparse(x1nm)
#         xlab <- x1nm
#     }
#     if (is.null(ylab)) {
#         #ylab = deparse(x2nm)
#         ylab <- x2nm
#     }
#     if (is.null(zlab)) {
#         if (surface == "mu") {
#             if (fml == "binomial") {
#                 zlab <- paste("Pr(", ynm, ")")
#             } else if (fml == "poisson" | fml == "gaussian" | fml == "Gamma") {
#                 zlab <- paste("Est mean of", ynm)
#             }
#         }
#         if (surface == "eta") {
#             if (fml == "binomial") {
#                 zlab <- paste("Est log odds ratio of", ynm)
#             }  else if (fml == "poisson" | fml == "Gamma") {
#                 zlab <- paste("Est log mean of", ynm)
#             } else if (fml == "gaussian") {
#                 zlab <- paste("Est mean of", ynm)
#             }
#         }
#     }
#     #if (is.null(zlim)) {
#     #	zlim <- range(xgmat, na.rm = TRUE)
#     #}
#     palette <- c("peachpuff", "lightblue", "limegreen", "grey", "wheat", "yellowgreen", "seagreen1", "palegreen", "azure", "whitesmoke")
#     if (!is.null(categ) & surface == "mu") {
#         #palette = c("peachpuff", "lightblue", "grey", "wheat", "yellowgreen", "plum", "limegreen", "paleturqoise", "azure", "whitesmoke")
#         kxgm <- length(xgmats)
#         if (is.null(col)) {
#             #if (kxgm == 2) {
#             #	col = c("peachpuff", "lightblue")
#             #} else if (kxgm == 3) {
#             #	col = c("peachpuff", "lightblue", "grey")
#             #} else if (kxgm > 3 & kxgm < 11) {
#             #	col = sample(palette, replace = FALSE)
#             if (random) {
#                 col <- topo.colors(kxgm)
#                 #col <- sample(palette, size = kxgm, replace = FALSE)
#                 #print (col)
#             } else {
#                 #print (kxgm)
#                 if (kxgm > 1 & kxgm < 11) {
#                     col <- palette[1:kxgm]
#                 } else {
#                     #integ <- floor(kxgm / 10)
#                     #rem <- kxgm %% 10
#                     #kint <- length(integ)
#                     #col = character(length = kxgm)
#                     #print (col)
#                     #col <- NULL
#                     #for (i in 1:kint) {
#                     #print (col[1 + (i - 1) * 10: i * 10])
#                     #print (palette)
#                     #col[1 + (i - 1) * 10: i * 10] = palette
#                     #	col <- c(col, palette)
#                     #}
#                     #print (col)
#                     #col[(kint * 10 + 1):kxgm] = palette[(kint * 10 + 1):kxgm]
#                     #col <- c(col, palette[1:rem])
#                     #print ((kint * 10 + 1):kxgm)
#                     #print (col)
#                     #print (integ)
#                     #new: use rainbow
#                     col <- topo.colors(kxgm)
#                 }
#             }
#         } else {
#             col0 <- col
#             if (col0 == "heat" | col0 == "topo" | col0 == "terrain" | col0 == "cm") {
#                 #col0 <- col
#                 ncol <- 100
#                 facets <- facetcols <- list()
#                 col <- list()
#                 for (i in 1:kxgm) {
#                     nr <- nrow(xgmats[[i]])
#                     nc <- ncol(xgmats[[i]])
#                     facets[[i]] <- (xgmats[[i]])[-1,-1] + (xgmats[[i]])[-1,-nc] + (xgmats[[i]])[-nr,-1] + (xgmats[[i]])[-nr,-nc]
#                     facetcols[[i]] <- cut(facets[[i]], ncol)
#                     #print (head(facetcols[[i]]))
#                     if (col0 == "heat") {
#                         col[[i]] <- (heat.colors(ncol))[facetcols[[i]]]
#                         #print (head(col[[i]]))
#                     } else if (col0 == "topo") {
#                         col[[i]] <- (topo.colors(ncol))[facetcols[[i]]]
#                     } else if (col0 == "terrain") {
#                         col[[i]] <- (terrain.colors(ncol))[facetcols[[i]]]
#                     } else {
#                         col[[i]] <- (cm.colors(ncol))[facetcols[[i]]]
#                     }
#                 }
#             } else if (length(col0) < kxgm) {
#                 rem <- kxgm - length(col0)
#                 nrem <- length(rem)
#                 rem_col <- palette[1:nrem]
#                 col <- c(col0, rem_col)
#                 #new:
#                 #nr <- nrow(xgmat)
#                 #nc <- ncol(xgmat)
#                 #ncol <- 100
#                 #facet <- xgmat[-1,-1] + xgmat[-1,-nc] + xgmat[-nr,-1] + xgmat[-nr,-nc]
#                 #facetcol <- cut[facet, ncol]
#                 #col <- topo.colors[facetcol]
#                 
#             } else if (length(col0) > kxgm) {
#                 col <- col0[1:kxgm]
#                 #print (paste("The first", kxgm, "colors are used!"))
#             }
#             if (random) {
#                 print ("User defined colors are used!")
#             }
#         }
#         #print (col[[1]][1:10])
#         #print (kxgm)
#         #new: set th for decr or incr
#         decrs = shapes[c(x1id, x2id)] %in% c(2, 6, 8, 10, 15, 16)
#         incrs = shapes[c(x1id, x2id)] %in% c(1, 5, 7, 9, 13, 14)
#         if (is.null(th) | !is.numeric(th)) {
#             ang = NULL
#             if (incrs[1] & incrs[2]) {
#                 if (is.null(ang)) {
#                     ang = -40
#                 }
#             } else if (decrs[1] & incrs[2]) {
#                 if (is.null(ang)) {
#                     ang = 40
#                 }
#             } else if (incrs[1] & decrs[2]) {
#                 if (is.null(ang)) {
#                     ang = 240
#                 }
#             } else if (decrs[1] & decrs[2]) {
#                 if (is.null(ang)) {
#                     ang = 140
#                 }
#             } else {ang = -37}
#         } else {ang = th}
#         if (is.null(ltheta) | !is.numeric(ltheta)) {
#             ltheta <- -135
#         }
#         #print (col[[1]][1:10])
#         for (i in 1:kxgm) {
#             #print (col[i])
#             #i = 1
#             #print (i)
#             #print (length(col))
#             xgmat <- xgmats[[i]]
#             #if (is.null(th) | !is.numeric(th)) {
#             #	th <- -40
#             #}
#             #if (is.null(ltheta) | !is.numeric(ltheta)) {
#             #	ltheta <- -135
#             #}
#             #persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, theta = th)
#             #new: avoid thick labs
#             box = TRUE
#             axes = TRUE
#             if (i > 1) {
#                 xlab = ylab = zlab = " "
#                 box = FALSE
#                 axes = FALSE
#             }
#             #print (length(col))
#             if (is.list(col)) {
#                 coli = unlist(col[[i]])
#                 #print (head(coli))
#             } else {coli = col[i]}
#             #print (head(coli))
#             #print (head(xgmat[[i]]))
#             #print (head(col[[1]]))
#             if (is.null(zlim)) {
#                 lwr = min(mins)
#                 upp = max(maxs)
#                 zlim0 = c(lwr - (upp-lwr)/3, upp + (upp-lwr)/3)
#             } else {
#                 zlim0 = zlim
#             }
#             persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab, col = coli, main = main, theta = ang, ltheta = ltheta, ticktype = ticktype, box = box, axes = axes,...)
#             par(new = TRUE)
#         }
#         par(new = FALSE)
#     } else {
#         if (is.null(col)) {
#             if (random) {
#                 col <- sample(palette, size = 1, replace = FALSE)
#             } else {
#                 #col <- "white"
#                 #col <- color[facetcol]
#                 nr <- nrow(xgmat)
#                 nc <- ncol(xgmat)
#                 ncol <- 100
#                 facet <- xgmat[-1,-1] + xgmat[-1,-nc] + xgmat[-nr,-1] + xgmat[-nr,-nc]
#                 #print (facet)
#                 facetcol <- cut(facet, ncol)
#                 col <- heat.colors(ncol)[facetcol]
#             }
#         } else {
#             #if (length(col) > 1) {
#             #	col <- col[1]
#             #	print ("The first color is used!")
#             #	col <- heat.colors(x_grid*y_grid)
#             #}
#             if (col == "heat" | col == "topo" | col == "terrain" | col == "cm") {
#                 nr <- nrow(xgmat)
#                 nc <- ncol(xgmat)
#                 ncol <- 100
#                 facet <- xgmat[-1,-1] + xgmat[-1,-nc] + xgmat[-nr,-1] + xgmat[-nr,-nc]
#                 facetcol <- cut(facet, ncol)
#                 if (col == "heat") {
#                     col <- heat.colors(ncol)[facetcol]
#                 } else if (col == "topo") {
#                     col <- topo.colors(ncol)[facetcol]
#                 } else if (col == "terrain") {
#                     col <- terrain.colors(ncol)[facetcol]
#                 } else {
#                     col <- cm.colors(ncol)[facetcol]
#                 }
#             } 
#             if (random) {
#                 print ("User defined color is used!")
#             } 
#         }
#         #if (is.null(th) | !is.numeric(th)) {
#         #	th <- -40
#         #}
#         #if (is.null(ltheta) | !is.numeric(ltheta)) {
#         #	ltheta <- -135
#         #}
#         #new: set th for decr or incr
#         decrs = shapes[c(x1id, x2id)] %in% c(2, 6, 8, 10, 15, 16)
#         incrs = shapes[c(x1id, x2id)] %in% c(1, 5, 7, 9, 13, 14)
#         if (is.null(th) | !is.numeric(th)) {
#             ang = NULL
#             if (incrs[1] & incrs[2]) {
#                 if (is.null(ang)) {
#                     ang = -40
#                 }
#             } else if (decrs[1] & incrs[2]) { 
#                 if (is.null(ang)) {
#                     ang = 40
#                 }
#             } else if (incrs[1] & decrs[2]) {
#                 if (is.null(ang)) {
#                     ang = 240
#                 }
#             } else if (decrs[1] & decrs[2]) {
#                 if (is.null(ang)) {
#                     ang = 140
#                 }
#             } else {ang = -37}
#         } else {ang = th}
#         if (is.null(ltheta) | !is.numeric(ltheta)) {
#             ltheta <- -135
#         }
#         if (is.null(zlim)) {
#             lwr = min(xgmat)
#             upp = max(xgmat)
#             zlim0 = c(lwr - (upp-lwr)/3, upp + (upp-lwr)/3)
#         } else {
#             zlim0 = zlim
#         }
#         persp(x1g, x2g, xgmat, xlim = xlim, ylim = ylim, zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab, col = col, main = main, theta = ang, ltheta = ltheta, ticktype = ticktype,...)
#         rslt = list(zlim = zlim0, xlab = xlab, ylab = ylab, zlab = zlab, theta = ang, ltheta = ltheta, col = col, cex.axis = .75, main = main, ticktype = ticktype, z_add = z_add, x3_add = x3_add)
#         invisible(rslt)
#     }
#     #print (col)
# }

#############################################################
#apply plotpersp on a predict.cgam or predict.cgamm object
#############################################################
plotpersp.cgamp = function(object, x1=NULL, x2=NULL, x1nm=NULL, x2nm=NULL, data=NULL, up = TRUE, main=NULL, cex.main=.8, xlab = NULL, ylab = NULL, zlab = NULL, zlim = NULL, th = NULL, ltheta = NULL, ticktype = "detailed",...) {
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
    persp(x1g, x2g, xgmat, zlim = res$zlim, xlab = "", ylab = "", zlab = "", theta = res$theta, ltheta = res$ltheta, cex.axis = res$cex.axis, main = main, cex.main = cex.main, ticktype = res$ticktype, col=mycol, box=FALSE, axes=FALSE)
    par(new=FALSE)
}



