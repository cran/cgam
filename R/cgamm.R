#######
#cgamm#
#######
cgamm = function(formula, nsim = 0, family = gaussian, cpar = 1.2, data = NULL, weights = NULL, sc_x = FALSE, sc_y = FALSE, bisect = TRUE, reml = TRUE) {
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
    mf[[1]] <- quote(lme4::lFormula)
    mf <- eval(mf, parent.frame(1L))$fr
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
#test:
    id <- mf[,nc]
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
    rslt <- list(muhat = ans$muhat, coefs = ans$coefs, bh = ans$bh, zcoefs = ans$zcoefs, vcoefs = ans$vcoefs, ahat = ans$ahat, sig2hat = ans$sig2hat, siga2hat = ans$siga2hat, thhat = ans$thhat, bigmat = ans$bigmat, id = id, szs = szs, shapes = shapes0, numknots = ans$numknots, knots = ans$knots, space = sps0, d0 = ans$np, xmat_add = xmat, xmat0 = ans$xmat2, knots0 = ans$knots2, numknots0 = ans$numknots2, sps0 = ans$sps2, ms0 = ans$ms2, etacomps = ans$etacomps, xnms_add = xnms, xid1 = ans$xid1, xid2 = ans$xid2, ynm = ynm, y = y, znms = znms, zmat = zmat, zid = zid, zid1 = zid1, zid2 = zid2, family = family, is_fac = is_fac, is_param = is_param, tms = mt, capm = ans$capm, capms = ans$capms, capk = ans$capk, capt = ans$capt, capu = ans$capu, mod.lmer = ans$mod.lmer, pv.siga2 = ans$pv.siga2, ci.siga2 = ans$ci.siga2, ci.th = ans$ci.th, ci.rho = ans$ci.rho, ci.sig2 = ans$ci.sig2, ones = ans$ones)
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
		if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk > 0) {
			bigmat = rbind(1:n*0 + 1, t(zmat), t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
			np = 1 + capk + sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13)  + capms
		} else if (sum(shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13) > 0 & capk == 0) {
			bigmat = rbind(1:n*0 + 1, t(xmat[, shapes > 2 & shapes < 5 | shapes > 10 & shapes < 13]), delta)
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
            ans = coneB(y, dsend, zsend)
            edf = ans$df
            face = ans$face
            bh = coef(ans)
            #test:
            if (round(bh[1:np],6) < 0) {
                face = c(1:np, face)
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
                if (round(bh[1:np],6) < 0) {
                    face = c(1:np, face)
                }
            } else {
                bh = solve(crossprod(gtil), t(gtil)) %*% ytil
                edf = nrow(bigmat)
                face = 1:edf
            }
            muhat = t(bigmat) %*% bh
			diff = mean((oldmu - muhat)^2)
			oldmu = muhat
            if (reml) {
                dd = t(bigmat[face, ,drop = FALSE])
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
#new: tests and cis
#F-test for siga2
    pv.siga2 = ranef.test(ecl, szs, n, ncl)
#c.i. for siga2
    ci1 = NULL
    if (type == 'b') {
        ci1 = ranef.ci(ecl, szs, n, ncl, level = 0.95)
    }
    ci.siga2 = ci1
#c.i. for th and rho
    ci.th = ci.rho = NULL
    if (bisect) {
        #print (edf)
        ci2 = ranef.cith(thhat, sig2hat, siga2hat, ahat, ecl, szs, n, ncl, level = 0.95, xms=xms, p=edf, reml=reml)
        ci.th = ci2$ci.th
        ci.rho = ci2$ci.rho
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
    if (capk > 0) {
        zcoefs = bh[2:(1+capk)]
    } else {zcoefs = NULL}
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
    rslt = list(muhat = muhat, coefs = coefs, bh = bh, zcoefs = zcoefs, vcoefs = vcoefs, ahat = ahat, sig2hat = sig2hat, siga2hat = siga2hat, thhat = thhat, bigmat = bigmat, np = np, knots = knotsuse, knots2 = knotsuse2, numknots = numknotsuse, numknots2 = numknotsuse2, ms = mslst, ms2 = mslst2, xmat2 = xmat2, xid1 = xid1, xid2 = xid2, capm = capm, capms = capms, capk = capk, capt = capt, capu = capu, etacomps = thvecs, mod.lmer = mod.lmer, pv.siga2 = pv.siga2, ci.siga2 = ci.siga2, ci.th = ci.th, ci.rho = ci.rho, ci.sig2 = ci.sig2, ones = ones)
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
#C.I. for sig2                                   #
#balanced and unbalanced                         #
##################################################
ranef.cisig2 = function(ecl, N, ncl, level = 0.95) {
    evec = unlist(ecl)
    ebars = sapply(ecl, mean)
    sse = 0
    #ssb = 0
    ebar = mean(evec)
    for(icl in 1:ncl) {
        #sz = szs[icl]
        eibar = ebars[[icl]]
        ei = ecl[[icl]]
        ssei = sum((ei - eibar)^2)
        sse = sse + ssei
        #ssbi = sz*(eibar - ebar)^2
        #ssb = ssb + ssbi
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


#################
#new coef method#
#################
coef.cgamm <- function(object,...) {
    ans <- object$coefs
    ans
}

#ranef.cgamm <- function(object,...) {
#    ans <- object$ahat
#    ans
#}

fixef.cgamm <- function(object,...) {
    ans <- object$bh
    ans
}


###############
#predict.cgamm#
#smooth only  #
#guasisan only#
###############
predict.cgamm = function(object, newData, interval = c("confidence", "none"), type = c("response", "link"), level = 0.95, ...) {
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
    newData[[idlb]] = rep(1, nrow(newData))
    m = model.frame(Terms, newData)
    #print (m)
    #new: don't include id in newdata
    rm_id = NULL
    for (i in 1:ncol(m)) {
        mi = m[,i]
        if (is.null(attributes(mi)$shape)) {
           rm_id = i
        }
    }
    newdata = m[, -rm_id, drop=F]
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
#different  from cgam: simpler
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
            balanced = FALSE
            if (length(unique(szs)) == 1) {
                balanced = TRUE
            }
            ncl = length(szs)
            edges = t(bigmat)
#new:
            nd = ncol(edges)
            #face = 1:nd
            face = NULL
            nloop = 1000
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
                sector = 1:nsec*0
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
                    sec = 1:m_acc*0
                    sec[cf>0] = 1
                    r = makebin(sec)+1
                    sector[r] = sector[r]+1
                }
            } else {
                prior.w = 1:n*0 + 1
                sector = 1:nsec*0
                for (iloop in 1:nloop) {
                    ysim = muhat + rnorm(n)*sig2hat^.5
                    ysim = ysim * sqrt(prior.w)
                    dsend = edges[, (np+1):(m_acc+np)]
                    zsend = edges[, 1:np]
                    ans = coneB(ysim, dsend, zsend, face = NULL)
                    #face = ans$face
                    cf = round(ans$coefs[(np+1):(m_acc+np)], 10)
                    sec = 1:m_acc*0
                    sec[cf > 0] = 1
                    r = makebin(sec) + 1
                    sector[r] = sector[r] + 1
                }
            }
            ### calculate the mixture hat matrix:
            bsec = matrix(0, nrow=nsec, ncol=2)
            bsec[,1] = 1:nsec
            bsec[,2] = sector / nsec
            keep = sector > 0
            bsec = bsec[keep,]
            ns = dim(bsec)[1]
            bsec[,2] = bsec[,2] / sum(bsec[,2])
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
            #only xmatpr and splpr have interpolation points
            #xmatpr = cbind(newv, t(splpr))
            #temp:
            #muhatpr = xmatpr %*% coefs[1:ncol(xmatpr), ,drop=FALSE]
            #new: C.I. level
            mult = qnorm((1 - level)/2, lower.tail=FALSE)
            hl = mult*sqrt(diag(xmatpr%*%acov%*%t(xmatpr)))
            ans = list(fit = muhatpr, lower = muhatpr - hl, upper = muhatpr + hl)
        }
    }
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

