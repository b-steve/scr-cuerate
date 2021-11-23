## Main fitting function.
## - Arguments:
##   - capthist:     Capture history matrix.
##   - ids:          A vector of IDs, where each element indicates which individual is associated with each capture history.
##   - traps:        Matrix of detector locations.
##   - mask:         Mask object with points for numerical integration over activity centre locations. Must include an attribute called "area", providing the area covered by a single point.
##   - detfn:        The detection function to use; either "hn" (halfnormal) or "hhn" (hazard halfnormal).
##   - start:        Start values for numerical maximisation. Order of parameters is D, lambda0/g0, sigma, lambda_c, sigma_t. 
##   - toa:          Time of arrival matrix in the same structure as the capthist object (optional).
##   - speed_sound:  The speed of sound in metres per second.
cuerate.scr.fit <- function(capthist, ids, traps, mask, detfn = NULL, start, ss = NULL, toa = NULL, speed_sound = 330, trace = FALSE){
    ## Indicator for whether or not signal strengths are used.
    use_ss <- !is.null(ss)
    ## Indicator for whether or not times of arrival are used.
    use_toa <- !is.null(toa)
    if (is.list(capthist)){
        multi.sess <- TRUE
    } else {
        multi.sess <- FALSE
        capthist <- list(capthist)
        ids <- list(ids)
        traps <- list(traps)
        mask <- list(mask)
        if (use_toa){
            toa <- list(toa)
        }
        if (use_ss){
            ss <- list(ss)
        }
    }
    n.sessions <- length(capthist)
    aMask <- maskDists <- toa_ssq <- vector("list", n.sessions)
    if (!use_toa){
        toa <- vector("list", n.sessions)
    }
    if (!use_ss){
        ss <- vector("list", n.sessions)
    }
    for (i in 1:n.sessions){
        ## ids must be 1:N animals.
        ids[[i]] = as.numeric(factor(ids[[i]]))
        ## Area of a single mask point.
        aMask[[i]] <- attr(mask[[i]], "area")
        
        ## Calculating distances between mask points and detectors.
        maskDists[[i]] <- eucdist(mask[[i]], traps[[i]])
        if (use_toa){
            ## Creating TOA sum of squares matrix.
            toa_ssq[[i]] <- make_toa_ssq(toa[[i]], eucdist(traps[[i]], mask[[i]]), speed_sound)
        } else {
            ## Dummy objects if not used.
            toa[[i]] <- toa_ssq[[i]] <- matrix(0, nrow = 1, ncol = 1)
        }
        if (!use_ss){
            ss[[i]] <- matrix(0, nrow = 1, ncol = 1)
        }
    }
    ## Indicator for detection function.
    if (is.null(detfn) & !use_ss){
        stop("A detection function must be selected.")
    }
    if (use_ss){
        if (!is.null(detfn)){
            warning("The choice of detection function is being ignored because signal strengths have been provided.")
        }
        detfn <- "ss"
        hn <- FALSE
    } else if (detfn == "hn"){
        hn <- TRUE
    } else if (detfn == "hhn"){
        hn <- FALSE
    } else {
        stop("The argument detfn must either be 'hn' or 'hhn'.")
    }
    ## Converting parameters to link scale.
    start.link <- numeric(6)
    names(start.link) <- c("D", "df1", "df2", "lambda", "sigma_toa", "sigma_ss")
    start.link["D"] <- log(start["D"])
    start.link["lambda"] <- log(start["lambda"])
    if (use_ss){
        names(start.link)[c(2, 3)] <- c("b0_ss", "b1_ss")
        start.link["b0_ss"] <- log(start["b0_ss"])
        start.link["b1_ss"] <- log(start["b1_ss"])
        start.link["sigma_ss"] <- log(start["sigma_ss"])
    } else if (hn){
        names(start.link)[c(2, 3)] <- c("g0", "sigma")
        start.link["g0"] <- qlogis(start["g0"])
        start.link["sigma"] <- log(start["sigma"])
    } else {
        names(start.link)[c(2, 3)] <- c("lambda0", "sigma")
        start.link["lambda0"] <- log(start["lambda0"])
        start.link["sigma"] <- log(start["sigma"])
    }
    if (use_toa){
        start.link["sigma_toa"] <- log(start["sigma_toa"])
    }
    start <- c(start[1:4], start["sigma_toa"][use_toa], start["sigma_ss"][use_ss])
    start.link <- c(start.link[1:4], start.link["sigma_toa"][use_toa], start.link["sigma_ss"][use_ss])
    n.pars <- length(start)
    par.names <- names(start)
    ## Fitting model.
    fit <- nlminb(start.link, scr.nll.cuerate.multi,
                  caps = capthist,
                  aMask = aMask,
                  maskDists = maskDists,
                  ID = ids,
                  toa = toa,
                  toa_ssq = toa_ssq,
                  use_toa = use_toa,
                  ss = ss,
                  use_ss = use_ss,
                  hn = hn,
                  trace = trace)
    ## Approximating Hessian.
    hess <- optimHess(fit$par, scr.nll.cuerate.multi,
                      caps = capthist,
                      aMask = aMask,
                      maskDists = maskDists,
                      ID = ids,
                      toa = toa,
                      toa_ssq = toa_ssq,
                      use_toa = use_toa,
                      ss = ss,
                      use_ss = use_ss,
                      hn = hn,
                      trace = trace)
    
    ## Calculating confidence intervals
    ## - Using the (sqrt of) diagonals of (-ve) Hessian obtained from optimHess
    ##    - i.e. Information matrix
    ## - Wald CIs calculated by sapply() loop
    ##    - Loops through each of fitted parameters and calculates lower/upper bounds
    ## Note: fitted pars must be on LINK scale
    ##     : if matrix is singular, none of the SEs or CIs are calculated (inherits/try statement)
    fittedPars = fit$par
    names(fittedPars) <- par.names
    if(inherits(try(solve(hess), silent = TRUE), "try-error")) {
        ## Hessian is singular
        warning("Warning: singular hessian")
        ## SE and Wald CIs not calculated
        se = NA
        waldCI = matrix(NA, nrow = length(fittedPars), ncol = 2)
        ## But columns still need to be returned (if/when simulations are run)
        cnames = c("Estimate", "SE", "Lower", "Upper")
    } else {
        ## Calculating basic Wald CIs
        se = sqrt(diag(solve(hess)))
        waldCI.link = t(sapply(1:length(fittedPars),
                               function(i) fittedPars[i] + (c(-1, 1) * (qnorm(0.975) * se[i]))))
        G.mult <- numeric(n.pars)
        waldCI <- 0*waldCI.link
        ## Back-transforming the confidence limits.
        for (i in par.names){
            if (i %in% c("D", "b0_ss", "b1_ss", "sigma_ss", "lambda0", "sigma", "lambda", "sigma_toa")){
                waldCI[par.names == i, ] <- exp(waldCI.link[par.names == i, ])
                G.mult[par.names == i] <- exp(fittedPars[i])
                fittedPars[i] <- exp(fittedPars[i])
            }
            if (i == "g0"){
                waldCI[par.names == i, ] <- plogis(waldCI.link[par.names == i, ])
                G.mult[par.names == i] <- dlogis(fittedPars[i])
                fittedPars[i] <- plogis(fittedPars[i])
            }
        }
        ## Using the delta method to get the standard errors
        ## - G = jacobian matrix of partial derivatives of back-transformed
        ##    - i.e. log(D) -> exp(D) -- deriv. --> exp(D)
        ##    - Note: 1st deriv of plogis (CDF) = dlogis (PDF)
        G = diag(length(fittedPars)) * G.mult
        se = sqrt(diag(G %*% solve(hess) %*% t(G)))
    }
    results = cbind(fittedPars, se, waldCI)
    dimnames(results) = list(par.names, c("Estimate", "SE", "Lower", "Upper"))
    ## Returning a list with everything.
    list(results = results, capthist = capthist, 
         mask = mask, aMask = aMask, maskDists = maskDists, 
         speed_sound = speed_sound, ids = ids,
         traps = traps, detfn = detfn, toa = toa, toa_ssq = toa_ssq, hess = hess)
   
}

scr.nll.cuerate.multi <- function(pars, caps, aMask, maskDists, ID, toa, toa_ssq,
                                  use_toa, ss, use_ss, hn, trace, is.multi){
    n.sessions <- length(caps)
    sess.nll <- numeric(n.sessions)
    for (i in 1:n.sessions){
        sess.nll[i] <- scr_nll_cuerate(pars, caps[[i]], aMask[[i]], maskDists[[i]],
                                       ID[[i]], toa[[i]], toa_ssq[[i]], use_toa, ss[[i]],
                                       use_ss, hn, trace)
    }
    sum(sess.nll)
}
