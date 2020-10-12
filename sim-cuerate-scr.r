library(Rcpp)
library(ascr)
library(spatstat)
## Sourcing functions.
source("fit-cuerate-scr.r")
sourceCpp("fit-cuerate-scr.cpp")

## Loading in fitted models.
load("fit-cuerate.RData")

set.seed(1234)
## Speed of sound (metres per second).
speed.sound <- 330
## Number of detectors.
n.traps <- nrow(traps)
## The x- and y-limits of the mask area.
xlim <- range(mask[, 1])
ylim <- range(mask[, 2])
## The area of a bounding box that covers the mask.
a.bb <- diff(xlim)*diff(ylim)

## Setting parameter values.
D <- fit.id$results[1, 1]
lambda0 <- fit.id$results[2, 1]
sigma <- fit.id$results[3, 1]
lambda.c <- fit.id$results[4, 1]
sigma.toa <- fit.id$results[5, 1]/1000
pars <- c(D = D, lambda0 = lambda0, sigma = sigma, lambda.c = lambda.c, sigma.toa = sigma.toa)
## Start values for scr.fit.
start <- c(D = D, g0 = lambda0, sigma = sigma,
           lambda_c = lambda.c, sigma_toa = sigma.toa*1000)
## Expected number of animals within the mask bounding box.
e.n.animals <- D/10000*a.bb

## Number of simulation iterations.
n.sims <- 1000
## Initialising matrix to store point estimates.
ests1 <- matrix(0, nrow = n.sims, ncol = 4)
ests2 <- matrix(0, nrow = n.sims, ncol = 5)
## Initialising array to store confidence intervals.
cis1 <- array(0, dim = c(n.sims, 4, 2))
cis2 <- array(0, dim = c(n.sims, 5, 2))
## Initialising number of detected animals.
n.detected <- numeric(n.sims)
## Loop for simulation.
for (i in 1:n.sims){
    ## Simulated number of animals in the bounding box.
    n.animals <- rpois(1, e.n.animals)
    ## Simulated locations of animals.
    animal.locs <- cbind(runif(n.animals, xlim[1], xlim[2]),
                         runif(n.animals, ylim[1], ylim[2]))
    ## Distance matrix for each animal-microphone pairing.
    animal.dists <- crossdist(animal.locs[, 1], animal.locs[, 2],
                       traps[, 1], traps[, 2])
    ## Detection probabilities for each animal-microphone pairing.
    det.p <- 1 - exp(-lambda0*exp(-animal.dists^2/(2*sigma^2)))
    ## Number of calls made by each animal.
    n.calls <- rpois(n.animals, lambda.c)
    ## Initialising vector of indicators for detection.
    animal.detected <- logical(n.animals)
    ## Simulating the capture history matrix for each individual.
    capt.list <- vector(mode = "list", length = n.animals)
    for (j in 1:n.animals){
        capt.list[[j]] <- matrix(0, nrow = n.calls[j], ncol = n.traps)
        for (k in seq_along(numeric(n.calls[j]))){
            capt.list[[j]][k, ] <- rbinom(n.traps, 1, det.p[j, ])
        }
        animal.detected[j] <- sum(capt.list[[j]]) > 0
    }
    n.detected[i] <- sum(animal.detected)
    ## Removing undetected animals.
    capt.list <- capt.list[animal.detected]
    animal.locs <- animal.locs[animal.detected, ]
    ## For each animal detected, removing undetected calls.
    for (j in 1:n.detected[i]){
        capt.list[[j]] <- capt.list[[j]][apply(capt.list[[j]], 1, sum) > 0, , drop = FALSE]
    }
    ## Combining submatrices into one big matrix.
    capt <- matrix(0, nrow = 0, ncol = n.traps)
    for (j in 1:n.detected[i]){
        capt <- rbind(capt, capt.list[[j]])
    }
    ## Number of detected calls by each animal.
    n.calls.detected <- sapply(capt.list, nrow)
    ## Determining animal IDs.
    ids <- rep(1:n.detected[i], times = sapply(capt.list, nrow))
    ## Creating a matrix of call locations.
    call.locs <- matrix(0, nrow = sum(n.calls.detected), ncol = 2)
    call.locs[, 1] <- rep(animal.locs[, 1], times = n.calls.detected)
    call.locs[, 2] <- rep(animal.locs[, 2], times = n.calls.detected)
    ## Distance matrix for each detected call-microphone pairing.
    call.dists <- crossdist(call.locs[, 1], call.locs[, 2],
                            traps[, 1], traps[, 2])
    ## Simulating times of arrival. Note that source time doesn't
    ## matter, see Borchers et al (2015), so here we assume all
    ## animals called at time zero.
    call.toa <- call.dists/speed.sound +
        matrix(rnorm(length(call.dists), 0, sigma.toa),
               ncol = ncol(call.dists), nrow = nrow(call.dists))
    ## Setting call.toa to zero for undetected calls.
    call.toa[capt == 0] <- 0
    ## Fitting model in ascr.
    capt.ascr <- list(bincapt = capt, toa = call.toa)
    fit1 <- fit.ascr(capt = capt.ascr, traps = traps, mask = mask, detfn = "hhn", trace = FALSE)
    ## Fitting new model.
    fit2 <- cuerate.scr.fit(capt, ids, traps, mask, detfn = "hhn",
                            toa = call.toa, start = start, trace = FALSE)
    ## Saving point estimates.
    ests1[i, ] <- coef(fit1)
    ests2[i, ] <- fit2$results[, 1]
    ## Saving confidence intervals.
    cis1[i, , 1] <- confint(fit1)[, 1]
    cis1[i, , 2] <- confint(fit1)[, 2]
    cis2[i, , 1] <- fit2$results[, 3]
    cis2[i, , 2] <- fit2$results[, 4]
    ## Plotting results as they come in.
    ests1.Dc <- ests1[1:i, 1]
    ests2.Dc <- ests2[1:i, 1]*ests2[1:i, 4]
    par(mfrow = c(1, 2))
    boxplot(ests1.Dc, ests2.Dc, ylim = c(0, max(c(ests1.Dc, ests2.Dc))), names = c("CD", "AD"))
    abline(h = D*lambda.c, lty = "dashed")
    ests2.D <- ests2[1:i, 1]
    boxplot(ests2.D, names = "new", ylim = c(0, max(ests2.D)))
    abline(h = D, lty = "dashed")
    cat(i, " ")
}

save.image("sim-results.RData")

