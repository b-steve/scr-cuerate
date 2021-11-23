library(ascr)
library(Rcpp)
## Loading in data.
load("lightfooti.RData")
## This .RData file contains the following objects:

## capt.both: A list with a component for each survey, comprising
##            capture histories and times-of-arrival for detected
##            calls.
## capt:      Capture histories for calls from both surveys combined
##            into a single matrix.
## ids.both:  A list with a component for each survey, comprising
##            individual identifiers for each call detected.
## ids:       Individual identifiers from both surveys combined into a
##            single vector.
## mask:      A suitable SCR mask object for analysis.
## toa:       Times-of-arrival for detected calls from both surveys
##            combined into a single matrix.

## Sourcing functions.
source("fit-cuerate-scr.r")
sourceCpp("fit-cuerate-scr.cpp")
start <- c(D = 600, lambda0 = 9, sigma = 2, lambda = 9, sigma_toa = 10)

fit.cuerate <- cuerate.scr.fit(capt, ids, traps, mask, detfn = "hhn",
                               start = start, toa = toa, trace = FALSE)
## Parameter estimates.
fit.cuerate$results
# [1] 716.981120   7.495377   2.213246   9.070865   1.039855 

## Note that lambda is the number of calls per individual per 30
## s. The following provides an estimate of the number of calls per
## individual per minute, as presented in the paper.
2*fit.cuerate$results[4, ]

## Note that density is the number of calling animals per hectare, but
## accumulated across two different surveys. This is halved to give
## the density reported in the paper.
fit.cuerate$results[1, ]/2

## Calculating the estimated ESA.
## Distances from each detector to each mask point.
dists <- distances(mask, traps)
## Call detection probabilities for each mask point at each detector.
call.probs <- 1 - exp(-fit.cuerate$results[2, 1]*exp(-dists^2/(2*fit.cuerate$results[3, 1]^2)))
## Individual detection probabilities.
ind.probs <- 1 - exp(-call.probs*fit.cuerate$results[4, 1])
## Probability of detecting at least one call from  an individual at each mask point.
ind.det <- 1 - apply(1 - ind.probs, 1, prod)
esa <- sum(ind.det*attr(mask, "area"))
## Note that this matches n/D, which holds for other constant-density SCR models.
length(unique(ids))/fit.cuerate$results[1, 1]

## Model fitting using ascr, assuming call locations are independent.
fit.call <- fit.ascr(capt = list(bincapt = capt, toa = toa),
                     traps = traps, mask = mask, detfn = "hhn", trace = FALSE)
## Note that D here is call density (calls per hectare per minute).
summary(fit.call)

#save(fit.cuerate, fit.call, capt, ids, traps, mask, toa, capt.both, ids.both, esa, file = "fit-cuerate.RData")




## Trying out multi-session model fits.
capt.multi <- ids.multi <- traps.multi <- mask.multi <- toa.multi <- vector("list", length = 2)
capt.multi[[1]] <- capt.both[[1]]$bincapt
capt.multi[[2]] <- capt.both[[2]]$bincapt
ids.multi[[1]] <- ids[ids < 100]
ids.multi[[2]] <- ids[ids >= 100]
traps.multi[[1]] <- traps
traps.multi[[2]] <- traps
mask.multi[[1]] <- mask
mask.multi[[2]] <- mask
toa.multi[[1]] <- capt.both[[1]]$toa
toa.multi[[2]] <- capt.both[[2]]$toa
## Multi-session model.
fit.cuerate.multi <- cuerate.scr.fit(capt.multi, ids.multi, traps.multi, mask.multi, detfn = "hhn",
                                     start = start, toa = toa.multi, trace = FALSE)
fit.cuerate.multi$results
