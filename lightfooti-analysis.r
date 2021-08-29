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
start <- c(600, 9, 2, 9, 10)

mask <- create.mask(traps, buffer = 15)

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

save(fit.cuerate, fit.call, capt, ids, traps, mask, toa, capt.both, ids.both, esa, file = "fit-cuerate.RData")

