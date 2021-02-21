#' Compute confidence around rate estimates
#' 
#' corHMM gives a single point estimate of rates, but this point estimate could be very uncertain. A traditional way to evaluate confidence is to vary each parameter until the log likelihood gets 2 units worse, while holding other parameters at their maximum likelihood estimates. That's fine and fast, but it can be misled by ridges. So, instead, we want all values that lead to a likelihood within two log likelihood units of the best. The range will be at least as wide as the univariate estimates but probably much larger. 
#' 
#' @param corhmm.object The result of a corHMM search
#' @param desired.delta How many log likelihood units to deviate from the optimal likelihood from
#' @param n.points How many points to use
#' @param verbose If TRUE, print details of the search to the screen
#' @param likelihood.precision How precisely to calculate likelihoood (less precise=faster)
#' @param ... Other arguments to pass into the likelihood function.
#' @export 
#' @examples 
#' 
#' @author Brian C. O'Meara
ComputeConfidenceIntervals <- function(corhmm.object, desired.delta = 2, n.points=5000, verbose=TRUE, likelihood.precision=0.01, ...) {
	best.lnl <- corhmm.object$loglik
	index.mat <- corhmm.object$index.mat
	raw.rates <- corhmm.object$solution
	corhmm.object$node.states <- "none" #we don't want to waste time on this

	par <- rep(NA,max(index.mat, na.rm=TRUE))
	for (i in seq_along(par)) {
		par[i] <- raw.rates[which(index.mat==i)]
	}
	par.best <- par

	compute_likelihood <- function(par, corhmm_object) {
		return(dev.corhmm(log(par), corhmm_object$phy, liks=42, Q=42, rate=42, root.p=corhmm.object$root.p, rate.cat=42, order.test=42, lewis.asc.bias=ifelse(any(grepl("lewis.asc.bias", names(corhmm.object))), corhmm.object$lewis.asc.bias, FALSE)))
	}

	# Univariate
	results<-data.frame(data.frame(matrix(nrow=1, ncol=1+length(par))))
	results[1,] <- c(best.lnl, par.best)
	for(par_index in seq_along(par)) {
		par <- par.best
		current.lnl <- best.lnl
		while(current.lnl > best.lnl - desired.delta) {
			par[par_index] <- par[par_index]*.95
			current.lnl <- compute_likelihood(par, corhmm_object)
			results[nrow(results)+1,] <- c(current.lnl, par)
		}
		current.lnl <- best.lnl
		while(current.lnl > best.lnl - desired.delta) {
			par[par_index] <- par[par_index]*1.05
			current.lnl <- compute_likelihood(par, corhmm_object)
			results[nrow(results)+1,] <- c(current.lnl, par)
		}
	}

	# Multivariate
	par <- par.best
	lower <- apply(results[,-1], 2, min, na.rm=TRUE)
	upper <- apply(results[,-1], 2, max, na.rm=TRUE)
	starting.row <- nrow(results)
	results <- rbind(results, matrix(nrow=n.points, ncol=1+length(par)))
    min.multipliers <- rep(1, length(par))
    max.multipliers <- rep(1, length(par))
    for (i in sequence(n.points)) {
        sim.points<-NA
        while(is.na(sim.points[1]) | !is.numeric(sim.points[1])) {
            sim.points<-GenerateValues(par, lower, upper, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
        }
        results[i+starting.row,] <- c(compute_likelihood(sim.points, corhmm.object), sim.points)
        # if(i>5 & restart.mode) {
        #     if((best.lnl - min(results[,1], na.rm=TRUE) > likelihood.precision ) & allow.restart) {
        #         results <- results[sequence(i+1),] #stop here and restart
        #         return(results)
        #     }
        # }
        if (i%%20==0) {
            for (j in sequence(length(par))) {
                returned.range <- range(results[which((results[,1]-min(results[,1], na.rm=TRUE))<desired.delta), j+1], na.rm=TRUE)
                total.range <- range(results[,j+1], na.rm=TRUE)
                width.ratio <- diff(returned.range)/diff(total.range)
                if(is.na(width.ratio)) {
                    width.ratio=1
                }
                if(width.ratio > 0.5) { #we are not sampling widely enough
                    min.multipliers[j] <- min.multipliers[j] * 0.9
                    max.multipliers[j] <- max.multipliers[j] * 1.1 #expand the range
                } else {
                    min.multipliers[j] <- 1
                    max.multipliers[j] <- 1
                }
            }
        }
        if (verbose && i%%100==0) {
            print(paste(i, "of", n.points, "done"))
        }
    }
    return(results)
}

#' Gets random numbers to try
GenerateValues <- function(par, lower, upper, max.tries=100, expand.prob=0, examined.max, examined.min) {
    if(is.null(lower)) {
        lower <- 0.1*par
    }
    if(is.null(upper)) {
        upper <- 10*par
    }
    pass=FALSE
    tries=0
    while(!pass && tries<=max.tries) {
        tries <- tries+1
        pass=TRUE
        new.vals <- rep(NA, length(par))
        for(i in sequence(length(par))) {
            examined.max[i]<-max(0.001, examined.max[i])
            new.vals.bounds <- sort(c(max(lower[i], 0.9*examined.min[i]), min(upper[i], 1.1*examined.max[i])), decreasing=FALSE)
            new.vals[i]<-stats::runif(1, min=ifelse(is.finite(new.vals.bounds[1]),new.vals.bounds[1], 0.000001) , max=ifelse(is.finite(new.vals.bounds[2]), new.vals.bounds[2], 10000))

            if(new.vals[i]<lower[i]) {
                pass=FALSE
            }
            if(new.vals[i]>upper[i]) {
                pass=FALSE
            }
        }
    }
    if(tries>max.tries) {
        return(NA)
    }
    names(new.vals) <- names(par)
    return(new.vals)
}