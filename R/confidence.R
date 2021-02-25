#' Compute confidence around rate estimates
#' 
#' corHMM gives a single point estimate of rates, but this point estimate could be very uncertain. A traditional way to evaluate confidence is to vary each parameter until the log likelihood gets 2 units worse, while holding other parameters at their maximum likelihood estimates. That's fine and fast, but it can be misled by ridges. So, instead, we want all values that lead to a likelihood within two log likelihood units of the best. The range will be at least as wide as the univariate estimates but probably much larger. 
#' 
#' @param corhmm.object The result of a corHMM search
#' @param desired.delta How many log likelihood units to deviate from the optimal likelihood from
#' @param n.points How many points to use
#' @param verbose If TRUE, print details of the search to the screen
#' @param good.only If TRUE, only return the ones within the desired delta
#' @param ... Other arguments to pass into the likelihood function.
#' @export 
#' @examples 
#' 
#' @author Brian C. O'Meara
ComputeConfidenceIntervals <- function(corhmm.object, desired.delta = 2, n.points=5000, verbose=TRUE, good.only=TRUE, ...) {
	best.lnl <- corhmm.object$loglik
	index.mat <- corhmm.object$index.mat
	raw.rates <- corhmm.object$solution
	corhmm.object$node.states <- "none" #we don't want to waste time on this

	par <- rep(NA,max(index.mat, na.rm=TRUE))
	for (i in seq_along(par)) {
		par[i] <- raw.rates[which(index.mat==i)]
	}
	par.best <- par


	


	# Univariate
	results<-data.frame(data.frame(matrix(nrow=1, ncol=1+length(par))))
	results[1,] <- c(best.lnl, par.best)
	for(par_index in seq_along(par)) {
		par <- par.best
		current.lnl <- best.lnl
		while(current.lnl > best.lnl - 5*desired.delta) {
			par[par_index] <- par[par_index]*.99
			current.lnl <- compute_lnlikelihood(par, corhmm.object)
			results[nrow(results)+1,] <- c(current.lnl, par)
			if(verbose & nrow(results)%%100==0) {
				print(tail(results,100))
			}
		}
		current.lnl <- best.lnl
		while(current.lnl > best.lnl - 5*desired.delta) {
			par[par_index] <- par[par_index]*1.01
			current.lnl <- compute_lnlikelihood(par, corhmm.object)
			results[nrow(results)+1,] <- c(current.lnl, par)
			if(verbose & nrow(results)%%100==0) {
				print(tail(results,100))
			}
		}
	}

	# Multivariate
	lower <- apply(results[,-1], 2, min, na.rm=TRUE)
	upper <- apply(results[,-1], 2, max, na.rm=TRUE)
	new.params <- GetGridPoints(lower, upper, n.points)
	new.lnl <- pbapply(new.params, 1, compute_lnlikelihood, corhmm.object)

	more_results <- cbind(matrix(new.lnl, ncol=1), new.params)
	colnames(more_results) <- colnames(results)


	results <- rbind(results, more_results)
	results.good.enough <- subset(results, results[,1]>max(results[,1])-desired.delta)
	if(verbose) {
		print(apply(results.good.enough[,-1], 2, range))
	}

	if(good.only) {
		results <- results.good.enough
	}
    return(results)
}


HessianConfidence <- function(corhmm.object) {
	index.mat <- corhmm.object$index.mat
	raw.rates <- corhmm.object$solution
	corhmm.object$node.states <- "none" #we don't want to waste time on this

	par <- rep(NA,max(index.mat, na.rm=TRUE))
	for (i in seq_along(par)) {
		par[i] <- raw.rates[which(index.mat==i)]
	}
	hess <- numDeriv::hessian(compute_lnlikelihood, par, "Richardson", method.args=list(), corhmm.object)
	fisher_info<-solve(-hess)
	prop_sigma<-sqrt(diag(fisher_info))
	return.matrix <- rbind(par-1.96*prop_sigma, par+1.96*prop_sigma)
	rownames(return.matrix) <- c("lower", "upper")
	return(return.matrix)
}


# simple function; returns log likelihood
compute_lnlikelihood <- function(par, corhmm.object) {

	corhmm.object$order.test <- FALSE
	corhmm.object$phy$node.label <- NULL
	nObs <- length(corHMM:::corProcessData(corhmm.object$data)$ObservedTraits)
	model.set.final <- corHMM:::rate.cat.set.corHMM.JDB(phy = corhmm.object$phy, data = corhmm.object$data, rate.cat = corhmm.object$rate.cat, ntraits = nObs, model = "ARD")
	rate.mat <- MK_3state$index.mat
	rate.mat[rate.mat == 0] <- NA
	rate <- rate.mat
	model.set.final$np <- max(rate, na.rm=TRUE)
	rate[is.na(rate)]=max(rate, na.rm=TRUE)+1
	model.set.final$rate <- rate
	model.set.final$index.matrix <- rate.mat
	model.set.final$Q <- matrix(0, dim(rate.mat)[1], dim(rate.mat)[2])
	## for precursor type models ##
	col.sums <- which(colSums(rate.mat, na.rm=TRUE) == 0)
	row.sums <- which(rowSums(rate.mat, na.rm=TRUE) == 0)
	drop.states <- col.sums[which(col.sums == row.sums)]
	if(length(drop.states > 0)){
	model.set.final$liks[,drop.states] <- 0
	}
	result <- corHMM:::dev.corhmm(
		p = log(par), 
		phy = corhmm.object$phy, 
		liks = model.set.final$liks, 
		Q = model.set.final$Q, 
		rate = model.set.final$rate, 
		root.p = corhmm.object$root.p, 
		rate.cat = corhmm.object$rate.cat, 
		order.test = corhmm.object$order.test, 
		lewis.asc.bias = ifelse(any(grepl("lewis.asc.bias", names(corhmm.object))), corhmm.object$lewis.asc.bias, FALSE)
	)

	#return(dev.corhmm(log(par), corhmm_object$phy, liks=42, Q=42, rate=42, root.p=corhmm.object$root.p, rate.cat=42, order.test=42, lewis.asc.bias=ifelse(any(grepl("lewis.asc.bias", names(corhmm.object))), corhmm.object$lewis.asc.bias, FALSE)))
	return(-result)
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

GetGridPoints <- function(lower, upper, n.points) {
	n.grains <- ceiling(n.points^(1/length(lower)))
	vector.list <- list()
	for (i in seq_along(lower)) {
		vector.list[[i]] <- seq(from=lower[i], to=upper[i], length.out=n.grains)
	}
	parameter.matrix <- expand.grid(vector.list)
	return(parameter.matrix)
}