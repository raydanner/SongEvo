#' Parameter optimization
#' 
#' This function follows par.sens to help users optimize values for imperfectly known parameters for SongEvo. The goals are to maximize accuracy and precision of model prediction.
#'
#' @name par.opt
#' @param sens.results The sens.results array from par.sens(), which includes summary.results from SongEvo() for a range of parameter values.  summary.results from SongEvo() includes population summary values for each time step (dimension 1) in each iteration (dimension 2) of the model.  Population summary values are contained in five additional dimensions: population size for each time step of each iteration (“sample.n”), the population mean and variance of the song feature studied (“trait.pop.mean” and “trait.pop.variance”), with associated lower (“lci”) and upper (“uci”) confidence intervals.  
#' @param ts The timestep (“ts”) at which to compare simulated trait values to target trait values (“target.data”).
#' @param target.data A set of trait values against which to compare simulated values. target.data may be measured (e.g. from a training population) or hypothetical.
#' @param par.range Range of parameter values over which to optimize.
#'
#' @return Three measurements of accuracy and one measure of precision.  Accuracy is quantified by three different approaches: i) the mean of absolute residuals of the predicted population mean values in relation to observed values (smaller absolute residuals indicate a more accurate model), ii) the difference between the bootstrapped mean of predicted population means and the mean of the observed values, and iii) the proportion of simulated population trait means that fall within confidence intervals of the observed data (a higher proportion indicates greater accuracy). Precision is measured with the residuals of the predicted population variance to the variance of observed values (smaller residuals indicate a more precise model).
#'
#' @example inst/examples/par.optExamples.R
#' @seealso [SongEvo::SongEvo()], [SongEvo::par.sens()], [SongEvo::mod.val()], [SongEvo::h.test()], 'browseVignettes("SongEvo")'
#' @export
#' @importFrom stats var
#' @importFrom boot boot boot.ci
par.opt <- function(sens.results, ts, target.data, par.range) {
  
  iteration <- dim(sens.results)[1]
  #Calculate residuals
res <- array(NA, dim=c(length(par.range), iteration, 2), dimnames=list(paste("par.val", par.range), paste("Iteration", seq(1:iteration), sep=" "), c("Residuals of mean", "Residuals of variance")))
	for(p in 1:length(par.range)){
		res[p,,1] <- abs(mean(target.data)-sens.results[, ts, "trait.pop.mean", p])
		res[p,,2] <- abs(var(target.data)-sens.results[, ts, "trait.pop.variance", p])
		}

  #Calculate i) distance between mean of predicted means and the observed population mean and ii) proportion of predicted data points within confidence limits of observed data.

boot_curr <- boot(target.data, statistic=sample.mean, R=100)
ci.curr <- boot.ci(boot_curr, conf=0.95, type="basic")
low <- ci.curr$basic[4]
high <- ci.curr$basic[5]
cont <- array(NA, dim=c(length(par.range), 2), dimnames=list(paste("par.val", par.range), c("Difference in means", "Proportion contained")))
	for(p in 1:length(par.range)){
		cont[,1] <- abs(mean(target.data)-apply(sens.results[, ts, "trait.pop.mean", ], 2, mean, na.rm=TRUE))
		cont[,2] <- apply(sens.results[, ts, "trait.pop.mean", ], 2, function(x)sum(x > low & x < high, na.rm=TRUE))/iteration
		}
list("Residuals"=res, "Target.match"=cont)
}
