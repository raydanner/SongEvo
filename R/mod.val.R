#' Model validation
#' 
#' This function allows users to assess the validity of the specified model by testing model performance with a population different from the population used to build the model. The user first runs SongEvo with initial trait values from the validation population. 
#'
#' @name mod.val
#' @param summary.results The summary.results array (i.e. a multi-dimensional table) from SongEvo(), which includes population summary values for each time step (dimension 1) in each iteration (dimension 2) of the model.  Population summary values are contained in five additional dimensions: population size for each time step of each iteration (“sample.n”), the population mean and variance of the song feature studied (“trait.pop.mean” and “trait.pop.variance”), with associated lower (“lci”) and upper (“uci”) confidence intervals.  
#' @param ts The timestep (“ts”) at which to compare simulated trait values to empirical trait values (“empir.trait”).
#' @param empir.trait Trait values from the validation population to compare to simulated results. May be measured (i.e. empirical) or hypothetical. 
#' 
#' @param empirical values from a specified timestep
#'
#' @return Three measurements of accuracy: i) the mean of absolute residuals of the predicted population mean values in relation to observed values (smaller absolute residuals indicate a more accurate model), ii) the difference between the bootstrapped mean of predicted population means and the mean of the observed values, and iii) the proportion of simulated population trait means that fall within confidence intervals of the observed data (a higher proportion indicates greater accuracy). Precision is measured with the residuals of the predicted population variance to the variance of observed values (smaller residuals indicate a more precise model). Users specify the timestep (“ts”) at which to compare simulated trait values to empirical trait values (“empir.trait”).
#'
#' @example SongEvo/inst/examples/mod.valExamples.R
#' @references
#' @export
mod.val <- function(summary.results, ts, target.data) {
  #Calculate residuals
val.res <- array(NA, dim=c(iteration, 2), dimnames=list(paste("Iteration", seq(1:iteration), sep=" "), c("Residuals of mean", "Residuals of variance")))
		val.res[,1] <- abs(mean(target.data)-summary.results[, ts, "trait.pop.mean"])
		val.res[,2] <- abs(var(target.data)-summary.results[, ts, "trait.pop.variance"])
		
  #Calculate i) distance between mean of predicted means and the observed population mean and ii) proportion of predicted data points within confidence limits of observed data.
library("boot")
sample.mean <- function(d, x) {
	mean(d[x])
}

boot_curr <- boot(target.data, statistic=sample.mean, R=100)
ci.curr <- boot.ci(boot_curr, conf=0.95, type="basic")
low <- ci.curr$basic[4]
high <- ci.curr$basic[5]
cont <- array(NA, dim=2, dimnames=list(c("Difference in means", "Proportion contained")))
		cont[1] <- abs(mean(target.data)-mean(summary.results[, ts, "trait.pop.mean"], na.rm=TRUE))
		cont[2] <- sum(summary.results[, ts, "trait.pop.mean"] > low & summary.results[, ts, "trait.pop.mean"] < high, na.rm=TRUE)/iteration
		
list("Residuals"=val.res, "Target.match"=cont)
}
