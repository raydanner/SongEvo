#' Test hypotheses about song evolution
#' 
#' Test if cultural traits evolve through specific mechanisms (e.g. drift or selection).
#' 
#' @name h.test
#' @param summary.results The summary.results array (i.e. a multi-dimensional table) from SongEvo(), which includes population summary values for each time step (dimension 1) in each iteration (dimension 2) of the model.  Population summary values are contained in five additional dimensions: population size for each time step of each iteration (“sample.n”), the population mean and variance of the song feature studied (“trait.pop.mean” and “trait.pop.variance”), with associated lower (“lci”) and upper (“uci”) confidence intervals.  
#' @param ts The timestep (“ts”) at which to compare simulated trait values to empirical trait values (“empir.trait”).
#' @param empir.trait Trait values from the test population to compare to simulated results. May be measured (i.e. empirical) or hypothetical. 
#' 
#' @return a list with two measures of accuracy: 1. The proportion of observed points that fall within the confidence intervals of the simulated data and the residuals between simulated and observed population trait means; 2. Precision is measured as the residuals between simulated and observed population trait variances.
#' 
#' @example inst/examples/h.testExamples.R
#' 
#' @references
#' @export
h.test <- function(summary.results, ts, target.data) {
 
  #Calculate residuals
h.test.res <- array(NA, dim=c(iteration, 2), dimnames=list(paste("Iteration", seq(1:iteration), sep=" "), c("Residuals of mean", "Residuals of variance")))
		h.test.res[,1] <- abs(mean(target.data)-summary.results[, ts, "trait.pop.mean"])
		h.test.res[,2] <- abs(var(target.data)-summary.results[, ts, "trait.pop.variance"])
		
  #Calculate proportion contained
quants <- quantile(summary.results[, ts, "trait.pop.mean"], probs=c(0.95, 0.05), R=600, na.rm=TRUE)
h.test.cont <- sum(target.data > quants[2] & target.data < quants[1], na.rm=TRUE)/length(target.data)

list("Residuals"=h.test.res, "Prop.contained"=h.test.cont)
}
