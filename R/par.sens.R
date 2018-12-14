#' Parameter sensitivity
#' 
#' This function allows testing the sensitivity of SongEvo to different parameter values.
#'
#' @name par.sens
#' @param parm The parameter for which to test sensitivity over one or more values. 
#' @param par.range List of ranges of parameter values over which to test sensitivity.
#' @param iteration The number of iterations that the model will run. 
#' @param steps The number of steps (e.g. years) per iteration.
#' @param mate.comp Female preference for mates. Currently specified as “Yes” or “No”. 
#' @param fixed_parms Named boolean vector identifying which parameters to keep fixed.
#' @param all Save data for all individuals? Options are TRUE or FALSE. 
#'
#'The function currently allows examination of only one parameter at a time and requires at least two iterations.
#'
#' @return array named sens.results. The sens.results array from par.sens(), which includes summary.results from SongEvo() for a range of parameter values.  summary.results from SongEvo() includes population summary values for each time step (dimension 1) in each iteration (dimension 2) of the model.  Population summary values are contained in five additional dimensions: population size for each time step of each iteration (“sample.n”), the population mean and variance of the song feature studied (“trait.pop.mean” and “trait.pop.variance”), with associated lower (“lci”) and upper (“uci”) confidence intervals.  
#'
#' @example inst/examples/par.sensExamples.R
#' @seealso [SongEvo::SongEvo()], [SongEvo::par.opt()], [SongEvo::mod.val()], [SongEvo::h.test()], 'browseVignettes("SongEvo")'
#' @export
#' @importFrom stats quantile

par.sens <- function(parm, par.range, iteration, steps, mate.comp, fixed_parms, all) {
	par.rangel <- length(par.range)
	sens.results <- array(NA, dim=c(iteration, steps, 5, par.rangel), dimnames = list(paste("iteration", seq(1:iteration)), 1:steps, c("sample.n", "trait.pop.mean", "trait.pop.variance", "lci", "uci"), paste("par.val", par.range)))
	for (p in 1:par.rangel) {
		if (parm=="terr.turnover") {
		terr.turnover=par.range[p]			
			}
		if (parm=="learning.error.sd") {
		learning.error.sd=par.range[p]			
			}
		if (parm=="mortality.a") {
		mortality.a=par.range[p]			
			}
		if (parm=="mortality.j") {
		mortality.j=par.range[p]			
			}
		if (parm=="lifespan") {
		lifespan=par.range[p]			
			}
		if (parm=="phys.lim.min") {
		phys.lim.min=par.range[p]			
			}
		if (parm=="phys.lim.max") {
		phys.lim.max=par.range[p]			
			}
		if (parm=="male.fledge.n") {
		male.fledge.n=par.range[p]			
			}
		if (parm=="disp.age") {
		disp.age=par.range[p]			
			}
		if (parm=="disp.distance.mean") {
		disp.distance.mean=par.range[p]			
			}
		if (parm=="disp.distance.sd") {
		disp.distance.sd=par.range[p]			
			}
		if (parm=="n.territories") {
		n.territories=par.range[p]			
			}
		print(paste(parm, "= ", par.range[p]))			
		z <- with(fixed_parms[!names(fixed_parms) %in% c(parm,"iteration","steps","mate.comp","all")],
		          SongEvo(init.inds, iteration, steps, timestep, terr.turnover, mate.comp, 
		                  learning.method, integrate.dist, learning.error.d, learning.error.sd, 
		                  mortality.a, mortality.j, lifespan, phys.lim.min, phys.lim.max, 
		                  male.fledge.n.mean, male.fledge.n.sd, male.fledge.n, 
		                  disp.age, disp.distance.mean, disp.distance.sd, n.territories, prin, all))
		#Add summary.results to sens.results array. 
		sens.results[ , , , p] <- z$summary.results
		}
if(iteration>=3){
	#Calculate quantiles of trait values across iterations within a parameter value.  Must have at least 3 iterations per parameter value.
	sens.results.diff <- array(NA, dim=c(length(par.range), steps), dimnames=list(paste("par.val", par.range), paste("Quantile diff", seq(1:steps))))
	for(p in 1:length(par.range)){
		quant.means <- apply(sens.results[, , "trait.pop.mean", p], MARGIN=2, quantile, probs=c(0.95, 0.05), R=600, na.rm=TRUE)
		sens.results.diff[p, ] <- (quant.means[1,]-quant.means[2,])
		}
	} else if(iteration<3){
	  sens.results.diff <- NA
	  }
list("sens.results"=sens.results, "sens.results.diff"=sens.results.diff)	
}
