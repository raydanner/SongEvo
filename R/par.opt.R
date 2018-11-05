par.opt <- function(sens.results, ts, target.data, par.range) {
  #Calculate residuals
res <- array(NA, dim=c(length(par.range), iteration, 2), dimnames=list(paste("par.val", par.range), paste("Iteration", seq(1:iteration), sep=" "), c("Residuals of mean", "Residuals of variance")))
	for(p in 1:length(par.range)){
		res[p,,1] <- abs(mean(target.data)-sens.results[, ts, "trait.pop.mean", p])
		res[p,,2] <- abs(var(target.data)-sens.results[, ts, "trait.pop.variance", p])
		}

  #Calculate i) distance between mean of predicted means and the observed population mean and ii) proportion of predicted data points within confidence limits of observed data.
library("boot")
sample.mean <- function(d, x) {
	mean(d[x])
}
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
