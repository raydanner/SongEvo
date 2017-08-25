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