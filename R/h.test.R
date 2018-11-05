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
