### See vignette for an example that uses all functions in SongEvo.

#### Prepare current song values
target.data <- subset(song.data, Population=="PRBO" & Year==2005)$Trill.FBW

#### Specify and call `par.opt()`

# Users specify the timestep (“ts”) at which to compare simulated trait values
# to target trait data (“target.data”) and save the results in an object (called
# `par.opt1` here).
ts <- years
par.opt1 <- par.opt(sens.results=par.sens1$sens.results, ts=ts, target.data=target.data, par.range=par.range)

# Examine results objects (residuals and target match).  
par.opt1$Residuals
par.opt1$Target.match

#### Plot results of `par.opt()`
#### Accuracy

  #1. Difference in means.
plot(par.range, par.opt1$Target.match[,1], type="l", 
     xlab="Parameter range", ylab="Difference in means (Hz)")

  #2. Plot proportion contained.
plot(par.range, par.opt1$Prop.contained, type="l", 
     xlab="Parameter range", ylab="Proportion contained")

  #3. Calculate and plot mean and quantiles of residuals of mean trait values.
res.mean.means <- apply(par.opt1$Residuals[, , 1], MARGIN=1, 
                        mean, na.rm=TRUE)
res.mean.quants <- apply (par.opt1$Residuals[, , 1], MARGIN=1, 
                          quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
plot(par.range, res.mean.means, col="orange", 
     ylim=range(par.opt1$Residuals[,,1], na.rm=TRUE), 
     type="b", 
     xlab="Parameter value (territory turnover rate)", 
     ylab="Residual of trait mean (trill bandwidth, Hz)")
points(par.range, res.mean.quants[1,], col="orange")
points(par.range, res.mean.quants[2,], col="orange")
lines(par.range, res.mean.quants[1,], col="orange", lty=2)
lines(par.range, res.mean.quants[2,], col="orange", lty=2)

#### Precision
#Calculate and plot mean and quantiles of residuals of variance of trait values
res.var.mean <- apply(par.opt1$Residuals[, , 2], MARGIN=1, 
                      mean, na.rm=TRUE)
res.var.quants <- apply (par.opt1$Residuals[, , 2], MARGIN=1, 
                         quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
plot(par.range, res.var.mean, col="purple", 
     ylim=range(par.opt1$Residuals[,,2], na.rm=TRUE),
     type="b", 
     xlab="Parameter value (territory turnover rate)", 
     ylab="Residual of trait variance (trill bandwidth, Hz)")
points(par.range, res.var.quants[1,], col="purple")
points(par.range, res.var.quants[2,], col="purple")
lines(par.range, res.var.quants[1,], col="purple", lty=2)
lines(par.range, res.var.quants[2,], col="purple", lty=2)

#### Visual inspection of accuracy and precision: plot trait values for range of parameters
par(mfcol=c(3,2))
par(mar=c(4.1, 4.1, 1, 1))
par(cex=1.2)
for(i in 1:length(par.range)){
plot(par.sens1$sens.results[ , , "trait.pop.mean", ], 
     xlab="Year", ylab="Bandwidth (Hz)", 
     xaxt="n", type="n", 
     xlim=c(-0.5, years), ylim=range(par.sens1$sens.results[ , , "trait.pop.mean", ], na.rm=TRUE))
	for(p in 1:iteration){
		lines(par.sens1$sens.results[p, , "trait.pop.mean", i], col="light gray")
		}
freq.mean <- apply(par.sens1$sens.results[, , "trait.pop.mean", i], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))#, tcl=-0.25, mgp=c(2,0.5,0))
#Plot 95% quantiles
quant.means <- apply (par.sens1$sens.results[, , "trait.pop.mean", i], MARGIN=2, 
                      quantile, probs=c(0.95, 0.05), R=600, na.rm=TRUE)
lines(quant.means[1,], col="blue", lty=2)
lines(quant.means[2,], col="blue", lty=2)
#plot mean and CI for historic songs.  
 #plot original song values
library("boot")
sample.mean <- function(d, x) {
	mean(d[x])
}
boot_hist <- boot(starting.trait, statistic=sample.mean, R=100)#, strata=mn.res$iteration)	
ci.hist <- boot.ci(boot_hist, conf=0.95, type="basic")
low <- ci.hist$basic[4]
high <- ci.hist$basic[5]
points(0, mean(starting.trait), pch=20, cex=0.6, col="black")
library("Hmisc")
errbar(x=0, y=mean(starting.trait), high, low, add=TRUE)
 
 #plot current song values
boot_curr <- boot(target.data, statistic=sample.mean, R=100)#, strata=mn.res$iteration)	
ci.curr <- boot.ci(boot_curr, conf=0.95, type="basic")
low <- ci.curr$basic[4]
high <- ci.curr$basic[5]
points(years, mean(target.data), pch=20, cex=0.6, col="black")
errbar(x=years, y=mean(target.data), high, low, add=TRUE)
  #plot panel title
text(x=3, y=max(par.sens1$sens.results[ , , "trait.pop.mean", ], na.rm=TRUE)-100, 
     labels=paste("Par = ", par.range[i], sep=""))  
}
