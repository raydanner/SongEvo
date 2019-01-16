##Vignette code from Readme

#Load functions from files, keeping source references for future debugging

library(SongEvo)
# Load the example model for WCSP 
# To explore the SongEvo package, we will use a database of songs from 
# Nuttall’s white-crowned sparrow (Zonotrichia leucophrys nuttalli) 
# recorded at three locations in 1969 and 2005.
# example data: song.data
data(song.data)

# Preselected global parameters: glo.parms
data(glo.parms)
# direct load into global envrionment
list2env(glo.parms, globalenv()) 

# Examine song data Data including: 
# the population name (Bear Valley, PRBO, or Schooner), 
# the year of song recording (1969 or 2005), and 
# the frequency bandwidth of the trill.
str(song.data)
# RStudio alternative: View(song.data)

# Use individuals from PRBO in 1969 
starting.trait.obs <- subset(song.data, Population %in% "PRBO" & Year %in% 1969)$Trill.FBW

# Generate additional individuals to reach number of territories for this vignette (40)
starting.trait.gen <- rnorm(n.territories-length(starting.trait.obs), 
                            mean=mean(starting.trait.obs), sd=sd(starting.trait.obs))
starting.trait <- c(starting.trait.obs, starting.trait.gen)

# Build data frame of territory id, age, and song trait value (in this case trill bandwidth)
init.inds <- data.frame(id = 1:n.territories, age = 2, trait = starting.trait,
                        row.names=paste0("A",formatC(1:n.territories,width=2,flag = "0")))

# Generate lat long coordinates for center? of each territory
init.inds$x1 <-  round(runif(n.territories, min=-122.481858, max=-122.447270), digits=8)
init.inds$y1 <-  round(runif(n.territories, min=37.787768, max=37.805645), digits=8)

# Quick aside, plot generated values
inds.plot <- init.inds
inds.plot$type= factor(c(rep("PRBO1969OBS",30),rep("PRBO1969GEN",10)))
library(lattice)
# Trait comparison, seems a bit oddly bunched
dotplot(type~trait,data=inds.plot)
densityplot(~trait,subset=type%in%"PRBO1969OBS",data=inds.plot)
densityplot(~trait,data=inds.plot)

# Scattering of territories:
xyplot(x1~y1,data=inds.plot,group=type)

# Running the SongEvo model:
#   Doing things a bit differently, working directly from the function call (instead of keeping global):
SongEvo1 <- SongEvo(init.inds = init.inds, 
                    iteration = 10, 
                    steps = 36,  # years / timestep
                    timestep = 1, 
                    n.territories = nrow(init.inds), # perhaps replace with subset argument?
                    terr.turnover = 0.5, 
                    learning.method = "integrate", 
                    integrate.dist = 0.1, 
                    learning.error.d = learning.error.d, # From globals
                    learning.error.sd = learning.error.sd, # From globals
                    mortality.a = mortality.a, # From globals
                    mortality.j = mortality.j, # From globals
                    lifespan = NA, # Not sure why
                    phys.lim.min = phys.lim.min, # From globals
                    phys.lim.max = phys.lim.max, # From globals
                    male.fledge.n.mean = male.fledge.n.mean, # From globals
                    male.fledge.n.sd = male.fledge.n.sd, # From globals
                    male.fledge.n = male.fledge.n, # From globals
                    disp.age = disp.age, # From globals
                    disp.distance.mean = disp.distance.mean, # From globals
                    disp.distance.sd = disp.distance.sd, # From globals
                    mate.comp = FALSE, 
                    prin = FALSE, 
                    all = TRUE)
str(SongEvo1)
#Save result for repeatable testing

save(SongEvo1, init.inds, 
     file=(most_recent_SongEvo1<-paste0("SongEvo1_std_",format(Sys.time(),format="%Y%m%d%H%M"),".RData")))


# Alternative is to load parameters into list and attach to result:
SongEvo1_alt=list(parms=glo.parms)
parms_1alt <- list(init.inds = init.inds, 
                   iteration = 10,
                   steps = 36,  # years / timestep
                   timestep = 1, 
                   n.territories = nrow(init.inds), # perhaps replace with subset argument?
                   terr.turnover = 0.5, 
                   learning.method = "integrate", 
                   integrate.dist = 0.1, 
                   lifespan = NA, # Not sure why
                   mate.comp = FALSE, 
                   prin = FALSE, 
                   all = TRUE)
SongEvo1_alt$parms[names(parms_1alt)]=parms_1alt
SongEvo1_alt$res=do.call(SongEvo,SongEvo1_alt$parms)
str(SongEvo1_alt)

#Save result for repeatable testing
save(SongEvo1_alt,
     file=(most_recent_SongEvo1_alt<-paste0("SongEvo1_alt_",format(Sys.time(),format="%Y%m%d%H%M"),".RData")))

# Compare calculation time from the current runs:
rbind(SongEvo1$time,SongEvo1_alt$res$time)

# Look at currently alive individuals:
head(SongEvo1$inds, 5)
head(SongEvo1_alt$res$inds, 5)

# Breakdown time series of population summaries?
results.iterations =list()
for(i_run in 1:10){
  ts(data=as.data.frame(SongEvo1$summary.results[i_run,,]),
     start = 1970, end=2005,frequency = 1 
     )->results.iterations[[i_run]]
}
# Look at time series of population trait value, with ci, for iteration 1
plot(results.iterations[[1]][,c(2,4,5)],plot.type="single")
results.vars=list()
for(i_var in c("sample.n","trait.pop.mean","trait.pop.variance","lci","uci")){
  
  ts(data=as.data.frame(t(SongEvo1$summary.results[,,i_var])),
     start = 1970, end=2005,frequency = 1 
  )->results.vars[[i_var]]
}

# Look at time series of population trait value, for all iterations
plot(results.vars$trait.pop.mean,plot.type="multiple")
#Over plot versions
plot(results.vars$trait.pop.mean,plot.type="single")
plot(results.vars$trait.pop.mean,plot.type="single",col=1:5,lty=rep(2:3,each=5))


## Similar to Figure S1
# Note that model only generates 10 iterations, unlike the paper's 100
# Oddity, note that population size is often > 40 despite the "n.territories" value
plot(SongEvo1$summary.results[1, , "sample.n"], xlab="Year", ylab="Abundance", type="n", xaxt="n", 
     ylim=c(0, max(SongEvo1$summary.results[, , "sample.n"], na.rm=TRUE)))
axis(side=1, at=seq(0, 40, by=5), labels=seq(1970, 2010, by=5))
for(p in 1:10){
  lines(SongEvo1$summary.results[p, , "sample.n"], col="light gray")
}
n.mean <- apply(SongEvo1$summary.results[, , "sample.n"], 2, mean, na.rm=TRUE)
lines(n.mean, col="red")

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "sample.n"], MARGIN=2, 
                      quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
lines(quant.means[1,], col="red", lty=2)
lines(quant.means[2,], col="red", lty=2)
library("Hmisc")

## Similar to Figure 2
# Note that model only generates 10 iterations, unlike the paper's 100
plot(SongEvo1$summary.results[1, , "trait.pop.mean"], xlab="Year", ylab="Bandwidth (Hz)", 
     xaxt="n", type="n", xlim=c(-0.5, 36), ylim=range(SongEvo1$summary.results[, , "trait.pop.mean"], na.rm=TRUE))
iteration=dim(SongEvo1$summary.results)[1]
for(p in 1:iteration){
  lines(SongEvo1$summary.results[p, , "trait.pop.mean"], col="light gray")
}
freq.mean <- apply(SongEvo1$summary.results[, , "trait.pop.mean"], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))#, tcl=-0.25, mgp=c(2,0.5,0))

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "trait.pop.mean"], MARGIN=2, 
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
errbar(x=0, y=mean(starting.trait), high, low, add=TRUE)
#text and arrows
text(x=5, y=2720, labels="Historical songs", pos=1)
arrows(x0=5, y0=2750, x1=0.4, y1=mean(starting.trait), length=0.1)

#plot variance for each iteration per year
## Similar to Figure S2
# Again, only 10 iterations
plot(SongEvo1$summary.results[1, , "trait.pop.variance"], 
     xlab="Year", ylab="Bandwidth Variance (Hz)", type="n", xaxt="n", 
     ylim=range(SongEvo1$summary.results[, , "trait.pop.variance"], na.rm=TRUE))
axis(side=1, at=seq(0, 40, by=5), labels=seq(1970, 2010, by=5))
for(p in 1:iteration){
  lines(SongEvo1$summary.results[p, , "trait.pop.variance"], col="light gray")
}
n.mean <- apply(SongEvo1$summary.results[, , "trait.pop.variance"], 2, mean, na.rm=TRUE)
lines(n.mean, col="green")

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "trait.pop.variance"], MARGIN=2, 
                      quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
lines(quant.means[1,], col="green", lty=2)
lines(quant.means[2,], col="green", lty=2)



## Similar to Figure 3
library("reshape2")
library("lattice")
library("sp")

all.inds1 <- subset(SongEvo1$all.inds, iteration==1)
w <- dcast(as.data.frame(all.inds1), id ~ timestep, value.var="trait", fill=0)
all.inds1w <- merge(all.inds1, w, by="id")
names(all.inds1w) <- c(names(all.inds1), paste("Ts", seq(1:36), sep=""))

rbPal <- colorRampPalette(c('blue','red')) #Create a function to generate a continuous color palette

# Plot maps, including a separate panel for each timestep (each of 36 years).
# Our example shows that individuals move across the landscape and that regional
# dialects evolve and move. The x-axis is longitude, the y-axis is latitude, and
# the color ramp indicates trill bandwidth in Hz.
spplot(all.inds1w[,-c(1:ncol(all.inds1))], as.table=TRUE, cuts=c(0, seq(from=1500, to=4500, by=10)), 
       ylab="", col.regions=c("transparent", rbPal(1000)), #cuts specifies that the first level (e.g. <1500) is transparent.
       colorkey=list(
         right=list(
           fun=draw.colorkey,
           args=list( 
             key=list(
               at=seq(1500, 4500, 10),
               col=rbPal(1000),
               labels=list(
                 at=c(1500, 2000, 2500, 3000, 3500, 4000, 4500),
                 labels=c("1500", "2000", "2500", "3000", "3500", "4000", "4500")
               )
             )
           )
         )
       )
)
# In addition, you can plot simpler multi-panel maps that do not take advantage of the spatial data class.

#Lattice plot (not as a spatial frame)
it1 <- subset(SongEvo1$all.inds, iteration==1)
rbPal <- colorRampPalette(c('blue','red')) #Create a function to generate a continuous color palette
it1$Col <- rbPal(10)[as.numeric(cut(it1$trait, breaks = 10))]
xyplot(it1$y1~it1$x1 | it1$timestep, groups=it1$trait, asp="iso", col=it1$Col, xlab="Longitude", ylab="Latitude")
#Alternate call:
it1$trait_group <- cut(it1$trait, breaks = 10,dig.lab=4)
xyplot(y1~x1 | timestep,data=it1@data,# groups=trait_group, 
       asp="iso", 
       col=it1$Col, xlab="Longitude", ylab="Latitude")

# With factoring:
xyplot(y1~x1 | equal.count(timestep),data=it1@data,# groups=trait_group, 
       asp="iso", col=it1$Col, xlab="Longitude", ylab="Latitude")
it1$trait_group <- cut(it1$trait, breaks = 4,dig.lab=4)
it1$ts_group<-cut(it1$timestep,breaks=4,dig.lab=2)
xyplot(y1~x1 | ts_group*trait_group,data=it1@data,# groups=trait_group, 
       asp="iso", col=it1$Col, xlab="Longitude", ylab="Latitude")


##Test model sensitivity with par.sens() This function allows testing the sensitivity of SongEvo to different parameter values.

###Specify and call par.sens() Here we test the sensitivity of the Acquire a
###Territory submodel to variation in territory turnover rates, ranging from
###0.8–1.2 times the published rate (40–60% of territories turned over). The
###call for the par.sens function has a format similar to SongEvo. The user
###specifies the parameter to test and the range of values for that parameter.
###The function currently allows examination of only one parameter at a time and
###requires at least two iterations.

parm <- "terr.turnover"
par.range = seq(from=0.4, to=0.6, by=0.05)
sens.results <- NULL
# Now we call the par.sens function with our specifications.
years=36
extra_parms <- list(init.inds = init.inds, 
                    timestep = 1, 
                    n.territories = nrow(init.inds), 
                    learning.method = "integrate", 
                    integrate.dist = 0.1, 
                    lifespan = NA, 
                    terr.turnover = 0.5, 
                    mate.comp = FALSE, 
                    prin = FALSE,
                    all = TRUE)
global_parms_key <- which(!names(glo.parms) %in% names(extra_parms))
extra_parms[names(glo.parms[global_parms_key])]=glo.parms[global_parms_key]
par.sens1 <- par.sens(parm = parm, par.range = par.range, 
                      iteration = iteration, steps = years, mate.comp = FALSE, 
                      fixed_parms=extra_parms[names(extra_parms)!=parm], all = TRUE)
###Examine par.sens results Examine results objects, which include two arrays:

# The first array, sens.results, contains the SongEvo model results for each parameter. It has the following dimensions:
  
dimnames(par.sens1$sens.results)
# The second array, sens.results.diff contains the quantile range of trait
# values across iterations within a parameter value. It has the following
# dimensions:
  
dimnames(par.sens1$sens.results.diff)
#plot of range in trait quantiles by year by parameter value
# Similar to Figure S3
plot(1:years, par.sens1$sens.results.diff[1,], 
     ylim=range(par.sens1$sens.results.diff, na.rm=TRUE), type="l", 
     ylab="Quantile range (Hz)", xlab="Year", col="transparent", xaxt="n")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))

#Make a continuous color ramp from gray to black
grbkPal <- colorRampPalette(c('gray','black'))

#Plot a line for each parameter value
for(i in 1:length(par.range)){
  lines(1:years, par.sens1$sens.results.diff[i,], type="l", col=grbkPal(length(par.range))[i])
}

#Plot values from published parameter values
lines(1:years, par.sens1$sens.results.diff[2,], type="l", col="black", lwd=4)

#Calculate and plot mean and quantiles
quant.mean <- apply(par.sens1$sens.results.diff, 2, mean, na.rm=TRUE)
lines(quant.mean, col="orange")

#Plot 95% quantiles (which are similar to credible intervals)
#95% quantiles of population means (narrower)
quant.means <- apply (par.sens1$sens.results.diff, MARGIN=2, quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
lines(quant.means[1,], col="orange", lty=2)
lines(quant.means[2,], col="orange", lty=2)

##Optimize parameter values with par.opt() This function follows par.sens to
##help users optimize values for imperfectly known parameters for SongEvo. The
##goals are to maximize accuracy and precision of model prediction. Accuracy is
##quantified by three different approaches: i) the mean of absolute residuals of
##the predicted population mean values in relation to target data (e.g. observed
##or hypothetical values (smaller absolute residuals indicate a more accurate
##model)), ii) the difference between the bootstrapped mean of predicted
##population means and the mean of the target data, and iii) the proportion of
##simulated population trait means that fall within (i.e. are "contained by")
##the confidence intervals of the target data (a higher proportion indicates
##greater accuracy). Precision is measured with the residuals of the predicted
##population variance to the variance of target data (smaller residuals indicate
##a more precise model).

###Prepare current song values

target.data <- subset(song.data, Population=="PRBO" & Year==2005)$Trill.FBW
###Specify and call par.opt() Users specify the timestep (“ts”) at which to
###compare simulated trait values to target trait data (“target.data”) and save
###the results in an object (called par.opt1 here).

ts <- years
par.opt1 <- par.opt(sens.results=par.sens1$sens.results, ts=ts, target.data=target.data, par.range=par.range)
# Examine results objects (residuals and target match).

par.opt1$Residuals
par.opt1$Target.match

###Plot results of par.opt() ####Accuracy

# Difference in means.
# Similar to figure S5
plot(par.range, par.opt1$Target.match[,1], type="l", xlab="Parameter range", ylab="Difference in means (Hz)")
# Plot proportion contained.
# Similar to figure ??
plot(par.range, par.opt1$Prop.contained, type="l", xlab="Parameter range", ylab="Proportion contained")
# Calculate and plot mean and quantiles of residuals of mean trait values.
# Similar to figure S6
res.mean.means <- apply(par.opt1$Residuals[, , 1], MARGIN=1, mean, na.rm=TRUE)
res.mean.quants <- apply (par.opt1$Residuals[, , 1], MARGIN=1, quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
plot(par.range, res.mean.means, col="orange", ylim=range(par.opt1$Residuals[,,1], na.rm=TRUE), type="b", 
     xlab="Parameter value (territory turnover rate)", ylab="Residual of trait mean (trill bandwidth, Hz)")
points(par.range, res.mean.quants[1,], col="orange")
points(par.range, res.mean.quants[2,], col="orange")
lines(par.range, res.mean.quants[1,], col="orange", lty=2)
lines(par.range, res.mean.quants[2,], col="orange", lty=2)
####Precision

#Calculate and plot mean and quantiles of residuals of variance of trait values
# Similar to figure S7
res.var.mean <- apply(par.opt1$Residuals[, , 2], MARGIN=1, mean, na.rm=TRUE)
res.var.quants <- apply (par.opt1$Residuals[, , 2], MARGIN=1, quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
plot(par.range, res.var.mean, col="purple", 
     ylim=range(par.opt1$Residuals[,,2], na.rm=TRUE), 
     type="b", 
     xlab="Parameter value (territory turnover rate)", ylab="Residual of trait variance (trill bandwidth, Hz)")
points(par.range, res.var.quants[1,], col="purple")
points(par.range, res.var.quants[2,], col="purple")
lines(par.range, res.var.quants[1,], col="purple", lty=2)
lines(par.range, res.var.quants[2,], col="purple", lty=2)
####Visual inspection of accuracy and precision: plot trait values for range of parameters
# Similar to Figure S4
par(mfcol=c(3,2))
par(mar=c(4.1, 4.1, 1, 1))
par(cex=1.2)
for(i in 1:length(par.range)){
  plot(par.sens1$sens.results[ , , "trait.pop.mean", ], xlab="Year", ylab="Bandwidth (Hz)", 
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
  library("boot")
  sample.mean <- function(d, x) {
    mean(d[x])
  }
  boot_curr <- boot(target.data, statistic=sample.mean, R=100)#, strata=mn.res$iteration)	
  ci.curr <- boot.ci(boot_curr, conf=0.95, type="basic")
  low <- ci.curr$basic[4]
  high <- ci.curr$basic[5]
  points(years, mean(target.data), pch=20, cex=0.6, col="black")
  library("Hmisc")
  errbar(x=years, y=mean(target.data), high, low, add=TRUE)
  
  #plot panel title
  text(x=3, y=max(par.sens1$sens.results[ , , "trait.pop.mean", ], na.rm=TRUE)-100, labels=paste("Par = ", par.range[i], sep=""))  
}

starting.trait <- subset(song.data, Population=="Schooner" & Year==1969)$Trill.FBW
starting.trait2 <- c(starting.trait, rnorm(n.territories-length(starting.trait), mean=mean(starting.trait), sd=sd(starting.trait)))

init.inds <- data.frame(id = seq(1:n.territories), age = 2, trait = starting.trait2)
init.inds$x1 <-  round(runif(n.territories, min=-122.481858, max=-122.447270), digits=8)
init.inds$y1 <-  round(runif(n.territories, min=37.787768, max=37.805645), digits=8)


##Model validation with mod.val()
# Specify and call SongEvo() with validation data
{
iteration <- 100
years <- 36
timestep <- 1
terr.turnover <- 0.5

# SongEvo2 <- SongEvo(init.inds = init.inds, iteration = iteration, steps = years, 
#                     timestep = timestep, n.territories = n.territories, terr.turnover = terr.turnover, 
#                     learning.method = learning.method, 
#                     integrate.dist = integrate.dist, learning.error.d = learning.error.d, 
#                     learning.error.sd = learning.error.sd, 
#                     mortality.a = mortality.a, mortality.j = mortality.j, 
#                     lifespan = lifespan, phys.lim.min = phys.lim.min, phys.lim.max = phys.lim.max, 
#                     male.fledge.n.mean = male.fledge.n.mean, male.fledge.n.sd = male.fledge.n.sd, 
#                     male.fledge.n = male.fledge.n, 
#                     disp.age = disp.age, disp.distance.mean = disp.distance.mean, 
#                     disp.distance.sd = disp.distance.sd, mate.comp = mate.comp, prin = prin, all)
SongEvo2_alt=list(parms=glo.parms)
parms_2alt <- list(init.inds = init.inds, 
                   iteration = 100,
                   steps = 36,  # years / timestep
                   timestep = 1, 
                   n.territories = 40,
                   terr.turnover = 0.5, 
                   learning.method = "integrate", 
                   integrate.dist = 0.1, 
                   lifespan = NA, 
                   mate.comp = FALSE, 
                   prin = FALSE, 
                   all = TRUE)
SongEvo2_alt$parms[names(parms_2alt)]=parms_2alt
SongEvo2_alt$res=do.call(SongEvo,SongEvo2_alt$parms)
str(SongEvo2_alt)

#Save result for repeatable testing
save(SongEvo2_alt,
     file=(most_recent_SongEvo2_alt<-paste0("SongEvo2_alt_",format(Sys.time(),format="%Y%m%d%H%M"),".RData")))

SongEvo2=SongEvo2_alt$res

# Specify and call mod.val

ts <- 36
target.data <- subset(song.data, Population=="Schooner" & Year==2005)$Trill.FBW
mod.val1 <- mod.val(summary.results=SongEvo2$summary.results, ts=ts, target.data=target.data)

# Plot results from mod.val()
# Similar to figure S8

#Clear out plotting parameters from par.opt plots
par(mfcol=c(1,1))
par(mar=c(5, 4, 4, 2) + 0.1)
par(cex=1)
plot(SongEvo2$summary.results[1, , "trait.pop.mean"], 
     xlab="Year", ylab="Bandwidth (Hz)", xaxt="n", type="n", xlim=c(-0.5, 36.5), ylim=range(SongEvo2$summary.results[, , "trait.pop.mean"], na.rm=TRUE))
for(p in 1:iteration){
  lines(SongEvo2$summary.results[p, , "trait.pop.mean"], col="light gray")
}
freq.mean <- apply(SongEvo2$summary.results[, , "trait.pop.mean"], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))#, tcl=-0.25, mgp=c(2,0.5,0))

#Plot 95% quantiles 
quant.means <- apply (SongEvo2$summary.results[, , "trait.pop.mean"], MARGIN=2, quantile, probs=c(0.95, 0.05), R=600, na.rm=TRUE)
lines(quant.means[1,], col="blue", lty=2)
lines(quant.means[2,], col="blue", lty=2)

#plot mean and CI for historic songs.  
#plot original song values
library("boot")
sample.mean <- function(d, x) {
  mean(d[x])
}
boot_hist <- boot(starting.trait, statistic=sample.mean, R=100)
ci.hist <- boot.ci(boot_hist, conf=0.95, type="basic")
low <- ci.hist$basic[4]
high <- ci.hist$basic[5]
points(0, mean(starting.trait), pch=20, cex=0.6, col="black")
library("Hmisc")
errbar(x=0, y=mean(starting.trait), high, low, add=TRUE)

#text and arrows
text(x=5, y=2720, labels="Historical songs", pos=1)
arrows(x0=5, y0=2750, x1=0.4, y1=mean(starting.trait), length=0.1)

#plot current song values
library("boot")
sample.mean <- function(d, x) {
  mean(d[x])
}
boot_curr <- boot(target.data, statistic=sample.mean, R=100)
ci.curr <- boot.ci(boot_curr, conf=0.95, type="basic")
low <- ci.curr$basic[4]
high <- ci.curr$basic[5]
points(years, mean(target.data), pch=20, cex=0.6, col="black")
library("Hmisc")
errbar(x=years, y=mean(target.data), high, low, add=TRUE)

#text and arrows
text(x=25, y=3100, labels="Current songs", pos=3)
arrows(x0=25, y0=3300, x1=36, y1=mean(target.data), length=0.1)
}
##Hypothesis testing with h.test() This function allows hypothesis testing with
##SongEvo. To test if measured songs from two time points evolved through
##mechanisms described in the model (e.g. drift or selection), users initialize
##the model with historical data, parameterize the model based on their
##understanding of the mechanisms, and test if subsequently observed or
##predicted data match the simulated data. The output data list includes two
##measures of accuracy: the proportion of observed points that fall within the
##confidence intervals of the simulated data and the residuals between simulated
##and observed population trait means. Precision is measured as the residuals
##between simulated and observed population trait variances. We tested the
##hypothesis that songs of Z. l. nuttalli in Bear Valley, CA evolved through
##cultural drift from 1969 to 2005.

# Prepare initial song data for Bear Valley.

starting.trait <- subset(song.data, Population=="Bear Valley" & Year==1969)$Trill.FBW
starting.trait2 <- c(starting.trait, rnorm(n.territories-length(starting.trait), 
                                           mean=mean(starting.trait), sd=sd(starting.trait)))

init.inds <- data.frame(id = seq(1:n.territories), age = 2, trait = starting.trait2)
init.inds$x1 <-  round(runif(n.territories, min=-122.481858, max=-122.447270), digits=8)
init.inds$y1 <-  round(runif(n.territories, min=37.787768, max=37.805645), digits=8)
# Specify and call SongEvo() with test data

# SongEvo3 <- SongEvo(init.inds = init.inds, iteration = iteration, steps =
# years,  timetep = timetep, n.territories = n.territories, terr.turnover =
# terr.turnover, learning.method = learning.method, integrate.dist =
# integrate.dist, learning.error.d = learning.error.d, learning.error.sd =
# learning.error.sd, mortality.a = mortality.a, mortality.j = mortality.j,
# lifespan = lifespan, phys.lim.min = phys.lim.min, phys.lim.max = phys.lim.max,
# male.fledge.n.mean = male.fledge.n.mean, male.fledge.n.sd = male.fledge.n.sd,
# male.fledge.n = male.fledge.n, disp.age = disp.age, disp.distance.mean =
# disp.distance.mean, disp.distance.sd = disp.distance.sd, mate.comp =
# mate.comp, prin = prin, all)
SongEvo3_alt=list(parms=glo.parms)
parms_3alt <- list(init.inds = init.inds, 
                   iteration = 10,
                   steps = 36,  # years / timestep
                   timestep = 1, 
                   n.territories = 40,
                   terr.turnover = 0.5, 
                   learning.method = "integrate", 
                   integrate.dist = 0.1, 
                   lifespan = NA, 
                   mate.comp = FALSE, 
                   prin = FALSE, 
                   all = TRUE)
SongEvo3_alt$parms[names(parms_3alt)]=parms_3alt
SongEvo3_alt$res=do.call(SongEvo,SongEvo3_alt$parms)
str(SongEvo3_alt)

#Save result for repeatable testing
save(SongEvo3_alt,
     file=(most_recent_SongEvo3_alt<-paste0("SongEvo3_alt_",format(Sys.time(),format="%Y%m%d%H%M"),".RData")))

SongEvo3=SongEvo3_alt$res


# Specify and call h.test()

target.data <- subset(song.data, Population=="Bear Valley" & Year==2005)$Trill.FBW
h.test1 <- h.test(summary.results=SongEvo3$summary.results, ts=ts, target.data=target.data)
# The output data list includes two measures of accuracy: the proportion of
# observed points that fall within the confidence intervals of the simulated
# data and the residuals between simulated and observed population trait means.
# Precision is measured as the residuals between simulated and observed
# population trait variances.

# Eighty percent of the observed data fell within the central 95% of the
# simulated values, providing support for the hypothesis that cultural drift as
# described in this model is sufficient to describe the evolution of trill
# frequency bandwidth in this population.
h.test1

# We can plot simulated data in relation to measured data.

# Similar Plot to Figure 4
plot(SongEvo3$summary.results[1, , "trait.pop.mean"], xlab="Year", ylab="Bandwidth (Hz)", xaxt="n", type="n", 
     xlim=c(-0.5, 35.5), ylim=range(SongEvo3$summary.results[, , "trait.pop.mean"], na.rm=TRUE))
for(p in 1:SongEvo3_alt$parms$iteration){
  lines(SongEvo3$summary.results[p, , "trait.pop.mean"], col="light gray")
}
freq.mean <- apply(SongEvo3$summary.results[, , "trait.pop.mean"], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))#, tcl=-0.25, mgp=c(2,0.5,0))

#Plot 95% quantiles (which are similar to credible intervals)
quant.means <- apply (SongEvo3$summary.results[, , "trait.pop.mean"], MARGIN=2, quantile, probs=c(0.95, 0.05), R=600, na.rm=TRUE)
lines(quant.means[1,], col="blue", lty=2)
lines(quant.means[2,], col="blue", lty=2)

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
points(rep(ts, length(target.data)), target.data)

library("boot")
sample.mean <- function(d, x) {
  mean(d[x])
}
boot_curr <- boot(target.data, statistic=sample.mean, R=100)#, strata=mn.res$iteration)	
ci.curr <- boot.ci(boot_curr, conf=0.95, type="basic")
low <- ci.curr$basic[4]
high <- ci.curr$basic[5]
points(years, mean(target.data), pch=20, cex=0.6, col="black")
library("Hmisc")
errbar(x=years, y=mean(target.data), high, low, add=TRUE)

#text and arrows
text(x=11, y=2850, labels="Historical songs", pos=1)
arrows(x0=5, y0=2750, x1=0.4, y1=mean(starting.trait), length=0.1)
text(x=25, y=2900, labels="Current songs", pos=1)
arrows(x0=25, y0=2920, x1=years, y1=mean(target.data), length=0.1)
