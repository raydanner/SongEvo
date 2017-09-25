##Vignette code from Readme

#Load functions from files, keeping source references for future debugging

SongEvo <- dget("SongEvo.R", TRUE)
par.sens <- dget("par.sens.R", TRUE)
par.opt <- dget("par.opt.R", TRUE)
mod.val <- dget("mod.val.R", TRUE)
h.test <- dget("h.test.R", TRUE)

# Load the example data: WCSP 
# To explore the SongEvo package, we will use a database of songs from 
# Nuttall’s white-crowned sparrow (Zonotrichia leucophrys nuttalli) 
# recorded at three locations in 1969 and 2005.

data("WCSP")

## discard malformed glo.parms entry "global.parms$learning.error.d"

glo.parms=glo.parms[ ! names(glo.parms) %in% "global.parms$learning.error.d" ]

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

plot(SongEvo1$summary.results[1, , "sample.n"], xlab="Year", ylab="Abundance", type="n", xaxt="n", ylim=c(0, max(SongEvo1$summary.results[, , "sample.n"], na.rm=TRUE)))
axis(side=1, at=seq(0, 40, by=5), labels=seq(1970, 2010, by=5))
for(p in 1:10){
  lines(SongEvo1$summary.results[p, , "sample.n"], col="light gray")
}
n.mean <- apply(SongEvo1$summary.results[, , "sample.n"], 2, mean, na.rm=TRUE)
lines(n.mean, col="red")

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "sample.n"], MARGIN=2, quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
lines(quant.means[1,], col="red", lty=2)
lines(quant.means[2,], col="red", lty=2)
library("Hmisc")

plot(SongEvo1$summary.results[1, , "trait.pop.mean"], xlab="Year", ylab="Bandwidth (Hz)", xaxt="n", type="n", xlim=c(-0.5, 36), ylim=c(min(SongEvo1$summary.results[, , "trait.pop.mean"], na.rm=TRUE), max(SongEvo1$summary.results[, , "trait.pop.mean"], na.rm=TRUE)))
iteration=dim(SongEvo1$summary.results)[1]
for(p in 1:iteration){
  lines(SongEvo1$summary.results[p, , "trait.pop.mean"], col="light gray")
}
freq.mean <- apply(SongEvo1$summary.results[, , "trait.pop.mean"], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))#, tcl=-0.25, mgp=c(2,0.5,0))

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "trait.pop.mean"], MARGIN=2, quantile, probs=c(0.95, 0.05), R=600, na.rm=TRUE)
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
plot(SongEvo1$summary.results[1, , "trait.pop.variance"], xlab="Year", ylab="Bandwidth Variance (Hz)", type="n", xaxt="n", ylim=c(min(SongEvo1$summary.results[, , "trait.pop.variance"], na.rm=TRUE), max(SongEvo1$summary.results[, , "trait.pop.variance"], na.rm=TRUE)))
axis(side=1, at=seq(0, 40, by=5), labels=seq(1970, 2010, by=5))
for(p in 1:iteration){
  lines(SongEvo1$summary.results[p, , "trait.pop.variance"], col="light gray")
}
n.mean <- apply(SongEvo1$summary.results[, , "trait.pop.variance"], 2, mean, na.rm=TRUE)
lines(n.mean, col="green")

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "trait.pop.variance"], MARGIN=2, quantile, probs=c(0.975, 0.025), R=600, na.rm=TRUE)
lines(quant.means[1,], col="green", lty=2)
lines(quant.means[2,], col="green", lty=2)

library("reshape2")
library("lattice")
library("sp")

all.inds1 <- subset(SongEvo1$all.inds, iteration==1)
w <- dcast(as.data.frame(all.inds1), id ~ timestep, value.var="trait", fill=0)
all.inds1w <- merge(all.inds1, w, by="id")
names(all.inds1w) <- c(names(all.inds1), paste("Ts", seq(1:36), sep=""))

rbPal <- colorRampPalette(c('blue','red')) #Create a function to generate a continuous color palette

# Plot maps, including a separate panel for each timestep (each of 36 years). Our example shows that individuals move across the landscape and that regional dialects evolve and move. The x-axis is longitude, the y-axis is latitude, and the color ramp indicates trill bandwidth in Hz.
spplot(all.inds1w[,-c(1:ncol(all.inds1))], as.table=TRUE, cuts=c(0, seq(from=1500, to=4500, by=10)), ylab="", col.regions=c("transparent", rbPal(1000)), #cuts specifies that the first level (e.g. <1500) is transparent.
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
