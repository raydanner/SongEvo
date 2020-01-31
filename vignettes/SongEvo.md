‘SongEvo’ package
================
Raymond M. Danner, Elizabeth P. Derryberry, Graham E. Derryberry,  
Julie E. Danner, Sara J. Evans, David A. Luther
2020-01-31

## Introduction

SongEvo simulates the cultural evolution of quantitative traits of bird
song. SongEvo is an individual- (agent-) based model. SongEvo is
spatially-explicit and can be parameterized with, and tested against,
measured song data. Functions are available for model implementation,
sensitivity analyses, parameter optimization, model validation, and
hypothesis testing.

### Overview of Functions

1.  `SongEvo` implements the model
2.  `par.sens` allows sensitivity analyses
3.  `par.opt` allows parameter optimization
4.  `mod.val` allows model validation
5.  `h.test` allows hypothesis testing

## Getting Started

### Load and attach SongEvo package

``` r
library(SongEvo)
```

### Functions

`SongEvo` implements the model `par.sens` allows sensitivity analyses
`par.opt` allows parameter optimization `mod.val` allows model
validation `h.test` allows hypothesis testing

## Examples

### EXAMPLE 1

### Load the example data: song.data and global parameters

To explore the SongEvo package, we will use a database of songs from
Nuttall’s white-crowned sparrow (*Zonotrichia leucophrys nuttalli*)
recorded at three locations in 1969 and 2005.

``` r
data("song.data")
```

Examine global parameters. Global parameters describe our understanding
of the system and may be measured or hypothesized. They are called
“global” because they are used by many many functions and subroutines
within functions. For descriptions of all adjustable parameters, see
`?song.data`.

``` r
data("glo.parms")
glo.parms$mortality.a.m <- glo.parms$mortality.a.f <- glo.parms$mortality.a
glo.parms$mortality.j.m <- glo.parms$mortality.j.f <- glo.parms$mortality.j
glo.parms$male.fledge.n.mean <- glo.parms$male.fledge.n.mean*2
glo.parms$male.fledge.n.sd <- glo.parms$male.fledge.n.sd*2
glo.parms <- glo.parms[!names(glo.parms) %in% c("mortality.a","mortality.j")]
str(glo.parms)
#> List of 17
#>  $ learning.error.d  : num 0
#>  $ learning.error.sd : num 430
#>  $ n.territories     : num 40
#>  $ lifespan          : num 2.08
#>  $ phys.lim.min      : num 1559
#>  $ phys.lim.max      : num 4364
#>  $ male.fledge.n.mean: num 2.7
#>  $ male.fledge.n.sd  : num 1
#>  $ disp.age          : num 2
#>  $ disp.distance.mean: num 110
#>  $ disp.distance.sd  : num 100
#>  $ terr.turnover     : num 0.5
#>  $ male.fledge.n     : num [1:40] 1 1 2 1 0 2 2 2 2 1 ...
#>  $ mortality.a.f     : num 0.468
#>  $ mortality.a.m     : num 0.468
#>  $ mortality.j.f     : num 0.5
#>  $ mortality.j.m     : num 0.5
```

Share global parameters with the global environment. We make these
parameters available in the global environment so that we can access
them with minimal code.

``` r
list2env(glo.parms, globalenv())
#> <environment: R_GlobalEnv>
```

#### Examine song data

Data include the population name (Bear Valley, PRBO, or Schooner), year
of song recording (1969 or 2005), and the frequency bandwidth of the
trill.

``` r
str(song.data)
#> 'data.frame':    89 obs. of  3 variables:
#>  $ Population: Factor w/ 3 levels "Bear Valley",..: 3 3 3 3 3 3 3 3 3 3 ...
#>  $ Year      : int  1969 1969 1969 1969 1969 1969 1969 1969 1969 1969 ...
#>  $ Trill.FBW : num  3261 2494 2806 2878 2758 ...
```

### Simulate bird song evolution with `SongEvo()`

#### Define initial individuals

In this example, we use songs from individual birds recorded in one
population (PRBO) in the year 1969, which we will call
`starting.trait`.

``` r
starting.trait <- subset(song.data, Population=="PRBO" & Year==1969)$Trill.FBW
```

We want a starting population of 40 individuals, so we generate
additional trait values to complement those from the existing 30
individuals. Then we create a data frame that includes a row for each
individual; we add identification numbers, ages, and geographical
coordinates for each
individual.

``` r
starting.trait2 <- c(starting.trait, rnorm(n.territories-length(starting.trait), 
    mean=mean(starting.trait), sd=sd(starting.trait)))
init.inds <- data.frame(id = seq(1:n.territories), age = 2, trait = starting.trait2)
init.inds$x1 <-  round(runif(n.territories, min=-122.481858, max=-122.447270), digits=8)
init.inds$y1 <-  round(runif(n.territories, min=37.787768, max=37.805645), digits=8)
```

#### Specify and call the SongEvo model

`SongEvo()` includes several settings, which we specify before running
the model. For this example, we run the model for 10 iterations, over 36
years (i.e. 1969–2005). When conducting research with `SongEvo()`, users
will want to increase the number iterations (e.g. to 100 or 1000). Each
timestep is one year in this model (i.e. individuals complete all
components of the model in 1 year). We specify territory turnover rate
here as an example of how to adjust parameter values. We could adjust
any other parameter value here also. The learning method specifies that
individuals integrate songs heard from adults within the specified
integration distance (intigrate.dist, in kilometers). In this example,
we do not includ a lifespan, so we assign it NA. In this example, we do
not model competition for mates, so specify it as FALSE. Last, specify
all as TRUE in order to save data for every single simulated individual
because we will use those data later for mapping. If we do not need data
for each individual, we set all to FALSE because the all.inds data.frame
becomes very large\!

``` r
iteration <- 10
years <- 36
timestep <- 1
terr.turnover <- 0.5
integrate.dist <- 0.1
lifespan <- NA
mate.comp <- FALSE
prin <- FALSE
all <- TRUE
```

Now we call SongEvo with our specifications and save it in an object
called
SongEvo1.

``` r
SongEvo1 <- SongEvo(init.inds = init.inds, females = 1.0, iteration = iteration, steps = years,  
    timestep = timestep, n.territories = n.territories, terr.turnover = terr.turnover, 
    integrate.dist = integrate.dist, 
    learning.error.d = learning.error.d, learning.error.sd = learning.error.sd, 
    mortality.a.m = mortality.a.m, mortality.a.f = mortality.a.f,
    mortality.j.m = mortality.j.m, mortality.j.f = mortality.j.f, lifespan = lifespan, 
    phys.lim.min = phys.lim.min, phys.lim.max = phys.lim.max, 
    male.fledge.n.mean = male.fledge.n.mean, male.fledge.n.sd = male.fledge.n.sd, male.fledge.n = male.fledge.n,
    disp.age = disp.age, disp.distance.mean = disp.distance.mean, disp.distance.sd = disp.distance.sd, 
    mate.comp = mate.comp, prin = prin, all = TRUE)
```

#### Examine results from SongEvo model

The model required the following time to run on your computer:

``` r
SongEvo1$time
#>    user  system elapsed 
#>  36.553   0.660  39.598
```

Three main objects hold data regarding the SongEvo model. Additional
objects are used temporarily within modules of the model.

First, currently alive individuals are stored in a data frame called
“inds.” Values within “inds” are updated throughout each of the
iterations of the model, and “inds” can be viewed after the model is
completed.

``` r
head(SongEvo1$inds, min(5,nrow(SongEvo1$inds)))
#>                 coordinates   id age    trait        x1       y1
#> M1586 (-122.4756, 37.79727) 1586   6 3284.121 -122.4756 37.79727
#> M1593 (-122.4815, 37.80091) 1593   6 3097.741 -122.4815 37.80091
#> M1596 (-122.4517, 37.80072) 1596   6 3362.291 -122.4517 37.80072
#> M1631 (-122.4718, 37.79526) 1631   6 2987.358 -122.4718 37.79526
#> M1638  (-122.482, 37.79895) 1638   5 3444.528 -122.4820 37.79895
#>       male.fledglings female.fledglings territory father sex fitness learn.dir
#> M1586               1                 0         1   1334   M       1         0
#> M1593               1                 3         1   1459   M       1         0
#> M1596               0                 0         0   1486   M       1         0
#> M1631               2                 2         1   1571   M       1         0
#> M1638               0                 0         0   1346   M       1         0
#>              x0       y0
#> M1586 -122.4755 37.79644
#> M1593 -122.4823 37.79847
#> M1596 -122.4525 37.80115
#> M1631 -122.4722 37.79508
#> M1638 -122.4823 37.79873
```

Second, an array (i.e. a multi-dimensional table) entitled
“summary.results” includes population summary values for each time
step (dimension 1) in each iteration (dimension 2) of the model.
Population summary values are contained in five additional dimensions:
population size for each time step of each iteration (“sample.n”), the
population mean and variance of the song feature studied
(“trait.pop.mean” and “trait.pop.variance”), with associated lower
(“lci”) and upper (“uci”) confidence intervals.

``` r
dimnames(SongEvo1$summary.results)
#> $iteration
#>  [1] "iteration 1"  "iteration 2"  "iteration 3"  "iteration 4"  "iteration 5" 
#>  [6] "iteration 6"  "iteration 7"  "iteration 8"  "iteration 9"  "iteration 10"
#> 
#> $step
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
#> [16] "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30"
#> [31] "31" "32" "33" "34" "35" "36"
#> 
#> $feature
#> [1] "sample.n"           "trait.pop.mean"     "trait.pop.variance"
#> [4] "lci"                "uci"
```

Third, individual values may optionally be concatenated and saved to one
data frame entitled “all.inds.” all.inds can become quite large, and is
therefore only recommended if additional data analyses are desired.

``` r
head(SongEvo1$all.inds,  min(5,nrow(SongEvo1$all.inds)))
#>                   coordinates id age  trait        x1       y1 male.fledglings
#> I1.T1.1 (-122.4636, 37.80077)  1   2 4004.8 -122.4636 37.80077               0
#> I1.T1.2  (-122.4812, 37.7882)  2   2 3765.0 -122.4812 37.78820               1
#> I1.T1.3 (-122.4552, 37.79986)  3   2 3237.4 -122.4552 37.79986               1
#> I1.T1.4 (-122.4709, 37.79962)  4   2 3621.1 -122.4709 37.79962               1
#> I1.T1.5 (-122.4473, 37.79041)  5   2 3285.4 -122.4473 37.79041               0
#>         female.fledglings territory father sex fitness learn.dir x0 y0 timestep
#> I1.T1.1                 1         1      0   M       1         0  0  0        1
#> I1.T1.2                 0         1      0   M       1         0  0  0        1
#> I1.T1.3                 1         1      0   M       1         0  0  0        1
#> I1.T1.4                 0         1      0   M       1         0  0  0        1
#> I1.T1.5                 0         1      0   M       1         0  0  0        1
#>         iteration
#> I1.T1.1         1
#> I1.T1.2         1
#> I1.T1.3         1
#> I1.T1.4         1
#> I1.T1.5         1
```

#### Simulated population size

We see that the simulated population size remains relatively stable over
the course of 36 years. This code uses the summary.results
array.

``` r
plot(SongEvo1$summary.results[1, , "sample.n"], xlab="Year", ylab="Abundance", type="n", 
    xaxt="n", ylim=c(0, max(SongEvo1$summary.results[, , "sample.n"], na.rm=TRUE)))
axis(side=1, at=seq(0, 40, by=5), labels=seq(1970, 2010, by=5))
    for(p in 1:iteration){
        lines(SongEvo1$summary.results[p, , "sample.n"], col="light gray")
        }
n.mean <- apply(SongEvo1$summary.results[, , "sample.n"], 2, mean, na.rm=TRUE)
lines(n.mean, col="red")

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "sample.n"], MARGIN=2, quantile, 
    probs=c(0.975, 0.025), R=600, na.rm=TRUE)
lines(quant.means[1,], col="red", lty=2)
lines(quant.means[2,], col="red", lty=2)
```

![](SongEvo_files/figure-gfm/Examine%20population%20size%20over%20time-1.png)<!-- -->

Load Hmisc package for plotting functions.

``` r
library("Hmisc")
```

#### Simulated trait values

We see that the mean trait values per iteration varied widely, though
mean trait values over all iterations remained relatively stable. This
code uses the summary.results
array.

``` r
plot(SongEvo1$summary.results[1, , "trait.pop.mean"], xlab="Year", ylab="Bandwidth (Hz)", 
    xaxt="n", type="n", xlim=c(-0.5, 36), 
    ylim=c(min(SongEvo1$summary.results[, , "trait.pop.mean"], na.rm=TRUE), 
    max(SongEvo1$summary.results[, , "trait.pop.mean"], na.rm=TRUE)))
    for(p in 1:iteration){
        lines(SongEvo1$summary.results[p, , "trait.pop.mean"], col="light gray")
        }
freq.mean <- apply(SongEvo1$summary.results[, , "trait.pop.mean"], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))#, tcl=-0.25, mgp=c(2,0.5,0))

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "trait.pop.mean"], MARGIN=2, quantile, 
    probs=c(0.95, 0.05), R=600, na.rm=TRUE)
lines(quant.means[1,], col="blue", lty=2)
lines(quant.means[2,], col="blue", lty=2)

#plot mean and CI for historic songs.  
 #plot original song values
library("boot")
#> 
#> Attaching package: 'boot'
#> The following object is masked from 'package:survival':
#> 
#>     aml
#> The following object is masked from 'package:lattice':
#> 
#>     melanoma
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
```

![](SongEvo_files/figure-gfm/Examine%20simulated%20trait%20evolution-1.png)<!-- -->

#### Trait variance

We see that variance for each iteration per year increased in the first
few years and then stabilized. This code uses the summary.results array.

``` r
 #plot variance for each iteration per year
plot(SongEvo1$summary.results[1, , "trait.pop.variance"], xlab="Year", 
    ylab="Bandwidth Variance (Hz)", type="n", xaxt="n", 
    ylim=c(min(SongEvo1$summary.results[, , "trait.pop.variance"], na.rm=TRUE), 
    max(SongEvo1$summary.results[, , "trait.pop.variance"], na.rm=TRUE)))
axis(side=1, at=seq(0, 40, by=5), labels=seq(1970, 2010, by=5))
    for(p in 1:iteration){
        lines(SongEvo1$summary.results[p, , "trait.pop.variance"], col="light gray")
        }
n.mean <- apply(SongEvo1$summary.results[, , "trait.pop.variance"], 2, mean, na.rm=TRUE)
lines(n.mean, col="green")

#Plot 95% quantiles
quant.means <- apply (SongEvo1$summary.results[, , "trait.pop.variance"], MARGIN=2, quantile, 
    probs=c(0.975, 0.025), R=600, na.rm=TRUE)
lines(quant.means[1,], col="green", lty=2)
lines(quant.means[2,], col="green", lty=2)
```

![](SongEvo_files/figure-gfm/Examine%20simulated%20trait%20variance-1.png)<!-- -->

#### Maps

The simulation results include geographical coordinates and are in a
standard spatial data format, thus allowing calculation of a wide
variety of spatial statistics.

Load packages for making maps.

``` r
library("sp")
library("reshape2")
library("lattice")
```

Convert data frame from long to wide format. This is necessary for
making a multi-panel plot.

``` r
all.inds1 <- subset(SongEvo1$all.inds, SongEvo1$all.inds$iteration==1)
w <- dcast(as.data.frame(all.inds1), id ~ timestep, value.var="trait", fill=0)
all.inds1w <- merge(all.inds1, w, by="id")
years.SongEvo1 <- (dim(w)[2]-1 )
names(all.inds1w@data)[-(1:length(all.inds1@data))] <-paste("Ts", 1:(dim(w)[2]-1 ), sep="")
```

Create a function to generate a continuous color palette–we will use the
palette in the next call to make color ramp to represent the trait
value.

``` r
rbPal <- colorRampPalette(c('blue','red')) #Create a function to generate a continuous color palette
```

Plot maps, including a separate panel for each timestep (each of 36
years). Our example shows that individuals move across the landscape and
that regional dialects evolve and move. The x-axis is longitude, the
y-axis is latitude, and the color ramp indicates trill bandwidth in Hz.

``` r
spplot(all.inds1w[,-c(1:ncol(all.inds1))], as.table=TRUE, 
    cuts=c(0, seq(from=1500, to=4500, by=10)), ylab="", 
    col.regions=c("transparent", rbPal(1000)), 
    #cuts specifies that the first level (e.g. <1500) is transparent.
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
```

![](SongEvo_files/figure-gfm/Plot%20map-1.png)<!-- -->

In addition, you can plot simpler multi-panel maps that do not take
advantage of the spatial data class.

``` r
 #Lattice plot (not as a spatial frame)
it1 <- subset(SongEvo1$all.inds, iteration==1)
rbPal <- colorRampPalette(c('blue','red')) #Create a function to generate a continuous color palette
it1$Col <- rbPal(10)[as.numeric(cut(it1$trait, breaks = 10))]
xyplot(it1$y1~it1$x1 | it1$timestep, groups=it1$trait, asp="iso", col=it1$Col, 
    xlab="Longitude", ylab="Latitude")
```

![](SongEvo_files/figure-gfm/Plot%20maps%20without%20using%20the%20spatial%20data%20class-1.png)<!-- -->

### Test model sensitivity with `par.sens()`

This function allows testing the sensitivity of SongEvo to different
parameter values.

#### Specify and call `par.sens()`

Here we test the sensitivity of the Acquire a Territory submodel to
variation in territory turnover rates, ranging from 0.8–1.2 times the
published rate (40–60% of territories turned over). The call for the
par.sens function has a format similar to SongEvo. The user specifies
the parameter to test and the range of values for that parameter. The
function currently allows examination of only one parameter at a time
and requires at least two iterations.

``` r
parm <- "terr.turnover"
par.range = seq(from=0.4, to=0.6, by=0.025)
sens.results <- NULL
```

Now we call the par.sens function with our specifications.

``` r
extra_parms <- list(init.inds = init.inds, 
                    females = 1,  # New in SongEvo v2
                    timestep = 1, 
                    n.territories = nrow(init.inds),
                    integrate.dist = 0.1,
                    lifespan = NA, 
                    terr.turnover = 0.5, 
                    mate.comp = FALSE, 
                    prin = FALSE,
                    all = TRUE,
                    # New in SongEvo v2
                    selectivity = 3,
                    content.bias = FALSE,
                    n.content.bias.loc = "all",
                    content.bias.loc = FALSE,
                    content.bias.loc.ranges = FALSE,
                    affected.traits = FALSE,
                    conformity.bias = FALSE,
                    prestige.bias=FALSE,
                    learn.m="default",
                    learn.f="default",
                    learning.error.d=0,
                    learning.error.sd=200)
global_parms_key <- which(!names(glo.parms) %in% names(extra_parms))
extra_parms[names(glo.parms[global_parms_key])]=glo.parms[global_parms_key]
par.sens1 <- par.sens(parm = parm, par.range = par.range, 
                      iteration = iteration, steps = years, mate.comp = FALSE, 
                      fixed_parms=extra_parms[names(extra_parms)!=parm], all = TRUE)
#> [1] "terr.turnover =  0.4"
#> [1] "terr.turnover =  0.425"
#> [1] "terr.turnover =  0.45"
#> [1] "terr.turnover =  0.475"
#> [1] "terr.turnover =  0.5"
#> [1] "terr.turnover =  0.525"
#> [1] "terr.turnover =  0.55"
#> [1] "terr.turnover =  0.575"
#> [1] "terr.turnover =  0.6"
```

#### Examine par.sens results

Examine results objects, which include two arrays:

The first array, `sens.results`, contains the SongEvo model results for
each parameter. It has the following dimensions:

``` r
dimnames(par.sens1$sens.results)
#> [[1]]
#>  [1] "iteration 1"  "iteration 2"  "iteration 3"  "iteration 4"  "iteration 5" 
#>  [6] "iteration 6"  "iteration 7"  "iteration 8"  "iteration 9"  "iteration 10"
#> 
#> [[2]]
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
#> [16] "16" "17" "18" "19" "20" "21" "22" "23" "24" "25" "26" "27" "28" "29" "30"
#> [31] "31" "32" "33" "34" "35" "36"
#> 
#> [[3]]
#> [1] "sample.n"           "trait.pop.mean"     "trait.pop.variance"
#> [4] "lci"                "uci"               
#> 
#> [[4]]
#> [1] "par.val 0.4"   "par.val 0.425" "par.val 0.45"  "par.val 0.475"
#> [5] "par.val 0.5"   "par.val 0.525" "par.val 0.55"  "par.val 0.575"
#> [9] "par.val 0.6"
```

The second array, `sens.results.diff` contains the quantile range of
trait values across iterations within a parameter value. It has the
following dimensions:

``` r
dimnames(par.sens1$sens.results.diff)
#> [[1]]
#> [1] "par.val 0.4"   "par.val 0.425" "par.val 0.45"  "par.val 0.475"
#> [5] "par.val 0.5"   "par.val 0.525" "par.val 0.55"  "par.val 0.575"
#> [9] "par.val 0.6"  
#> 
#> [[2]]
#>  [1] "Quantile diff 1"  "Quantile diff 2"  "Quantile diff 3"  "Quantile diff 4" 
#>  [5] "Quantile diff 5"  "Quantile diff 6"  "Quantile diff 7"  "Quantile diff 8" 
#>  [9] "Quantile diff 9"  "Quantile diff 10" "Quantile diff 11" "Quantile diff 12"
#> [13] "Quantile diff 13" "Quantile diff 14" "Quantile diff 15" "Quantile diff 16"
#> [17] "Quantile diff 17" "Quantile diff 18" "Quantile diff 19" "Quantile diff 20"
#> [21] "Quantile diff 21" "Quantile diff 22" "Quantile diff 23" "Quantile diff 24"
#> [25] "Quantile diff 25" "Quantile diff 26" "Quantile diff 27" "Quantile diff 28"
#> [29] "Quantile diff 29" "Quantile diff 30" "Quantile diff 31" "Quantile diff 32"
#> [33] "Quantile diff 33" "Quantile diff 34" "Quantile diff 35" "Quantile diff 36"
```

To assess sensitivity of SongEvo to a range of parameter values, plot
the range in trait quantiles per year by the parameter value. We see
that territory turnover values of 0.4–0.6 provided means and quantile
ranges of trill bandwidths that are similar to those obtained with the
published estimate of 0.5, indicating that the Acquire a Territory
submodel is robust to realistic variation in those parameter values.

In the figure, solid gray and black lines show the quantile range of
song frequency per year over all iterations as parameterized with the
published territory turnover rate (0.5; thick black line) and a range of
values from 0.4 to 0.6 (in steps of 0.05, light to dark gray). Orange
lines show the mean and 2.5th and 97.5th quantiles of all quantile
ranges.

``` r
 #plot of range in trait quantiles by year by parameter value
plot(1:years, par.sens1$sens.results.diff[1,], ylim=c(min(par.sens1$sens.results.diff, 
    na.rm=TRUE), max(par.sens1$sens.results.diff, na.rm=TRUE)), type="l", 
    ylab="Quantile range (Hz)", xlab="Year", col="transparent", xaxt="n")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))

  #Make a continuous color ramp from gray to black
grbkPal <- colorRampPalette(c('gray','black'))
  
  #Plot a line for each parameter value
for(i in 1:length(par.range)){
lines(1:years, par.sens1$sens.results.diff[i,], type="l", 
    col=grbkPal(length(par.range))[i])
}

  #Plot values from published parameter values
lines(1:years, par.sens1$sens.results.diff[2,], type="l", col="black", lwd=4)

  #Calculate and plot mean and quantiles
quant.mean <- apply(par.sens1$sens.results.diff, 2, mean, na.rm=TRUE)
lines(quant.mean, col="orange")

#Plot 95% quantiles (which are similar to credible intervals)
  #95% quantiles of population means (narrower)
quant.means <- apply (par.sens1$sens.results.diff, MARGIN=2, quantile, 
    probs=c(0.975, 0.025), R=600, na.rm=TRUE)
lines(quant.means[1,], col="orange", lty=2)
lines(quant.means[2,], col="orange", lty=2)
```

![](SongEvo_files/figure-gfm/plot%20of%20range%20in%20trait%20quantiles%20by%20year%20by%20parameter%20value-1.png)<!-- -->

### Optimize parameter values with `par.opt()`

This function follows par.sens to help users optimize values for
imperfectly known parameters for SongEvo. The goals are to maximize
accuracy and precision of model prediction. Accuracy is quantified by
three different approaches: i) the mean of absolute residuals of the
predicted population mean values in relation to target data
(e.g. observed or hypothetical values (smaller absolute residuals
indicate a more accurate model)), ii) the difference between the
bootstrapped mean of predicted population means and the mean of the
target data, and iii) the proportion of simulated population trait means
that fall within (i.e. are “contained by”) the confidence intervals of
the target data (a higher proportion indicates greater accuracy).
Precision is measured with the residuals of the predicted population
variance to the variance of target data (smaller residuals indicate a
more precise
model).

#### Prepare current song values

``` r
target.data <- subset(song.data, Population=="PRBO" & Year==2005)$Trill.FBW
```

#### Specify and call `par.opt()`

Users specify the timestep (“ts”) at which to compare simulated trait
values to target trait data (“target.data”) and save the results in an
object (called `par.opt1` here).

``` r
ts <- years
par.opt1 <- par.opt(sens.results=par.sens1$sens.results, ts=ts, 
    target.data=target.data, par.range=par.range)
```

Examine results objects (residuals and target match).

``` r
par.opt1$Residuals
#> , , Residuals of mean
#> 
#>               Iteration 1 Iteration 2 Iteration 3 Iteration 4 Iteration 5
#> par.val 0.4      375.1796    172.2582   486.93384    293.5727   299.89124
#> par.val 0.425    393.2643    484.5063   434.94376    250.6562   521.20525
#> par.val 0.45     208.7989    447.4654   815.17267    328.1507   540.09510
#> par.val 0.475    461.3758    342.8075   260.04323    388.3182   127.79368
#> par.val 0.5      375.1772    355.9300    56.42830    284.2852   392.95899
#> par.val 0.525    372.0067    526.4238     5.27347    232.8296    78.82136
#> par.val 0.55     296.9348    300.1850   532.16051    209.0245    17.16042
#> par.val 0.575    110.4704    307.8186   385.18824    135.2806   171.29309
#> par.val 0.6      437.8608    602.9041   139.47967    331.1494   183.35176
#>               Iteration 6 Iteration 7 Iteration 8 Iteration 9 Iteration 10
#> par.val 0.4      363.5481    264.4108   268.74419   318.51976    514.46377
#> par.val 0.425    213.6552    481.6147   140.15586    29.05265    160.90980
#> par.val 0.45     667.2729    249.2244   218.55530   480.92863    590.86455
#> par.val 0.475    344.1920    378.4391    35.62631   258.13527    179.25004
#> par.val 0.5      143.8095    253.1561   346.94532   185.10014    265.70876
#> par.val 0.525    250.9261    403.2667   544.05379   337.30352     84.29053
#> par.val 0.55     192.9230    464.1890   199.02752   364.67413    458.63284
#> par.val 0.575    248.5390    568.7111   489.62999   337.19055    387.88652
#> par.val 0.6      225.5181    106.8162   399.71736   520.10971    481.25120
#> 
#> , , Residuals of variance
#> 
#>               Iteration 1 Iteration 2 Iteration 3 Iteration 4 Iteration 5
#> par.val 0.4      3856.590  10839.0036    858.4907    5480.707    205.0313
#> par.val 0.425    6974.166   8348.9359  13508.9644    7973.211  17586.8458
#> par.val 0.45     2906.645  15223.2451    572.3942   24314.537   9449.4499
#> par.val 0.475    3117.883    323.9328   1128.4665    2849.572   8725.9833
#> par.val 0.5      6628.449   1630.8279   6447.6606   28137.731   4922.0719
#> par.val 0.525   14715.638  13755.8743  15204.6969   14808.396  12393.7876
#> par.val 0.55     3895.445   6389.2736  12818.3057    8390.959    303.3174
#> par.val 0.575    2049.725   8417.3492    854.6628   22361.477   8541.6948
#> par.val 0.6      6686.687  14260.9655   7648.5817   10428.393   3074.9688
#>               Iteration 6 Iteration 7 Iteration 8 Iteration 9 Iteration 10
#> par.val 0.4      1666.563   9958.3390   30733.739    8992.862     9430.897
#> par.val 0.425    4928.805    126.4336   12577.322    2795.613    20217.278
#> par.val 0.45    12770.417  15509.8828   23157.496    8153.749     4657.414
#> par.val 0.475   16489.731   1812.4683    8195.665    3936.264     1834.785
#> par.val 0.5      1740.912  10386.9964    6686.988   11536.405    18165.757
#> par.val 0.525   12142.004   2333.7526    2334.686    1349.201     3863.080
#> par.val 0.55    12469.707   2541.5544    3792.644    1233.921    23819.710
#> par.val 0.575    6421.378  12710.1271   31545.377    7594.087     2867.104
#> par.val 0.6     16332.996    523.5722    6339.977    3384.769    10154.483
par.opt1$Target.match
#>               Difference in means Proportion contained
#> par.val 0.4              335.7522                  0.0
#> par.val 0.425            305.1859                  0.1
#> par.val 0.45             454.6529                  0.0
#> par.val 0.475            270.4729                  0.1
#> par.val 0.5              265.9500                  0.1
#> par.val 0.525            265.6067                  0.1
#> par.val 0.55             303.4912                  0.1
#> par.val 0.575            314.2008                  0.0
#> par.val 0.6              342.8158                  0.0
```

#### Plot results of `par.opt()`

#### Accuracy of `par.opt()`

1.  Difference in
means.

<!-- end list -->

``` r
plot(par.range, par.opt1$Target.match[,1], type="l", xlab="Parameter range", 
    ylab="Difference in means (Hz)")
```

![](SongEvo_files/figure-gfm/difference%20in%20means-1.png)<!-- -->

2.  Plot proportion
contained.

<!-- end list -->

``` r
plot(par.range, par.opt1$Prop.contained, type="l", xlab="Parameter range", 
    ylab="Proportion contained")
```

![](SongEvo_files/figure-gfm/proportion%20contained-1.png)<!-- -->

3.  Calculate and plot mean and quantiles of residuals of mean trait
    values.

<!-- end list -->

``` r
res.mean.means <- apply(par.opt1$Residuals[, , 1], MARGIN=1, mean, na.rm=TRUE)
res.mean.quants <- apply (par.opt1$Residuals[, , 1], MARGIN=1, quantile, 
    probs=c(0.975, 0.025), R=600, na.rm=TRUE)
plot(par.range, res.mean.means, col="orange", ylim=c(min(par.opt1$Residuals[,,1], 
    na.rm=TRUE), max(par.opt1$Residuals[,,1], na.rm=TRUE)), type="b", 
    xlab="Parameter value (territory turnover rate)", 
    ylab="Residual of trait mean (trill bandwidth, Hz)")
points(par.range, res.mean.quants[1,], col="orange")
points(par.range, res.mean.quants[2,], col="orange")
lines(par.range, res.mean.quants[1,], col="orange", lty=2)
lines(par.range, res.mean.quants[2,], col="orange", lty=2)
```

![](SongEvo_files/figure-gfm/residuals%20of%20the%20mean-1.png)<!-- -->

#### Precision of `par.opt()`

``` r
#Calculate and plot mean and quantiles of residuals of variance of trait values
res.var.mean <- apply(par.opt1$Residuals[, , 2], MARGIN=1, mean, na.rm=TRUE)
res.var.quants <- apply (par.opt1$Residuals[, , 2], MARGIN=1, quantile, 
    probs=c(0.975, 0.025), R=600, na.rm=TRUE)
plot(par.range, res.var.mean, col="purple", 
    ylim=c(min(par.opt1$Residuals[,,2], na.rm=TRUE), 
    max(par.opt1$Residuals[,,2], na.rm=TRUE)), type="b", 
    xlab="Parameter value (territory turnover rate)", 
    ylab="Residual of trait variance (trill bandwidth, Hz)")
points(par.range, res.var.quants[1,], col="purple")
points(par.range, res.var.quants[2,], col="purple")
lines(par.range, res.var.quants[1,], col="purple", lty=2)
lines(par.range, res.var.quants[2,], col="purple", lty=2)
```

![](SongEvo_files/figure-gfm/residuals%20of%20the%20variance-1.png)<!-- -->

#### Visual inspection of accuracy and precision of `par.opt()`: plot trait values for range of parameters

``` r
par(mfcol=c(3,2),
    mar=c(2.1, 2.1, 0.1, 0.1),
    cex=0.8)
for(i in 1:length(par.range)){
plot(par.sens1$sens.results[ , , "trait.pop.mean", ], xlab="Year", ylab="Bandwidth (Hz)", 
    xaxt="n", type="n", xlim=c(-0.5, years), 
    ylim=c(min(par.sens1$sens.results[ , , "trait.pop.mean", ], na.rm=TRUE), 
    max(par.sens1$sens.results[ , , "trait.pop.mean", ], na.rm=TRUE)))
    for(p in 1:iteration){
        lines(par.sens1$sens.results[p, , "trait.pop.mean", i], col="light gray")
        }
freq.mean <- apply(par.sens1$sens.results[, , "trait.pop.mean", i], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))

#Plot 95% quantiles
quant.means <- apply (par.sens1$sens.results[, , "trait.pop.mean", i], MARGIN=2, quantile, 
    probs=c(0.95, 0.05), R=600, na.rm=TRUE)
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
text(x=3, y=max(par.sens1$sens.results[ , , "trait.pop.mean", ], na.rm=TRUE)-100, 
    labels=paste("Par = ", par.range[i], sep=""))  
}
```

![](SongEvo_files/figure-gfm/visual%20inspection%20of%20simulated%20data-1.png)<!-- -->![](SongEvo_files/figure-gfm/visual%20inspection%20of%20simulated%20data-2.png)<!-- -->

### Model validation with `mod.val()`

This function allows users to assess the validity of the specified model
by testing model performance with a population different from the
population used to build the model. The user first runs SongEvo with
initial trait values from the validation population. `mod.val()` uses
the summary.results array from SongEvo, along with target values from a
specified timestep, to calculate the same three measures of accuracy and
one measure of precision that are calculated in par.opt.

We parameterized SongEvo with initial song data from Schooner Bay, CA in
1969, and then compared simulated data to target (i.e. observed) data in
2005.

Prepare initial song data for Schooner
Bay.

``` r
starting.trait <- subset(song.data, Population=="Schooner" & Year==1969)$Trill.FBW
starting.trait2 <- c(starting.trait, rnorm(n.territories-length(starting.trait), 
    mean=mean(starting.trait), sd=sd(starting.trait)))

init.inds <- data.frame(id = seq(1:n.territories), age = 2, trait = starting.trait2)
init.inds$x1 <-  round(runif(n.territories, min=-122.481858, max=-122.447270), digits=8)
init.inds$y1 <-  round(runif(n.territories, min=37.787768, max=37.805645), digits=8)
```

Specify and call SongEvo() with validation data

``` r
iteration <- 10
years <- 36
timestep <- 1
terr.turnover <- 0.5

SongEvo2 <- SongEvo(init.inds = init.inds, females = 1.0, iteration = iteration, steps = years,  
    timestep = timestep, n.territories = n.territories, terr.turnover = terr.turnover, 
    integrate.dist = integrate.dist, 
    learning.error.d = learning.error.d, learning.error.sd = learning.error.sd, 
    mortality.a.m = mortality.a.m, mortality.a.f = mortality.a.f,
    mortality.j.m = mortality.j.m, mortality.j.f = mortality.j.f, lifespan = lifespan, 
    phys.lim.min = phys.lim.min, phys.lim.max = phys.lim.max, 
    male.fledge.n.mean = male.fledge.n.mean, male.fledge.n.sd = male.fledge.n.sd, male.fledge.n = male.fledge.n,
    disp.age = disp.age, disp.distance.mean = disp.distance.mean, disp.distance.sd = disp.distance.sd, 
    mate.comp = mate.comp, prin = prin, all = TRUE)
```

Specify and call `mod.val()`

``` r
ts <- 36
target.data <- subset(song.data, Population=="Schooner" & Year==2005)$Trill.FBW
mod.val1 <- mod.val(summary.results=SongEvo2$summary.results, ts=ts, 
    target.data=target.data)
```

Plot results from
`mod.val()`

``` r
plot(SongEvo2$summary.results[1, , "trait.pop.mean"], xlab="Year", ylab="Bandwidth (Hz)", 
    xaxt="n", type="n", xlim=c(-0.5, 36.5), 
    ylim=c(min(SongEvo2$summary.results[, , "trait.pop.mean"], na.rm=TRUE), 
    max(SongEvo2$summary.results[, , "trait.pop.mean"], na.rm=TRUE)))
    for(p in 1:iteration){
        lines(SongEvo2$summary.results[p, , "trait.pop.mean"], col="light gray")
        }
freq.mean <- apply(SongEvo2$summary.results[, , "trait.pop.mean"], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))

#Plot 95% quantiles 
quant.means <- apply (SongEvo2$summary.results[, , "trait.pop.mean"], MARGIN=2, quantile, 
    probs=c(0.95, 0.05), R=600, na.rm=TRUE)
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
```

![](SongEvo_files/figure-gfm/Plot%20results-1.png)<!-- -->

The model did reasonably well predicting trait evolution in the
validation population, suggesting that it is valid for our purposes: the
mean bandwidth was `abs(mean(target.data)-freq.mean)`Hz from the
observed values, ~21% of predicted population means fell within the 95%
confidence intervals of the observed data, and residuals of means (~545
Hz) and variances (~415181 Hz) were similar to those produced by the
training data set.

### Hypothesis testing with `h.test()`

This function allows hypothesis testing with SongEvo. To test if
measured songs from two time points evolved through mechanisms described
in the model (e.g. drift or selection), users initialize the model with
historical data, parameterize the model based on their understanding of
the mechanisms, and test if subsequently observed or predicted data
match the simulated data. The output data list includes two measures of
accuracy: the proportion of observed points that fall within the
confidence intervals of the simulated data and the residuals between
simulated and observed population trait means. Precision is measured as
the residuals between simulated and observed population trait variances.
We tested the hypothesis that songs of *Z. l. nuttalli* in Bear Valley,
CA evolved through cultural drift from 1969 to 2005.

Prepare initial song data for Bear
Valley.

``` r
starting.trait <- subset(song.data, Population=="Bear Valley" & Year==1969)$Trill.FBW
starting.trait2 <- c(starting.trait, rnorm(n.territories-length(starting.trait), 
    mean=mean(starting.trait), sd=sd(starting.trait)))

init.inds <- data.frame(id = seq(1:n.territories), age = 2, trait = starting.trait2)
init.inds$x1 <-  round(runif(n.territories, min=-122.481858, max=-122.447270), digits=8)
init.inds$y1 <-  round(runif(n.territories, min=37.787768, max=37.805645), digits=8)
```

Specify and call SongEvo() with test
data

``` r
SongEvo3 <- SongEvo(init.inds = init.inds, females = 1.0, iteration = iteration, steps = years,  
    timestep = timestep, n.territories = n.territories, terr.turnover = terr.turnover, 
    integrate.dist = integrate.dist, 
    learning.error.d = learning.error.d, learning.error.sd = learning.error.sd, 
    mortality.a.m = mortality.a.m, mortality.a.f = mortality.a.f,
    mortality.j.m = mortality.j.m, mortality.j.f = mortality.j.f, lifespan = lifespan, 
    phys.lim.min = phys.lim.min, phys.lim.max = phys.lim.max, 
    male.fledge.n.mean = male.fledge.n.mean, male.fledge.n.sd = male.fledge.n.sd, male.fledge.n = male.fledge.n,
    disp.age = disp.age, disp.distance.mean = disp.distance.mean, disp.distance.sd = disp.distance.sd, 
    mate.comp = mate.comp, prin = prin, all = TRUE)
```

Specify and call
`h.test()`

``` r
target.data <- subset(song.data, Population=="Bear Valley" & Year==2005)$Trill.FBW
h.test1 <- h.test(summary.results=SongEvo3$summary.results, ts=ts, 
    target.data=target.data)
```

The output data list includes two measures of accuracy: the proportion
of observed points that fall within the confidence intervals of the
simulated data and the residuals between simulated and observed
population trait means. Precision is measured as the residuals between
simulated and observed population trait variances.

Eighty percent of the observed data fell within the central 95% of the
simulated values, providing support for the hypothesis that cultural
drift as described in this model is sufficient to describe the evolution
of trill frequency bandwidth in this population.

``` r
h.test1
#> $Residuals
#>              Residuals of mean Residuals of variance
#> Iteration 1          702.12246            25858.5280
#> Iteration 2          571.31453           137168.0396
#> Iteration 3          143.36055            10476.1041
#> Iteration 4          215.02870            35821.5017
#> Iteration 5          578.96351            13080.5858
#> Iteration 6           22.45852             7269.0398
#> Iteration 7          356.71645              821.7307
#> Iteration 8          855.96356            19854.4851
#> Iteration 9          719.00214            20768.7905
#> Iteration 10          84.65474            27935.5851
#> 
#> $Prop.contained
#> [1] 0.4
```

We can plot simulated data in relation to measured data.

``` r
#Plot
plot(SongEvo3$summary.results[1, , "trait.pop.mean"], xlab="Year", ylab="Bandwidth (Hz)", 
    xaxt="n", type="n", xlim=c(-0.5, 35.5), 
    ylim=c(min(SongEvo3$summary.results[, , "trait.pop.mean"], na.rm=TRUE), 
    max(SongEvo3$summary.results[, , "trait.pop.mean"], na.rm=TRUE)))
    for(p in 1:iteration){
        lines(SongEvo3$summary.results[p, , "trait.pop.mean"], col="light gray")
        }
freq.mean <- apply(SongEvo3$summary.results[, , "trait.pop.mean"], 2, mean, na.rm=TRUE)
lines(freq.mean, col="blue")
axis(side=1, at=seq(0, 35, by=5), labels=seq(1970, 2005, by=5))#, tcl=-0.25, mgp=c(2,0.5,0))

#Plot 95% quantiles (which are similar to credible intervals)
quant.means <- apply (SongEvo3$summary.results[, , "trait.pop.mean"], MARGIN=2, quantile, 
    probs=c(0.95, 0.05), R=600, na.rm=TRUE)
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
```

![](SongEvo_files/figure-gfm/Plot%20song%20simulation%20data-1.png)<!-- -->
