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
                    timetep = 1, # misspelled?
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

# Alternative is to load parameters into list and attach to result:
SongEvo1_alt=list(parms=glo.parms)
parms_1alt <- list(init.inds = init.inds, 
                   iteration = 10,
                   steps = 36,  # years / timestep
                   timetep = 1, # misspelled?
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
