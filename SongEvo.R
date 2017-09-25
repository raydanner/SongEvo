SongEvo <- function(init.inds, iteration, steps, timetep, terr.turnover, mate.comp, learning.method, integrate.dist, learning.error.d, learning.error.sd, mortality.a, mortality.j, lifespan, phys.lim.min, phys.lim.max, male.fledge.n.mean, male.fledge.n.sd, male.fledge.n, disp.age, disp.distance.mean, disp.distance.sd, n.territories, prin, all) {

ptm <- proc.time()

#===============================================
# (1) Create data structures
#===============================================
#steps <- 5

summary.results <- array(NA, dim=c(iteration, steps, 5), dimnames = list(paste("iteration", seq(1:iteration)), 1:steps, c("sample.n", "trait.pop.mean", "trait.pop.variance", "lci", "uci"))) 
trait.results <- NULL

#Further prepare initial individuals 
init.inds$male.fledglings <- male.fledge.n[1:n.territories]
init.inds$territory <- rep(1, n.territories)
init.inds$father <- 0
init.inds$x0 <- 0
init.inds$y0 <- 0
init.inds$x <- init.inds$x1
init.inds$y <- init.inds$y1
library("sp")
coordinates(init.inds) = ~x+y 
proj4string(init.inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 

all.inds <- inds0 <- init.inds

#Add timestep and iteration
all.inds$timestep <- 0
all.inds$iteration <- 0
all.inds <- all.inds[which(is.na(all.inds$trait)),] #This removes all data from all.inds.  
maxid <- 0

#Template for one individual
new.bird <- data.frame(id = 0, age = 1, 
trait = 0,
male.fledglings = 0,
territory = 0,
father=0,
x = 0, #current Longitude
y = 0, #current Latitude
x0 = 0, #starting Longitude
y0 = 0, #starting Latitude
x1 = 0, #ending Longitude
y1 = 0) #ending Latitude

#===============================================
# (2) life methods of the individuals
#===============================================

hatch <- function(inds){
	ninds <- length(inds$age)
	newinds <- NULL
	for (i in 1: ninds) {
		if (inds$male.fledglings[i] > 0){
			for (j in 1:inds$male.fledglings[i]) {
				new.bird2 <- new.bird
				new.bird2$id <- maxid + NROW(newinds) + 1
				new.bird2$father <- inds$id[i]
				new.bird2$x <- new.bird2$x0 <- inds$x1[i] 
				new.bird2$y <- new.bird2$y0 <- inds$y1[i]
				newinds <- rbind(newinds, new.bird2)
				}	
		inds$male.fledglings[i] <- 0 #This will be repopulated later in life (in Compete for mates, and potentially 			Reproduce).  
			}
		}
	library("sp")
	coordinates(newinds) = ~x+y 
	proj4string(newinds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
	inds <- rbind(inds, newinds)
	}

learn <- function(inds){
	ninds <- length(inds$age)

	#1. Young learn from their fathers:
	if (learning.method=="father") {
	for (i in 1: ninds) {
		if (inds$age[i]==1) {
		tutor <- inds$father[i]
		inds$trait[i] <- inds[inds$id==tutor, ]$trait + rnorm(1, mean=learning.error.d, sd=learning.error.sd)  
		}
	}
	}

	#2. Young learn by integrating songs from the neighborhood within a specified distance:
	if (learning.method=="integrate") {
	for (i in 1: ninds) {
		if (inds$age[i]==1) {
			singing.inds <- subset(inds, age>1)
			tutors <- which(spDistsN1(pts=singing.inds, pt=inds[i,], longlat=TRUE) <= integrate.dist) #distance in km
			inds$trait[i] <- mean(inds$trait[tutors]) + rnorm(1, mean=learning.error.d, sd=learning.error.sd)  			
		}
	}
	}

	#restrict learned song values such that they cannot exceed range of physical possibility: 
	inds$trait[inds$trait < phys.lim.min] <- phys.lim.min
	inds$trait[inds$trait > phys.lim.max] <- phys.lim.max
	inds <- inds #Not sure why this is needed.  But, if it or print(inds) is not included, inds=NULL.
}

die <- function(inds) {
	j <- subset(inds, age == 1)
	j2 <- subset(j, runif(j$age) > mortality.j)
	a <- subset(inds, age > 1)
	if (!is.na(lifespan)) {
	a2 <- subset(a, runif(a$age) > mortality.a & a$age <= lifespan)
		} 
	if (is.na(lifespan)) {
	a2 <- subset(a, runif(a$age) > mortality.a)
	}
	if ((NROW(a2) > 0) | (NROW(j2) > 0)) { #if/else statement to 	control for special case when no adults and no juveniles survive.
	inds <- rbind(a2, j2)
	} else 
	(inds <- NULL)
}

grow <- function(inds){
	inds$age <- inds$age + timestep
inds <- inds
}

disperse <- function(inds){
	ninds <- length(inds$age)
	inds <- as.data.frame(inds)
	for (i in 1: ninds) {
		if (inds$age[i] == disp.age){
		disp.dist <- rnorm(1, mean=disp.distance.mean, sd=disp.distance.sd)  
		disp.direction <- runif(1, min=0, max=360) #Assume 0 is due North.
		library("geosphere")
		inds[i, c("x", "y")] <- inds[i, c("x1", "y1")] <- destPoint(inds[i, c("x0", "y0")], disp.direction, disp.dist)
    }
  }
	coordinates(inds) = ~x+y 
	proj4string(inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
inds <- inds
}

compete.for.territories <- function(inds){
	ninds <- length(inds$age)
	#Allows a flexible population size, but caps breeding territories at a specific number. 
	without <- which(inds$territory==0)
	with <- which(inds$territory==1)
	nwithout <- length(without)
	nwith <- length(with)
	nopen <- length(n.territories-nwith)

	#Bring available territories up to number based on specified turnover rate.  Some have already opened because of death in the Die submodel.
	if (nwith > n.territories-terr.turnover*n.territories) {		
		ran <- sample(with, size=round(nwith-(n.territories-terr.turnover*n.territories)), replace=FALSE) 
		inds[c(ran), "territory"] <- 0		
	} else

	with <- which(inds$territory==1)
	nwith <- length(with)

 #Give open territories to birds that previously lacked territories. There are two scenarios:
 #First scenario: There are more open territories than specified by the turnover rate because of mortality in the die submodel:
	if (nwith < n.territories-terr.turnover*n.territories) {		
		if (nwithout > 40-nwith) { #open spaces	
		ran2 <- sample(without, size=n.territories-nwith, replace=FALSE)
		inds[c(ran2), "territory"] <- 1
		} 			
		else if (nwithout <= n.territories-nwith) { #open spaces	
		ran3 <- without #all get territories!			
		inds[c(ran3), "territory"] <- 1
		} 

	}

 #Second scenario: The exact turnover is available:
	if (nwith == n.territories-terr.turnover*n.territories) {		
		if (nwithout > n.territories-nwith) { 
		ran4 <- sample(without, size=n.territories-nwith, replace=FALSE)
		inds[c(ran4), "territory"] <- 1
		} 			
		else if (nwithout <= n.territories-nwith) { 	
		ran5 <- without #all get territories!			
		inds[c(ran5), "territory"] <- 1
		} 
	}
inds <- inds
}

compete.for.mates <- function(inds){
	ninds <- length(inds$age)
	prev.songs <- subset(inds, age > 1)$trait
	avg.song <- mean(prev.songs)
	for (i in 1: ninds) {
	 if (mate.comp){
		if (abs(inds$trait[i]-avg.song) < (2*sd(prev.songs))){ 
			if (inds$territory[i]==1){
			inds$male.fledglings[i] <- round(rnorm(1, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
		 	inds$male.fledglings[inds$male.fledglings < 0] <- 0						 }
		 							 }
		  } 
	 if (!mate.comp){
		if (inds$territory[i]==1){
			inds$male.fledglings[i] <- round(rnorm(1, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
			inds$male.fledglings[inds$male.fledglings < 0] <- 0
			}
		  }
	}
inds <- inds
}

reproduce <- NULL #allows mated individuals to have differential reproductive success, depending on chance, environmental variables, etc. Not currently implemented

#===============================================
# (3) life loop
#===============================================

for (b in 1:iteration){

inds <- init.inds
maxid <- max(inds$id) #store max id, which determines new id numbers in Hatch

for (k in 1:steps) {

if (prin==TRUE){
	print(paste("iteration ", b, ", timestep ", k, ", n males ", NROW(inds), ", n territorial males ", length(which(inds$territory==1)), sep=""))}

	if (NROW(inds) >= 3){

	#Data frames used for recording individual values at the start of each time step:
timestep.inds <- inds
timestep.inds$timestep <- k
timestep.inds$iteration <- b
all.inds <- rbind(all.inds, timestep.inds)
	
inds <- hatch(inds)

maxid <- max(inds$id)

inds <- learn(inds)

inds <- die(inds)

if (NROW(inds) >= 3){

inds <- grow(inds)

inds <- disperse(inds)

inds <- compete.for.territories(inds)

inds <- compete.for.mates(inds)

#Calculate summary values
summary.results[b, , "sample.n"][k] <- length(inds$age)
summary.results[b, , "trait.pop.mean"][k] <- mean(subset(inds, age==2)$trait)
summary.results[b, , "trait.pop.variance"][k] <- var(subset(inds, age==2)$trait)

library("boot")
sample.mean <- function(d, x) {
	mean(d[x])
	}

boot_obj <- boot(inds$trait, statistic=sample.mean, R=100)#, strata=mn.res$iteration)	
ci.res <- boot.ci(boot_obj, conf=0.95, type="basic")
summary.results[b, , "lci"][k] <- ci.res$basic[4]
summary.results[b, , "uci"][k] <- ci.res$basic[5]
	}
	}
	}
	} 
if (all==TRUE){
z <- list("summary.results"=summary.results, "inds"=inds, "all.inds"=all.inds, "time"=proc.time()-ptm)}
if (all!=TRUE){
z <- list("summary.results"=summary.results, "inds"=inds, "time"=proc.time()-ptm)}
z
}
#End of I. SongEvo function
#########################