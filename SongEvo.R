
library("boot")
sample.mean <- function(d, x) {
  mean(d[x])
}

library("geosphere")
library("sp")

SongEvo <- function(init.inds, 
                    iteration, 
                    steps,
                    timestep,
                    terr.turnover,
                    mate.comp,
                    learning.method,
                    integrate.dist,
                    learning.error.d,
                    learning.error.sd,
                    mortality.a,
                    mortality.j,
                    lifespan,
                    phys.lim.min,
                    phys.lim.max,
                    male.fledge.n.mean,
                    male.fledge.n.sd,
                    male.fledge.n,
                    disp.age,
                    disp.distance.mean,
                    disp.distance.sd,
                    n.territories,
                    prin,
                    all){

ptm <- proc.time()
# m_pnt(2)
summary.results <- array(NA, 
                         dim=c(iteration, steps, 5), 
                         dimnames = list(iteration=paste("iteration", seq(1:iteration)), 
                                         step=1:steps, 
                                         feature=c("sample.n", "trait.pop.mean", "trait.pop.variance", "lci", "uci"))) 
trait.results <- NULL

#Further prepare initial individuals 
init.inds$male.fledglings <- male.fledge.n[1:n.territories]
init.inds$territory <- rep(1, n.territories)
init.inds$father <- 0
init.inds$x0 <- 0
init.inds$y0 <- 0
init.inds$x <- init.inds$x1
init.inds$y <- init.inds$y1
coordinates(init.inds) = ~x+y 
proj4string(init.inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 

all.inds <- inds0 <- init.inds

#Add timestep and iteration
all.inds$timestep <- 0
all.inds$iteration <- 0
all.inds0 <- all.inds[which(is.na(all.inds$trait)),] #This removes all data from all.inds.  
all.inds<-NULL
step_list=list()
length(step_list)<-steps
names(step_list)<-paste0("T",1:steps)
inds.all_list<-lapply(1:iteration,function(x) step_list)
names(inds.all_list)<-paste0("I",1:iteration)

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

hatch <- function(inds){
	ninds <- nrow(inds)
	n.hatch<- sum(inds$male.fledglings)
	if(n.hatch==0) return(inds)
	newinds <- new.bird[rep(1,n.hatch),]
	newinds$id <- maxid +(1:n.hatch)
	row.names(newinds)<-as.character(newinds$id)
	map.key<-do.call(c,mapply(rep,1:ninds,inds$male.fledglings, SIMPLIFY = FALSE))
	newinds$father <- inds$id[map.key]
	newinds$x <- newinds$x0 <- inds$x1[map.key] 
	newinds$y <- newinds$y0 <- inds$y1[map.key]
	# coordinates(newinds) = ~x+y 
	# proj4string(newinds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
	# inds$male.fledglings <- 0
	# rbind(inds, newinds)
	newinds
	}

if (learning.method=="father") {
  #1. Young learn from their fathers:
  learn <- function(children,tutors){
    
    # child=which(inds$age==1)
    children$trait=sapply(children$father, function(x) tutors[tutors$id==x, ]$trait ) + 
      rnorm(nrow(children), mean=learning.error.d, sd=learning.error.sd) 
    #restrict learned song values such that they cannot exceed range of physical possibility: 
    children$trait[children$trait < phys.lim.min] <- phys.lim.min
    children$trait[children$trait > phys.lim.max] <- phys.lim.max
    return(children)}
} else if (learning.method=="integrate") {
  #2. Young learn by integrating songs from the neighborhood within a specified distance:
    learn <- function(children,tutors){
      
      map.key<-do.call(c,mapply(rep,1:nrow(tutors),tutors$male.fledglings, SIMPLIFY = FALSE))
      # stopifnot(children$father==tutors$id[map.key])
      # child=which(inds$age==1)
      # singing.inds <- subset(inds, age>1)
      key=spDists(tutors,tutors,longlat = TRUE)[map.key,] <= integrate.dist
      children$trait=
        key%*%tutors$trait /rowSums(key) +
        # sapply(1:nrow(children), function(x) {
        # key=spDistsN1(pts=singing.inds, pt=inds[x,], longlat=TRUE) <= integrate.dist
        # mean(tutors[key[,x], ]$trait) }) + 
        rnorm(nrow(children), mean=learning.error.d, sd=learning.error.sd) 
      #restrict learned song values such that they cannot exceed range of physical possibility: 
      children$trait[children$trait < phys.lim.min] <- phys.lim.min
      children$trait[children$trait > phys.lim.max] <- phys.lim.max
      return(children)
    } }

if (is.na(lifespan)) {
  die <- function(inds) {
    ninds<-nrow(inds)
    subset(inds,  runif(ninds) > ifelse(age > 1,
                                        rep(mortality.a,ninds),
                                        rep(mortality.j,ninds)))
  }
} else {
  die <- function(inds) {
    ninds<-nrow(inds)
    subset(inds,  age <= lifespan & 
             runif(ninds) > ifelse(age > 1,
                                   rep(mortality.a,ninds),
                                   rep(mortality.j,ninds)))
  }
}
grow <- function(inds){
	inds$age <- inds$age + timestep
inds <- inds
}

disperse <- function(inds){
	ninds <- sum(key<-(inds$age == disp.age))
	# inds <- as.data.frame(inds)
	# for (i in 1: ninds) {
	# 	if (inds$age[i] == disp.age){
		disp.dist <- rnorm(ninds, mean=disp.distance.mean, sd=disp.distance.sd)  
		disp.direction <- runif(ninds, min=0, max=360) #Assume 0 is due North.
		inds[key, c("x", "y")] <- 
		  inds[key, c("x1", "y1")] <- 
		  destPoint(inds[key, c("x0", "y0")], disp.direction, disp.dist)
  #   }
  # }
	coordinates(inds) = ~x+y 
	proj4string(inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
inds 
}

compete.for.territories <- function(inds){
	ninds <- length(inds$age)
	#Allows a flexible population size, but caps breeding territories at a specific number. 
	without <- which(inds$territory==0)
	with <- which(inds$territory==1)
	nwithout <- length(without)
	nwith <- length(with)
	nopen <- n.territories-nwith

	#Bring available territories up to number based on specified turnover rate.  Some have already opened because of death in the Die submodel.
	if (terr.turnover*n.territories > nopen) {		
		ran <- sample(with, size=round(terr.turnover*n.territories - nopen), replace=FALSE) 
		inds[c(ran), "territory"] <- 0		
	} 

	with <- which(inds$territory==1)
	nwith <- length(with)
	nopen <- n.territories-nwith

 #Give open territories to birds that previously lacked territories. There are two scenarios:
 #First scenario: There are more open territories than specified by the turnover rate because of mortality in the die submodel:
	# if (nwith < n.territories-terr.turnover*n.territories) {		
		if (nwithout > nopen) { #open spaces	
		inds[sample(without, size=nopen, replace=FALSE), "territory"] <- 1
		} else { #open spaces	
		#ran3 <- without #all get territories!			
		inds[without, "territory"] <- 1
		} 

	# }

#  #Second scenario: The exact turnover is available:
# 	if (nwith == n.territories-terr.turnover*n.territories) {		
# 		if (nwithout > n.territories-nwith) { 
# 		ran4 <- sample(without, size=n.territories-nwith, replace=FALSE)
# 		inds[c(ran4), "territory"] <- 1
# 		} 			
# 		else if (nwithout <= n.territories-nwith) { 	
# 		ran5 <- without #all get territories!			
# 		inds[c(ran5), "territory"] <- 1
# 		} 
# 	}
inds 
}

compete.for.mates <- function(inds){
	ninds <- length(inds$age)
	prev.songs <- subset(inds, age > 1)$trait
	# avg.song <- mean(prev.songs)
	comp.res<- (!mate.comp | abs(inds$trait-mean(prev.songs)) < (2*sd(prev.songs)))
	inds$male.fledglings <- (inds$territory==1) * comp.res * round(rnorm(ninds, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
	# for (i in 1: ninds) {
	#  if (mate.comp){
	# 	if (abs(inds$trait[i]-avg.song) < (2*sd(prev.songs))){ 
	# 		if (inds$territory[i]==1){
	# 		inds$male.fledglings[i] <- round(rnorm(1, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
	# 	 	inds$male.fledglings[inds$male.fledglings < 0] <- 0						 }
	# 	 							 }
	# 	  } 
	#  if (!mate.comp){
	# 	if (inds$territory[i]==1){
	# 		inds$male.fledglings[i] <- round(rnorm(1, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
	# 		inds$male.fledglings[inds$male.fledglings < 0] <- 0
	# 		}
	# 	  }
	# }
	inds$male.fledglings[inds$male.fledglings < 0] <- 0	
inds 
}

reproduce <- NULL #allows mated individuals to have differential reproductive success, depending on chance, environmental variables, etc. Not currently implemented

#===============================================
# (3) life loop
#===============================================

for (b in 1:iteration){

inds <- init.inds
maxid <- max(inds$id) #store max id, which determines new id numbers in Hatch

for (k in 1:steps) {
  # m_pnt(3)
  
if (prin==TRUE){
	print(paste("iteration ", b, ", timestep ", k, ", n males ", NROW(inds), ", n territorial males ", length(which(inds$territory==1)), sep=""))}

	if (NROW(inds) >= 3){

	#Data frames used for recording individual values at the start of each time step:
timestep.inds <- inds
timestep.inds$timestep <- k
timestep.inds$iteration <- b
# all.inds <- rbind(all.inds, timestep.inds)
	inds.all_list[[b]][[k]]<-timestep.inds
chicks <- hatch(inds)
# m_pnt(4)

maxid <- max(chicks$id)

chicks <- learn(chicks,inds)
# m_pnt(5)
inds <- die(inds)
chicks <- die(chicks)
# m_pnt(6)

if (NROW(inds)+NROW(chicks) >= 3){

inds <- grow(inds)
chicks <- grow(chicks)
# m_pnt(7)
chicks<-disperse(chicks)
inds <- rbind(inds,chicks)
# m_pnt(8)
inds <- compete.for.territories(inds)
# m_pnt(9)
inds <- compete.for.mates(inds)
# m_pnt(10)

#Calculate summary values
summary.results[b, , "sample.n"][k] <- length(inds$age)
summary.results[b, , "trait.pop.mean"][k] <- mean(subset(inds, age==2)$trait)
summary.results[b, , "trait.pop.variance"][k] <- var(subset(inds, age==2)$trait)

boot_obj <- boot(inds$trait, statistic=sample.mean, R=100)#, strata=mn.res$iteration)	
ci.res <- boot.ci(boot_obj, conf=0.95, type="basic")
summary.results[b, , "lci"][k] <- ci.res$basic[4]
summary.results[b, , "uci"][k] <- ci.res$basic[5]
	}
	}
	}
	} 
if (all==TRUE){
  all.inds=rbind(init=all.inds0,  do.call(rbind,do.call(c,inds.all_list)))
z <- list("summary.results"=summary.results, "inds"=inds, "all.inds"=all.inds, "time"=proc.time()-ptm)}
if (all!=TRUE){
z <- list("summary.results"=summary.results, "inds"=inds, "time"=proc.time()-ptm)}
# m_pnt(11)

z
}
#End of I. SongEvo function
#########################