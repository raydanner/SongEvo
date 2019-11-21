#' Model bird song evolution
#'
#' This function simulates bird song evolution. Submodels are performed once per time step, and include fledging from the nest, song learning, ageing and death, dispersal, competition for territories, mate attraction, and reproduction.
#'
#' @name SongEvo
#' @param init.inds Initial population data. A data frame that includes columns for “id,” “age,” “trait,” “x1” (longitude) and “y1” (latitude). 
#' @param females Either sex ratio expressed in terms of females to males (e.g. 2 would represent 2 females : 1 male), or female dataframe.
#' @param iteration The number of iterations that the model will run. 
#' @param steps The number of steps per iteration.
#' @param timestep The length of time that passes in each step. For annually breeding species, timestep = 1 year.
#' @param terr.turnover The proportion of territories that change ownership during a step.
#' @param mate.comp Female preference for mates. Options are TRUE, FALSE, or the name of a user defined function. User functions must paramatize male id number and return number of total offspring for that step.
#' @param selectivity The ability of females to choose males according to their preference. Expressed in number of standard deviations from her preferred trait value.
#' @param content.bias Reduces the heritability of individuals with affected traits. Options are FALSE or a vector of lowest and highest fitness reductions, bias strength (e.g. c(.9,.45)).
#' @param n.content.bias.loc Number of content bias locations, options are either an integer, 'random' (1:5 locations), or 'all' where content bias is delocalized.
#' @param content.bias.loc Location centers of content bias effects, options are a user defined dataframe with x and y positions, 'random' (random points on spatial plane), or FALSE.
#' @param content.bias.loc.ranges Radius of content bias effects at each location. Options are a user defined vector, 'random' (0.5:10), or FALSE.
#' @param affected.traits Either a dataframe containing max and min affected traits at each content bias center, 'random' (phys.lim.min:phys.lim.max), or FALSE.
#' @param conformity.bias If an individual learns from their ('father'), all males within a specified radius ('integrate'), or FALSE if no conformity bias.
#' @param integrate.dist Distance over which tutor values are integrated for learning.
#' @param prestige.bias Learners will preferentially learn from males with more offspring. Options are a user defined vector of fitness reduction (e.g. c(.95,.25)), or FALSE.
#' @param learn.m How males aquire a trait. Options are 'default', where males take the weighted average of tutors' traits according to their fitness values, or the name of a user defined function that takes id as a parameter and returns traits.
#' @param learn.f How females aquire a trait. Options are 'default', where females take the weighted average of tutors traits according to their fitness values, or the name of a user defined function that takes id as a parameter and returns traits.
#' @param learning.error.d Direction of learning error (measured in trait units, e.g. Hz).
#' @param learning.error.sd The standard deviation of imitation error.  
#' @param mortality.a.m Annual mortality of adult males.
#' @param mortality.a.f Annual mortality of adult females.
#' @param mortality.j.m Annual mortality of juvenile males.
#' @param mortality.j.f Annual mortality of juvenile females.
#' @param lifespan Maximum age for individuals; any number is accepted. “NA” causes SongEvo to disregard lifespan and sets population size based on mortality rates alone.
#' @param phys.lim.min The minimum physical limit of trait production.
#' @param phys.lim.max The maximum physical limit of trait production.
#' @param male.fledge.n.mean The mean number of offspring produced per time step per individual breeding male. Includes only offspring raised in that breeding male’s nest (i.e. it does not account for extra-pair offspring in other nests).
#' @param male.fledge.n.sd Standard deviation of the number of male fledglings.
#' @param disp.age The age at which individual males disperse from their birth location.
#' @param disp.distance.mean The distance that individual males disperse (meters).
#' @param disp.distance.sd The standard deviation of dispersal distance.
#' @param n.territories The maximum number of potential territories in the population. This number is fixed for all iterations.
#' @param prin Print summary values after each timestep has completed? Options are TRUE or FALSE. 
#' @param all Save data for all individuals? Options are TRUE or FALSE. 
#' 
#' @return three objects. First, currently alive individuals are stored in a data frame called “inds.”  Values within “inds” are updated throughout each of the iterations of the model, and “inds” can be viewed after the model is completed.  Second, an array (i.e. a multi-dimensional table) entitled “summary.results” includes population summary values for each time step (dimension 1) in each iteration (dimension 2) of the model.  Population summary values are contained in five additional dimensions: population size for each time step of each iteration (“sample.n”), the population mean and variance of the song feature studied (“trait.pop.mean” and “trait.pop.variance”), with associated lower (“lci”) and upper (“uci”) confidence intervals.  Third, individual values may optionally be concatenated and saved to one data frame entitled “all.inds.”  all.inds can become quite large, and is therefore only recommended if additional data analyses are desired. 
#' 
#' @example inst/examples/SongEvoExamples.R
#' 
#' @seealso \code{\link{par.sens}}, \code{\link{par.opt}}, \code{\link{mod.val}}, \code{\link{h.test}} 
#' 
#'
#' @importFrom boot boot boot.ci
#' @import lattice
#' @import sp
#' @import geosphere
#' @importFrom stats runif rnorm sd
#' @export
SongEvo <- function(init.inds, 
                        females,
                        iteration, 
                        steps,
                        timestep,
                        terr.turnover,
                        mate.comp,
                        selectivity,
                        content.bias = FALSE, 
                        n.content.bias.loc = 'all', 
                        content.bias.loc = FALSE, 
                        content.bias.loc.ranges = FALSE, 
                        affected.traits = FALSE,
                        conformity.bias = FALSE, 
                        integrate.dist = 0.1, 
                        prestige.bias = FALSE, 
                        learn.m = 'default', 
                        learn.f = 'default', 
                        learning.error.d = 0,
                        learning.error.sd = 200,
                        mortality.a.m = 0.5,
                        mortality.a.f = mortality.a.m,
                        mortality.j.m = 0.5,
                        mortality.j.f = mortality.j.m,
                        lifespan = 2,
                        phys.lim.min = 1000,
                        phys.lim.max = 8000,
                        male.fledge.n.mean = 2,
                        male.fledge.n.sd = 0.5,
                        disp.age = 1,
                        disp.distance.mean = 100,
                        disp.distance.sd = 50,
                        n.territories = 4 * nrow(init.inds),
                        prin = FALSE,
                        all = FALSE){
  
  ptm <- proc.time()
  # m_pnt(2)
  summary.results <- array(NA, 
                           dim=c(iteration, steps, 5), 
                           dimnames = list(iteration=paste("iteration", seq(1:iteration)), 
                                           step=1:steps, 
                                           feature=c("sample.n", "trait.pop.mean", "trait.pop.variance", "lci", "uci"))) 
  trait.results <- NULL
  
  #Further prepare initial individuals 
  init.inds$male.fledglings <- as.integer(rnorm(nrow(init.inds), male.fledge.n.mean, male.fledge.n.sd))
  init.inds$female.fledglings <- as.integer(rnorm(nrow(init.inds), male.fledge.n.mean, male.fledge.n.sd))
  init.inds$territory <- rep(1, n.territories)
  init.inds$father <- 0
  init.inds$sex <- 'M'
  init.inds$fitness <- 1
  init.inds$learn.dir <- 0
  init.inds$x0 <- 0
  init.inds$y0 <- 0
  init.inds$x <- init.inds$x1
  init.inds$y <- init.inds$y1
  # coordinates(init.inds) = ~x+y 
  # proj4string(init.inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
  
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
  
  maxid.m <- 0
  maxid.f <- 0
  
  #create females with the same traits as males or modify an existing dataframe
  create.females <- function(inds, females){
    if (typeof(females)=="double"){
      n.fem <- as.integer(females*nrow(inds))
      f.trait <- rnorm(n.fem, mean = mean(inds$trait), sd = sd(inds$trait))
      init.fem <- data.frame(id = seq(1:n.fem), age = 2, trait = f.trait)
      init.fem$x1 <- round(runif(n.fem, min=-122.481858, max=-122.447270), digits=8)
      init.fem$y1 <- round(runif(n.fem, min=37.787768, max=37.805645), digits=8)
      init.fem$father <- 0
      init.fem$sex <- 'F'
      init.fem$x0 <- 0
      init.fem$y0 <- 0
      init.fem$x <- init.fem$x1
      init.fem$y <- init.fem$y1
      return(init.fem)
    } else if (typeof(females)=="list"){

      females$father <- 0
      females$sex <- 'F'
      females$x0 <- 0
      females$y0 <- 0
      females$x <- females$x1
      females$y <- females$y1
      return(females)
      
    }
    else{
      print("user input error females must be either a ratio (f:m) or a dataframe")
    }
  }
  
  #Template for one individual
  new.bird <- data.frame(id = 0, age = 1, 
                         trait = 0,
                         male.fledglings = 0,
                         female.fledglings = 0,
                         territory = 0,
                         father=0,
                         sex='F',
                         x = 0, #current Longitude
                         y = 0, #current Latitude
                         x0 = 0, #starting Longitude
                         y0 = 0, #starting Latitude
                         x1 = 0, #ending Longitude
                         y1 = 0 #ending Latitude
                         )

  hatch <- function(inds){
    ninds <- nrow(inds)
    n.hatch<- sum(inds$male.fledglings, inds$female.fledglings)
    if(n.hatch==0) {
      newinds <- new.bird[-c(1),]
      print(paste('no hatch row n: ', NROW(newinds)))
      #print('no hatch df')
      #print(newinds)
      return(newinds)
      
    }
    else {
      fathers.m <- subset(inds, male.fledglings != 0)
      id.m <- maxid.m
      newinds <- new.bird
      for (i in 1:nrow(fathers.m)){
        for (n in 1:fathers.m$male.fledglings[i]){
          id.m <- as.numeric(id.m)+1
          father <- fathers.m$id[i]
          x.loc <- fathers.m$x1[i]
          y.loc <- fathers.m$y1[i]
          newinds <- rbind(newinds, c(id.m, 1, 0, 0, 0, 0, father, 'M', x.loc, y.loc, x.loc, y.loc, 0, 0))
          newinds$sex <- 'M'
        }
      }
      fathers.f <- subset(inds, female.fledglings !=0)
      id.f <- maxid.f
      for (i in 1:nrow(fathers.f)){
        for (n in 1:fathers.f$female.fledglings[i]){
          id.f <- as.numeric(id.f)+1
          father <- fathers.f$id[i]
          x.loc <- fathers.f$x1[i]
          y.loc <- fathers.f$y1[i]
          newinds <- rbind(newinds, c(id.f, 1, 0, 0, 0, 0, father, 'F', x.loc, y.loc, x.loc, y.loc, 0, 0))
        }
      }
      newinds <- newinds[c(-1),]
      newinds$fitness <- 1
      newinds$learn.dir <- 0
      return(newinds)
      'newinds <- new.bird[rep(1,n.hatch),]
      newinds$id <- maxid +(1:n.hatch)
      row.names(newinds)<-as.character(newinds$id)
      map.key<-do.call(c,mapply(rep,1:ninds,inds$male.fledglings, SIMPLIFY = FALSE))
      newinds$father <- inds$id[map.key]
      newinds$x <- newinds$x0 <- inds$x1[map.key] 
      newinds$y <- newinds$y0 <- inds$y1[map.key]
      return(newinds)
      # coordinates(newinds) = ~x+y 
      # proj4string(newinds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
      # inds$male.fledglings <- 0
      # rbind(inds, newinds)'
    }
  }
  
  #create content bias location spots and associated affected trait values
  if (typeof(content.bias) == 'double' & n.content.bias.loc != 'all'){
    if (n.content.bias.loc=='random'){
      n.content.bias.loc <- round(runif(1, min=1, max=5),0)
      print(paste('number of content bias locations:', n.content.bias.loc))
    }
    if (content.bias.loc=='random'){
      noise.loc.x <- round(runif(n.content.bias.loc, min=-122.481858, max=-122.447270), digits=8)
      noise.loc.y <- round(runif(n.content.bias.loc, min= 37.787768, max= 37.805645), digits=8)
      
      content.bias.loc <- data.frame(x = noise.loc.x, y = noise.loc.y)
      
      print('content bias locations:')
      print(content.bias.loc)
    }
    if (content.bias.loc.ranges=='random'){
      content.bias.loc.ranges <- runif(n.content.bias.loc, min = .5, max = 5)
      print(paste('content bias ranges:', content.bias.loc.ranges))
    }
    if (affected.traits=='random'){
      max <- runif(n.content.bias.loc, min = phys.lim.min, max = phys.lim.max)
      freq.range <- runif(n.content.bias.loc, min = 100, max = 1500)
      min <- max - freq.range
      traits.affected = data.frame(max=max,min=min)
      print('affected traits:')
      print(traits.affected)
      
    }
  } else if (typeof(content.bias) == 'double' & n.content.bias.loc == 'all'){
    if (affected.traits=='random'){
      max <- runif(1, min = phys.lim.min, max = phys.lim.max)
      freq.range <- runif(1, min = 100, max = 1500)
      min <- max-freq.range
      affected.traits <- c(max, min)
      print(paste('affected traits:', affected.traits))
    }
  }
  
  #create dataframe to access later
  if (n.content.bias.loc != 'all'){
    content.bias.info <- data.frame(id = seq(1,n.content.bias.loc))
    content.bias.info <- cbind(content.bias.info, content.bias.loc)
    content.bias.info$range <- content.bias.loc.ranges
    content.bias.info <- cbind(content.bias.info, traits.affected)
  } else if (n.content.bias.loc == 'all') {
    content.bias.info <- data.frame(max.freq <- max(affected.traits), min.freq <- min(affected.traits))
  }
  
  #individuals in content bias zones with affected traits have their fitness and traits adjusted  
  content <- function(males,noise.loc,ranges,affected.freqs, bias.strength){
    prop.range <- round(seq(0,1,by = .01), digits = 2)
    bias.range <- seq(min(bias.strength),max(bias.strength), length.out = 101)
    #data frame that creates a bigger fitness effect the closer the trait is to the middle of the affected trait range
    fitness.effects <- data.frame(proportion = prop.range, bias = bias.range)
    diffusion.range <- rev(seq(.001,1,length.out = 101))
    #data frame that creates a bigger fitness effect the closer an individual is to the center of the noise location
    sound.effects <- data.frame(proportion = prop.range, diffusion = diffusion.range)
    males$learn.dir <- 0
    #for content biases with discrete locations
    if (n.content.bias.loc != 'all') {
      for (i in 1:nrow(content.bias.info)){
        disturbance <- c(content.bias.info$x[i], content.bias.info$y[i])
        range <- as.numeric(content.bias.info$range[i])
        max.freq <- as.numeric(content.bias.info$max[i])
        min.freq <- as.numeric(content.bias.info$min[i])
        males$disturbance.dist <- spDistsN1(pts = as.matrix(males[,c('x','y')]), pt = disturbance, longlat = TRUE)
        #males within disturbance range and have traits within the affected range are affected
        affected.inds <- subset(males, disturbance.dist <= range & trait < max.freq & trait > min.freq)
        if (nrow(affected.inds)==0){
          #if no individuals are affected do nothing
        } else{
          middle.point <- (min.freq+max.freq)/2
          half.range <- max.freq-middle.point
          for (x in 1:nrow(affected.inds)){
            ind.trait <- as.numeric(affected.inds$trait[x])
            ind.id <- affected.inds$id[x]
            #if affected trait is within 50 of the affected ranges max or min, individuals trait is moved just out of range and their fitness is not affected
            if (ind.trait > (max.freq-50)){
              if (phys.lim.max >= (max.freq+1)){
                males$trait[males$id==ind.id] <- max.freq+1
              }
              
            }
            else if (ind.trait < (min.freq+50)){
              if (phys.lim.min <= (min.freq-1)){
                males$trait[males$id==ind.id] <- min.freq-1
              }
            }
            else{
              #individuals fitness will be affected, individuals who learn from them will have a directional learning error
              direction <- ind.trait-middle.point
              if (direction > 0){
                males$learn.dir[males$id==ind.id] <- males$learn.dir[males$id==ind.id] + 50
                prop <- round(direction/middle.point,digits = 2)
              }
              else if (direction < 0){
                males$learn.dir[males$id==ind.id] <- males$learn.dir[males$id==ind.id] - 50
                prop <- round(abs(direction)/middle.point,digits = 2)
              }
              range.prop <- round(affected.inds$disturbance.dist[x]/range, digits = 2)
              sound.pressure <- sound.effects$diffusion[sound.effects$proportion==range.prop]
              bias.adj <- fitness.effects$bias[fitness.effects$proportion==prop]
              fitness.adj <- bias.adj/sound.pressure
              males$fitness[males$id==ind.id] <- males$fitness[males$id==ind.id] * fitness.adj
            }
          }
        }
        
      }
      males <- subset(males,select = -c(disturbance.dist))
    } else if (n.content.bias.loc=='all'){
      #all males with traits within the afected range will be affected in the same way
      max.freq <- max(affected.freqs)
      min.freq <- min(affected.freqs)
      middle.point <- mean(affected.freqs)
      half.range <- max.freq-middle.point
      affected.inds <- affected.inds <- subset(males, trait < max.freq & trait > min.freq)
      if (nrow(affected.inds)==0){
        
      } else{
        for (x in 1:nrow(affected.inds)){
          ind.trait <- as.numeric(affected.inds$trait[x])
          ind.id <- as.numeric(affected.inds$id[x])
          if (ind.trait > (max.freq-50)){
            if (phys.lim.max >= (max.freq+1)){
              males$trait[males$id==ind.id] <- max.freq+1
              
            }
            
          }
          else if (ind.trait < (min.freq+50)){
            if (phys.lim.min <= (min.freq-1)){
              males$trait[males$id==ind.id] <- min.freq-1
            }
          }
          else{
            direction <- ind.trait-middle.point
            if (direction > 0){
              males$learn.dir[males$id==ind.id] <- males$learn.dir[males$id==ind.id] + 50
              prop <- round(direction/middle.point,digits = 2)
            }
            else if (direction < 0){
              males$learn.dir[males$id==ind.id] <- males$learn.dir[males$id==ind.id] - 50
              prop <- round(abs(direction)/middle.point,digits = 2)
            }
            bias.adj <- fitness.effects$bias[fitness.effects$proportion==prop]
            males$fitness[males$id==ind.id] <- males$fitness[males$id==ind.id] * bias.adj
          }
        }
      }
    }
    return(males)
  }
  
  
  conformity <- function(children, tutors, conformity.bias){
    if (conformity.bias =='father'){
      #print(paste('n chicks: ', NROW(children)))
      #print(paste('father: ', children$father, 'trait: ',tutors[tutors$id==children$father, ]$trait ))
      #print('fathers')
      #print(subset(tutors, male.fledglings != 0)[,1:8])
      #print('children')
      #print(children)
      # child=which(inds$age==1)
      fathers <- list()
      for (i in 1:nrow(children))
        fathers[[i]] <- tutors$id[tutors$id==children$father[i]]
      return(fathers)
      'children$trait=sapply(children$father, function(x) tutors[tutors$id==x, ]$trait ) + 
          rnorm(nrow(children), mean=learning.error.d, sd=learning.error.sd) 
        #restrict learned song values such that they cannot exceed range of physical possibility: 
        children$trait[children$trait < phys.lim.min] <- phys.lim.min
        children$trait[children$trait > phys.lim.max] <- phys.lim.max
        return(children)'
    } else if (conformity.bias=="integrate") {
      #2. Young learn by integrating songs from the neighborhood within a specified distance:
      if (children$sex=='M'){
        map.key<-do.call(c,mapply(rep,1:nrow(tutors),tutors$male.fledglings, SIMPLIFY = FALSE))
      } else if (children$sex=='F'){
        #print(paste('m chicks: ', sum(tutors$male.fledglings[tutors$male.fledglings != 0]), ' m key: ', length(map.key.m)))
        map.key<-do.call(c,mapply(rep,1:nrow(tutors),tutors$female.fledglings, SIMPLIFY = FALSE))
      }
      #print(paste('f chicks: ', sum(tutors$female.fledglings[tutors$female.fledglings != 0]), ' f key: ', length(map.key.f)))
      #map.key <- c(map.key.m,map.key.f)
      #print(length(map.key))
      # stopifnot(children$father==tutors$id[map.key])
      # child=which(inds$age==1)
      # singing.inds <- subset(inds, inds$age>1)
      key=spDists(as.matrix(tutors[,c("x","y")]),longlat = TRUE)[map.key,] <= integrate.dist
      rownames(key) <- children$id
      colnames(key) <- tutors$id
      tutors.df <- as.data.frame(key)
      tutors.in.range <- list()
      for (i in 1:nrow(children)){
        child.id <- i
        tutors.near <- as.vector(colnames(tutors.df[i,])[tutors.df[i,]==TRUE])
        tutors.in.range[[child.id]] <- tutors.near
      }
      return(tutors.in.range)
    }}
  
  #tutors within range are ranked by how many children they had in the last step. More successful individuals have a higher fitness than those with few or no offspring
  prestige <- function(tutors, bias.strength){
    tutors$children <- rowSums(tutors[,c('male.fledglings','female.fledglings')])
    max.n.children <- max(tutors$children)
    child.range <- seq(0, max.n.children, by = 1)
    bias.range <- seq(min(bias.strength),max(bias.strength), length.out = length(child.range))
    fitness.effects <- data.frame(children = child.range, bias = bias.range)
    for (i in 1:nrow(tutors)){
      fitness.adj <- fitness.effects$bias[fitness.effects$children == tutors[i,]$children]
      tutors$fitness[i] <- tutors$fitness[i] * fitness.adj
    }
    tutors <- subset(tutors, select = -c(children))
    return(tutors)
  }
  
  #conformity.bias <- conformity.bias
  
  learn <- function(tutors, children, conformity.bias, integrate.dist, prestige.bias.strength){
    if (typeof(conformity.bias)=='character'){ #with conformity bias
      tutors.close <- conformity(children = children, tutors = tutors, conformity.bias = conformity.bias)
      #print(tutors)
      #print(tutors[[3]])
      chick.num <- min(as.numeric(children$id))
      for (chick in tutors.close){
        tutors.id <- chick
        tutors.df <- tutors[1,]
        for (i in tutors.id){
          tutor.add <- tutors[tutors$id==i,]
          tutors.df <- rbind(tutors.df,tutor.add)

        }
        tutors.df <- tutors.df[c(-1),]
        if (typeof(prestige.bias.strength)=='double'){ #with prestige bias
          tutors.df <- prestige(tutors = tutors.df, bias.strength = prestige.bias.strength)

        }
        total.fitness <- sum(tutors.df$fitness)
        trait.wighted.avg <- 0
        weighted.learn.dir <- 0
        trait <- as.numeric(tutors.df$trait)
        fitness <- as.numeric(tutors.df$fitness)
        total.fitness <- sum(fitness)
        fitness.ratio <- fitness/total.fitness
        trait.wighted.avg <- trait*fitness.ratio
        chick.trait <- sum(trait.wighted.avg)
        learn.dir <- as.numeric(tutors.df$learn.dir)

        adj.learn.dir <- learn.dir*fitness.ratio
        chick.learn.dir <- sum(adj.learn.dir)
        new.learn.dir <- learning.error.d + chick.learn.dir
        real.trait <- chick.trait + rnorm(1, mean = new.learn.dir, sd = learning.error.sd)
        children$trait[children$id==chick.num] <- real.trait
        
        chick.num <- chick.num+1
      }
      children$trait[children$trait < phys.lim.min] <- phys.lim.min
      children$trait[children$trait > phys.lim.max] <- phys.lim.max
      return(children)
    } else { #without conformity bias
      tutors.df <- tutors
      if (typeof(prestige.bias.strength)=='double'){ #with prestige bias
        tutors.df <- prestige(tutors = tutors.df, bias.strength = prestige.bias.strength)
      }
      total.fitness <- sum(tutors.df$fitness)
      trait.wighted.avg <- 0
      weighted.learn.dir <- 0
      trait <- tutors.df$trait
      fitness <- tutors.df$fitness
      total.fitness <- sum(fitness)
      fitness.ratio <- fitness/total.fitness
      trait.wighted.avg <- trait*fitness.ratio
      chick.trait <- sum(trait.wighted.avg)
      learn.dir <- tutors.df$learn.dir
      adj.learn.dir <- learn.dir*fitness.ratio
      chick.learn.dir <- sum(adj.learn.dir)
      new.learn.dir <- learning.error.d + chick.learn.dir
      
      children$trait <- chick.trait + rnorm(1, mean = new.learn.dir, sd = learning.error.sd)
      
      chick.num <- chick.num+1
    }
    children$trait[children$trait < phys.lim.min] <- phys.lim.min
    children$trait[children$trait > phys.lim.max] <- phys.lim.max
    if (children$sex=='M'){
      id.start <- maxid.m+1
      id.end <- maxid.m+nrow(children)
      ids <- seq(id.start,id.end, by = 1)
      children$id <- ids
    } else if (children$sex=='F'){
      id.start <- maxid.f+1
      id.end <- maxid.f+nrow(children)
      ids <- seq(id.start,id.end, by = 1)
      children$id <- ids
    }
    return(children)
  }
  
  if (is.na(lifespan)) {
    die <- function(inds) {
      ninds<-nrow(inds)
      subset(inds,  runif(ninds) > ifelse(inds$age > 1 & inds$sex=='M', rep(mortality.a.m,ninds),
                                          ifelse(inds$age > 1 & inds$sex=='F', rep(mortality.a.f,ninds),
                                                 ifelse(inds$age == 1 & inds$sex=='M', rep(mortality.j.m,ninds),
                                                        rep(mortality.j.f,ninds)))))
    }
  } else {
    die <- function(inds) {
      ninds<-nrow(inds)
      subset(inds,  inds$age <= lifespan & 
               runif(ninds) > ifelse(inds$age > 1 & inds$sex=='M', rep(mortality.a.m,ninds),
                                     ifelse(inds$age > 1 & inds$sex=='F', rep(mortality.a.f,ninds),
                                            ifelse(inds$age == 1 & inds$sex=='M', rep(mortality.j.m,ninds),
                                                   rep(mortality.j.f,ninds)))))
    }
  }
  grow <- function(inds){
    inds$age <- as.integer(inds$age) + timestep
    inds <- inds
  }
  
  disperse <- function(inds){
    ninds <- nrow(inds)
    # inds <- as.data.frame(inds)
    # for (i in 1: ninds) {
    # 	if (inds$age[i] == disp.age){
    disp.dist <- rnorm(ninds, mean=disp.distance.mean, sd=disp.distance.sd) 
    disp.direction <- runif(ninds, min=0, max=360) #Assume 0 is due North.
    coords <- inds[, c("x", "y")]
    coords <- sapply(coords, as.numeric)
    inds[, c("x", "y")] <- 
      inds[, c("x1", "y1")] <- 
      destPoint(coords, disp.direction, disp.dist)
    #   }
    # }
    # coordinates(inds) = ~x+y 
    # proj4string(inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
    inds 
  }
  
  compete.for.territories <- function(inds){
    # (assuming turnover rate of 0.5 and 40 available territories)
    # 1. Set A may contain at most 20 birds (N_T*(1-t_rate))
    # 2. If set (A U B) has less than 20 birds then set B has no birds
    # 3. Set (A U C) may contain at most 40 birds (N_t)
    # 5. Sets A, B, C and D are mutually exclusive (ie every bird is in exactly on of these four sets)
    # 6. Set C can contain at most 20 birds 
    # 7. Set C can only contain first year birds
    
    ninds <- length(inds$age)
    #Allows a flexible population size, but caps breeding territories at a specific number. 
    without <- which(inds$territory==0 & (inds$age >= disp.age) &  (inds$age < disp.age+timestep))
    with <- which(inds$territory==1)
    nwithout <- length(without)
    nwith <- length(with)
    nopen <- n.territories-nwith
    
    #Bring available territories up to number based on specified turnover rate.  Some have already opened because of death in the Die submodel.
    if (terr.turnover*n.territories > nopen) {		
      ran <- sample(with, size=round(terr.turnover*n.territories - nopen), replace=FALSE) 
      inds[c(ran), "territory"] <- 0	
      with <- which(inds$territory==1)
      nwith <- length(with)
      nopen <- n.territories-nwith	
    } 
    nopen=min(c(nopen,round(terr.turnover*n.territories)))
    
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
  
  compete.for.mates <- function(inds, fems, selectivity){
    inds$male.fledglings = 0
    inds$female.fledglings = 0
    if (mate.comp) {
      for (i in 1:nrow(fems)){
        fem.loc.x <- fems$x[i]
        fem.loc.y <- fems$y[i]
        fem.loc <- c(fem.loc.x, fem.loc.y)
        #males who can mate with females, those with territories and are near her
        competitors <- subset(inds, territory==1)
        competitors$dist <- spDistsN1(pts = as.matrix(competitors[,c("x","y")]), pt = fem.loc, longlat = TRUE)
        competitors.near <- subset(competitors, dist <= integrate.dist*10)
        if (nrow(competitors.near)==0){ #if no potential mates female will disperse up to 10 times to find a mate
          ndisp <- 0
          while (nrow(competitors.near)==0 & ndisp < 10) {
            move <- disperse(fems[i,])
            competitors$dist <- spDistsN1(pts = as.matrix(competitors[,c("x","y")]), pt = c(move$x[1],move$y[1]), longlat = TRUE)
            competitors.near <- subset(competitors, dist <= integrate.dist)
            ndisp <- ndisp + 1
          #print(fems[i],)
          }
        }
        #print(competitors.near)
        #sd dev of all potential mates a female can hear
        comp.traits.sd <- sd(competitors.near$trait)
        winner.id <- NULL
        #if there is only one male, she mates
        if (is.na(comp.traits.sd)){
          if (nrow(competitors.near) > 0){
          winner.id <- competitors.near$id[1]
          prev.fledge.m <- inds[inds$id == winner.id,]$male.fledglings
          prev.fledge.f <- inds[inds$id == winner.id,]$female.fledglings
          n.children <- round(rnorm(1, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
          if (n.children < 0) {
            n.children = 0
          }
          sex.chick <- round(runif(n.children,0,1))
          for (i in sex.chick) {
            if (i==1){
              inds[inds$id == winner.id,]$male.fledglings <- prev.fledge.m + 1
            }
            else {
              inds[inds$id == winner.id,]$female.fledglings <- prev.fledge.f + 1
            }
          }
          }
          #if more potential mates those that fall within her selective ability have an equal chance of mating
        } else if (!is.na(comp.traits.sd)){
          selective.level <- comp.traits.sd*selectivity
          fem.pref <-as.numeric(fems$trait[i]) #females prefer the traits they learned
          competitors.near$trait.var <- abs(as.numeric(competitors.near$trait)-fem.pref)
          close.competitors <- subset(competitors.near, trait.var <= selective.level)
          if (nrow(close.competitors)==0){ #if none are within her selective ability, she'll mate with the closest to her preference
            winner.id <- competitors.near$id[competitors.near$trait.var == min(competitors.near$trait.var)]
            prev.fledge.m <- inds[inds$id == winner.id,]$male.fledglings
            prev.fledge.f <- inds[inds$id == winner.id,]$female.fledglings
            n.children <- round(rnorm(1, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
            if (n.children < 0) {
              n.children = 0
            }
            sex.chick <- round(runif(n.children,0,1))
            for (i in sex.chick) {
              if (i==1){
                inds[inds$id == winner.id,]$male.fledglings <- prev.fledge.m + 1
              }
              else {
                inds[inds$id == winner.id,]$female.fledglings <- prev.fledge.f + 1
              }
            }
          }
          else { #if there are competitors that fall within her selective ability she choses one at random
            winner <- sample_n(close.competitors, 1)
            winner.id <- winner$id
            prev.fledge.m <- inds[inds$id == winner.id,]$male.fledglings
            prev.fledge.f <- inds[inds$id == winner.id,]$female.fledglings
            n.children <- round(rnorm(1, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
            if (n.children < 0) {
              n.children = 0
            }
            sex.chick <- round(runif(n.children,0,1))
            for (i in sex.chick) {
              if (i==1){
                inds[inds$id == winner.id,]$male.fledglings <- prev.fledge.m + 1
              }
              else {
                inds[inds$id == winner.id,]$female.fledglings <- prev.fledge.f + 1
              }
            }
          }
        }
        
      }
    } else if (!mate.comp) { #without mate competition, all males with territories have children
      n.children <- round(rnorm(nrow(inds[inds$territory==1,]), mean=male.fledge.n.mean, sd=male.fledge.n.sd))
      male.children <- c()
      female.children <- c()
      for (i in 1:length(n.children)){
        sex.chick <- round(runif(n.children[i],0,1))
        males <- 0
        females <- 0
        for (n in sex.chick) {
          if (n==1){
            males <- males+1
          } else {
            females <- females+1
          }
        }
        male.children[i] <- males
        female.children[i] <- females
      }
      
      inds$male.fledglings <- (inds$territory==1) * male.children
      inds$female.fledglings <- (inds$territory==1) * female.children
    }
    #competitors <- spDists(as.matrix(inds[,c("x","y")]),longlat = TRUE) <= integrate.dist
    #print(competitors)
    #ninds <- length(inds$age)
    #prev.songs <- subset(inds, inds$age > 1)$trait
    #mated.males <- subset(inds, abs(inds$trait-mean(prev.songs)) < (selectivity*sd(prev.songs)) & inds$territory==1)
    #print(paste('n mated males: ', nrow(mated.males))) 
    #comp.res<- (!mate.comp | abs(inds$trait-mean(prev.songs)) < (selectivity*sd(prev.songs)))
    #inds$male.fledglings <- (inds$territory==1) * comp.res * round(rnorm(ninds, mean=male.fledge.n.mean, sd=male.fledge.n.sd))
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
    inds$female.fledglings[inds$female.fledglings < 0] <- 0	
    inds 
  }
  ################################ LIFE LOOP ##################################
  ratio.nums <- c()
  for (b in 1:iteration){
  inds <- init.inds
  fems <- create.females(inds = init.inds, females = females)
  maxid.m <- max(inds$id) #store max id, which determines new id numbers in Hatch
  maxid.f <- max(fems$id)
  
  
  k <- 1
  while (k <= steps & nrow(inds) >= 1 & nrow(fems) >=1){
    
    if (prin==TRUE){
      sex.ratio <- nrow(fems)/nrow(inds)
      print(paste("iteration ", b, ", timestep ", k, ", n males ", NROW(inds),", n females ", nrow(fems), ", n territorial males ", length(which(inds$territory==1)), ', sex ratio: ', sex.ratio, sep=""))}
    if (NROW(inds) >= 3 & NROW(fems) >=1){
      timestep.inds <- inds
      timestep.inds$timestep <- k
      timestep.inds$iteration <- b
      # all.inds <- rbind(all.inds, timestep.inds)
      inds.all_list[[b]][[k]]<-timestep.inds
      inds$fitness <- 1
      chicks <- hatch(inds = inds)
      chicks.m <- subset(chicks, sex=='M')
      chicks.f <- subset(chicks, sex=='F', select = -c(male.fledglings,female.fledglings, territory, fitness, learn.dir))
      if (typeof(content.bias)=='double'){
        inds <- content(inds, content.bias.loc, content.bias.loc.ranges,affected.traits, content.bias)
      }
    
      if (learn.m=='default'){
        chicks.m <- learn(tutors = inds, children = chicks.m, conformity.bias = conformity.bias, integrate.dist = integrate.dist, prestige.bias.strength = prestige.bias)
      } else { #user function
        ids <- chicks.m$id
        chicks.m$trait <- sapply(ids, learn.m)
      }
      if (learn.f == 'default'){
        chicks.f <- learn(tutors = inds, children = chicks.f, conformity.bias = conformity.bias, integrate.dist = integrate.dist, prestige.bias.strength = prestige.bias)
      } else{ #user function
        ids <- chicks.f$id
        chicks.f$trait <- sapply(ids,learn.f)
      }
      maxid.m <- max(chicks.m$id)
      maxid.f <- max(chicks.f$id)
      inds <- die(inds)
      fems <- die(fems)
      chicks.m <- die(chicks.m)
      chicks.f <- die(chicks.f)
      
      if (NROW(inds)+NROW(chicks.m) >= 3 & nrow(fems)+nrow(chicks.f) >= 3){
        inds <- grow(inds)
        fems <- grow(fems)
        chicks.m <-grow(chicks.m)
        chicks.f <- grow(chicks.f)
        fems <- rbind(fems, chicks.f)
        fems <- disperse(fems)
        chicks.m <- disperse(chicks.m)
        inds <- rbind(inds, chicks.m)
        inds <- compete.for.territories(inds)
        if (typeof(mate.comp)=='logical'){
          inds <- compete.for.mates(inds, fems, selectivity)
        } else{ #user function
          ids <- as.numeric(inds$id)
          fledglings <- sapply(ids, mate.comp)
          for (i in 1:length(fledglings)){
            if (fledglings[i] < 0){
              fledglings[i] <- 0
            }
            sex.chick <- round(runif(fledglings, 0,1))
            males <- 0
            females <- 0
            for (n in sex.chick){
              if (n==1){
                males <- males+1
              } else {
                females <- females+1
              }
            }
            inds$male.fledglings[i] <- males
            inds$female.fledglings[i] <- females
          }
        }
        
        
        
        
        summary.results[b, , "sample.n"][k] <- length(inds$age)
        summary.results[b, , "trait.pop.mean"][k] <- mean(as.numeric(subset(inds, inds$age==2)$trait))
        summary.results[b, , "trait.pop.variance"][k] <- var(as.numeric(subset(inds, inds$age==2)$trait))

        
        boot_obj <- boot(as.numeric(inds$trait), statistic=sample.mean, R=100)#, strata=mn.res$iteration)	

        ci.res <- boot.ci(boot_obj, conf=0.95, type="basic")

        summary.results[b, , "lci"][k] <- ci.res$basic[4]
        summary.results[b, , "uci"][k] <- ci.res$basic[5]

      }
    }
  k <- k+1}
  }

  if (all==TRUE){
    all.inds=rbind(init=all.inds0,  do.call(rbind,do.call(c,inds.all_list)))
    # coordinates(all.inds) = ~x+y 
    # proj4string(all.inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
    all.inds=fast.coords.frame(all.inds)
    # coordinates(inds) = ~x+y 
    # proj4string(inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
    inds=fast.coords.frame(inds)
    z <- list("summary.results"=summary.results, "inds"=inds, "all.inds"=all.inds, "content.bias.info"=content.bias.info, "time"=proc.time()-ptm)
  }else if (all=="sparse"){
    # all.inds0=fast.coords.frame(all.inds0)
    inds=fast.coords.frame(inds)
    all.inds<-rapply(inds.all_list,fast.coords.frame,classes = "data.frame",how = "replace")
    # all.inds=rbind(init=all.inds0,  do.call(rbind,do.call(c,inds.all_list)))
    # coordinates(all.inds0) = ~x+y 
    # proj4string(all.inds0) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
    # coordinates(inds) = ~x+y 
    # proj4string(inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
    # 
    z <- list(summary.results=summary.results, inds.last=inds, inds.init=all.inds0, inds.slices=all.inds, "content.bias.info"=content.bias.info, time=proc.time()-ptm)
  }else{
    inds=fast.coords.frame(inds)
    # coordinates(inds) = ~x+y 
    # proj4string(inds) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84") 
    z <- list("summary.results"=summary.results, "inds"=inds, "content.bias.info"=content.bias.info, "time"=proc.time()-ptm)
  }
  
  # m_pnt(11)
  
  z
}

#End of I. SongEvo function
#########################

fast.coords.frame <- function(data.src,x.col="x",y.col="y"){
  coor.cols=c(which(colnames(data.src)%in% x.col),
              which(colnames(data.src)%in% y.col));
  SpatialPointsDataFrame(coords = data.src[,coor.cols],
                         data = data.src[,-coor.cols],
                         coords.nrs = coor.cols,
                         match.ID = FALSE,
                         proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  ) 
}
sample.mean <- function(d, x) {
  mean(d[x])
}
