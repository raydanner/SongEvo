# Benchmarking / stress test code

source("SongEvo.R")
reg.count=40
fledge_counts=c(1, 1, 2, 1, 0, 2, 2, 2, 2, 1, 1, 2, 1, 
                1, 1, 1, 1, 1, 0, 2, 1, 2, 1, 2, 2, 2, 
                1, 1, 1, 2, 1, 1, 
                1, 2, 2, 1, 1, 2, 1, 1)

default.parms <- list(init.inds = NULL, iteration = 10, steps = NULL,  # years / timestep
                      timestep = 1, #steps in years by default
                      n.territories = reg.count, terr.turnover = 0.5, 
                      learning.method = "integrate", integrate.dist = 0.1, 
                      learning.error.d = 0, learning.error.sd = 430, 
                      mortality.a = 0.468, mortality.j = 0.5, lifespan = NA, 
                      phys.lim.min = 1559, phys.lim.max = 4364, 
                      male.fledge.n.mean = 1.35, male.fledge.n.sd = 0.5,
                      male.fledge.n = fledge_counts, disp.age = 2, 
                      disp.distance.mean = 110, disp.distance.sd = 100,
                      mate.comp = FALSE, prin = FALSE, all = FALSE)


ptm_laps <- do.call(rbind,structure(list(NULL,proc.time()),.Names=c(NA,"start")))
for(r_count in 1:10 ){
  gc()
  
  # Make starting population for prbo in 1969
  local({
    data("WCSP")
    starting.trait.obs <- subset(song.data, Population %in% "PRBO" & Year %in% 1969)$Trill.FBW
    
    # Generate additional individuals to reach number of territories for this vignette (40)
    starting.trait.gen <- rnorm(reg.count-length(starting.trait.obs), 
                                mean=mean(starting.trait.obs), sd=sd(starting.trait.obs))
    starting.trait <- c(starting.trait.obs, starting.trait.gen)
    
    # Build data frame of territory id, age, and song trait value (in this case trill bandwidth)
    init.inds <- data.frame(id = 1:reg.count, age = 2, trait = starting.trait,
                            row.names=paste0("A",formatC(1:reg.count, width=2, flag = "0")))
    
    # Generate lat long coordinates for center? of each territory
    init.inds$x1 <-  round(runif(reg.count, min=-122.481858, max=-122.447270), digits=8)
    init.inds$y1 <-  round(runif(reg.count, min=37.787768, max=37.805645), digits=8)
    init.inds
  }) -> trill.prbo.1969
  
  # Make starting population for Schooner in 1969
  local({
    data("WCSP")
    starting.trait.obs <- subset(song.data, Population %in% "Schooner" & Year %in% 1969)$Trill.FBW
    
    # Generate additional individuals to reach number of territories for this vignette (40)
    starting.trait.gen <- rnorm(reg.count-length(starting.trait.obs), 
                                mean=mean(starting.trait.obs), sd=sd(starting.trait.obs))
    starting.trait <- c(starting.trait.obs, starting.trait.gen)
    
    # Build data frame of territory id, age, and song trait value (in this case trill bandwidth)
    init.inds <- data.frame(id = 1:reg.count, age = 2, trait = starting.trait,
                            row.names=paste0("A",formatC(1:reg.count, width=2, flag = "0")))
    
    # Generate lat long coordinates for center? of each territory
    init.inds$x1 <-  round(runif(reg.count, min=-122.481858, max=-122.447270), digits=8)
    init.inds$y1 <-  round(runif(reg.count, min=37.787768, max=37.805645), digits=8)
    init.inds
  }) -> trill.sb.1969
  
  # Make starting population for Bear Valley in 1969
  local({
    data("WCSP")
    starting.trait.obs <- subset(song.data, Population %in% "Bear Valley" & Year %in% 1969)$Trill.FBW
    
    # Generate additional individuals to reach number of territories for this vignette (40)
    starting.trait.gen <- rnorm(reg.count-length(starting.trait.obs), 
                                mean=mean(starting.trait.obs), sd=sd(starting.trait.obs))
    starting.trait <- c(starting.trait.obs, starting.trait.gen)
    
    # Build data frame of territory id, age, and song trait value (in this case trill bandwidth)
    init.inds <- data.frame(id = 1:reg.count, age = 2, trait = starting.trait,
                            row.names=paste0("A",formatC(1:reg.count, width=2, flag = "0")))
    
    # Generate lat long coordinates for center? of each territory
    init.inds$x1 <-  round(runif(reg.count, min=-122.481858, max=-122.447270), digits=8)
    init.inds$y1 <-  round(runif(reg.count, min=37.787768, max=37.805645), digits=8)
    init.inds
  }) -> trill.bear.1969
  # SongEvo args for prbo 1969
  trill.prbo.test=list(parms=default.parms)
  trill.prbo.test$parms$init.inds = trill.prbo.1969 
  trill.prbo.test$parms$steps = 2005-1969
  
  # SongEvo args for Schooner Bay 1969
  trill.sb.test=list(parms=default.parms)
  trill.sb.test$parms$init.inds = trill.sb.1969 
  trill.sb.test$parms$steps = 2005-1969
  
  # SongEvo args for Bear Valley 1969
  trill.bear.test=list(parms=default.parms)
  trill.bear.test$parms$init.inds = trill.bear.1969 
  trill.bear.test$parms$steps = 2005-1969
  ptm_laps <- do.call(rbind,structure(list(ptm_laps,proc.time()),
                          .Names=c(NA,paste0("init",formatC(r_count,width=3,flag="0")) )))
  
  
  trill.prbo.test$res=do.call(SongEvo,trill.prbo.test$parms)
  ptm_laps <- do.call(rbind,structure(list(ptm_laps,proc.time()),
                                      .Names=c(NA,paste0("prbo",formatC(r_count,width=3,flag="0")) )))
  trill.sb.test$res=  do.call(SongEvo,trill.sb.test$parms)
  ptm_laps <- do.call(rbind,structure(list(ptm_laps,proc.time()),
                                      .Names=c(NA,paste0("sbay",formatC(r_count,width=3,flag="0")) )))
  trill.bear.test$res=do.call(SongEvo,trill.bear.test$parms)
  ptm_laps <- do.call(rbind,structure(list(ptm_laps,proc.time()),
                                      .Names=c(NA,paste0("bear",formatC(r_count,width=3,flag="0")) )))
}

ptm_laps <- do.call(rbind,structure(list(ptm_laps,proc.time()),
                                    .Names=c(NA,"end" )))
print(ptm_laps)

write.csv(ptm_laps,
          file=(most_recent_benchmark<-
                  paste0("profile_SongEvo_3x10x10_",
                         format(Sys.time(),format="%Y%m%d%H%M"),".csv")))
