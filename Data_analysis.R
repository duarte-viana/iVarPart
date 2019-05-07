###########################################
# Supporting Information, Appendix S3
# "Partitioning environment and space in species-by-site matrices: a comparison of methods for community ecology and macroecology"
# Duarte S. Viana, Petr Keil, Alienor Jeliazkov
###########################################


# Code for analysing the data according to the different statistical methods
# See README.txt file
# This is the code for the analysis of one dataset corresponding to the scenario identified as "task.id",
# which corresponds to a given data simulation

# ------------------------------------------------------------------------------


# load functions
source("Methods_functions.R")
source("R2D2_functions.R")

# Load data 
# (list with each element consisting of the object "lout" produced in "Data_simulation.R")
load("data.Rdata") 
simn <- sim[[task.id]]
mat.sp <- simn$mat.sp[,which(apply(simn$mat.sp,2,sum)>0)] # exclude possibly extinct species (due to cropping)
mat.env <- simn$mat.env
mat.spa <- simn$mat.spa
mat.xy <- simn$mat.xy
mem.listw <- simn$mem.listw

# Gather methods names
meths <- c("R2.CCA.raw", "R2.CCA.hel", "R2.RDA.hel", "R2.dbRDA.raw", "R2.MRM", 
           "R2.GLM.quasiPois","R2.hmsc.Pois", "R2.hmsc.overPois", "R2.GAM",
           "R2.BRT", "R2.UniRndForest", "R2.MVRndForest", "R2.MVRegTree")

# Number of differrent R2 indices calculated
n.r2 <- 6

# Prepare output object
res.env <- array(NA, c(1, length(meths), n.r2))
res.shared <- array(NA, c(1, length(meths), n.r2))
res.spa <- array(NA, c(1, length(meths), n.r2))
res1 <- array(NA, c(1, length(meths), n.r2))
res2 <- array(NA, c(1, length(meths), n.r2))
res3 <- array(NA, c(1, length(meths), n.r2))

# Loop to apply each method to given simulated dataset
for(m in 1:length(meths)){
  
  # MSR
  nresamp <- 99
  emsr <- msr(mat.env,mem.listw,nrepet=nresamp,simplify=FALSE)
  
  try({
  
  # Constrained ordination
  # varpart for methods using MEM
  if(m %in% 1:4){
    nsites <- nrow(mat.sp)
    nenv <- ncol(mat.env)
    nspa <- ncol(mat.spa)
    # Calculate the R2 for each model
    r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
    r2.env.adj <- 1-(((nsites-1)/(nsites-nenv-1))*(1-r2.env.raw))
    r2.spa.raw <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
    r2.spa.adj <- 1-(((nsites-1)/(nsites-nspa-1))*(1-r2.spa.raw))
    r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
    r2.env.spa.adj <- 1-(((nsites-1)/(nsites-(nenv+nspa)-1))*(1-r2.env.spa.raw))
    
    # resampled R2 
    msr.env <- matrix(nrow=nresamp,ncol=n.r2)
    msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
    for(i in 1:nresamp){
      R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]])))
      msr.env[i,] <- R2e
      R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa)))
      msr.env.spa[i,] <- R2s
    }
    
    # corrected variation partitioning (Clappe et al. 2018)
    vp.R2avr <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
    vp.R2com <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
    vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    vp.D2avr <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    vp.D2com <- VarPartClap(r2.env.raw[5],r2.spa.raw[5],r2.spa.adj[5],r2.env.spa.raw[5],msr.env[,5],msr.env.spa[,5])
    vp.R2log <- VarPartClap(r2.env.raw[6],r2.spa.raw[6],r2.spa.adj[6],r2.env.spa.raw[6],msr.env[,6],msr.env.spa[,6])
  }
  
  # Distance based regression
  if(m %in% c(5)){
    # Calculate the R2 for each model
    r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
    r2.env.adj <- r2.env.raw
    r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
    r2.spa.adj <- r2.spa.raw
    r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
    r2.env.spa.adj <- r2.env.spa.raw
    
    # resampled R2 
    msr.env <- matrix(nrow=nresamp,ncol=n.r2)
    msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
    for(i in 1:nresamp){
      R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]])))
      msr.env[i,] <- R2e
      R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy)))
      msr.env.spa[i,] <- R2s
    }
    
    # corrected variation partitioning (Clappe et al. 2018)
    vp.R2avr <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
    vp.R2com <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
    vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    vp.D2avr <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    vp.D2com <- VarPartClap(r2.env.raw[5],r2.spa.raw[5],r2.spa.adj[5],r2.env.spa.raw[5],msr.env[,5],msr.env.spa[,5])
    vp.R2log <- VarPartClap(r2.env.raw[6],r2.spa.raw[6],r2.spa.adj[6],r2.env.spa.raw[6],msr.env[,6],msr.env.spa[,6])
  }
  
  # GLM quasiPoisson
  # R2 is adjusted 
  if(m %in% 6){
    nsites <- nrow(mat.sp)
    # Calculate the R2 for each model
    r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
    r2.env.raw <- r2.env[[1]]
    np <- r2.env[[2]]
    r2.env.adj <- 1-(((nsites-1)/(nsites-np-1))*(1-r2.env.raw))
    r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
    r2.spa.raw <- r2.spa[[1]]
    np <- r2.spa[[2]]
    r2.spa.adj <- 1-(((nsites-1)/(nsites-np-1))*(1-r2.spa.raw))
    r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
    r2.env.spa.raw <- r2.env.spa[[1]]
    np <- r2.env.spa[[2]]
    r2.env.spa.adj <- 1-(((nsites-1)/(nsites-np-1))*(1-r2.env.spa.raw))
    
    # resampled R2 
    msr.env <- matrix(nrow=nresamp,ncol=n.r2)
    msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
    for(i in 1:nresamp){
      R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]])))[[1]]
      msr.env[i,] <- R2e
      R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa)))[[1]]
      msr.env.spa[i,] <- R2s
    }
    
    # corrected variation partitioning (Clappe et al. 2018)
    vp.R2avr <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
    vp.R2com <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
    vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    vp.D2avr <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    vp.D2com <- VarPartClap(r2.env.raw[5],r2.spa.raw[5],r2.spa.adj[5],r2.env.spa.raw[5],msr.env[,5],msr.env.spa[,5])
    vp.R2log <- VarPartClap(r2.env.raw[6],r2.spa.raw[6],r2.spa.adj[6],r2.env.spa.raw[6],msr.env[,6],msr.env.spa[,6])
  }
  
  # GLMM
  if(m %in% 7:8){
    # Calculate the R2 for each model
    r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
    r2.env.adj <- r2.env.raw
    r2.spa.raw <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
    r2.spa.adj <- r2.spa.raw
    r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
    r2.env.spa.adj <- r2.env.spa.raw
      
    # resampled R2 
    msr.env <- matrix(nrow=nresamp,ncol=n.r2)
    msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
    for(i in 1:nresamp){
      R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]])))
      msr.env[i,] <- R2e
      R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa)))
      msr.env.spa[i,] <- R2s
    }
      
    # corrected variation partitioning (Clappe et al. 2018)
    vp.R2avr <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
    vp.R2com <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
    vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    vp.D2avr <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    vp.D2com <- VarPartClap(r2.env.raw[5],r2.spa.raw[5],r2.spa.adj[5],r2.env.spa.raw[5],msr.env[,5],msr.env.spa[,5])
    vp.R2log <- VarPartClap(r2.env.raw[6],r2.spa.raw[6],r2.spa.adj[6],r2.env.spa.raw[6],msr.env[,6],msr.env.spa[,6])
  }
  
  # GAM
  if(m %in% c(9)){
    nsites <- nrow(mat.sp)
    nenv <- 4
    n.coef.xy <- 0 # as s in GAM penalises overfitting
    # Calculate the R2 for each model
    r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
    r2.env.adj <- 1-(((nsites-1)/(nsites-nenv-1))*(1-r2.env.raw))
    r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
    r2.spa.adj <- 1-(((nsites-1)/(nsites-n.coef.xy-1))*(1-r2.spa.raw))
    r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
    r2.env.spa.adj <- 1-(((nsites-1)/(nsites-(nenv+n.coef.xy)-1))*(1-r2.env.spa.raw))
    
    # resampled R2 
    msr.env <- matrix(nrow=nresamp,ncol=n.r2)
    msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
    for(i in 1:nresamp){
      R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]])))
      msr.env[i,] <- R2e
      R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy)))
      msr.env.spa[i,] <- R2s
    }
    
    # corrected variation partitioning (Clappe et al. 2018)
    vp.R2avr <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
    vp.R2com <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
    vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    vp.D2avr <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    vp.D2com <- VarPartClap(r2.env.raw[5],r2.spa.raw[5],r2.spa.adj[5],r2.env.spa.raw[5],msr.env[,5],msr.env.spa[,5])
    vp.R2log <- VarPartClap(r2.env.raw[6],r2.spa.raw[6],r2.spa.adj[6],r2.env.spa.raw[6],msr.env[,6],msr.env.spa[,6])
  }
  
  # Tree-based methods
  # varpart for methods using the xy coordinates
  if(m %in% 10:13){
    # Calculate the R2 for each model
    r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
    r2.env.adj <- r2.env.raw
    r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
    r2.spa.adj <- r2.spa.raw
    r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
    r2.env.spa.adj <- r2.env.spa.raw
    
    # resampled R2 
    msr.env <- matrix(nrow=nresamp,ncol=n.r2)
    msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
    for(i in 1:nresamp){
      R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]])))
      msr.env[i,] <- R2e
      R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy)))
      msr.env.spa[i,] <- R2s
    }
    
    # corrected variation partitioning (Clappe et al. 2018)
    vp.R2avr <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.raw[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
    vp.R2com <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.raw[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
    vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.raw[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    vp.D2avr <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.raw[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    vp.D2com <- VarPartClap(r2.env.raw[5],r2.spa.raw[5],r2.spa.raw[5],r2.env.spa.raw[5],msr.env[,5],msr.env.spa[,5])
    vp.R2log <- VarPartClap(r2.env.raw[6],r2.spa.raw[6],r2.spa.adj[6],r2.env.spa.raw[6],msr.env[,6],msr.env.spa[,6])
  }
  
  # Store results
  vp.msr <- rbind(vp.R2avr,vp.R2com,vp.R2mv,vp.D2avr,vp.D2com,vp.R2log)
  for(n in 1:n.r2){
    res.env[1,m,n] <- vp.msr[n,1]
    res.shared[1,m,n] <- vp.msr[n,2]
    res.spa[1,m,n] <- vp.msr[n,3]
    res1[1,m,n] <- r2.env.adj[n]
    res2[1,m,n] <- r2.spa.adj[n]
    res3[1,m,n] <- r2.env.spa.adj[n]
  }
  })
}

# Output object
lout <- list(res.env, res.shared, res.spa, res1, res2, res3)





