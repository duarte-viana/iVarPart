###########################################
# Supporting Information
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


########################################################################################
########################################################################################

# Abundance data

# Load data 
# (list with each element consisting of the object "lout" produced in "Data_simulation.R")
load("Abundance_data.Rdata") 
simn <- sim[[task.id]]
mat.sp <- simn$mat.sp[,which(apply(simn$mat.sp,2,sum)>0)] # exclude possibly extinct species (due to cropping)
mat.env <- simn$mat.env
mat.spa <- simn$mat.spa
mat.xy <- simn$mat.xy
mem.listw <- simn$mem.listw

# MSR
nresamp <- 99
emsr <- msr(mat.env,mem.listw,nrepet=nresamp,simplify=FALSE)

# Gather methods names
meths <- c("R2.CCA", "R2.dbRDA", "R2.RDA", "R2.MRM", 
           "R2.GLM", "R2.hmsc", "R2.hmsc", "R2.GAM",
           "R2.BRT", "R2.BRT", "R2.BRT", "R2.BRT", "R2.UniRndForest", 
           "R2.MVRndForest", "R2.MVRndForest", "R2.MVRegTree", "R2.MVRegTree")

# Number of differrent R2 indices calculated
n.r2 <- 4

# Prepare output object
res.env <- array(NA, c(1, length(meths), n.r2))
res.shared <- array(NA, c(1, length(meths), n.r2))
res.spa <- array(NA, c(1, length(meths), n.r2))
res1 <- array(NA, c(1, length(meths), n.r2))
res2 <- array(NA, c(1, length(meths), n.r2))
res3 <- array(NA, c(1, length(meths), n.r2))

# Loop to apply each method to given simulated dataset
for(m in 1:length(meths)){
  
  try({
    
    # Constrained ordination
    # varpart for methods using MEM
    if(m %in% 1:2){
      nsites <- nrow(mat.sp)
      nenv <- ncol(mat.env)
      nspa <- ncol(mat.spa)
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
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    
    if(m %in% 3){
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
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    # Distance based regression
    if(m %in% c(4)){
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
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    # GLM quasiPoisson
    # R2 is adjusted 
    if(m %in% 5){
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
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    # GLMM
    if(m %in% 6){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="poisson"))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),family="poisson"))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),family="poisson"))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),family="poisson"))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa),family="poisson"))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    if(m %in% 7){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="overPoisson"))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),family="overPoisson"))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),family="overPoisson"))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),family="overPoisson"))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa),family="overPoisson"))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    # GAM
    if(m %in% c(8)){
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
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    # Tree-based methods
    # varpart for methods using the xy coordinates
    if(m %in% c(9,13,14,16)){
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
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    if(m %in% c(10)){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),inter.depth=1))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),inter.depth=1))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),inter.depth=1))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),inter.depth=1))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),inter.depth=1))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    if(m %in% c(11)){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),CV=TRUE,inter.depth=1))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),CV=TRUE,inter.depth=1))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),CV=TRUE,inter.depth=1))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),CV=TRUE,inter.depth=1))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),CV=TRUE,inter.depth=1))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    if(m %in% c(12,15,17)){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),CV=TRUE))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),CV=TRUE))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),CV=TRUE))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),CV=TRUE))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),CV=TRUE))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2cla <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2log <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2mv <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
      vp.R2fad <- VarPartClap(r2.env.raw[4],r2.spa.raw[4],r2.spa.adj[4],r2.env.spa.raw[4],msr.env[,4],msr.env.spa[,4])
    }
    
    # Store results
    vp.msr <- rbind(vp.R2cla,vp.R2log,vp.R2mv,vp.R2fad)
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


########################################################################################
########################################################################################

# Binary occurrence data

# Load data 
# (list with each element consisting of the object "lout" produced in "Data_simulation.R")
load("Binary_data.Rdata")
simn <- sim[[task.id]]
mat.sp <- simn$mat.sp[,which(apply(simn$mat.sp,2,sum)>0)] # exclude possibly extinct species (due to cropping)
mat.sp <- mat.sp[,which(apply(mat.sp,2,sum)<=(nrow(mat.sp)-1))] # exclude species present everywhere
mat.env <- simn$mat.env
mat.spa <- simn$mat.spa
mat.xy <- simn$mat.xy
mem.listw <- simn$mem.listw

# MSR
nresamp <- 99
emsr <- msr(mat.env,mem.listw,nrepet=nresamp,simplify=FALSE)

# Gather methods names
meths <- c("R2.CCA", "R2.dbRDA", "R2.RDA", "R2.MRM", 
           "R2.GLM", "R2.hmsc", "R2.GAM",
           "R2.BRT", "R2.BRT", "R2.BRT", "R2.BRT", "R2.UniRndForest", 
           "R2.MVRndForest", "R2.MVRndForest", "R2.MVRegTree", "R2.MVRegTree")

# Number of differrent R2 indices calculated
n.r2 <- 3

# Prepare output object
res.env <- array(NA, c(1, length(meths), n.r2))
res.shared <- array(NA, c(1, length(meths), n.r2))
res.spa <- array(NA, c(1, length(meths), n.r2))
res1 <- array(NA, c(1, length(meths), n.r2))
res2 <- array(NA, c(1, length(meths), n.r2))
res3 <- array(NA, c(1, length(meths), n.r2))

# Loop to apply each method to given simulated dataset
for(m in 1:length(meths)){
  
  try({
    
    # Constrained ordination
    # varpart for methods using MEM
    if(m %in% 1:2){
      nsites <- nrow(mat.sp)
      nenv <- ncol(mat.env)
      nspa <- ncol(mat.spa)
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),binary=TRUE))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),binary=TRUE))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),binary=TRUE))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa),binary=TRUE))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    if(m %in% 3){
      nsites <- nrow(mat.sp)
      nenv <- ncol(mat.env)
      nspa <- ncol(mat.spa)
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
      r2.env.adj <- 1-(((nsites-1)/(nsites-nenv-1))*(1-r2.env.raw))
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),binary=TRUE))
      r2.spa.adj <- 1-(((nsites-1)/(nsites-nspa-1))*(1-r2.spa.raw))
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),binary=TRUE))
      r2.env.spa.adj <- 1-(((nsites-1)/(nsites-(nenv+nspa)-1))*(1-r2.env.spa.raw))
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),binary=TRUE))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa),binary=TRUE))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    # Distance based regression
    if(m %in% c(4)){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),binary=TRUE))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),binary=TRUE))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),binary=TRUE))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),binary=TRUE))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    # GLM binomial
    # R2 is adjusted 
    if(m %in% 5){
      nsites <- nrow(mat.sp)
      # Calculate the R2 for each model
      r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="binomial"))
      r2.env.raw <- r2.env[[1]]
      np <- r2.env[[2]]
      r2.env.adj <- 1-(((nsites-1)/(nsites-np-1))*(1-r2.env.raw))
      r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),family="binomial"))
      r2.spa.raw <- r2.spa[[1]]
      np <- r2.spa[[2]]
      r2.spa.adj <- 1-(((nsites-1)/(nsites-np-1))*(1-r2.spa.raw))
      r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),family="binomial"))
      r2.env.spa.raw <- r2.env.spa[[1]]
      np <- r2.env.spa[[2]]
      r2.env.spa.adj <- 1-(((nsites-1)/(nsites-np-1))*(1-r2.env.spa.raw))
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),family="binomial"))[[1]]
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa),family="binomial"))[[1]]
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    # GLMM
    if(m %in% 6){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="probit"))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),family="probit"))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),family="probit"))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),family="probit"))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],mem=mat.spa),family="probit"))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    # GAM
    if(m %in% c(7)){
      nsites <- nrow(mat.sp)
      nenv <- 4
      n.coef.xy <- 0 # as s in GAM penalises overfitting
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="binomial"))
      r2.env.adj <- 1-(((nsites-1)/(nsites-nenv-1))*(1-r2.env.raw))
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="binomial"))
      r2.spa.adj <- 1-(((nsites-1)/(nsites-n.coef.xy-1))*(1-r2.spa.raw))
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="binomial"))
      r2.env.spa.adj <- 1-(((nsites-1)/(nsites-(nenv+n.coef.xy)-1))*(1-r2.env.spa.raw))
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),family="binomial"))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),family="binomial"))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    # Tree-based methods
    # varpart for methods using the xy coordinates
    if(m %in% c(8,12,13,15)){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),binary=TRUE))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),binary=TRUE))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),binary=TRUE))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),binary=TRUE))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    if(m %in% c(9)){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),inter.depth=1,binary=TRUE))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),inter.depth=1,binary=TRUE))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),inter.depth=1,binary=TRUE))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),inter.depth=1,binary=TRUE))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),inter.depth=1,binary=TRUE))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    if(m %in% c(10)){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),CV=TRUE,inter.depth=1,binary=TRUE))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),CV=TRUE,inter.depth=1,binary=TRUE))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),CV=TRUE,inter.depth=1,binary=TRUE))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),CV=TRUE,inter.depth=1,binary=TRUE))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),CV=TRUE,inter.depth=1,binary=TRUE))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    if(m %in% c(11,14,16)){
      # Calculate the R2 for each model
      r2.env.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env),CV=TRUE,binary=TRUE))
      r2.env.adj <- r2.env.raw
      r2.spa.raw <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),CV=TRUE,binary=TRUE))
      r2.spa.adj <- r2.spa.raw
      r2.env.spa.raw <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),CV=TRUE,binary=TRUE))
      r2.env.spa.adj <- r2.env.spa.raw
      
      # resampled R2 
      msr.env <- matrix(nrow=nresamp,ncol=n.r2)
      msr.env.spa <- matrix(nrow=nresamp,ncol=n.r2)
      for(i in 1:nresamp){
        R2e <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]]),CV=TRUE,binary=TRUE))
        msr.env[i,] <- R2e
        R2s <- do.call(meths[m],list(mat.sp,list(env=emsr[[i]],xy=mat.xy),CV=TRUE,binary=TRUE))
        msr.env.spa[i,] <- R2s
      }
      
      # corrected variation partitioning (Clappe et al. 2018)
      vp.R2fad <- VarPartClap(r2.env.raw[1],r2.spa.raw[1],r2.spa.adj[1],r2.env.spa.raw[1],msr.env[,1],msr.env.spa[,1])
      vp.R2efr <- VarPartClap(r2.env.raw[2],r2.spa.raw[2],r2.spa.adj[2],r2.env.spa.raw[2],msr.env[,2],msr.env.spa[,2])
      vp.R2tju <- VarPartClap(r2.env.raw[3],r2.spa.raw[3],r2.spa.adj[3],r2.env.spa.raw[3],msr.env[,3],msr.env.spa[,3])
    }
    
    # Store results
    vp.msr <- rbind(vp.R2fad,vp.R2efr,vp.R2tju)
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
