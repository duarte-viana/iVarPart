###########################################
# Supporting Information
# "Partitioning environment and space in site-by-species matrices: a comparison of methods for community ecology and macroecology"
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

# Load parameters
load("pars_VP.Rdata")

# Load data 
# (list with each element consisting of the object "lout" produced in "Data_simulation.R")
load("sim_data_vp.Rdata")
simn <- sim[[task.id]]
mat.sp <- simn$mat.sp
mat.sp <- mat.sp[,which(apply(mat.sp,2,sum)>0)] # exclude possibly extinct species
mat.env <- simn$mat.env
mat.spa <- simn$mat.spa
mat.xy <- simn$mat.xy


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# normal data

if(pars$data.type[task.id]=="normal"){
  
  # Gather methods names
  meths <- c("R2.CCA", "R2.dbRDA", "R2.RDA", "R2.MRM", 
             "R2.GLM", "R2.GAM",
             "R2.BRT", "R2.BRT", "R2.UniRndForest", "R2.UniRndForest", 
             "R2.MVRndForest", "R2.MVRegTree")
  
  # Number of differrent R2 indices calculated
  n.r2 <- 4
  
  # Prepare output object
  res.env <- array(NA, c(1, length(meths), n.r2, 2))
  res.shared <- array(NA, c(1, length(meths), n.r2, 2))
  res.spa <- array(NA, c(1, length(meths), n.r2, 2))
  res1 <- array(NA, c(1, length(meths), n.r2, 2))
  res2 <- array(NA, c(1, length(meths), n.r2, 2))
  res3 <- array(NA, c(1, length(meths), n.r2, 2))
  
  # Loop to apply each method to given simulated dataset
  for(m in 1:length(meths)){
    
    try({
      
      # CCA and dbRDA
      if(m %in% 1:2){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
      }
      
      # RDA
      if(m %in% 3 && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
      }
      
      # RDA with quadratic effect of X (env)
      if(m %in% 3  && pars$resp[task.id]=="gaussian"){
        mat.sp.log <- log(mat.sp+1)
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp.log,list(env=mat.env),env.eff="quadratic"))
        r2.spa <- do.call(meths[m],list(mat.sp.log,list(mem=mat.spa),env.eff="quadratic"))
        r2.env.spa <- do.call(meths[m],list(mat.sp.log,list(env=mat.env,mem=mat.spa),env.eff="quadratic"))
      }
      
      # Distance based regression
      if(m %in% c(4)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
      }
      
      # GLM quasiPoisson
      # R2 is adjusted 
      if(m %in% 5  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="gaussian"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),env.eff="linear",family="gaussian"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),env.eff="linear",family="gaussian"))
      }
      
      if(m %in% 5  && pars$resp[task.id]=="gaussian"){
        mat.sp.log <- log(mat.sp+1)
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp.log,list(env=mat.env),family="gaussian"))
        r2.spa <- do.call(meths[m],list(mat.sp.log,list(mem=mat.spa),family="gaussian"))
        r2.env.spa <- do.call(meths[m],list(mat.sp.log,list(env=mat.env,mem=mat.spa),family="gaussian"))
      }
      
      # GAM
      if(m %in% 6  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="gaussian"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="gaussian"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="gaussian"))
      }
      
      if(m %in% 6  && pars$resp[task.id]=="gaussian"){
        mat.sp.log <- log(mat.sp+1)
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp.log,list(env=mat.env),family="gaussian"))
        r2.spa <- do.call(meths[m],list(mat.sp.log,list(xy=mat.xy),family="gaussian"))
        r2.env.spa <- do.call(meths[m],list(mat.sp.log,list(env=mat.env,xy=mat.xy),family="gaussian"))
      }
      
      # Tree-based methods
      # varpart for methods using the xy coordinates
      if(m %in% c(7,9,11,12)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
      }
      
      if(m %in% c(8)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),inter.depth=1))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),inter.depth=1))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),inter.depth=1))
      }
      
      if(m %in% c(10)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),sampsize=0.2))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),sampsize=0.2))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),sampsize=0.2))
      }
      
      # variation partitioning 
      vp.R2cla <- VarPart(r2.env[[1]][1],r2.spa[[1]][1],r2.env.spa[[1]][1])
      vp.R2log <- VarPart(r2.env[[1]][2],r2.spa[[1]][2],r2.env.spa[[1]][2])
      vp.R2mv <- VarPart(r2.env[[1]][3],r2.spa[[1]][3],r2.env.spa[[1]][3])
      vp.R2fad <- VarPart(r2.env[[1]][4],r2.spa[[1]][4],r2.env.spa[[1]][4])
      
      vp.R2cla.cv <- VarPart(r2.env[[2]][1],r2.spa[[2]][1],r2.env.spa[[2]][1])
      vp.R2log.cv <- VarPart(r2.env[[2]][2],r2.spa[[2]][2],r2.env.spa[[2]][2])
      vp.R2mv.cv <- VarPart(r2.env[[2]][3],r2.spa[[2]][3],r2.env.spa[[2]][3])
      vp.R2fad.cv <- VarPart(r2.env[[2]][4],r2.spa[[2]][4],r2.env.spa[[2]][4])
      
      # Store results
      vp.msr <- rbind(vp.R2cla,vp.R2log,vp.R2mv,vp.R2fad)
      vp.msr.cv <- rbind(vp.R2cla.cv,vp.R2log.cv,vp.R2mv.cv,vp.R2fad.cv)
      for(n in 1:n.r2){
        res.env[1,m,n,1] <- vp.msr[n,1]
        res.shared[1,m,n,1] <- vp.msr[n,2]
        res.spa[1,m,n,1] <- vp.msr[n,3]
        res1[1,m,n,1] <- r2.env[[1]][n]
        res2[1,m,n,1] <- r2.spa[[1]][n]
        res3[1,m,n,1] <- r2.env.spa[[1]][n]
        # CV
        res.env[1,m,n,2] <- vp.msr.cv[n,1]
        res.shared[1,m,n,2] <- vp.msr.cv[n,2]
        res.spa[1,m,n,2] <- vp.msr.cv[n,3]
        res1[1,m,n,2] <- r2.env[[2]][n]
        res2[1,m,n,2] <- r2.spa[[2]][n]
        res3[1,m,n,2] <- r2.env.spa[[2]][n]
      }
    })
  }
}



#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# count data

if(pars$data.type[task.id]=="counts"){
  
  # Gather methods names
  meths <- c("R2.CCA", "R2.dbRDA", "R2.RDA", "R2.MRM", 
             "R2.GLM", "R2.GAM",
             "R2.BRT", "R2.BRT", "R2.UniRndForest", "R2.UniRndForest", 
             "R2.MVRndForest", "R2.MVRegTree")
  
  # Number of differrent R2 indices calculated
  n.r2 <- 4
  
  # Prepare output object
  res.env <- array(NA, c(1, length(meths), n.r2, 2))
  res.shared <- array(NA, c(1, length(meths), n.r2, 2))
  res.spa <- array(NA, c(1, length(meths), n.r2, 2))
  res1 <- array(NA, c(1, length(meths), n.r2, 2))
  res2 <- array(NA, c(1, length(meths), n.r2, 2))
  res3 <- array(NA, c(1, length(meths), n.r2, 2))
  
  # Loop to apply each method to given simulated dataset
  for(m in 1:length(meths)){
    
    try({
      
      # CCA and dbRDA
      if(m %in% 1:2){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
      }
      
      # RDA
      if(m %in% 3 && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
      }
      
      # RDA with quadratic effect of X (env)
      if(m %in% 3  && pars$resp[task.id]=="gaussian"){
        mat.sp.log <- log(mat.sp+1)
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp.log,list(env=mat.env),env.eff="quadratic"))
        r2.spa <- do.call(meths[m],list(mat.sp.log,list(mem=mat.spa),env.eff="quadratic"))
        r2.env.spa <- do.call(meths[m],list(mat.sp.log,list(env=mat.env,mem=mat.spa),env.eff="quadratic"))
      }
      
      # Distance based regression
      if(m %in% c(4)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
      }
      
      # GLM quasiPoisson
      # R2 is adjusted 
      if(m %in% 5  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="quasipoisson"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),env.eff="linear",family="quasipoisson"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),env.eff="linear",family="quasipoisson"))
      }
      
      if(m %in% 5  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="quasipoisson"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),family="quasipoisson"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),family="quasipoisson"))
      }
      
      # GAM
      if(m %in% 6  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="poisson"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="poisson"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="poisson"))
      }
      
      if(m %in% 6  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="poisson"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="poisson"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="poisson"))
      }
      
      # Tree-based methods
      # varpart for methods using the xy coordinates
      if(m %in% c(7,9,11,12)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
      }
      
      if(m %in% c(8)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),inter.depth=1))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),inter.depth=1))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),inter.depth=1))
      }
      
      if(m %in% c(10)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),sampsize=0.2))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),sampsize=0.2))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),sampsize=0.2))
      }
      
      # variation partitioning 
      vp.R2cla <- VarPart(r2.env[[1]][1],r2.spa[[1]][1],r2.env.spa[[1]][1])
      vp.R2log <- VarPart(r2.env[[1]][2],r2.spa[[1]][2],r2.env.spa[[1]][2])
      vp.R2mv <- VarPart(r2.env[[1]][3],r2.spa[[1]][3],r2.env.spa[[1]][3])
      vp.R2fad <- VarPart(r2.env[[1]][4],r2.spa[[1]][4],r2.env.spa[[1]][4])
      
      vp.R2cla.cv <- VarPart(r2.env[[2]][1],r2.spa[[2]][1],r2.env.spa[[2]][1])
      vp.R2log.cv <- VarPart(r2.env[[2]][2],r2.spa[[2]][2],r2.env.spa[[2]][2])
      vp.R2mv.cv <- VarPart(r2.env[[2]][3],r2.spa[[2]][3],r2.env.spa[[2]][3])
      vp.R2fad.cv <- VarPart(r2.env[[2]][4],r2.spa[[2]][4],r2.env.spa[[2]][4])
      
      # Store results
      vp.msr <- rbind(vp.R2cla,vp.R2log,vp.R2mv,vp.R2fad)
      vp.msr.cv <- rbind(vp.R2cla.cv,vp.R2log.cv,vp.R2mv.cv,vp.R2fad.cv)
      for(n in 1:n.r2){
        res.env[1,m,n,1] <- vp.msr[n,1]
        res.shared[1,m,n,1] <- vp.msr[n,2]
        res.spa[1,m,n,1] <- vp.msr[n,3]
        res1[1,m,n,1] <- r2.env[[1]][n]
        res2[1,m,n,1] <- r2.spa[[1]][n]
        res3[1,m,n,1] <- r2.env.spa[[1]][n]
        # CV
        res.env[1,m,n,2] <- vp.msr.cv[n,1]
        res.shared[1,m,n,2] <- vp.msr.cv[n,2]
        res.spa[1,m,n,2] <- vp.msr.cv[n,3]
        res1[1,m,n,2] <- r2.env[[2]][n]
        res2[1,m,n,2] <- r2.spa[[2]][n]
        res3[1,m,n,2] <- r2.env.spa[[2]][n]
      }
    })
  }
}


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# binary data

if(pars$data.type[task.id]=="binary"){
  
  mat.sp <- mat.sp[,which(apply(mat.sp,2,sum)<=(nrow(mat.sp)-1))] # exclude species present everywhere
  
  # Gather methods names
  meths <- c("R2.CCA", "R2.dbRDA", "R2.RDA", "R2.MRM", 
             "R2.GLM", "R2.GAM",
             "R2.BRT", "R2.BRT", "R2.UniRndForest", "R2.UniRndForest", 
             "R2.MVRndForest", "R2.MVRegTree")
  
  # Number of differrent R2 indices calculated
  n.r2 <- 3
  
  # Prepare output object
  res.env <- array(NA, c(1, length(meths), n.r2, 2))
  res.shared <- array(NA, c(1, length(meths), n.r2, 2))
  res.spa <- array(NA, c(1, length(meths), n.r2, 2))
  res1 <- array(NA, c(1, length(meths), n.r2, 2))
  res2 <- array(NA, c(1, length(meths), n.r2, 2))
  res3 <- array(NA, c(1, length(meths), n.r2, 2))
  
  # Loop to apply each method to given simulated dataset
  for(m in 1:length(meths)){
    
    try({
      
      # Constrained ordination
      # varpart for methods using MEM
      if(m %in% 1:2){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),binary=TRUE))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),binary=TRUE))
      }
      
      # RDA
      if(m %in% 3 && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),binary=TRUE))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),binary=TRUE))
      }
      
      # RDA with quadratic effect of X (env)
      if(m %in% 3  && pars$resp[task.id]=="gaussian"){
        mat.sp.log <- log(mat.sp+1)
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp.log,list(env=mat.env),env.eff="quadratic",binary=TRUE))
        r2.spa <- do.call(meths[m],list(mat.sp.log,list(mem=mat.spa),env.eff="quadratic",binary=TRUE))
        r2.env.spa <- do.call(meths[m],list(mat.sp.log,list(env=mat.env,mem=mat.spa),env.eff="quadratic",binary=TRUE))
      }
      
      # Distance based regression
      if(m %in% c(4)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),binary=TRUE))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),binary=TRUE))
      }
      
      # GLM
      # R2 is adjusted 
      if(m %in% 5  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="binomial"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),env.eff="linear",family="binomial"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),env.eff="linear",family="binomial"))
      }
      
      if(m %in% 5  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="binomial"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),family="binomial"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),family="binomial"))
      }
      
      # GAM
      if(m %in% 6  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="binomial"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="binomial"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="binomial"))
      }
      
      if(m %in% 6  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="binomial"))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="binomial"))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="binomial"))
      }
      
      # Tree-based methods
      # varpart for methods using the xy coordinates
      if(m %in% c(7,9,11,12)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),binary=TRUE))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),binary=TRUE))
      }
      
      if(m %in% c(8)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),inter.depth=1,binary=TRUE))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),inter.depth=1,binary=TRUE))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),inter.depth=1,binary=TRUE))
      }
      
      if(m %in% c(10)){
        # Calculate the R2 for each model
        r2.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),sampsize=0.2,binary=TRUE))
        r2.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),sampsize=0.2,binary=TRUE))
        r2.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),sampsize=0.2,binary=TRUE))
      }
      
      # variation partitioning
      vp.R2fad <- VarPart(r2.env[[1]][1],r2.spa[[1]][1],r2.env.spa[[1]][1])
      vp.R2efr <- VarPart(r2.env[[1]][2],r2.spa[[1]][2],r2.env.spa[[1]][2])
      vp.R2tju <- VarPart(r2.env[[1]][3],r2.spa[[1]][3],r2.env.spa[[1]][3])
      
      vp.R2fad.cv <- VarPart(r2.env[[2]][1],r2.spa[[2]][1],r2.env.spa[[2]][1])
      vp.R2efr.cv <- VarPart(r2.env[[2]][2],r2.spa[[2]][2],r2.env.spa[[2]][2])
      vp.R2tju.cv <- VarPart(r2.env[[2]][3],r2.spa[[2]][3],r2.env.spa[[2]][3])
      
      # Store results
      vp.msr <- rbind(vp.R2fad,vp.R2efr,vp.R2tju)
      vp.msr.cv <- rbind(vp.R2fad.cv,vp.R2efr.cv,vp.R2tju.cv)
      for(n in 1:n.r2){
        res.env[1,m,n,1] <- vp.msr[n,1]
        res.shared[1,m,n,1] <- vp.msr[n,2]
        res.spa[1,m,n,1] <- vp.msr[n,3]
        res1[1,m,n,1] <- r2.env[[1]][n]
        res2[1,m,n,1] <- r2.spa[[1]][n]
        res3[1,m,n,1] <- r2.env.spa[[1]][n]
        # CV
        res.env[1,m,n,2] <- vp.msr.cv[n,1]
        res.shared[1,m,n,2] <- vp.msr.cv[n,2]
        res.spa[1,m,n,2] <- vp.msr.cv[n,3]
        res1[1,m,n,2] <- r2.env[[2]][n]
        res2[1,m,n,2] <- r2.spa[[2]][n]
        res3[1,m,n,2] <- r2.env.spa[[2]][n]
      }
    })
  }
}


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# Output object
lout <- list(res.env, res.shared, res.spa, res1, res2, res3)


