##############
## iVarPart ##
##############

# setwd("~/Documents/Papers/iDiv/Paper_process_inference/EVE")

# Passing arguments to R from command lines
args = commandArgs(trailingOnly=TRUE)
output.file <- args[1]

# try to get SGE_TASK_ID from submit script, otherwise fall back to 1
task.id = as.integer(Sys.getenv("SGE_TASK_ID", "1"))

# load functions
source("Methods_functions.R")

# Load parameters
load("pars.Rdata")

# Load data 
# (list with each element consisting of the object "lout" produced in "Data_simulation.R")
load("sim_data.Rdata")
simn <- sim[[task.id]]
mat.sp <- simn$mat.sp
mat.sp <- mat.sp[,which(apply(mat.sp,2,sum)>0)] # exclude possibly extinct species
mat.env <- simn$mat.env
mat.spa <- simn$mat.spa
mat.xy <- simn$mat.xy

# Gather methods names
meths <- c("RDA", "MRM", 
           "GLM", "GAM", "GAM",
           "BRT", "BRT", "BRT", "UniRndForest", "UniRndForest", "UniRndForest", 
           "MVRndForest", "MVRndForest", "MVRndForest", "MVRegTree", "MVRegTree", "MVRegTree")

# Output results
lout <- list()

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# normal data

if(pars$data.type[task.id]=="normal"){
  
  # Loop to apply each method to given simulated dataset
  for(m in 1:length(meths)){
    
    try({
      
      # RDA
      if(m %in% 1 && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
      }
      # RDA with quadratic effect of X (env)
      if(m %in% 1  && pars$resp[task.id]=="gaussian"){
        mat.sp.log <- log(mat.sp+1)
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp.log,list(env=mat.env),env.eff="quadratic"))
        preds.spa <- do.call(meths[m],list(mat.sp.log,list(mem=mat.spa),env.eff="quadratic"))
        preds.env.spa <- do.call(meths[m],list(mat.sp.log,list(env=mat.env,mem=mat.spa),env.eff="quadratic"))
      }
      
      # Distance based regression (MRM)
      if(m %in% c(2)){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
      }
      
      # GLM 
      if(m %in% 3  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="gaussian"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),env.eff="linear",family="gaussian"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),env.eff="linear",family="gaussian"))
      }
      if(m %in% 3  && pars$resp[task.id]=="gaussian"){
        mat.sp.log <- log(mat.sp+1)
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp.log,list(env=mat.env),family="gaussian"))
        preds.spa <- do.call(meths[m],list(mat.sp.log,list(mem=mat.spa),family="gaussian"))
        preds.env.spa <- do.call(meths[m],list(mat.sp.log,list(env=mat.env,mem=mat.spa),family="gaussian"))
      }
      
      # GAM
      if(m %in% 4  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="gaussian"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="gaussian"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="gaussian"))
      }
      if(m %in% 4  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="gaussian",k.env=3))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="gaussian",k.env=3))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="gaussian",k.env=3))
      }
      if(m %in% 5  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="gaussian",k.spa=10,fx=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="gaussian",k.spa=10,fx=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="gaussian",k.spa=10,fx=TRUE))
      }
      if(m %in% 5  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="gaussian",k.env=3,k.spa=10,fx=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="gaussian",k.env=3,k.spa=10,fx=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="gaussian",k.env=3,k.spa=10,fx=TRUE))
      }
      
      # Tree-based methods
      
      if(m %in% c(6,9,12,15)){
        # BRT, UniRndForest, MVRndForest, MVRegTree with default settings
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
      }
      
      if(m %in% c(7)){
        # BRT
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),shrink=0.1))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),shrink=0.1))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),shrink=0.1))
      }
      
      if(m %in% c(8)){
        # BRT with MEMs
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),inter.depth=1))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa),inter.depth=1))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa),inter.depth=1))
      }
      
      if(m %in% c(10,13)){
        # UniRndForest, MVRndForest 
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),sampsize=0.2))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),sampsize=0.2))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),sampsize=0.2))
      }
      
      if(m %in% c(11,14)){
        # UniRndForest, MVRndForest with MEMs
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa)))
      }
      
      if(m %in% c(16)){
        # MVRegTree
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),CV=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),CV=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),CV=TRUE))
      }
      
      if(m %in% c(17)){
        # MVRegTree with MEMs
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa)))
      }
      
      # Store results
      preds.env <- as.matrix(preds.env)
      preds.spa <- as.matrix(preds.spa)
      preds.env.spa <- as.matrix(preds.env.spa)
      res <- array(NA, c(nrow(preds.env), ncol(preds.env), 3))
      res[,,1] <- preds.env
      res[,,2] <- preds.spa
      res[,,3] <- preds.env.spa
      lout[[m]] <- res
    })
  }
}



#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# count data

if(pars$data.type[task.id]=="counts"){
  
  # Loop to apply each method to given simulated dataset
  for(m in 1:length(meths)){
    
    try({
      
      # RDA
      if(m %in% 1 && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa)))
      }
      # RDA with quadratic effect of X (env)
      if(m %in% 1  && pars$resp[task.id]=="gaussian"){
        mat.sp.log <- log(mat.sp+1)
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp.log,list(env=mat.env),env.eff="quadratic"))
        preds.spa <- do.call(meths[m],list(mat.sp.log,list(mem=mat.spa),env.eff="quadratic"))
        preds.env.spa <- do.call(meths[m],list(mat.sp.log,list(env=mat.env,mem=mat.spa),env.eff="quadratic"))
      }
      
      # Distance based regression (MRM)
      if(m %in% c(2)){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
      }
      
      # GLM
      if(m %in% 3  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="quasipoisson"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),env.eff="linear",family="quasipoisson"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),env.eff="linear",family="quasipoisson"))
      }
      if(m %in% 3  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="quasipoisson"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),family="quasipoisson"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),family="quasipoisson"))
      }
      
      # GAM
      if(m %in% 4  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="poisson"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="poisson"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="poisson"))
      }
      if(m %in% 4  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="poisson",k.env=3))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="poisson",k.env=3))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="poisson",k.env=3))
      }
      if(m %in% 5  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="poisson",k.spa=10,fx=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="poisson",k.spa=10,fx=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="poisson",k.spa=10,fx=TRUE))
      }
      if(m %in% 5  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="poisson",k.env=3,k.spa=10,fx=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="poisson",k.env=3,k.spa=10,fx=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="poisson",k.env=3,k.spa=10,fx=TRUE))
      }
      
      # Tree-based methods
      
      if(m %in% c(6)){
        # BRT
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),distr="poisson"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),distr="poisson"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),distr="poisson"))
      }
      
      if(m %in% c(7)){
        # BRT
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),shrink=0.1,distr="poisson"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),shrink=0.1,distr="poisson"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),shrink=0.1,distr="poisson"))
      }
      
      if(m %in% c(8)){
        # BRT with MEMs
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),inter.depth=1,distr="poisson"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa),inter.depth=1,distr="poisson"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa),inter.depth=1,distr="poisson"))
      }
      
      if(m %in% c(9,12,15)){
        # UniRndForest, MVRndForest, MVRegTree with default settings
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy)))
      }
      
      if(m %in% c(10,13)){
        # UniRndForest, MVRndForest
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),sampsize=0.2))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),sampsize=0.2))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),sampsize=0.2))
      }
      
      if(m %in% c(11,14)){
        # UniRndForest, MVRndForest with MEMs
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa)))
      }
      
      if(m %in% c(16)){
        # MVRegTree
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),CV=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),CV=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),CV=TRUE))
      }
      
      if(m %in% c(17)){
        # MVRegTree with MEMs
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env)))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa)))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa)))
      }
      
      # Store results
      preds.env <- as.matrix(preds.env)
      preds.spa <- as.matrix(preds.spa)
      preds.env.spa <- as.matrix(preds.env.spa)
      res <- array(NA, c(nrow(preds.env), ncol(preds.env), 3))
      res[,,1] <- preds.env
      res[,,2] <- preds.spa
      res[,,3] <- preds.env.spa
      lout[[m]] <- res
    })
  }
}


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# binary data

if(pars$data.type[task.id]=="binary"){
  
  mat.sp <- mat.sp[,which(apply(mat.sp,2,sum)<=(nrow(mat.sp)-1))] # exclude species present everywhere
  
  # Loop to apply each method to given simulated dataset
  for(m in 1:length(meths)){
    
    try({
      
      # RDA
      if(m %in% 1 && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),binary=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),binary=TRUE))
      }
      # RDA with quadratic effect of X (env)
      if(m %in% 1  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="quadratic",binary=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),env.eff="quadratic",binary=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),env.eff="quadratic",binary=TRUE))
      }
      
      # Distance based regression (MRM)
      if(m %in% c(2)){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),binary=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),binary=TRUE))
      }
      
      # GLM
      # R2 is adjusted 
      if(m %in% 3  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="binomial"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),env.eff="linear",family="binomial"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),env.eff="linear",family="binomial"))
      }
      if(m %in% 3  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="binomial"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(mem=mat.spa),family="binomial"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,mem=mat.spa),family="binomial"))
      }
      
      # GAM
      if(m %in% 4  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="binomial"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="binomial"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="binomial"))
      }
      if(m %in% 4  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="binomial",k.env=3))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="binomial",k.env=3))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="binomial",k.env=3))
      }
      if(m %in% 5  && pars$resp[task.id]=="linear"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),env.eff="linear",family="binomial",k.spa=10,fx=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),env.eff="linear",family="binomial",k.spa=10,fx=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),env.eff="linear",family="binomial",k.spa=10,fx=TRUE))
      }
      if(m %in% 5  && pars$resp[task.id]=="gaussian"){
        # Calculate the R2 for each model
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),family="binomial",k.env=3,k.spa=10,fx=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),family="binomial",k.env=3,k.spa=10,fx=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),family="binomial",k.env=3,k.spa=10,fx=TRUE))
      }
      
      # Tree-based methods
      # methods using the xy coordinates
      
      if(m %in% c(6)){
        # BRT
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),distr="bernoulli"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),distr="bernoulli"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),distr="bernoulli"))
      }
      
      if(m %in% c(7)){
        # BRT
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),shrink=0.1,distr="bernoulli"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),shrink=0.1,distr="bernoulli"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),shrink=0.1,distr="bernoulli"))
      }
      
      if(m %in% c(8)){
        # BRT
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),inter.depth=1,distr="bernoulli"))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa),inter.depth=1,distr="bernoulli"))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa),inter.depth=1,distr="bernoulli"))
      }
      
      if(m %in% c(9,12,15)){
        # UniRndForest, MVRndForest, MVRegTree with default settings
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),binary=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),binary=TRUE))
      }
      
      if(m %in% c(10,13)){
        # UniRndForest, MVRndForest
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),sampsize=0.2,binary=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),sampsize=0.2,binary=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),sampsize=0.2,binary=TRUE))
      }
      
      if(m %in% c(11,14)){
        # UniRndForest, MVRndForest with MEMs
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa),binary=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa),binary=TRUE))
      }
      
      if(m %in% c(16)){
        # MVRegTree
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),CV=TRUE,binary=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.xy),CV=TRUE,binary=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.xy),CV=TRUE,binary=TRUE))
      }
      
      if(m %in% c(17)){
        # MVRegTree with MEMs
        preds.env <- do.call(meths[m],list(mat.sp,list(env=mat.env),binary=TRUE))
        preds.spa <- do.call(meths[m],list(mat.sp,list(xy=mat.spa),binary=TRUE))
        preds.env.spa <- do.call(meths[m],list(mat.sp,list(env=mat.env,xy=mat.spa),binary=TRUE))
      }
      
      # Store results
      preds.env <- as.matrix(preds.env)
      preds.spa <- as.matrix(preds.spa)
      preds.env.spa <- as.matrix(preds.env.spa)
      res <- array(NA, c(nrow(preds.env), ncol(preds.env), 3))
      res[,,1] <- preds.env
      res[,,2] <- preds.spa
      res[,,3] <- preds.env.spa
      lout[[m]] <- res
    })
  }
}


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

# Save output to a .Rdata file
save(lout,file=output.file)




