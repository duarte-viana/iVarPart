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


# load libraries and functions
library(adespatial)
library(mgcv)
library(randomForest)
library(gbm)

# Load parameters
load("pars1.Rdata")

# Load data 
load("sim1.Rdata")

# output object
preds.out <- list()

# Extract data
simn <- sim[[task.id]]
Y <- simn$Y
X <- simn$X
if(is.vector(X)) X <- data.frame(X=X)

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

if(pars$effect[task.id]=="env"){
  
  # LM
  if(pars$data.type[task.id]=="normal"){
    rm(m)
    glm.preds <- rep(NA, length(Y))
    glm.preds.cv <- rep(NA, length(Y))
    try({
      if(pars$resp[task.id]=="linear") m <- lm(Y ~ ., data=X)
      if(pars$resp[task.id]=="gaussian") m <- lm(Y ~ poly(X,2), data=X)
      glm.preds <- predict(m)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        if(pars$resp[task.id]=="linear") m <- lm(Y[-i] ~ ., data=X[-i,,drop=F])
        if(pars$resp[task.id]=="gaussian") m <- lm(Y[-i] ~ poly(X,2), data=X[-i,,drop=F])
        glm.preds.cv[i] <- predict(m, newdata=X[i,,drop=F])
      }
    })
  }
  
  # GLM Poisson
  if(pars$data.type[task.id]=="counts"){
    rm(m)
    glm.preds <- rep(NA, length(Y))
    glm.preds.cv <- rep(NA, length(Y))
    try({
      if(pars$resp[task.id]=="linear") m <- glm(Y ~ ., data=X, family="poisson")
      if(pars$resp[task.id]=="gaussian") m <- glm(Y ~ poly(X,2), data=X, family="poisson")
      glm.preds <- predict(m, type="response")
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        if(pars$resp[task.id]=="linear") m <- glm(Y[-i] ~ ., data=X[-i,,drop=F], family="poisson")
        if(pars$resp[task.id]=="gaussian") m <- glm(Y[-i] ~ poly(X,2), data=X[-i,,drop=F], family="poisson")
        glm.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response")
      }
    })
  }
  
  # GLM Binomial
  if(pars$data.type[task.id]=="binary"){
    rm(m)
    glm.preds <- rep(NA, length(Y))
    glm.preds.cv <- rep(NA, length(Y))
    try({
      if(pars$resp[task.id]=="linear") m <- glm(Y ~ ., data=X, family="binomial")
      if(pars$resp[task.id]=="gaussian") m <- glm(Y ~ poly(X,2), data=X, family="binomial")
      glm.preds <- predict(m, type="response")
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        if(pars$resp[task.id]=="linear") m <- glm(Y[-i] ~ ., data=X[-i,,drop=F], family="binomial")
        if(pars$resp[task.id]=="gaussian") m <- glm(Y[-i] ~ poly(X,2), data=X[-i,,drop=F], family="binomial")
        glm.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response")
      }
    })
  }
  
  # BRT Gaussian
  if(pars$data.type[task.id]=="normal"){
    rm(m)
    brt.preds <- rep(NA, length(Y))
    brt.preds.cv <- rep(NA, length(Y))
    try({
      m <- gbm(Y ~ X, data=X, distribution = "gaussian", n.minobsinnode = 5,
               shrinkage = 0.01, interaction.depth = 1, n.trees = 1000, cv.folds = 0)
      brt.preds <- predict(m, type="response", n.trees=1000)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gbm(Y[-i] ~ X, data=X[-i,,drop=F], distribution = "gaussian", n.minobsinnode = 5,
                 shrinkage = 0.01, interaction.depth = 1, n.trees = 1000, cv.folds = 0)
        brt.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response", n.trees=1000)
      }
    })
  }
  
  # BRT Poisson
  if(pars$data.type[task.id]=="counts"){
    rm(m)
    brt.preds <- rep(NA, length(Y))
    brt.preds.cv <- rep(NA, length(Y))
    try({
      m <- gbm(Y ~ X, data=X, distribution = "poisson", n.minobsinnode = 5,
               shrinkage = 0.01, interaction.depth = 1, n.trees = 1000, cv.folds = 0)
      brt.preds <- predict(m, type="response", n.trees=1000)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gbm(Y[-i] ~ X, data=X[-i,,drop=F], distribution = "poisson", n.minobsinnode = 5,
                 shrinkage = 0.01, interaction.depth = 1, n.trees = 1000, cv.folds = 0)
        brt.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response", n.trees=1000)
      }
    })
  }
  
  # BRT Binomial
  if(pars$data.type[task.id]=="binary"){
    rm(m)
    brt.preds <- rep(NA, length(Y))
    brt.preds.cv <- rep(NA, length(Y))
    try({
      m <- gbm(Y ~ X, data=X, distribution = "bernoulli", n.minobsinnode = 5,
               shrinkage = 0.01, interaction.depth = 1, n.trees = 1000, cv.folds = 0)
      brt.preds <- predict(m, type="response", n.trees=1000)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gbm(Y[-i] ~ X, data=X[-i,,drop=F], distribution = "bernoulli", n.minobsinnode = 5,
                 shrinkage = 0.01, interaction.depth = 1, n.trees = 1000, cv.folds = 0)
        brt.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response", n.trees=1000)
      }
    })
  }
  
  # Random Forest for normal and count data
  if(pars$data.type[task.id] %in%  c("normal","counts")){
    rm(m)
    rf.preds <- rep(NA, length(Y))
    rf.preds.cv <- rep(NA, length(Y))
    try({
      m <- randomForest(Y ~ X, data=X, ntree = 500, sampsize = length(Y)/5)
      rf.preds <- predict(m)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- randomForest(Y[-i] ~ X, data=X[-i,,drop=F], ntree = 500, sampsize = length(Y)/5)
        rf.preds.cv[i] <- predict(m, newdata=X[i,,drop=F])
      }
    })
  }
  
  # Random Forest for binary data
  if(pars$data.type[task.id]=="binary"){
    rm(m)
    rf.preds <- rep(NA, length(Y))
    rf.preds.cv <- rep(NA, length(Y))
    try({
      m <- randomForest(factor(Y) ~ X, data=X, ntree = 500, sampsize = length(Y)/5)
      rf.preds <- predict(m, type="prob")[,2]
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- randomForest(factor(Y[-i]) ~ X, data=X[-i,,drop=F], ntree = 500, sampsize = length(Y)/5)
        rf.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="prob")[,2]
      }
    })
  }
  # output 
  preds.out[[task.id]] <- list(glm.preds, glm.preds.cv, brt.preds, brt.preds.cv, rf.preds, rf.preds.cv)
}

#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------

if(pars$effect[task.id]=="spa"){
  
  mem <- dbmem(X,MEM.autocor = "positive",silent=TRUE) # Calculate MEMs
  Xmem <- as.data.frame(mem)
  
  # LM (using MEMs)
  if(pars$data.type[task.id]=="normal"){
    rm(m)
    glm.preds <- rep(NA, length(Y))
    glm.preds.cv <- rep(NA, length(Y))
    try({
      m <- lm(Y ~ ., data=Xmem)
      glm.preds <- predict(m)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- lm(Y[-i] ~ ., data=Xmem[-i,])
        glm.preds.cv[i] <- predict(m, newdata=Xmem[i,,drop=F])
      }
    })
  }
  
  # GLM Poisson (using MEMs)
  if(pars$data.type[task.id]=="counts"){
    rm(m)
    glm.preds <- rep(NA, length(Y))
    glm.preds.cv <- rep(NA, length(Y))
    try({
      m <- glm(Y ~ ., data=Xmem, family="poisson")
      glm.preds <- predict(m, type="response")
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- glm(Y[-i] ~ ., data=Xmem[-i,], family="poisson")
        glm.preds.cv[i] <- predict(m, newdata=Xmem[i,,drop=F], type="response")
      }
    })
  }
  
  # GLM Binomial (using MEMs)
  if(pars$data.type[task.id]=="binary"){
    rm(m)
    glm.preds <- rep(NA, length(Y))
    glm.preds.cv <- rep(NA, length(Y))
    try({
      m <- glm(Y ~ ., data=Xmem, family="binomial")
      glm.preds <- predict(m, type="response")
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- glm(Y[-i] ~ ., data=Xmem[-i,], family="binomial")
        glm.preds.cv[i] <- predict(m, newdata=Xmem[i,,drop=F], type="response")
      }
    })
  }
  
  # GAM Gaussian
  if(pars$data.type[task.id]=="normal"){
    rm(m)
    gam.preds <- rep(NA, length(Y))
    gam.preds.cv <- rep(NA, length(Y))
    try({
      m <- gam(Y ~ s(x,y,k=-1), data=X, family="gaussian", method="GCV.Cp")
      gam.preds <- predict(m, type="response")
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gam(Y[-i] ~ s(x,y,k=-1), data=X[-i,], family="gaussian", method="GCV.Cp")
        gam.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response")
      }
    })
  }
  
  # GAM Poisson
  if(pars$data.type[task.id]=="counts"){
    rm(m)
    gam.preds <- rep(NA, length(Y))
    gam.preds.cv <- rep(NA, length(Y))
    try({
      m <- gam(Y ~ s(x,y,k=-1), data=X, family="poisson", method="GCV.Cp")
      gam.preds <- predict(m, type="response")
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gam(Y[-i] ~ s(x,y,k=-1), data=X[-i,], family="poisson", method="GCV.Cp")
        gam.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response")
      }
    })
  }
  
  # GAM Binomial
  if(pars$data.type[task.id]=="binary"){
    rm(m)
    gam.preds <- rep(NA, length(Y))
    gam.preds.cv <- rep(NA, length(Y))
    try({
      m <- gam(Y ~ s(x,y,k=-1), data=X, family="binomial", method="GCV.Cp")
      gam.preds <- predict(m, type="response")
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gam(Y[-i] ~ s(x,y,k=-1), data=X[-i,], family="binomial", method="GCV.Cp")
        gam.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response")
      }
    })
  }
  
  # BRT Gaussian
  if(pars$data.type[task.id]=="normal"){
    rm(m)
    brt.preds <- rep(NA, length(Y))
    brt.preds.cv <- rep(NA, length(Y))
    try({
      m <- gbm(Y ~ ., data=X, distribution = "gaussian", n.minobsinnode = 5,
               shrinkage = 0.01, interaction.depth = 2, n.trees = 1000, cv.folds = 0)
      brt.preds <- predict(m, type="response", n.trees=1000)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gbm(Y[-i] ~ ., data=X[-i,], distribution = "gaussian", n.minobsinnode = 5,
                 shrinkage = 0.01, interaction.depth = 2, n.trees = 1000, cv.folds = 0)
        brt.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response", n.trees=1000)
      }
    })
  }
  
  # BRT Poisson
  if(pars$data.type[task.id]=="counts"){
    rm(m)
    brt.preds <- rep(NA, length(Y))
    brt.preds.cv <- rep(NA, length(Y))
    try({
      m <- gbm(Y ~ ., data=X, distribution = "poisson", n.minobsinnode = 5,
               shrinkage = 0.01, interaction.depth = 2, n.trees = 1000, cv.folds = 0)
      brt.preds <- predict(m, type="response", n.trees=1000)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gbm(Y[-i] ~ ., data=X[-i,], distribution = "poisson", n.minobsinnode = 5,
                 shrinkage = 0.01, interaction.depth = 2, n.trees = 1000, cv.folds = 0)
        brt.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response", n.trees=1000)
      }
    })
  }
  
  # BRT Binomial
  if(pars$data.type[task.id]=="binary"){
    rm(m)
    brt.preds <- rep(NA, length(Y))
    brt.preds.cv <- rep(NA, length(Y))
    try({
      m <- gbm(Y ~ ., data=X, distribution = "bernoulli", n.minobsinnode = 5,
               shrinkage = 0.01, interaction.depth = 2, n.trees = 1000, cv.folds = 0)
      brt.preds <- predict(m, type="response", n.trees=1000)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- gbm(Y[-i] ~ ., data=X[-i,], distribution = "bernoulli", n.minobsinnode = 5,
                 shrinkage = 0.01, interaction.depth = 2, n.trees = 1000, cv.folds = 0)
        brt.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="response", n.trees=1000)
      }
    })
  }
  
  # Random Forest for normal and count data
  if(pars$data.type[task.id] %in%  c("normal","counts")){
    rm(m)
    rf.preds <- rep(NA, length(Y))
    rf.preds.cv <- rep(NA, length(Y))
    try({
      m <- randomForest(Y ~ x + y, data=X, ntree = 500)
      rf.preds <- predict(m)
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- randomForest(Y[-i] ~ x + y, data=X[-i,], ntree = 500)
        rf.preds.cv[i] <- predict(m, newdata=X[i,,drop=F])
      }
    })
  }
  
  # Random Forest for binary data
  if(pars$data.type[task.id]=="binary"){
    rm(m)
    rf.preds <- rep(NA, length(Y))
    rf.preds.cv <- rep(NA, length(Y))
    try({
      m <- randomForest(factor(Y) ~ x + y, data=X, ntree = 500)
      rf.preds <- predict(m, type="prob")[,2]
      # LOO cross-valiadation
      for(i in 1:length(Y)){
        rm(m)
        m <- randomForest(factor(Y[-i]) ~ x + y, data=X[-i,], ntree = 500)
        rf.preds.cv[i] <- predict(m, newdata=X[i,,drop=F], type="prob")[,2]
      }
    })
  }
  
  # output 
  preds.out[[task.id]] <- list(glm.preds, glm.preds.cv, gam.preds, gam.preds.cv, 
                               brt.preds, brt.preds.cv, rf.preds, rf.preds.cv)
}

