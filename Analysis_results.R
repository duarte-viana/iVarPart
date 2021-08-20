###########################################
# Supporting Information
# "Disentangling spatial and environmental effects: flexible methods for community ecology and macroecology"
# Duarte S. Viana, Petr Keil, Alienor Jeliazkov
###########################################

# Performance of the different methods

library(vegan)
library(ltm)

# Load parameters
load("pars.Rdata")
# Load simulated data 
load("sim_data.Rdata")
# Load predictions from the different methods
load("Model_preds.Rdata")

meths <- c("RDA", "RDA", "MRM", 
           "GLM", "GAM", "GAM",
           "BRT", "BRT", "BRT", "UniRndForest", "UniRndForest", "UniRndForest", 
           "MVRndForest", "MVRndForest", "MVRndForest", "MVRegTree", "MVRegTree", "MVRegTree")

mnam <- c("RDA-raw", "RDA-Hel", "MRM", 
          "GLM", "GAM-kdef", "GAM-k10",
          "BRT-lr0.01", "BRT-lr0.1", "BRT-MEM", "UniRF-SS0", "UniRF-SS20", "UniRF-MEM",
          "MVRF-SS0", "MVRF-SS20", "MVRF-MEM", "MVRT-noCV", "MVRT-CV", "MVRT-MEM")

meth.family <- character(length(mnam))
meth.family[1:2] <- "brown3"
meth.family[3] <- "deepskyblue3"
meth.family[4] <- "violet"
meth.family[5:6] <- "darkseagreen"
meth.family[7:18] <- "goldenrod1"
names(meth.family) <- mnam

#######################################################################################
#######################################################################################

# Exercise 1: model fitting

res.true <- array(NA, c(nrow(pars), length(meths), 3))
res.obs <- array(NA, c(nrow(pars), length(meths), 3))
res.sim <- matrix(nrow=nrow(pars), ncol=length(meths))
for(i in 1:nrow(pars)){
  preds0 <- lpreds[[i]]
  Ytrue0 <- sim[[i]]$mat.sp.true
  Yobs0 <- sim[[i]]$mat.sp
  sp.in <- which(apply(Yobs0,2,sum)>0)
  Yobs0 <- Yobs0[,sp.in]
  Ytrue0 <- Ytrue0[,sp.in]
  if(pars$data.type[i]=="binary"){
    sp.in <- which(apply(Yobs0,2,sum)<=(nrow(Ytrue0)-1)) # Exclude species present everywhere
    Yobs0 <- Yobs0[,sp.in]
    Ytrue0 <- Ytrue0[,sp.in]
  }
  
  for(m in 1:length(meths)){
    preds <- preds0[[m]]
    Yobs <- Yobs0
    Ytrue <- Ytrue0
    if(!is.null(preds)){
      if(meths[m]=="MRM"){
        if(pars$data.type[i]=="binary"){binary=TRUE} else{binary=FALSE}
        Yobs <- as.matrix(as.numeric(vegdist(Yobs0, method="bray", binary = binary)))
        Ytrue <- as.matrix(as.numeric(vegdist(Ytrue0, method="bray", binary = binary)))
        wna <- which(is.na(Yobs))
        if(length(wna)>0){
          Yobs <- Yobs[-wna,,drop=FALSE]
          Ytrue <- Ytrue[-wna,,drop=FALSE]
        }
      }
      if(meths[m]=="RDA"){
        if((pars$data.type[i]=="normal"  || pars$data.type[i]=="counts") && pars$resp[i]=="gaussian"){
          Ytrue <- as.matrix(scale(log(Ytrue0+1)))
          Yobs <- as.matrix(scale(log(Yobs0+1)))
        } 
        # Ytrue <- as.matrix(decostand(as.matrix(Ytrue), "hellinger"))
        # Ytrue <- as.matrix(scale(Ytrue,scale=FALSE))
        # Yobs <- as.matrix(decostand(as.matrix(Yobs), "hellinger"))
        # Yobs <- as.matrix(scale(Yobs,scale=FALSE))
      }
      if(meths[m]=="GLM" && pars$data.type[i]=="normal" && pars$resp[i]=="gaussian"){
        Ytrue <- log(Ytrue0+1)
        Yobs <- log(Yobs0+1)
      }
      
      cor.true1 <- c()
      cor.true2 <- c()
      cor.true3 <- c()
      cor.obs1 <- c()
      cor.obs2 <- c()
      cor.obs3 <- c()
      cor.sim <- c()
      for(j in 1:ncol(as.matrix(preds[,,3]))){
        cor.true1 <- c(cor.true1, cor(preds[,j,1], Ytrue[,j], method="pearson"))
        cor.true2 <- c(cor.true2, cor(preds[,j,2], Ytrue[,j], method="pearson"))
        cor.true3 <- c(cor.true3, cor(preds[,j,3], Ytrue[,j], method="pearson"))
        if(pars$data.type[i]!="binary"){
          cor.obs1 <- c(cor.obs1, cor(preds[,j,1], Yobs[,j], method="pearson"))
          cor.obs2 <- c(cor.obs2, cor(preds[,j,2], Yobs[,j], method="pearson"))
          cor.obs3 <- c(cor.obs3, cor(preds[,j,3], Yobs[,j], method="pearson"))
          cor.sim <- c(cor.sim, cor(Ytrue[,j], Yobs[,j], method="pearson"))
        } 
        if(pars$data.type[i]=="binary" && meths[m] %in% c("RDA","MRM")){
          cor.obs1 <- c(cor.obs1, cor(preds[,j,1], Yobs[,j], method="pearson"))
          cor.obs2 <- c(cor.obs2, cor(preds[,j,2], Yobs[,j], method="pearson"))
          cor.obs3 <- c(cor.obs3, cor(preds[,j,3], Yobs[,j], method="pearson"))
          cor.sim <- c(cor.sim, cor(Ytrue[,j], Yobs[,j], method="pearson"))
        }
        if(pars$data.type[i]=="binary" && !(meths[m] %in% c("RDA","MRM"))){
          cor.obs1 <- c(cor.obs1, biserial.cor(preds[,j,1], Yobs[,j], level=2))
          cor.obs2 <- c(cor.obs2, biserial.cor(preds[,j,2], Yobs[,j], level=2))
          cor.obs3 <- c(cor.obs3, biserial.cor(preds[,j,3], Yobs[,j], level=2))
          cor.sim <- c(cor.sim, biserial.cor(Ytrue[,j], Yobs[,j], level=2))
        }
      }
      res.true[i,m,1] <- mean(cor.true1,na.rm=T)
      res.true[i,m,2] <- mean(cor.true2,na.rm=T)
      res.true[i,m,3] <- mean(cor.true3,na.rm=T)
      res.obs[i,m,1] <- mean(cor.obs1,na.rm=T)
      res.obs[i,m,2] <- mean(cor.obs2,na.rm=T)
      res.obs[i,m,3] <- mean(cor.obs3,na.rm=T)
      res.sim[i,m] <- mean(cor.sim,na.rm=T)
    }
  }
}

# RMSE is not fair because RDA, GLM (normal, log-trasformed), and MRM are on different scales


#####################################################################################
#####################################################################################

# Exercise 2: variation partitioning
source("R2D2_functions.R")

# Variation partitioning
vp.sim <- array(NA, c(nrow(pars), length(meths), 3)) 
vp.est <- array(NA, c(nrow(pars), length(meths), 3)) # for (pseudo) R2
vp.est2 <- array(NA, c(nrow(pars), length(meths), 3)) # for correlation-based R2

for(i in 1:nrow(pars)){
  preds0 <- lpreds[[i]]
  Ytrue0 <- sim[[i]]$mat.sp.true
  Yobs0 <- sim[[i]]$mat.sp
  sp.in <- which(apply(Yobs0,2,sum)>0)
  Yobs0 <- Yobs0[,sp.in]
  Ytrue0 <- Ytrue0[,sp.in]
  Xtrue0 <- sim[[i]]$X
  Xtrue0 <- Xtrue0[,sp.in]
  Wtrue0 <- sim[[i]]$W
  Wtrue0 <- Wtrue0[,sp.in]
  if(pars$data.type[i]=="binary"){
    sp.in <- which(apply(Yobs0,2,sum)<=(nrow(Yobs0)-1)) # Exclude species present everywhere
    Yobs0 <- Yobs0[,sp.in]
    Ytrue0 <- Ytrue0[,sp.in]
    Xtrue0 <- Xtrue0[,sp.in]
    Wtrue0 <- Wtrue0[,sp.in]
    Xtrue0 <- sapply(1:ncol(Xtrue0), function(x) 1/(1+exp(-Xtrue0[,x]))) # pass through an inv-logit function
    Wtrue0 <- sapply(1:ncol(Wtrue0), function(x) 1/(1+exp(-Wtrue0[,x]))) # pass through an inv-logit function
  }
  
  # calculate number of model parameters
  if(pars$resp[i]=="linear") PX <- ncol(sim[[i]]$mat.env)
  if(pars$resp[i]=="gaussian") PX <- ncol(sim[[i]]$mat.env)*2
  PW <- ncol(sim[[i]]$mat.spa)
  
  for(m in 1:length(meths)){
    preds <- preds0[[m]]
    Yobs <- Yobs0
    Ytrue <- Ytrue0
    Xtrue <- Xtrue0
    Wtrue <- Wtrue0
    if(!is.null(preds)){
      #try({
      if(meths[m]=="MRM"){
        if(pars$data.type[i]=="binary"){binary=TRUE} else{binary=FALSE}
        Yobs.tr <- na.exclude(as.matrix(as.numeric(vegdist(Yobs0, method="bray", binary = binary))))
      }
      if(mnam[m]=="RDA-raw"){
        if((pars$data.type[i]=="normal"  || pars$data.type[i]=="counts") && pars$resp[i]=="gaussian"){
          Yobs.tr <- as.matrix(scale(log(Yobs0+1)))
          # Ytrue <- log(Ytrue0+1)
          # Xtrue <- log(Xtrue0+1)
          # Wtrue <- log(Wtrue0+1)
        } 
      }
      if(mnam[m]=="RDA-Hel"){
        if((pars$data.type[i]=="normal"  || pars$data.type[i]=="counts") && pars$resp[i]=="gaussian"){
          Yobs.tr <- log(Yobs0+1)
          Yobs.tr <- as.matrix(decostand(as.matrix(Yobs.tr), "hellinger"))
          Yobs.tr <- as.matrix(scale(Yobs.tr,scale=FALSE))
          # Ytrue <- log(Ytrue0+1)
          # Xtrue <- log(Xtrue0+1)
          # Wtrue <- log(Wtrue0+1)
        } else{
          Yobs.tr <- as.matrix(decostand(as.matrix(Yobs0), "hellinger"))
          Yobs.tr <- as.matrix(scale(Yobs.tr,scale=FALSE))
        }
      }
      if(meths[m]=="GLM" && pars$data.type[i]=="normal" && pars$resp[i]=="gaussian"){
        Yobs.tr <- log(Yobs0+1)
        # Ytrue <- log(Ytrue0+1)
        # Xtrue <- log(Xtrue0+1)
        # Wtrue <- log(Wtrue0+1)
      }
      
      if(pars$data.type[i]=="normal"){
        r2.X.sim <- cor2(Yobs, Xtrue, method="pearson")
        r2.W.sim <- cor2(Yobs, Wtrue, method="pearson")
        r2.XW.sim <- cor2(Yobs, Ytrue, method="pearson")
        # r2.X.sim <- R2.classic(Yobs, Xtrue)
        # r2.W.sim <- R2.classic(Yobs, Wtrue)
        # r2.XW.sim <- R2.classic(Yobs, Ytrue)
        if(meths[m]=="MRM"){
          r2.X.est <- R2.classic(Yobs.tr, preds[,,1])
          r2.W.est <- R2.classic(Yobs.tr, preds[,,2])
          r2.XW.est <- R2.classic(Yobs.tr, preds[,,3])
          r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="pearson")
          r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="pearson")
          r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="pearson")
        } 
        if(mnam[m]=="RDA-raw"){
          if((pars$data.type[i]=="normal"  || pars$data.type[i]=="counts") && pars$resp[i]=="gaussian"){
            r2.X.est <- R2.classic(Yobs.tr, preds[,,1], adjust=TRUE, N=nrow(Yobs.tr), P=PX)
            r2.W.est <- R2.classic(Yobs.tr, preds[,,2], adjust=TRUE, N=nrow(Yobs.tr), P=PW)
            r2.XW.est <- R2.classic(Yobs.tr, preds[,,3], adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
            r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX)
            r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PW)
            r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
          } else{
            r2.X.est <- R2.classic(Yobs, preds[,,1], adjust=TRUE, N=nrow(Yobs), P=PX)
            r2.W.est <- R2.classic(Yobs, preds[,,2], adjust=TRUE, N=nrow(Yobs), P=PW)
            r2.XW.est <- R2.classic(Yobs, preds[,,3], adjust=TRUE, N=nrow(Yobs), P=PX+PW)
            r2.X.est2 <- cor2(Yobs, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PX)
            r2.W.est2 <- cor2(Yobs, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PW)
            r2.XW.est2 <- cor2(Yobs, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PX+PW)
          }
        }
        if(mnam[m]=="RDA-Hel"){
          r2.X.est <- R2.classic(Yobs.tr, preds[,,1], adjust=TRUE, N=nrow(Yobs.tr), P=PX)
          r2.W.est <- R2.classic(Yobs.tr, preds[,,2], adjust=TRUE, N=nrow(Yobs.tr), P=PW)
          r2.XW.est <- R2.classic(Yobs.tr, preds[,,3], adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
          r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX)
          r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PW)
          r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
        }
        if(meths[m]=="GLM"){
          if(pars$data.type[i]=="normal" && pars$resp[i]=="gaussian"){
            r2.X.est <- R2.classic(Yobs.tr, preds[,,1], adjust=TRUE, N=nrow(Yobs.tr), P=PX)
            r2.W.est <- R2.classic(Yobs.tr, preds[,,2], adjust=TRUE, N=nrow(Yobs.tr), P=PW)
            r2.XW.est <- R2.classic(Yobs.tr, preds[,,3], adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
            r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX)
            r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PW)
            r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
          } else{
            r2.X.est <- R2.classic(Yobs, preds[,,1], adjust=TRUE, N=nrow(Yobs), P=PX)
            r2.W.est <- R2.classic(Yobs, preds[,,2], adjust=TRUE, N=nrow(Yobs), P=PW)
            r2.XW.est <- R2.classic(Yobs, preds[,,3], adjust=TRUE, N=nrow(Yobs), P=PX+PW)
            r2.X.est2 <- cor2(Yobs, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PX)
            r2.W.est2 <- cor2(Yobs, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PW)
            r2.XW.est2 <- cor2(Yobs, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PX+PW)
          }
        }
        if(meths[m] %in% c("GAM","BRT","UniRndForest","MVRndForest","MVRegTree")){
          r2.X.est <- R2.classic(Yobs, preds[,,1])
          r2.W.est <- R2.classic(Yobs, preds[,,2])
          r2.XW.est <- R2.classic(Yobs, preds[,,3])
          r2.X.est2 <- cor2(Yobs, preds[,,1], method="pearson")
          r2.W.est2 <- cor2(Yobs, preds[,,2], method="pearson")
          r2.XW.est2 <- cor2(Yobs, preds[,,3], method="pearson")
        }
      }
      
      if(pars$data.type[i]=="counts"){
        r2.X.sim <- cor2(Yobs, Xtrue, method="spearman")
        r2.W.sim <- cor2(Yobs, Wtrue, method="spearman")
        r2.XW.sim <- cor2(Yobs, Ytrue, method="spearman")
        if(meths[m]=="MRM"){
          # r2.X.sim <- R2.classic(Yobs, Xtrue)
          # r2.W.sim <- R2.classic(Yobs, Wtrue)
          # r2.XW.sim <- R2.classic(Yobs, Ytrue)
          r2.X.est <- R2.classic(Yobs.tr, preds[,,1])
          r2.W.est <- R2.classic(Yobs.tr, preds[,,2])
          r2.XW.est <- R2.classic(Yobs.tr, preds[,,3])
          r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="spearman")
          r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="spearman")
          r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="spearman")
        }
        if(mnam[m]=="RDA-raw"){
          # r2.X.sim <- R2.classic(Yobs, Xtrue)
          # r2.W.sim <- R2.classic(Yobs, Wtrue)
          # r2.XW.sim <- R2.classic(Yobs, Ytrue)
          if((pars$data.type[i]=="normal"  || pars$data.type[i]=="counts") && pars$resp[i]=="gaussian"){
            r2.X.est <- R2.classic(Yobs.tr, preds[,,1], adjust=TRUE, N=nrow(Yobs.tr), P=PX)
            r2.W.est <- R2.classic(Yobs.tr, preds[,,2], adjust=TRUE, N=nrow(Yobs.tr), P=PW)
            r2.XW.est <- R2.classic(Yobs.tr, preds[,,3], adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
            r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX)
            r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PW)
            r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
          } else{
            r2.X.est <- R2.classic(Yobs, preds[,,1], adjust=TRUE, N=nrow(Yobs), P=PX)
            r2.W.est <- R2.classic(Yobs, preds[,,2], adjust=TRUE, N=nrow(Yobs), P=PW)
            r2.XW.est <- R2.classic(Yobs, preds[,,3], adjust=TRUE, N=nrow(Yobs), P=PX+PW)
            r2.X.est2 <- cor2(Yobs, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PX)
            r2.W.est2 <- cor2(Yobs, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PW)
            r2.XW.est2 <- cor2(Yobs, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PX+PW)
          }
        }
        if(mnam[m]=="RDA-Hel"){
          # r2.X.sim <- R2.classic(Yobs, Xtrue)
          # r2.W.sim <- R2.classic(Yobs, Wtrue)
          # r2.XW.sim <- R2.classic(Yobs, Ytrue)
          r2.X.est <- R2.classic(Yobs.tr, preds[,,1], adjust=TRUE, N=nrow(Yobs.tr), P=PX)
          r2.W.est <- R2.classic(Yobs.tr, preds[,,2], adjust=TRUE, N=nrow(Yobs.tr), P=PW)
          r2.XW.est <- R2.classic(Yobs.tr, preds[,,3], adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
          r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="spearman", adjust=TRUE, N=nrow(Yobs.tr), P=PX)
          r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="spearman", adjust=TRUE, N=nrow(Yobs.tr), P=PW)
          r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="spearman", adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
        }
        if(meths[m]=="GLM"){
          # r2.X.sim <- D2.Poisson(Yobs, Xtrue)
          # r2.W.sim <- D2.Poisson(Yobs, Wtrue)
          # r2.XW.sim <- D2.Poisson(Yobs, Ytrue)
          r2.X.est <- D2.Poisson(Yobs, preds[,,1], adjust=TRUE, N=nrow(Yobs), P=PX)
          r2.W.est <- D2.Poisson(Yobs, preds[,,2], adjust=TRUE, N=nrow(Yobs), P=PW)
          r2.XW.est <- D2.Poisson(Yobs, preds[,,3], adjust=TRUE, N=nrow(Yobs), P=PX+PW)
          r2.X.est2 <- cor2(Yobs, preds[,,1], method="spearman", adjust=TRUE, N=nrow(Yobs), P=PX)
          r2.W.est2 <- cor2(Yobs, preds[,,2], method="spearman", adjust=TRUE, N=nrow(Yobs), P=PW)
          r2.XW.est2 <- cor2(Yobs, preds[,,3], method="spearman", adjust=TRUE, N=nrow(Yobs), P=PX+PW)
        }
        if(meths[m] %in% c("UniRndForest","MVRndForest","MVRegTree")){
          # r2.X.sim <- R2.classic(Yobs, Xtrue)
          # r2.W.sim <- R2.classic(Yobs, Wtrue)
          # r2.XW.sim <- R2.classic(Yobs, Ytrue)
          r2.X.est <- R2.classic(Yobs, preds[,,1])
          r2.W.est <- R2.classic(Yobs, preds[,,2])
          r2.XW.est <- R2.classic(Yobs, preds[,,3])
          r2.X.est2 <- cor2(Yobs, preds[,,1], method="spearman")
          r2.W.est2 <- cor2(Yobs, preds[,,2], method="spearman")
          r2.XW.est2 <- cor2(Yobs, preds[,,3], method="spearman")
        } 
        if(meths[m] %in% c("GAM","BRT")){
          # r2.X.sim <- D2.Poisson(Yobs, Xtrue)
          # r2.W.sim <- D2.Poisson(Yobs, Wtrue)
          # r2.XW.sim <- D2.Poisson(Yobs, Ytrue)
          r2.X.est <- D2.Poisson(Yobs, preds[,,1])
          r2.W.est <- D2.Poisson(Yobs, preds[,,2])
          r2.XW.est <- D2.Poisson(Yobs, preds[,,3])
          r2.X.est2 <- cor2(Yobs, preds[,,1], method="spearman")
          r2.W.est2 <- cor2(Yobs, preds[,,2], method="spearman")
          r2.XW.est2 <- cor2(Yobs, preds[,,3], method="spearman")
        }
      }
      
      if(pars$data.type[i]=="binary"){
        r2.X.sim <- cor2(Yobs, Xtrue, method="binary")
        r2.W.sim <- cor2(Yobs, Wtrue, method="binary")
        r2.XW.sim <- cor2(Yobs, Ytrue, method="binary")
        if(meths[m]=="MRM"){
          # r2.X.sim <- R2.classic(Yobs, Xtrue)
          # r2.W.sim <- R2.classic(Yobs, Wtrue)
          # r2.XW.sim <- R2.classic(Yobs, Ytrue)
          r2.X.est <- R2.classic(Yobs.tr, preds[,,1])
          r2.W.est <- R2.classic(Yobs.tr, preds[,,2])
          r2.XW.est <- R2.classic(Yobs.tr, preds[,,3])
          r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="pearson")
          r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="pearson")
          r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="pearson")
        }
        if(mnam[m]=="RDA-raw"){
          # r2.X.sim <- R2.classic(Yobs, Xtrue)
          # r2.W.sim <- R2.classic(Yobs, Wtrue)
          # r2.XW.sim <- R2.classic(Yobs, Ytrue)
          r2.X.est <- R2.classic(Yobs, preds[,,1], adjust=TRUE, N=nrow(Yobs), P=PX)
          r2.W.est <- R2.classic(Yobs, preds[,,2], adjust=TRUE, N=nrow(Yobs), P=PW)
          r2.XW.est <- R2.classic(Yobs, preds[,,3], adjust=TRUE, N=nrow(Yobs), P=PX+PW)
          r2.X.est2 <- cor2(Yobs, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PX)
          r2.W.est2 <- cor2(Yobs, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PW)
          r2.XW.est2 <- cor2(Yobs, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs), P=PX+PW)
        }
        if(mnam[m]=="RDA-Hel"){
          # r2.X.sim <- R2.classic(Yobs, Xtrue)
          # r2.W.sim <- R2.classic(Yobs, Wtrue)
          # r2.XW.sim <- R2.classic(Yobs, Ytrue)
          r2.X.est <- R2.classic(Yobs.tr, preds[,,1], adjust=TRUE, N=nrow(Yobs.tr), P=PX)
          r2.W.est <- R2.classic(Yobs.tr, preds[,,2], adjust=TRUE, N=nrow(Yobs.tr), P=PW)
          r2.XW.est <- R2.classic(Yobs.tr, preds[,,3], adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
          r2.X.est2 <- cor2(Yobs.tr, preds[,,1], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX)
          r2.W.est2 <- cor2(Yobs.tr, preds[,,2], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PW)
          r2.XW.est2 <- cor2(Yobs.tr, preds[,,3], method="pearson", adjust=TRUE, N=nrow(Yobs.tr), P=PX+PW)
        }
        if(meths[m]=="GLM"){
          # r2.X.sim <- R2.Tjur(Yobs, Xtrue)
          # r2.W.sim <- R2.Tjur(Yobs, Wtrue)
          # r2.XW.sim <- R2.Tjur(Yobs, Ytrue)
          r2.X.est <- R2.Tjur(Yobs, preds[,,1], adjust=TRUE, N=nrow(Yobs), P=PX)
          r2.W.est <- R2.Tjur(Yobs, preds[,,2], adjust=TRUE, N=nrow(Yobs), P=PW)
          r2.XW.est <- R2.Tjur(Yobs, preds[,,3], adjust=TRUE, N=nrow(Yobs), P=PX+PW)
          r2.X.est2 <- cor2(Yobs, preds[,,1], method="binary", adjust=TRUE, N=nrow(Yobs), P=PX)
          r2.W.est2 <- cor2(Yobs, preds[,,2], method="binary", adjust=TRUE, N=nrow(Yobs), P=PW)
          r2.XW.est2 <- cor2(Yobs, preds[,,3], method="binary", adjust=TRUE, N=nrow(Yobs), P=PX+PW)
        }
        if(meths[m] %in% c("GAM","BRT","UniRndForest","MVRndForest","MVRegTree")){
          # r2.X.sim <- R2.Tjur(Yobs, Xtrue)
          # r2.W.sim <- R2.Tjur(Yobs, Wtrue)
          # r2.XW.sim <- R2.Tjur(Yobs, Ytrue)
          r2.X.est <- R2.Tjur(Yobs, preds[,,1])
          r2.W.est <- R2.Tjur(Yobs, preds[,,2])
          r2.XW.est <- R2.Tjur(Yobs, preds[,,3])
          r2.X.est2 <- cor2(Yobs, preds[,,1], method="binary")
          r2.W.est2 <- cor2(Yobs, preds[,,2], method="binary")
          r2.XW.est2 <- cor2(Yobs, preds[,,3], method="binary")
        }
      }
      # Variation partitioning
      vp.sim[i,m,1] <- r2.X.sim
      vp.sim[i,m,2] <- r2.XW.sim-r2.X.sim
      vp.sim[i,m,3] <- r2.XW.sim
      vp.est[i,m,1] <- r2.X.est
      vp.est[i,m,2] <- r2.XW.est-r2.X.est
      vp.est[i,m,3] <- r2.XW.est
      vp.est2[i,m,1] <- r2.X.est2
      vp.est2[i,m,2] <- r2.XW.est2-r2.X.est2
      vp.est2[i,m,3] <- r2.XW.est2
      # Full variation partitioning
      # vpsim <- VarPart(r2.X.sim, r2.W.sim, r2.XW.sim)
      # ab.sim <- vpsim[1]+vpsim[2]
      # c.sim <- vpsim[3]
      # vp.sim[i,m,1] <- ab.sim
      # vp.sim[i,m,2] <- c.sim
      # vp.sim[i,m,3] <- r2.XW.sim
      # vpest <- VarPart(r2.X.est, r2.W.est, r2.XW.est)
      # ab.est <- vpest[1]+vpest[2]
      # c.est <- vpest[3]
      # vp.est[i,m,1] <- ab.est
      # vp.est[i,m,2] <- c.est
      # vp.est[i,m,3] <- r2.XW.est
      # vpest2 <- VarPart(r2.X.est2, r2.W.est2, r2.XW.est2)
      # ab.est2 <- vpest2[1]+vpest2[2]
      # c.est2 <- vpest2[3]
      # vp.est2[i,m,1] <- ab.est2
      # vp.est2[i,m,2] <- c.est2
      # vp.est2[i,m,3] <- r2.XW.est2
      #})
    }
  }
}

vp.sim[vp.sim<0] <- 0
vp.est[vp.est<0] <- 0
vp.est2[vp.est2<0] <- 0





