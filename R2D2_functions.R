###########################################
# Supporting Information
# "Partitioning environment and space in species-by-site matrices: a comparison of methods for community ecology and macroecology"
# Duarte S. Viana, Petr Keil, Alienor Jeliazkov
###########################################


# Functions to calculate R2 (explained squared error) and explained (Poisson) deviance,
# as well as to perform variation partitioning
# See README.txt file
# Here we present functions that can be used to to measure "explained variation"of species abundances.
# The ARGUMENTS are always the same (except for the variation partitioning function):
# Y.obs - a matrix or data.frame of OBSERVED species abundances,
#         where species are as columns and sites are rows
# Y.pred - a matrix or data.frame of PREDICTED species abundances,
#         where species are columns and sites are rows


# ------------------------------------------------------------------------------
# 0. Supporting function: Poisson deviance
# Arguments:
# Y - the observations
# mu - the predicted means (of Poisson distribution)
Poisson.deviance <- function(Y, mu)
{
  2*sum(ifelse(Y == 0, 
               -(Y-mu), # when the observation is 0
               Y*log(Y/mu) - (Y-mu))) # else
}


# ------------------------------------------------------------------------------
# 1. Explained deviance calculated for each species, then averaged
D2.Poisson <- function(Y.obs, Y.pred)
{
  Y <- as.matrix(Y.obs)
  mu <- as.matrix(Y.pred) # we assume that the predicted values are the Poisson mu
  
  resid.devs <- null.devs <- numeric(ncol(Y))
  for(i in 1:ncol(Y))
  {
    resid.devs[i] <- Poisson.deviance(Y[,i], mu[,i])
    null.devs[i]  <- Poisson.deviance(Y[,i], mean(Y[,i]))
  }
  
  D2s <- 1 - resid.devs/null.devs
  D2s[D2s < 0] <- 0 # change negative values to 0
  D2s[is.nan(D2s)] <-0
  D2 <- mean(D2s)
  
  return(D2)
}


# ------------------------------------------------------------------------------
# 2. R2 function for multivariate responses using the "trace" statistic
R2.multivar <- function(Y.obs, Y.pred)
{
  Y <- as.matrix(Y.obs)
  pred.mat <- as.matrix(Y.pred)
  
  S <- cov(pred.mat)
  eigenS <- eigen(S)
  kc <- length(which(eigenS$values>0.00000001))
  ev <- eigenS$values[1:kc]
  trace <- sum(diag(cov(as.matrix(Y))))
  r2 <- sum(ev/trace)
  return(r2)
  
}


# ------------------------------------------------------------------------------
# 3. R2 classic function averaged over species
R2.classic <- function(Y.obs, Y.pred)
{
  Y <- as.matrix(Y.obs)
  pred.mat <- as.matrix(Y.pred)
  r2s <- numeric(ncol(Y))
  for(i in 1:ncol(Y)){
    null.sum.sq <- sum((Y[,i] - mean(Y[,i]))^2)
    resid.sum.sq <- sum((Y[,i] - pred.mat[,i])^2)
    r2s[i] <- 1 - resid.sum.sq/null.sum.sq
  }
  
  r2s[r2s < 0] <- 0 # change negative values to 0
  r2s[is.nan(r2s)] <- 0
  r2 <- mean(r2s)
  return(r2)
}


# ------------------------------------------------------------------------------
# 4. log-R2 function averaged over species
R2.log <- function(Y.obs, Y.pred)
{
  Y <- as.matrix(Y.obs)
  pred.mat <- as.matrix(Y.pred)
  r2s <- numeric(ncol(Y))
  for(i in 1:ncol(Y)){
    null.sum.sq <- sum((log(Y[,i]+1) - log(mean(Y[,i])+1))^2)
    resid.sum.sq <- sum((log(Y[,i]+1) - log(pred.mat[,i]+1))^2)
    r2s[i] <- 1 - resid.sum.sq/null.sum.sq
  }
  
  r2s[r2s < 0] <- 0 # change negative values to 0
  r2s[is.nan(r2s)] <-0
  r2 <- mean(r2s)
  return(r2)
}



# ------------------------------------------------------------------------------
# Function that calculates all of the functions above

R2D2 <- function(Y.obs, Y.pred)
{
  res <-  c(R2.classic = R2.classic(Y.obs, Y.pred), 
            R2.log = R2.log(Y.obs, Y.pred),
            R2.multivar  = R2.multivar(Y.obs, Y.pred), 
            R2.McFadden = D2.Poisson(Y.obs, Y.pred))
  res[res < 0] <- 0 # change negative values to 0
  return(res)
}


# ------------------------------------------------------------------------------
# Function to perform variation partitioning (Clappe et al. 2018)
# Code adapted from the adespatial package
# Dray's github for the MSR (Clappe et al.)
# https://github.com/sdray/adespatial/blob/master/R/msr.varipart.R
VarPartClap <- function(ab.ini, bc.ini, bc.ini.adj, abc.ini, msr.ab, msr.abc){
  a.ini <- abc.ini-bc.ini
  ab.adj <- 1 - (1 - ab.ini) / (1-mean(msr.ab))
  a.adj <- 1 - (1 - a.ini) / (1 - mean(msr.abc - bc.ini))
  b.adj <- ab.adj - a.adj
  c.adj <- bc.ini.adj - b.adj
  d.adj <- 1 - (a.adj + b.adj + c.adj)
  R2.adj.msr <- c(a.adj, b.adj, c.adj, d.adj)
  return(R2.adj.msr)
}



################################################################################
# Measures for BINARY outcomes and probability predictions

# ------------------------------------------------------------------------------
# MCFADDEN PSEUDO R2 (gives identical result to explained deviance)
# ------------------------------------------------------------------------------

R2.McFadden <- function(Y.obs, Y.pred)
{
  L.base <- sum(dbinom(x=Y.obs, size = 1, prob = sum(Y.obs/length(Y.obs)), log=TRUE))
  L.full <- sum(dbinom(x=Y.obs, size = 1, prob = Y.pred, log=TRUE))
  1 - (L.full/L.base)
}

# ------------------------------------------------------------------------------
# EFFRON PSEUDO R2
#  - This implementation is taken from  the PseudoR2 function from package DescTools
# ------------------------------------------------------------------------------

R2.Efron <- function(Y.obs, Y.pred)
{
  1 - (sum((Y.obs - Y.pred)^2))/(sum((Y.obs - mean(Y.obs))^2))
}

# ------------------------------------------------------------------------------
# TJUR PSEUDO R2
# This implementation is taken from  the PseudoR2 function from package DescTools
# ------------------------------------------------------------------------------

R2.Tjur <- function(Y.obs, Y.pred)
{
  unname(diff(tapply(Y.pred, Y.obs, 
                     mean, na.rm = TRUE)))
}

# ------------------------------------------------------------------------------
# PROPORTION OF EXPLAINED BERNOULLI DEVIANCE (identical results to McFadden pseudo R2)
# ------------------------------------------------------------------------------

D2.binom <- function(Y.obs, Y.pred) 
{
  null.D = -2*sum(dbinom(x=Y.obs, size = 1, prob = sum(Y.obs/length(Y.obs)), log=TRUE))
  pred.D = -2*sum(dbinom(x=Y.obs, size = 1, prob = Y.pred, log=TRUE))
  1 - pred.D/null.D
}

# ------------------------------------------------------------------------------
# Function that calculates the functions above for all species and averages them
# across the species

R2D2.binom <- function(Y.obs, Y.pred)
{
  Y.obs <- as.matrix(Y.obs)
  Y.pred <- as.matrix(Y.pred) 
  
  res <- list()
  for(i in 1:ncol(Y.obs))
  {
    res[[i]] <- c(R2.McFadden = D2.binom(Y.obs[,i], Y.pred[,i]), 
                  R2.Efron    = R2.Efron(Y.obs[,i], Y.pred[,i]), 
                  R2.Tjur     = R2.Tjur(Y.obs[,i], Y.pred[,i]))
  }
  res <- do.call(rbind, res)
  res[res < 0] <- 0 # change negative values to 0
  res[is.nan(res)] <- 0

  res <- colMeans(res, na.rm = TRUE)
  return(res)
}

# Testing:
# Y.pred <- R2.MVRegTree(Y = mat.sp, X = mat.env, binary = TRUE)
# Y.obs <- mat.sp
# Y.obs[Y.obs > 1] <- 1
# plot(Y.pred, Y.obs); abline(a=0, b=1)
# R2D2.binom(Y.obs, Y.pred)


# ------------------------------------------------------------------------------

# testing of the binomial measures
# 1. simulated data
#      x <- sort(rnorm(30))
#      a <- 0
#      b <- 5
#      p <- exp(a + b*x)/(1 + exp(a + b*x))
#      plot(x, p, type = "l")
#      y <- rbinom(n=length(x), size= 1, prob = p)
#      points(x, y)
#      m1 <- glm(y ~ x, family="binomial")
#      y.pred <- predict(m1, type = "response")
#      lines(x, y.pred, col="red")

# 2. test the new functions against the implementations in other packages
#
#      PseudoR2(m1, which="all")
#      
#      summary(m1)
#      D2.binom(y, y.pred)
#      
#      PseudoR2(m1, which="Effron")
#      R2.Efron(y, y.pred)
#      
#      PseudoR2(m1, which="Tjur")
#      R2.Tjur(y, y.pred)
#      
#      PseudoR2(m1, which="McFadden")
#      R2.McFadden(y, y.pred)
