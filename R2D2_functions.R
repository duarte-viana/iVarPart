###########################################
# Supporting Information, Appendix S3
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
# 1. Explained deviance that uses a global mean
D2.glob.Poisson <- function(Y.obs, Y.pred)
{
  Y <- as.matrix(Y.obs)
  mu <- as.matrix(Y.pred) # we assume that the predicted values are the Poisson mu
  
  resid.dev <- Poisson.deviance(Y, mu)
  null.dev  <- Poisson.deviance(Y, mean(Y))
  D2 <- 1 - resid.dev/null.dev
  return(D2)
}


# ------------------------------------------------------------------------------
# 2. Explained deviance calculated for each species, then averaged
D2.avg.Poisson <- function(Y.obs, Y.pred)
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
# 3. R2 function for multivariate responses using the "trace" statistic
R2.glob.multivar <- function(Y.obs, Y.pred)
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
# 4. R2 function applying the "ordinary" R2 logic
R2.glob.classic <- function(Y.obs, Y.pred)
{
  Y <- as.matrix(Y.obs)
  pred.mat <- as.matrix(Y.pred)
  null.sum.sq <- sum((Y - mean(Y))^2)
  resid.sum.sq <- sum((Y - pred.mat)^2)
  r2 <- 1 - resid.sum.sq/null.sum.sq
  return(r2)
}


# ------------------------------------------------------------------------------
# 5. R2 function averaged over species
R2.avg.classic <- function(Y.obs, Y.pred)
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
# 6. log-R2 function averaged over species
R2.avg.log <- function(Y.obs, Y.pred)
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
  res <-  c(R2avr = R2.avg.classic(Y.obs, Y.pred), 
            R2com = R2.glob.classic(Y.obs, Y.pred), 
            R2mv  = R2.glob.multivar(Y.obs, Y.pred), 
            D2avr = D2.avg.Poisson(Y.obs, Y.pred), 
            D2com = D2.glob.Poisson(Y.obs, Y.pred),
            R2log = R2.avg.log(Y.obs, Y.pred))
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









