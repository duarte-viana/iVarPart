###########################################
# Supporting Information
# "Disentangling spatial and environmental effects: flexible methods for community ecology and macroecology"
# Duarte S. Viana, Petr Keil, Alienor Jeliazkov
###########################################

# Code for simulating the data
# See README.txt file
# This is the code for one simnulation corresponding to the scenario identified as "task.id",
# which corresponds to a given row number of the table of parameters

# ------------------------------------------------------------------------------

# Passing arguments to R from command lines
args = commandArgs(trailingOnly=TRUE)
output.file <- args[1]

# try to get SLURM_ARRAY_TASK_ID  from submit script, otherwise fall back to 1
task.id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

# Load libraries
library(sp)
library(gstat)
library(adespatial)
library(plotrix)
library(ltm)

# Load parameters
load("pars.Rdata")

# Function to simulate variables with a desired correlation rho
complement <- function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}

# Number of response variables (e.g. species) to simulate
R <- 20

# Define type of data ("normal", "counts" or "binary")
dat <- pars$data.type[task.id] 

# Define lattice
grid.side <- sqrt(pars$Nsamples[task.id])
grid.size<-grid.side^2
xy <- expand.grid(0:(grid.side-1), 0:(grid.side-1))
names(xy) <- c("x","y")
gridded(xy) = ~x+y
coords<-coordinates(xy)
coords2<-as.data.frame(coords)

# Spatial structure of X
# E.g. autocorrelated environment
g.dummy.e <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                   model=vgm(psill=0.025,model="Exp",range=grid.side*pars$rangeX[task.id]), nmax=20)
# make 1 simulation based on the stat object
yy.e <- predict(g.dummy.e, newdata=xy, nsim=1)
# Transform variable from normal to uniform
ecum.e<-ecdf(yy.e@data[,1])
ecum.yy.e<-ecum.e(yy.e@data[,1])
yy.e@data[,1]<-ecum.yy.e
yy.data.e<-yy.e@data[,1]
E<-yy.data.e


# Spatial structure of W
g.dummy.s <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                   model=vgm(psill=0.025,model=pars$spa.mod[task.id],range=grid.side*pars$rangeW[task.id]), nmax=20)
# make 1 simulation based on the stat object
yy.s <- predict(g.dummy.s, newdata=xy, nsim=1)
# Transform variable from normal to uniform
ecum.e<-ecdf(yy.s@data[,1])
ecum.yy.s<-ecum.e(yy.s@data[,1])
yy.s@data[,1]<-ecum.yy.s
yy.data.s<-yy.s@data[,1]
S<-data.frame(S=yy.data.s)

# Generate X (Gaussian response to E)

if(pars$resp[task.id]=="linear"){
  opt <- rnorm(R,10,2)
  X <- matrix(ncol=R,nrow=grid.size)
  for(i in 1:R){
    respi <- opt[i] * E
    X[,i] <- rescale(respi, c(0,20)) # rescale to have the same range as in Gaussian responses
  } 
}

if(pars$resp[task.id]=="gaussian"){
  opt <- seq(0.05,0.95,length.out=R)
  X <- matrix(ncol=R,nrow=grid.size)
  for(i in 1:R) X[,i] <- 20*exp((-(opt[i]-E)^2)/pars$sd.niche[task.id])
}

# Generate W
W <- as.matrix(S[,rep(1,R)])*20
# Force W to have the same frequency distribution as X
#W <- X
#for(i in 1:R) W[,i] <- X[match(rank(S@data[,i],ties.method="random"),rank(X[,i],ties.method="random")),i]




# Calculate species abundance - build species matrix
betaX <- pars$betaX[task.id]
betaW <- 1-betaX

if(dat=="normal"){
  Y.hat <- as.matrix(betaX*X + betaW*W)
  Y <- as.matrix(Y.hat + rnorm(grid.size, 0, 0.1*20))
  for(i in 1:R){
    minY <- abs(min(Y[,i]))
    Y[,i] <- Y[,i] + minY # make sure abundance is always non-negative
    Y.hat[,i] <- Y.hat[,i] + minY # scale the fitted accordingly
  } 
}

if(dat=="counts"){
  Y.hat <- betaX*X + betaW*W
  Y <- sapply(1:R, function(x) rpois(nrow(Y.hat),Y.hat[,x]))
}

if(dat=="binary"){
  X <- scale(X,scale=F)
  W <- scale(W,scale=F)
  Y.hat0 <- betaX*X + betaW*W
  Y.hat <- sapply(1:R, function(x) 1/(1+exp(-Y.hat0[,x]))) # pass through an inv-logit function
  Y <- sapply(1:R, function(x) rbinom(nrow(Y.hat),1,Y.hat[,x]))
}

# Add random X ("environmental") covariates?
# These are generated to be orthogonal to the response
E <- data.frame(E=E)
NXrand <- pars$NXrand[task.id]
if(NXrand>0){
  for(i in 1:NXrand){
    for(j in 1:R){
      Erand <- complement(Y.hat[,j], rho=0)
      Erand <- rescale(Erand, c(0,1))
      E[,i+1] <- Erand
    }
  }
}

# Undersampling
if(pars$underN[task.id]=="yes"){
  wr <- sample(1:nrow(Y),50)
  Y <- Y[wr,]
  Y.hat <- Y.hat[wr,]
  E <- E[wr,]
  X <- X[wr,]
  W <- W[wr,]
  coords2 <- coords2[wr,]
}

# MEM
xy.pcnm <- dbmem(coords2,MEM.autocor = "positive",silent=TRUE)
all.pcnm <- as.data.frame(xy.pcnm)

# Store relevant data
lout <- list(mat.sp=Y, mat.sp.true=Y.hat, X=X, W=W, mat.env=E, mat.xy=coords2, mat.spa=all.pcnm)

# the "lout" object was stored in a list called "sim", containing all simulations 

# Save output to a .Rdata file
save(lout,file=output.file)
