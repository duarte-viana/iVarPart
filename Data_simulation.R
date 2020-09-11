###########################################
# Supporting Information
# "Partitioning environment and space in site-by-species matrices: a comparison of methods for community ecology and macroecology"
# Duarte S. Viana, Petr Keil, Alienor Jeliazkov
###########################################

# Code for simulating the data
# See README.txt file
# This is the code for one simnulation corresponding to the scenario identified as "task.id",
# which corresponds to a given row number of the table of parameters

#####################################################################
#####################################################################

# Data simulation for Exercise 1

# Load libraries
library(sp)
library(gstat)
library(plotrix)

# Load parameters
load("pars1.Rdata")

# Define parameters
dat <- pars$data.type[task.id] 
resp <- pars$resp[task.id]
N <- pars$N[task.id]

# Simulation of environmental effect
#===================================

if(pars$effect[task.id]=="env"){
  # predictor
  X <- runif(N,0,1)
  
  # response
  if(resp=="linear") Yfit <- 20*X
  if(resp=="gaussian") Yfit <- 20*exp((-(0.5-X)^2)/pars$sd.gaussian[task.id])
}


# Simulation of spatial effect
#===================================

if(pars$effect[task.id]=="spa"){
  # Define lattice
  grid.side <- sqrt(pars$N[task.id])
  grid.size<-grid.side^2
  xy <- expand.grid(0:(grid.side-1), 0:(grid.side-1))
  names(xy) <- c("x","y")
  gridded(xy) = ~x+y
  coords<-coordinates(xy)
  coords2<-as.data.frame(coords)
  
  # response with autocorrelated spatial structure
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                   model=vgm(psill=0.025,model="Exp",range=pars$range[task.id]), nmax=20)
  Yfit <- unlist(predict(g.dummy, newdata=xy, nsim=1)@data)
  Yfit <- rescale(Yfit, c(0,20))
  
  # predictors
  X <- coords2
}

#============================

# Add error term to response (depending on the type of error distribution and/or magnitude)
if(pars$effect[task.id]=="env"){
  if(dat=="normal"){
    Y <- Yfit + rnorm(N,0,20*pars$error[task.id])
    minY <- abs(min(Y))
    if(resp=="linear"){
      Y <- Y + minY # avoid negative values and log
      Yfit <- Yfit + minY # rescale to match Y
    }
    if(resp=="gaussian"){
      Y <- log(Y + minY + 1) # avoid negative values and log
      Yfit <- log(Yfit + minY + 1) # rescale to match Y
    }
  } 
}
if(pars$effect[task.id]=="spa" && dat=="normal") Y <- Yfit + rnorm(N,0,20*pars$error[task.id])
if(dat=="counts") Y <- rpois(N,Yfit)
if(dat=="binary"){
  Y.centred <- scale(Yfit,scale=F)
  Yfit <- 1/(1+exp(-Y.centred)) # pass through an inv-logit function
  Y <- rbinom(N, 1, Yfit)
}

# Store relevant data
sim[[task.id]] <- list(Yfit=Yfit, Y=Y, X=X)


#####################################################################
#####################################################################

# Data simulation for Exercise 2

# Load libraries
library(sp)
library(gstat)
library(adespatial)
library(ltm)

# Load parameters
load("pars_VP.Rdata")

# Number of species to simulate
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

# Spatial structure of the environment
# Autocorrelated environment
g.dummy.e <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                   model=vgm(psill=0.025,model="Exp",range=pars$range[task.id]), nmax=20)
# make 1 simulation based on the stat object
yy.e <- predict(g.dummy.e, newdata=xy, nsim=1)
# Transform environmental variable from normal to uniform
ecum.e<-ecdf(yy.e@data[,1])
ecum.yy.e<-ecum.e(yy.e@data[,1])
yy.e@data[,1]<-ecum.yy.e
yy.data.e<-yy.e@data[,1]
E<-yy.data.e

# Spatial structure of W
g.dummy.s <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                   model=vgm(psill=0.025,model="Exp",range=pars$range[task.id]), nmax=20)
S <- predict(g.dummy.s, newdata=xy, nsim=R)

# Generate X (Gaussian response to E)

if(pars$resp[task.id]=="linear"){
  opt <- runif(R,0,20)
  X <- matrix(ncol=R,nrow=grid.size)
  for(i in 1:R) X[,i] <- opt[i] * E
}

if(pars$resp[task.id]=="gaussian"){
  opt <- seq(0.05,0.95,length.out=R)
  X <- matrix(ncol=R,nrow=grid.size)
  for(i in 1:R) X[,i] <- 20*exp((-(opt[i]-E)^2)/pars$sd.niche[task.id])
}

# Generate W
W <- X
for(i in 1:R) W[,i] <- X[match(rank(S@data[,i],ties.method="random"),rank(X[,i],ties.method="random")),i]


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


# Calculate true fractions of explained variation
# Before random sampling
ab.sp <- c()
bc.sp <- c()
for(i in 1:ncol(Y)){
  ab.sp[i] <- cor(X[,i],Y.hat[,i])^2
  bc.sp[i] <- cor(W[,i],Y.hat[,i])^2
}
ab <- mean(ab.sp)
bc <- mean(bc.sp)
abc <- 1
a <- abc-bc
b <- ab+bc-abc
c <- abc-ab
d <- 1-abc
# if(a<0) a <- 0
# if(b<0) b <- 0
# if(c<0) c <- 0
vp.true.determ <- c(a,b,c,ab,bc,abc)


# After random sampling
if(dat=="normal" | dat=="counts"){
  ab.sp <- c()
  bc.sp <- c()
  abc.sp <- c()
  for(i in 1:ncol(Y)){
    ab.sp[i] <- cor(X[,i],Y[,i])^2
    bc.sp[i] <- cor(W[,i],Y[,i])^2
    abc.sp[i] <- cor(Y.hat[,i],Y[,i])^2
  }
  ab <- mean(ab.sp)
  bc <- mean(bc.sp)
  abc <- mean(abc.sp)
  a <- abc-bc
  b <- ab+bc-abc
  c <- abc-ab
  d <- 1-abc
  # if(a<0) a <- 0
  # if(b<0) b <- 0
  # if(c<0) c <- 0
  vp.true <- c(a,b,c,ab,bc,abc)
}

if(dat=="binary"){
  ab.sp <- c()
  bc.sp <- c()
  abc.sp <- c()
  for(i in 1:ncol(Y)){
    ab.sp[i] <- biserial.cor(X[,i],Y[,i],level=2)^2
    bc.sp[i] <- biserial.cor(W[,i],Y[,i],level=2)^2
    abc.sp[i] <- biserial.cor(Y.hat[,i],Y[,i],level=2)^2
  }
  ab <- mean(ab.sp)
  bc <- mean(bc.sp)
  abc <- mean(abc.sp)
  a <- abc-bc
  b <- ab+bc-abc
  c <- abc-ab
  d <- 1-abc
  # if(a<0) a <- 0
  # if(b<0) b <- 0
  # if(c<0) c <- 0
  vp.true <- c(a,b,c,ab,bc,abc)
}


# Undersampling
if(pars$underN[task.id]=="yes"){
  wr <- sample(1:nrow(Y),50)
  Y <- Y[wr,]
  Y.hat <- Y.hat[wr,]
  E <- E[wr]
  X <- X[wr,]
  W <- W[wr,]
  coords2 <- coords2[wr,]
}

# MEM
xy.pcnm <- dbmem(coords2,MEM.autocor = "positive",silent=TRUE)
all.pcnm <- as.data.frame(xy.pcnm)

# Store relevant data
sim[[task.id]] <- list(mat.sp=Y, mat.sp.true=Y.hat, X=X, W=W, mat.env=data.frame(E=E),
             mat.xy=coords2, mat.spa=all.pcnm, vp.true.determ=vp.true.determ, vp.true=vp.true)


