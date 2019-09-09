###########################################
# Supporting Information
# "Partitioning environment and space in species-by-site matrices: a comparison of methods for community ecology and macroecology"
# Duarte S. Viana, Petr Keil, Alienor Jeliazkov
###########################################

# Code for simulating the data
# See README.txt file
# This is the code for one simnulation corresponding to the scenario identified as "task.id",
# which corresponds to a given row number of the table of parameters

# ------------------------------------------------------------------------------

# Load libraries
library(sp)
library(gstat)
library(rgeos)
library(adespatial)
library(plotrix)

# Load parameters
load("pars_GB.Rdata")

# Define type of data ("counts" or "binary")
dat <- "binary" 

# Define lattice
edge <- floor(pars$Nsamples[task.id]/2)
grid.side <- pars$Nsamples[task.id] + edge*2
grid.size<-grid.side^2
xy <- expand.grid(0:(grid.side-1), 0:(grid.side-1))
names(xy) <- c("x","y")
gridded(xy) = ~x+y
coords<-coordinates(xy)
coords2<-as.data.frame(coords)

# Autocorrelated environment
type <- "mosaic"
if(type=="mosaic") g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                                    model=vgm(psill=0.025,model="Exp",range=pars$range[task.id]), nmax=20)
# trend in x and y
if(type=="gradient") g.dummy <- gstat(formula=z~1+x+y, locations=~x+y, dummy=T, 
                                      beta=c(1,0.01,0.01), model=vgm(psill=0.025, range=pars$range[task.id], model='Exp'), nmax=20)
# make 1 simulation based on the stat object
yy <- predict(g.dummy, newdata=xy, nsim=1)
# Transform environmental variable from normal to uniform
ecum1<-ecdf(yy@data[,1])
ecum.yy<-ecum1(yy@data[,1])
yy@data[,1]<-ecum.yy
yy.data<-yy@data[,1]
env1<-yy.data

# Randomly distributed environment
yy2<-SpatialPixelsDataFrame(coordinates(xy),as.data.frame(runif(grid.size,0,1)))
env2 <- yy2@data[,1]

# Distribute species - calculate abundance for each species in each cell according to environment
# Niche component
R <- pars$R[task.id]
opt <- seq(0.05,0.95,length.out=R)
mat1 <- matrix(ncol=R,nrow=grid.size)
for(i in 1:R){
  fe1 <- exp((-(opt[i]-env1)^2)/pars$sd.niche[task.id])
  fe2 <- exp((-(opt[i]-env2)^2)/pars$sd.niche[task.id])
  mat1[,i] <- 0.75*fe1 + 0.25*fe2
}

# Dispersal component
matd <- matrix(ncol=R,nrow=grid.size)
for(i in 1:R) matd[,i] <- rlnorm(nrow(matd),1,0.5)
neighb <- edge
mat2 <- matd
for(i in 1:grid.size){
  # Calculate the dispersal contribution form a neighbourhood of neighb-units distance  
  wini <- 1:grid.size
  wini <- c(wini[i],wini[-i])
  cmat <- rbind(coords[i,],coords[-i,])
  dists <- as.matrix(dist(cmat))[,1]
  wnei <- which(dists<=neighb)[-1]
  jnei <- wini[wnei]
  dists <- dists[wnei]
  for(j in 1:R){
    mat2[i,j] <- sum(matd[jnei,j]*(exp((dists*-1)/0.5)))
  }
}

# crop the grid and data to avoid edge effects
crop<-which(coords2$x>edge&coords2$x<=(grid.side-edge)&coords2$y>edge&coords2$y<=(grid.side-edge))
xyf <- xy[crop,]
xyd <- coords2[crop,]
env1 <- env1[crop]
env2 <- env2[crop]
mat.env <- cbind(env1,env2)
mat.env0 <- mat1[crop,]
mat.spa <- mat2[crop,]

# Scale components
mat.spa <- mat.spa*sum(mat.env0)/sum(mat.spa)
mat.spa <- scale(mat.spa)
mat.env0 <- scale(mat.env0)

# Calculate species abundance - build species matrix
ws <- pars$ws[task.id]
we <- 1-ws

if(dat=="counts"){
  mat.spa <- rescale(mat.spa,c(0,20))
  mat.env0 <- rescale(mat.env0,c(0,20))
  mean.abund <- as.vector(ws*mat.spa + we*mat.env0)
  final.abund <- rpois(length(mean.abund),mean.abund)
  mat.sp <- matrix(final.abund,nrow=nrow(mat.env0),ncol=ncol(mat.env0))
}

if(dat=="binary"){
  mean.abund <- as.vector(ws*mat.spa + we*mat.env0)
  pr <- 1/(1+exp(-mean.abund)) # pass through an inv-logit function
  final.pr <- rbinom(length(mean.abund),1,pr)
  mat.sp <- matrix(final.pr,nrow=nrow(mat.env0),ncol=ncol(mat.env0))
}

# PCNM
xy.pcnm <- dbmem(xyd,MEM.autocor = "positive",silent=TRUE)
all.pcnm <- as.data.frame(xy.pcnm)

# Store relevant data
lout <- list(mat.sp=mat.sp,pred.env=mat.env0,pred.spa=mat.spa,mat.env=mat.env,mat.xy=xyd,mat.spa=all.pcnm, 
             mem.listw=attr(xy.pcnm, "listw"))
# the "lout" object was stored in a list called "sim", containing all simulations 
 

