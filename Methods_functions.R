###########################################
# Supporting Information
# "Disentangling spatial and environmental effects: flexible methods for community ecology and macroecology"
# Duarte S. Viana, Petr Keil, Alienor Jeliazkov
###########################################

# Functions for the different statistical methods
# See README.txt file
# In the following functions,
# Y is a species-by-sites matrix of abundance (columns are species and rows are sites)
# X is a list with one or more matrices of predictor variables (numeric variables): environment, x-y coordinates, MEMs



##############################################################################
# Constrained ordination
##############################################################################

# INSTALLATION OF THE PACKAGES

#install.packages(c("vegan", "stats"))

# ------------------------------------------------------------------------------
# RDA
# ------------------------------------------------------------------------------

# On Hellinger-transformed data
RDA <- function(Y, X, binary = FALSE, env.eff = "linear"){
  require(vegan)
  require(Matrix)
  
  if(ncol(Y)==1){
    Y.hel <- Y
  } else{
    Y.hel <- as.matrix(decostand(as.matrix(Y), "hellinger"))
  }
  
  env <- X$env
  if(!is.null(env) && env.eff=="quadratic"){
    n.env <- ncol(X$env)
    env <- as.matrix(sparse.model.matrix(as.formula(paste("~poly(X$env[,",n.env,"],2)", collapse = " + "))))[,-1]
    X$env <- env
  }
  
  Y.stand <- as.matrix(scale(Y.hel,scale=FALSE))
  X0 <- as.data.frame(do.call("cbind", X))
  Y.pred <- Y; Y.pred[] <- NA
  for(i in 1:ncol(Y)){
    lmi <- lm(Y.stand[,i]~., data=X0)
    Y.pred[,i] <- predict(lmi)
  }
  
  return(Y.pred)
}



# ------------------------------------------------------------------------------
# dbRDA
# ------------------------------------------------------------------------------

## On raw abundance data
dbRDA <- function(Y, X, binary = FALSE){
  require(vegan)
  pcoa <- cmdscale(vegdist(Y, method = "bray", binary = binary), 
                   k=(NROW(Y) - 1), list.=TRUE, eig=TRUE, add=TRUE)
  Y.pcoa <- pcoa$points
  
  X0 <- as.data.frame(do.call("cbind", X))
  Y.pred <- Y.pcoa; Y.pred[] <- NA
  for(i in 1:ncol(Y.pcoa)){
    lmi <- lm(Y.pcoa[,i]~., data=X0)
    Y.pred[,i] <- predict(lmi)
  }
  
  return(Y.pred)
}



##############################################################################
# Distance-based regression
##############################################################################

# ------------------------------------------------------------------------------
# MRM
# ------------------------------------------------------------------------------

## On raw abundance data
MRM <- function(Y, X, binary = FALSE){
  require(vegan)
  Ydist <- as.numeric(vegdist(Y, method="bray", binary = binary))
  #Ydist <- as.numeric(dist(Y))
  if(length(X)==1) mod <- lm(Ydist ~ as.numeric(dist(X[[1]])))
  if(length(X)==2) mod <- lm(Ydist ~ as.numeric(dist(X$env)) + as.numeric(dist(X$xy)))
  Y.pred <- predict(mod)
  return(Y.pred)
}



##############################################################################
# Generalised linear models
##############################################################################

# INSTALLATION OF THE PACKAGES

# for installing the HMSC package, visit https://github.com/guiblanchet/HMSC
# install.packages("Matrix")

# ------------------------------------------------------------------------------
# LM or GLM
# ------------------------------------------------------------------------------

GLM <- function(Y, X, family = "quasipoisson", env.eff = "quadratic"){
  require(Matrix)
  mem <- X$mem
  env <- X$env
  if(!is.null(env)) n.env <- 1:ncol(env)
  
  if(!is.null(env)){
    if(env.eff=="linear"){
      env <- as.matrix(sparse.model.matrix(as.formula(paste("~X$env[,",n.env,"]", collapse = " + "))))[,-1]
    }
    if(env.eff=="quadratic"){
      env <- as.matrix(sparse.model.matrix(as.formula(paste("~poly(X$env[,",n.env,"],2)", collapse = " + "))))[,-1]
    }
    X$env <- env
  }
  
  Y.pred <- Y; Y.pred[] <- NA
  X0 <- as.data.frame(do.call("cbind",X))
  names(X0) <- 1:ncol(X0)
  for(i in 1:ncol(Y)){
    # Estimate GLM
    glmi <- glm(Y[,i]~., data=X0, family=family)
    Y.pred[,i] <- predict(glmi,type="response")
  }
  
  return(Y.pred)
}



# ------------------------------------------------------------------------------
# GENERALIZED ADDITIVE MODELS (GAM) FROM PACKAGE 'mgcv'
# ------------------------------------------------------------------------------

# This is one particular implementation of the GAM approach -- GAMS and mgcv are such a 
# flexible toolboxes, that many types of GAMs can be specified.
# The GAM model is fitted to each species separately.


GAM <- function(Y, X, family = "poisson", env.eff = "splines", k.env = 3, k.spa = -1, fx = FALSE)
{
  require(mgcv) 
  
  if(is.list(X)) names(X) <- NULL
  Y <- data.frame(Y)
  X0 <- data.frame(X)
  xy.indcs <- names(X) %in% c("x","y")
  
  # if only x or only y coordinate is provided
  if(sum(xy.indcs) == 1) return("Error, two coordinates must be provided.")
  
  operator <- " + "
  
  # COMPOSING THE GAM MODEL FORMULA ------------
  
  # if only environmental data are provided
  if(sum(xy.indcs) == 0)
  {
    if(env.eff=="linear"){ 
      env.part <- paste(names(X0), collapse = operator)
    } 
    if(env.eff=="quadratic"){ 
      env.part <- paste("poly(", names(X0), ", 2)", collapse = operator)
    } 
    if(env.eff=="splines"){ 
      env.part <- paste("s(", names(X0), ", k = ", k.env, ")", collapse = operator)
    } 
    
    Formula  <- as.formula(paste("Y.i ~", env.part)) 
  }
  
  # if only xy data are provided
  if(sum(xy.indcs) == ncol(X0))
  {
    Formula <- as.formula(paste("Y.i ~ s(x, y, k = ", k.spa, ", fx = ", fx, ")", sep=""))  
  }
  
  # if xy and environmental data are provided
  if(sum(xy.indcs) == 2 && ncol(X0) > 2)
  {
    env.indcs <- xy.indcs == FALSE
    env.names <-  names(X0)[env.indcs]
    
    if(env.eff=="linear"){ 
      env.part <- paste(env.names, collapse = operator)
    } 
    if(env.eff=="quadratic"){
      env.part <- paste("poly(", env.names, ", 2)", collapse = operator)
    }
    if(env.eff=="splines"){
      env.part <- paste("s(", env.names, ", k = ", k.env, ")", collapse = operator)
    }
    
    Formula <- as.formula(paste(paste("Y.i ~ s(x, y, k=", k.spa, ", fx = ", fx, ")", sep=""), env.part, sep= " + ")) 
  }
  
  # if only environmental data are provided
  if(sum(xy.indcs) == 0)
  {
    if(env.eff=="linear"){ 
      env.part <- paste(names(X0), collapse = operator)
    } else{
      env.part <- paste("poly(", names(X0), ", 2)", collapse = operator)
    } 
    Formula  <- as.formula(paste("y.i ~", env.part)) 
  }
  
  # FITTING THE GAM MODELS -----------------
  
  # the model fitting loop -- doing GAM for each species
  Y.pred <- Y; Y.pred[] <- NA # empty container for predictions
  for(i in 1:ncol(Y)){
    y.i <- Y[,i]
    XY <- data.frame(y.i = y.i, X0)
    # Fit the GAM model
    m.i <- gam(Formula, data = XY, family=family, method="REML")
    Y.pred[,i] <- predict(m.i, type="response")
  }
  
  return(Y.pred)
}



##############################################################################
# Machine learning, tree-based methods
##############################################################################

# INSTALLATION OF THE PACKAGES

# devtools::install_github("cran/mvpart")
# install.packages("randomForestSRC")
# install.packages("randomForest")
# install.packages("MultivariateRandomForest")
# install.packages("mvtboost") 

# ------------------------------------------------------------------------------
# MULTIVARIATE REGRESSION TREE from package 'mvpart'
# ------------------------------------------------------------------------------

# Multivariate regression tree is growed on the data, up to a maximum size
# of 20 splits. At each split the tree predictive performance is assessed
# by 5-fold crossvalidation. The best tree complexity is then chosen based on 
# the crossvalidation error, and the CP (complexity) value corresponding to the
# lowes crossvalidation error is used to prune the tree. The pruned tree is then
# used to generate the predictions.

MVRegTree <- function(Y, X, binary = FALSE, CV = FALSE)
{
  require(mvpart)
  
  if(binary)  Y[Y > 1] <- 1
  
  # convert the list to a data.frame
  XY <- data.frame(Y, X)
  
  # indices of the responses (Y) and the predictors (X) in the XY data.frame
  Y.indcs <- 1:ncol(Y)
  X.indcs <- (ncol(Y)+1):ncol(XY)
  
  if(CV)
  {
    # fit the multivariate tree 
    mvtree <- mvpart(data.matrix(XY[, Y.indcs]) ~ data.matrix(XY[, X.indcs]), 
                     xval = 5, # 5-fold crossvalidation 
                     data = XY,
                     minbucket = 5,  # minimum required N per terminal node
                     minauto = FALSE,
                     xv = "1se", ### should be "1se" instead of "lse",
                     #method = method,
                     plot.add=FALSE)
    
    # the tree size that gives the minimal crossvalidation error (xerror)
    best.tree <- which.min(data.frame(mvtree$cptable)$xerror)
    
    # prune the tree, using the best.tree criterion
    mvtree <- prune(mvtree, cp=mvtree$cptable[best.tree, "CP"])
    Y.pred <- predict(mvtree, type="matrix") # modif AJ 30.08.18
  }
  
  else{
    mvtree <- mvpart(data.matrix(XY[, Y.indcs]) ~ data.matrix(XY[, X.indcs]), 
                     data = XY,
                     minbucket= 5,  # minimum required N per terminal node
                     minauto = FALSE,
                     xv = "none", 
                     plot.add=FALSE)
    Y.pred <- predict(mvtree, type="matrix") # modif AJ 30.08.18
  }
  
  return(Y.pred)
}


# ------------------------------------------------------------------------------
# UNIVARIATE RANDOM FORESTS from package 'randomForest'
# ------------------------------------------------------------------------------

# A classical univariate randomForest model is fitted to each species separately,
# using the X predictors. 

UniRndForest <- function(Y, X, binary = FALSE, sampsize="default")
{
  
  X <- data.frame(X)
  Y <- Y.pred <- data.frame(Y)
  Y.pred[] <- 0 # empty container for predictions
  
  if(binary) {
    Y[Y > 1] <- 1
    Y.orig <- Y
    for(i in 1:ncol(Y)) Y[,i] <- as.factor(Y[,i])
  }
  
  if(sampsize=="default") sampsize <- nrow(Y)-1
  if(sampsize==0.2) sampsize <- ceiling(nrow(Y)/5)
  
  # calculate random forest for each species separately
  for(spec.i in 1:ncol(Y)){
    rf.i <- randomForest::randomForest(x = X, 
                                       y = Y[,spec.i], 
                                       ntree = 500,
                                       #nodesize = 5,
                                       sampsize = sampsize)
    if(binary){
      Y.pred[,spec.i] <- predict(rf.i, type = "prob")[,2]
    }
    else Y.pred[,spec.i] <- predict(rf.i)
  }
  
  return(Y.pred)
}


# ------------------------------------------------------------------------------
# MULTIVARIATE RANDOM FOREST from package 'randomForestSRC'
# ------------------------------------------------------------------------------

# A multivariate random forest (MRF) is fitted here. The critical parameter is the
# 'nodesize', which strongly influences the resulting fit. Thus, a series of
# MRFs is fitted, with increasing nodesize parameter value. The best value is
# determined by the out-of-bag R2. The MRF with this nodesize is then fitted to
# the data, and used to make the predictions. 

MVRndForest <- function(Y, X, binary = FALSE, CV = FALSE, sampsize = "default")
{
  require(randomForestSRC)
  Y.orig <- Y
  
  if(binary) {
    Y[Y > 1] <- 1
    Y.orig <- Y
    Y <- data.frame(Y)
    for(i in 1:ncol(Y)) Y[,i] <- as.factor(Y[,i])
  }
  
  if(sampsize=="default") sampsize <- nrow(Y)-1
  if(sampsize==0.2) sampsize <- ceiling(nrow(Y)/5)
  
  XY <- data.frame(Y, X)
  
  # create model formula
  spec.names <- names(data.frame(Y))
  fmla <- paste("Multivar(",
                paste(spec.names, collapse=","),
                ") ~ .")
  fmla <- as.formula(fmla)
  
  if(CV) {
    # calculate OOB error for random forests with different node sizes
    # in order to choose the optimal one
    settings <- list(nodesizes = c(2:15))
    
    for(i in 1:length(settings$nodesizes)){
      # fit the multivariate random forest model    
      rf <- rfsrc(fmla,
                  nodesize=settings$nodesizes[i],
                  sampsize = sampsize,
                  data=XY,
                  tree.err=TRUE,
                  statistics =TRUE)
      # calculate the mean OOB squared error
      Y.pred <- get.mv.predicted(rf, oob=TRUE)
      if(binary){
        Y.pred <- Y.pred[,grep(x = colnames(Y.pred), pattern=".1", fixed=TRUE)]
      }
      
      settings$OOB.err[i] <- mean((Y.orig - Y.pred)^2) 
    }
    # get the nodesize with the lowest OOB error
    best.nodesize <- settings$nodesize[which.min(settings$OOB.err)]
  }
  else best.nodesize = 5
  
  # fit the random forest with the lowest OOB error
  rf.best <- rfsrc(fmla, nodesize=best.nodesize, sampsize = sampsize, data=XY)
  if(binary){
    Y.pred <- predict(rf.best, newdata=XY[,-(1:ncol(Y)),drop=F])
    Y.pred <- do.call("cbind",lapply(Y.pred$classOutput,function(x) x$predicted[,2]))
  }
  else{
    Y.pred <- predict(rf.best, newdata=XY[,-(1:ncol(Y)),drop=F])
    Y.pred <- do.call("cbind",lapply(Y.pred$regrOutput,function(x) x$predicted))
  } 
  
  return(Y.pred)
}


# ------------------------------------------------------------------------------
# MULTIVARIATE BOOSTED REGRESSION TREES from package 'mvtboost'
# ------------------------------------------------------------------------------

# An BRT model is fitted here. First, the code identifies if the data
# are count data (positive integer, Poisson family), or if they are of 
# other type (Gaussian family is chosen). The BRT model is then fitted, with
# a suite of given parameters. In this case, we fit only regression stumps 
# (interaction depth = 1). Five-fold crossvalidation is then used to determine
# the optimal complexity of BRT for each species. This optimally complex BRT
# is then used to make the predictions.

BRT <- function(Y, X, distr = "gaussian", CV = FALSE, inter.depth = 2, shrink = 0.01)
{
  require(mvtboost)
  
  X <- data.frame(X)
  
  # parameters of the trees
  max.trees <- 1000
  if(CV)
  {
    cv.folds = 5  
  }
  else cv.folds = 0
  
  # fit the boosted trees
  mvbrt <- mvtb(Y = Y, # responses
                X = X, # predictors
                distribution = distr, # since we have count data
                shrinkage = shrink, # a.k.a. learning rate, or weight of each tree
                interaction.depth=inter.depth,
                n.trees=max.trees,
                iter.details = TRUE,
                cv.folds = cv.folds) # 5-fold cross-validation on the training set
  
  # How many regression trees should be used for the predictions?
  if(CV){
    best.n.trees <-  mvbrt$best.trees$best.cv 
    #message(best.n.trees)
  }
  else { best.n.trees = 1000  }
  
  # Helper function (inverse logit transformation)
  inverse.logit <- function(x) exp(x)/(1+exp(x))
  
  # predicted values, using the best number of trees
  if(distr == "bernoulli") Y.pred <- inverse.logit(predict(mvbrt, newdata = X, n.trees = best.n.trees) )
  if(distr == "poisson") Y.pred <- exp( predict(mvbrt, newdata = X, n.trees = best.n.trees) )
  if(distr == "gaussian") Y.pred <-  predict(mvbrt, newdata = X, n.trees = best.n.trees) 
  
  return(Y.pred)
}


