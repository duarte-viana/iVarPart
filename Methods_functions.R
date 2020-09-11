###########################################
# Supporting Information, Appendix S2
# "Partitioning environment and space in site-by-species matrices: a comparison of methods for community ecology and macroecology"
# Duarte S. Viana, Petr Keil, Alienor Jeliazkov
###########################################

# Functions for the different statistical methods
# See README.txt file
# In the following functions,
# Y is a species-by-sites matrix of abundance (columns are species and rows are sites)
# X is a list with one or more matrices of predictor variables (numeric variables): environment, x-y coordinates, MEMs

# Load dependencies
source("R2D2_functions.R")



##############################################################################
# Constrained ordination
##############################################################################

# INSTALLATION OF THE PACKAGES

#install.packages(c("vegan", "stats"))

# ------------------------------------------------------------------------------
# RDA or univariate linear model
# ------------------------------------------------------------------------------

# On Hellinger-transformed data
R2.RDA <- function(Y, X, binary = FALSE, env.eff = "linear"){
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
  
  # Cross-validation (LOO)
  Y.pred.cv <- Y; Y.pred.cv[] <- NA
  for(r in 1:nrow(Y)){
    Yr <- Y.stand[-r,]
    X0r <- X0[-r,,drop=F]
    for(i in 1:ncol(Y)){
      lmi <- lm(Yr[,i]~., data=X0r)
      Y.pred.cv[r,i] <- predict(lmi, newdata=X0[r,,drop=F])
    }
  }
  
  nsites <- nrow(Y)
  nterms <- ncol(X0)
  
  if(binary){
    r2 <- R2D2(Y.stand,Y.pred)["R2.multivar"]
    r2.adj <- 1-(((nsites-1)/(nsites-nterms-1))*(1-r2))
    r2.cv <- R2D2(Y.stand,Y.pred.cv)["R2.multivar"]
    r2s.adj <- c(R2.McFadden=r2.adj,R2.Efron=NA,R2.Tjur=NA) # NOT McFadden!
    r2s.cv <- c(R2.McFadden=r2.cv,R2.Efron=NA,R2.Tjur=NA) # NOT McFadden!
  }
  else{
    r2s <- R2D2(Y.stand,Y.pred)
    r2s.adj <- 1-(((nsites-1)/(nsites-nterms-1))*(1-r2s))
    r2s.cv <- R2D2(Y.stand,Y.pred.cv)
    r2s[c("R2.McFadden","R2.log")] <- c(NA,NA)
    r2s.cv[c("R2.McFadden","R2.log")] <- c(NA,NA)
  }
  
  res <- list(r2s.adj,r2s.cv)
  return(res)
}


# ------------------------------------------------------------------------------
# CCA
# ------------------------------------------------------------------------------

## On raw abundance data
R2.CCA <- function(Y, X, binary = FALSE){
  require(vegan)
  mod <- cca(Y, do.call("cbind", X), scannf = F) 
  R2mv <- RsquareAdj(mod)$adj.r.squared
  
  if(binary){
    r2s.adj <- c(R2.McFadden=R2mv,R2.Efron=NA,R2.Tjur=NA) # NOT McFadden!
  }
  else{
    r2s.adj <- c(R2.classic=NA, R2.log=NA, R2.multivar=R2mv, R2.McFadden=NA)
  }
  r2s.cv <- c(R2.classic=NA, R2.log=NA, R2.multivar=NA, R2.McFadden=NA)
  
  res <- list(r2s.adj,r2s.cv)
  return(res)
}


# ------------------------------------------------------------------------------
# dbRDA
# ------------------------------------------------------------------------------

## On raw abundance data
R2.dbRDA <- function(Y, X, binary = FALSE){
  require(vegan)
  BC <- vegdist(Y, method="bray", binary=binary)
  mod <- capscale(BC ~.,data=as.data.frame(do.call("cbind", X)))
  R2mv <- RsquareAdj(mod)$adj.r.squared
  
  if(binary){
    r2s.adj <- c(R2.McFadden=R2mv,R2.Efron=NA,R2.Tjur=NA) # NOT McFadden!
  }
  else{
    r2s.adj <- c(R2.classic=NA, R2.log=NA, R2.multivar=R2mv, R2.McFadden=NA)
  }
  r2s.cv <- c(R2.classic=NA, R2.log=NA, R2.multivar=NA, R2.McFadden=NA)
  
  res <- list(r2s.adj,r2s.cv)
  return(res)
}



##############################################################################
# Disance-based regression
##############################################################################

# ------------------------------------------------------------------------------
# MRM
# ------------------------------------------------------------------------------

## On raw abundance data
R2.MRM <- function(Y, X, binary = FALSE){
  require(vegan)
  BC <- as.numeric(vegdist(Y, method="bray", binary = binary))
  if(length(X)==1) mod <- lm(BC ~ as.numeric(dist(X[[1]])))
  if(length(X)==2) mod <- lm(BC ~ as.numeric(dist(X$env)) + as.numeric(dist(X$xy)))
  Y.pred <- predict(mod)
  
  if(binary){
    r2 <- R2D2(BC,Y.pred)["R2.classic"]
    r2s <- c(R2.McFadden=r2,R2.Efron=NA,R2.Tjur=NA) # NOT McFadden!
  }
  else{
    r2s <- R2D2(BC,Y.pred)
    r2s[c("R2.multivar", "R2.McFadden")] <- c(NA,NA)
  }
  r2s.cv <- c(R2.classic=NA, R2.log=NA, R2.multivar=NA, R2.McFadden=NA)
  
  res <- list(r2s,r2s.cv)
  return(res)
}



##############################################################################
# Generalised linear models
##############################################################################

# INSTALLATION OF THE PACKAGES

# for installing the HMSC package, visit https://github.com/guiblanchet/HMSC
# install.packages("Matrix")

# ------------------------------------------------------------------------------
# GLM
# ------------------------------------------------------------------------------

R2.GLM <- function(Y, X, family = "quasipoisson", env.eff = "quadratic"){
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
  
  # Cross-validation (LOO)
  Y.pred.cv <- Y; Y.pred.cv[] <- NA
  for(r in 1:nrow(Y)){
    for(i in 1:ncol(Y)){
      # Estimate GLM
      glmi <- glm(Y[-r,i]~., data=X0[-r,,drop=F], family=family)
      Y.pred.cv[r,i] <- predict(glmi, newdata=X0[r,,drop=F], type="response")
    }
  }
  
  nsites <- nrow(Y)
  nterms <- ncol(X0)
  
  if(family %in% c("quasipoisson","gaussian")){
    r2s <- R2D2(Y,Y.pred)
    r2s.adj <- 1-(((nsites-1)/(nsites-nterms-1))*(1-r2s))
    r2s.cv <- R2D2(Y,Y.pred.cv)
  } 
  if(family=="binomial"){
    r2s <- R2D2.binom(Y,Y.pred)
    r2s.adj <- 1-(((nsites-1)/(nsites-nterms-1))*(1-r2s))
    r2s.cv <- R2D2.binom(Y,Y.pred.cv)
  } 
  res <- list(r2s.adj,r2s.cv)
  return(res)
}



# ------------------------------------------------------------------------------
# GENERALIZED ADDITIVE MODELS (GAM) FROM PACKAGE 'mgcv'
# ------------------------------------------------------------------------------

# This is one particular implementation of the GAM approach -- GAMS and mgcv are such a 
# flexible toolboxes, that many types of GAMs can be specified.
# The GAM model is fitted to each species separately.


R2.GAM <- function(Y, X, family = "poisson", interaction = FALSE, env.eff = "quadratic")
{
  require(mgcv) 
  
  if(is.list(X)) names(X) <- NULL
  Y <- data.frame(Y)
  X0 <- data.frame(X)
  xy.indcs <- names(X) %in% c("x","y")
  
  # number of spline basis functions:
  k = -1
  
  # if only x or only y coordinate is provided
  if(sum(xy.indcs) == 1) return("Error, two coordinates must be provided.")
  
  operator <- " + "
  
  # COMPOSING THE GAM MODEL FORMULA ------------
  
  # if only xy data are provided
  if(sum(xy.indcs) == ncol(X0))
  {
    Formula <- as.formula(paste("y.i ~ s(x, y, k = ", k, ")", sep=""))  
  }
  
  # if xy and environmental data are provided
  if(sum(xy.indcs) == 2 && ncol(X0) > 2)
  {
    env.indcs <- xy.indcs == FALSE
    env.names <-  names(X)[env.indcs]
    
    if(env.eff=="linear"){ 
      env.part <- paste(env.names, collapse = operator)
    } else{
      env.part <- paste("poly(", env.names, ", 2)", collapse = operator)
    } 
    
    Formula <- as.formula(paste(paste("y.i ~ s(x, y, k=", k, ")", sep=""), env.part, sep= " + ")) 
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
  preds <- Y; preds[] <- NA # empty container for predictions
  for(i in 1:ncol(Y)){
    y.i <- Y[,i]
    XY <- data.frame(y.i = y.i, X0)
    # Fit the GAM model
    m.i <- gam(Formula, data = XY, family=family, method="GCV.Cp")
    preds[,i] <- predict(m.i, type="response")
  }
  
  # Cross-validation (LOO)
  preds.cv <- Y; preds.cv[] <- NA # empty container for predictions
  for(r in 1:nrow(Y)){
    for(i in 1:ncol(Y)){
      y.i <- Y[-r,i]
      XYr <- data.frame(y.i = y.i, X0[-r,,drop=F])
      # Fit the GAM model
      m.i <- gam(Formula, data = XYr, family=family, method="GCV.Cp")
      preds.cv[r,i] <- predict(m.i, newdata=XY[r,-1,drop=F], type="response")
    }
  }
  
  nsites <- nrow(Y)
  if(!is.null(X$env) && env.eff=="linear") nenv <- 1
  if(!is.null(X$env) && env.eff=="quadratic") nenv <- 2
  
  # ASSESSING MODEL FIT (VARIATION AND DEVIANCE EXPLAINED)
  if(family %in% c("poisson","gaussian")){
    r2s <- R2D2(Y.obs = as.matrix(Y), Y.pred =  as.matrix(preds))
    if(!is.null(X$env)) r2s <- 1-(((nsites-1)/(nsites-nenv-1))*(1-r2s)) # only adjust according to the env terms
    r2s.cv <- R2D2(Y.obs = as.matrix(Y), Y.pred =  as.matrix(preds.cv))
  } 
  if(family=="binomial"){
    r2s <- R2D2.binom(Y.obs = as.matrix(Y), Y.pred =  as.matrix(preds))
    if(!is.null(X$env)) r2s <- 1-(((nsites-1)/(nsites-nenv-1))*(1-r2s))
    r2s.cv <- R2D2.binom(Y.obs = as.matrix(Y), Y.pred =  as.matrix(preds.cv))
  } 
  
  res <- list(r2s, r2s.cv)
  return(res)
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

R2.MVRegTree <- function(Y, X, binary = FALSE, CV = FALSE)
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
    preds <- predict(mvtree, type="matrix") # modif AJ 30.08.18
  }
  
  else{
    mvtree <- mvpart(data.matrix(XY[, Y.indcs]) ~ data.matrix(XY[, X.indcs]), 
                     data = XY,
                     minbucket= 5,  # minimum required N per terminal node
                     minauto = FALSE,
                     xv = "none", 
                     plot.add=FALSE)
    preds <- predict(mvtree, type="matrix") # modif AJ 30.08.18
    
    # Cross-validation (LOO)
    preds.cv <- Y; preds.cv[] <- NA # empty container for predictions
    for(r in 1:nrow(XY)){
      XYr <- XY[-r,]
      if(length(X.indcs)==1){
        mvtree <- mvpart(data.matrix(XYr[, Y.indcs]) ~ E, 
                         data = XYr,
                         minbucket= 5,  # minimum required N per terminal node
                         minauto = FALSE,
                         xv = "none", 
                         plot.add=FALSE)
      }
      if(length(X.indcs)==2){
        mvtree <- mvpart(data.matrix(XYr[, Y.indcs]) ~ xy.x + xy.y, 
                         data = XYr,
                         minbucket= 5,  # minimum required N per terminal node
                         minauto = FALSE,
                         xv = "none", 
                         plot.add=FALSE)
      }
      if(length(X.indcs)==3){
        mvtree <- mvpart(data.matrix(XYr[, Y.indcs]) ~ E + xy.x + xy.y, 
                         data = XYr,
                         minbucket= 5,  # minimum required N per terminal node
                         minauto = FALSE,
                         xv = "none", 
                         plot.add=FALSE)
      }
      preds.cv[r,] <- predict(mvtree, newdata=XY[c(r,r), X.indcs, drop=F], type="matrix")[1,]
    }
    
  }
  
  
  if(binary){
    fit <- R2D2.binom(Y.obs = Y, Y.pred = preds)
    fit.cv <- R2D2.binom(Y.obs = Y, Y.pred = preds.cv)
  }
  else{
    fit <- R2D2(Y.obs = Y, Y.pred = preds)
    fit.cv <- R2D2(Y.obs = Y, Y.pred = preds.cv)
  } 
  
  res <- list(fit, fit.cv)
  return(res)
}

# R2.MVRegTree(Y = mat.sp, X = mat.env, binary = TRUE, CV=TRUE)


# ------------------------------------------------------------------------------
# UNIVARIATE RANDOM FORESTS from package 'randomForest'
# ------------------------------------------------------------------------------

# A classical univariate randomForest model is fitted to each species separately,
# using the X predictors. 

R2.UniRndForest <- function(Y, X, binary = FALSE, sampsize="default")
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
  
  # Cross-validation (LOO)
  Y.pred.cv <- Y; Y.pred.cv[] <- NA # empty container for predictions
  for(r in 1:nrow(Y)){
    for(spec.i in 1:ncol(Y)){
      rf.i <- randomForest::randomForest(x = X[-r,,drop=F], 
                                         y = Y[-r,spec.i], 
                                         ntree = 500,
                                         #nodesize = 5,
                                         sampsize = sampsize)
      if(binary){
        Y.pred.cv[r,spec.i] <- predict(rf.i, newdata=X[r,,drop=F], type = "prob")[,2]
      }
      else Y.pred.cv[r,spec.i] <- predict(rf.i, newdata=X[r,,drop=F])
    }
  }
  
  
  if(binary){
    fit <- R2D2.binom(Y.obs = Y.orig, Y.pred = Y.pred)
    fit.cv <- R2D2.binom(Y.obs = Y.orig, Y.pred = Y.pred.cv)
  }
  else{
    fit <- R2D2(Y.obs = Y, Y.pred = Y.pred)
    fit.cv <- R2D2(Y.obs = Y, Y.pred = Y.pred.cv)
  } 
  
  res <- list(fit, fit.cv)
  return(res)
}

# R2.UniRndForest(Y = mat.sp, X = mat.env, binary = TRUE)



# ------------------------------------------------------------------------------
# MULTIVARIATE RANDOM FOREST from package 'randomForestSRC'
# ------------------------------------------------------------------------------

# A multivariate random forest (MRF) is fitted here. The critical parameter is the
# 'nodesize', which strongly influences the resulting fit. Thus, a series of
# MRFs is fitted, with increasing nodesize parameter value. The best value is
# determined by the out-of-bag R2. The MRF with this nodesize is then fitted to
# the data, and used to make the predictions. 

R2.MVRndForest <- function(Y, X, binary = FALSE, CV = FALSE)
{
  require(randomForestSRC)
  Y.orig <- Y
  
  if(binary) {
    Y[Y > 1] <- 1
    Y.orig <- Y
    Y <- data.frame(Y)
    for(i in 1:ncol(Y)) Y[,i] <- as.factor(Y[,i])
  }
  
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
  rf.best <- rfsrc(fmla, nodesize=best.nodesize, data=XY)
  if(binary){
    Y.pred <- predict(rf.best, newdata=XY[,-(1:ncol(Y)),drop=F])
    Y.pred <- do.call("cbind",lapply(Y.pred$classOutput,function(x) x$predicted[,2]))
  }
  else{
    Y.pred <- predict(rf.best, newdata=XY[,-(1:ncol(Y)),drop=F])
    Y.pred <- do.call("cbind",lapply(Y.pred$regrOutput,function(x) x$predicted))
  } 
  
  # Cross-validation (LOO)
  Y.pred.cv <- Y; Y.pred.cv[] <- NA # empty container for predictions
  for(r in 1:nrow(XY)){
    rf.best <- rfsrc(fmla, nodesize=best.nodesize, data=XY[-r,])
    if(binary){
      Y.pred.cv0 <- predict(rf.best, newdata=XY[r,-(1:ncol(Y)),drop=F])
      Y.pred.cv[r,] <- do.call("cbind",lapply(Y.pred.cv0$classOutput,function(x) x$predicted[,2]))
    }
    else{
      Y.pred.cv0 <- predict(rf.best, newdata=XY[r,-(1:ncol(Y)),drop=F])
      Y.pred.cv[r,] <- do.call("cbind",lapply(Y.pred.cv0$regrOutput,function(x) x$predicted))
    } 
  }
  
  # predictions
  if(binary){
    fit <- R2D2.binom(Y.obs = Y.orig, Y.pred = Y.pred)
    fit.cv <- R2D2.binom(Y.obs = Y.orig, Y.pred = Y.pred.cv)
  }
  else {
    fit <- R2D2(Y.obs = Y, Y.pred = Y.pred)
    fit.cv <- R2D2(Y.obs = Y, Y.pred = Y.pred.cv)
  }
  
  res <- list(fit, fit.cv)
  return(res)
}

# R2.MVRndForest(Y = mat.sp, X = mat.env, binary = TRUE, CV = TRUE)

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

R2.BRT <- function(Y, X, binary = FALSE, CV = FALSE, inter.depth = 3)
{
  require(mvtboost)
  
  X <- data.frame(X)
  
  # determine the appropriate distribution
  if(binary) 
  {
    distr = "bernoulli"
    Y.orig <- Y
    Y[Y > 1] <- 1
  }
  else
  {
    is.count <- sum(Y%%1 == 0) == length(Y)
    if(is.count){ distr = "poisson" }
    if(is.count == FALSE){ distr = "gaussian" }
  }
  
  # parameters of the trees
  shrink <- 0.01 # learning rate
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
  if(distr == "bernoulli") preds <- inverse.logit(predict(mvbrt, newdata = X, n.trees = best.n.trees) )
  if(distr == "poisson") preds <- exp( predict(mvbrt, newdata = X, n.trees = best.n.trees) )
  if(distr == "gaussian") preds <-  predict(mvbrt, newdata = X, n.trees = best.n.trees) 
  
  good.spec <- is.nan(colSums(preds)) == FALSE # some species have NaNs
  # remove the species for which NaNs are produced
  preds <- preds[, good.spec]
  Y.orig <- Y[, good.spec]
  
  # Cross-validation (LOO)
  preds.cv <- Y; preds.cv[] <- NA # empty container for predictions
  for(r in 1:nrow(Y)){
    mvbrt <- mvtb(Y = Y[-r,], # responses
                  X = X[-r,,drop=F], # predictors
                  distribution = distr, # since we have count data
                  shrinkage = shrink, # a.k.a. learning rate, or weight of each tree
                  interaction.depth=inter.depth,
                  n.trees=max.trees,
                  iter.details = TRUE,
                  cv.folds = 0) # 5-fold cross-validation on the training set
    
    # predicted values, using the best number of trees
    if(distr == "bernoulli") preds.cv[r,] <- inverse.logit(predict(mvbrt, newdata = X[r,,drop=F], n.trees = best.n.trees) )
    if(distr == "poisson") preds.cv[r,] <- exp( predict(mvbrt, newdata = X[r,,drop=F], n.trees = best.n.trees) )
    if(distr == "gaussian") preds.cv[r,] <-  predict(mvbrt, newdata = X[r,,drop=F], n.trees = best.n.trees) 
  }
  good.spec.cv <- is.nan(colSums(preds.cv)) == FALSE # some species have NaNs
  # remove the species for which NaNs are produced
  preds.cv <- preds.cv[, good.spec.cv]
  Y.cv <- Y[, good.spec.cv]
  
  # calculate R2s
  if(binary){
    fit <- R2D2.binom(Y.obs = Y.orig, Y.pred = preds)
    fit.cv <- R2D2.binom(Y.obs = Y.cv, Y.pred = preds.cv)
  }
  else{
    fit <- R2D2(Y = Y.orig, Y.pred = preds)
    fit.cv <- R2D2(Y = Y.cv, Y.pred = preds.cv)
  }
  
  res <- list(fit, fit.cv)
  return(res)
}


# R2.BRT(Y = mat.sp, X = mat.env, binary = TRUE, inter.depth = 1, CV = FALSE)

