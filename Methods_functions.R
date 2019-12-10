###########################################
# Supporting Information
# "Partitioning environment and space in species-by-site matrices: a comparison of methods for community ecology and macroecology"
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

#install.packages(c("vegan", "ecodist", "stats"))

# ------------------------------------------------------------------------------
# RDA
# ------------------------------------------------------------------------------

# On Hellinger-transformed abundance data
R2.RDA <- function(Y, X, binary = FALSE){
  require(vegan)
  Y.hel <- as.matrix(decostand(as.matrix(Y), "hellinger"))
  
  Y.stand <- scale(Y.hel,scale=FALSE)
  pred <- list()
  for(i in 1:ncol(Y)){
    lmi <- lm(Y.stand[,i]~., data=as.data.frame(do.call("cbind", X)))
    pred[[i]] <- predict(lmi)
  }
  Y.pred <- as.matrix(do.call("cbind",pred))
  
  
  if(binary){
    r2 <- R2D2(Y.stand,Y.pred)["R2.multivar"]
    r2s <- c(R2.McFadden=r2,R2.Efron=NA,R2.Tjur=NA) # NOT McFadden!
  }
  else{
    r2s <- R2D2(Y.stand,Y.pred)
    r2s[c("R2.McFadden","R2.log")] <- c(NA,NA)
  }
  return(r2s)
}


# ------------------------------------------------------------------------------
# CCA
# ------------------------------------------------------------------------------

R2.CCA <- function(Y, X, binary = FALSE){
  require(vegan)
  mod <- cca(Y, do.call("cbind", X), scannf = F) 
  R2mv <- RsquareAdj(mod)$adj.r.squared
  
  if(binary){
    r2s <- c(R2.McFadden=R2mv,R2.Efron=NA,R2.Tjur=NA) # NOT McFadden!
  }
  else{
    r2s <- c(R2.classic=NA, R2.log=NA, R2.multivar=R2mv, R2.McFadden=NA)
  }
  return(r2s)
}


# ------------------------------------------------------------------------------
# dbRDA
# ------------------------------------------------------------------------------

# Using Bray-Curtis distance

R2.dbRDA <- function(Y, X, binary = FALSE){
  require(vegan)
  BC <- vegdist(Y, method="bray", binary=binary)
  mod <- capscale(BC ~.,data=as.data.frame(do.call("cbind", X)))
  R2mv <- RsquareAdj(mod)$adj.r.squared
  
  if(binary){
    r2s <- c(R2.McFadden=R2mv,R2.Efron=NA,R2.Tjur=NA) # NOT McFadden!
  }
  else{
    r2s <- c(R2.classic=NA, R2.log=NA, R2.multivar=R2mv, R2.McFadden=NA)
  }
  return(r2s)
}



##############################################################################
# Disance-based regression
##############################################################################

# ------------------------------------------------------------------------------
# MRM
# ------------------------------------------------------------------------------

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
  return(r2s)
}



##############################################################################
# Generalised linear models
##############################################################################

# INSTALLATION OF THE PACKAGES

# for installing the HMSC package, visit https://github.com/guiblanchet/HMSC
# install.packages("Matrix")

# ------------------------------------------------------------------------------
# GLM quasi-Poisson 
# ------------------------------------------------------------------------------

R2.GLM <- function(Y, X, family = "quasipoisson", interaction = FALSE){
  require(Matrix)
  mem <- X$mem
  env <- X$env
  
  Y.pred <- Y
  
  wmem <- c()
  for(i in 1:ncol(Y)){
    if(!is.null(mem)){
      # selection of MEMs
      pvals <- numeric(ncol(mem))
      for(j in 1:ncol(mem)){
        glmj <- glm(Y[,i]~mem[,j], family=family)
        pvals[j] <- coef(summary(glmj))[2,4]
      }
      wmem <- c(wmem,which(pvals<0.05))
    }
    X$mem <- mem[,unique(wmem)]
  }
  
  dat <- as.data.frame(do.call("cbind",X))
  
  for(i in 1:ncol(Y)){
    # Estimate GLM
    if(ncol(dat)>0){
      if(!is.null(env)){
        ei <- 1:ncol(env)
        if(interaction==FALSE){
          env <- as.matrix(sparse.model.matrix(as.formula(paste("~poly(X$env[,",ei,"],2)", collapse = " + "))))[,-1]
        }
        if(interaction==TRUE){
          env <- as.matrix(sparse.model.matrix(as.formula(paste("~poly(X$env[,",ei,"],2)", collapse = " * "))))[,-1]
        }
        X$env <- env
      }
      Xi <- as.data.frame(do.call("cbind",X))
      names(Xi) <- 1:ncol(Xi)
      glmi <- glm(Y[,i]~., data=Xi, family=family)
      Y.pred[,i] <- predict(glmi,type="response")
    }
    
    if(ncol(dat)==0){
      glmi <- glm(Y[,i]~1, family=family)
      Y.pred[,i] <- predict(glmi,type="response")
    }
  }
  
  if(family=="quasipoisson") r2s <- R2D2(Y,Y.pred)
  if(family=="binomial") r2s <- R2D2.binom(Y,Y.pred)
  datf <- as.data.frame(do.call("cbind",X))
  np <- ncol(datf)
  res <- list(r2s,np)
  return(res)
}


# ------------------------------------------------------------------------------
# GLMM Poisson (using HMSC)
# ------------------------------------------------------------------------------

R2.hmsc <- function(Y, X, family = "poisson", interaction = FALSE){
  
  require(HMSC)
  require(Matrix)
  
  ## format data
  if(interaction==FALSE && !is.null(X$env)){
    ei <- 1:ncol(X$env)
    env <- as.matrix(sparse.model.matrix(as.formula(paste("~poly(X$env[,",ei,"],2)", collapse = " + "))))[,-1]
    X$env <- env
  }
  if(interaction==TRUE && !is.null(X$env)){
    ei <- 1:ncol(X$env)
    env <- as.matrix(sparse.model.matrix(as.formula(paste("~poly(X$env[,",ei,"],2)", collapse = " * "))))[,-1]
    X$env <- env
  }
  formdata <- as.HMSCdata(Y = Y, X = do.call("cbind",X),
                          scaleX = FALSE, interceptX = TRUE, interceptTr = FALSE)
  # use a Poisson distribution for abundance data
  m1 <- hmsc(formdata, family = family, niter = 20000, nburn = 1000, thin = 20, verbose = FALSE)
  Y.pred <- predict(m1)
  
  if(family=="poisson"|family=="overPoisson") r2s <- R2D2(Y,Y.pred)
  if(family=="probit") r2s <- R2D2.binom(Y,Y.pred)
  return(r2s)
}


# ------------------------------------------------------------------------------
# GENERALIZED ADDITIVE MODELS (GAM) FROM PACKAGE 'mgcv'
# ------------------------------------------------------------------------------

# This is one particular implementation of the GAM approach -- GAMS and mgcv are such a 
# flexible toolboxes, that many types of GAMs can be specified.
# The GAM model is fitted to each species separately.


R2.GAM <- function(Y, X, family = "poisson", interaction = FALSE)
{
  require(mgcv) 
  
  if(is.list(X)) names(X) <- NULL
  Y <- data.frame(Y)
  X <- data.frame(X)
  xy.indcs <- names(X) %in% c("x","y")
  
  # number of spline basis functions:
  k = -1
  
  # if only x or only y coordinate is provided
  if(sum(xy.indcs) == 1) return("Error, two coordinates must be provided.")
  
  # interaction or no interaction?
  if(interaction == TRUE) { operator <- " * " }
  if(interaction == FALSE){ operator <- " + " }
  
  # COMPOSING THE GAM MODEL FORMULA ------------
  
  # if only xy data are provided
  if(sum(xy.indcs) == ncol(X))
  {
    Formula <- as.formula(paste("y.i ~ s(x, y, k = ", k, ")", sep=""))  
  }
  
  # if xy and environmental data are provided
  if(sum(xy.indcs) == 2 && ncol(X) > 2)
  {
    env.indcs <- xy.indcs == FALSE
    env.names <-  names(X)[env.indcs]
    
    env.part <- paste("poly(", env.names, ", 2)", collapse = operator)
    
    # if we also want the interaction term between the environmental effects
    Formula <- as.formula(paste(paste("y.i ~ s(x, y, k=", k, ")", sep=""), 
                                env.part, sep= " + ")) 
  }
  
  # if only environmental data are provided
  if(sum(xy.indcs) == 0)
  {
    env.part <- paste("poly(", names(X), ", 2)", collapse = operator)
    Formula  <- as.formula(paste("y.i ~", env.part)) 
  }
  
  # FITTING THE GAM MODELS -----------------
  
  preds <- Y; preds[] <- NA # empty container for predictions
  
  # the model fitting loop -- doing GAM for each species
  for(i in 1:ncol(Y))
  {
    y.i <- Y[,i]
    XY <- data.frame(y.i = y.i, X)
    
    # Fit the GAM model
    
    m.i <- gam(Formula, data = XY, family=family, method="GCV.Cp")
    preds[,i] <- predict(m.i, type="response")
  }
  
  # ASSESSING MODEL FIT (VARIATION AND DEVIANCE EXPLAINED)
  if(family=="poisson") fit <- R2D2(Y.obs = as.matrix(Y), Y.pred =  as.matrix(preds))
  if(family=="binomial") fit <- R2D2.binom(Y.obs = as.matrix(Y), Y.pred =  as.matrix(preds))
  
  return(fit)
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
  
  else
  {
    mvtree <- mvpart(data.matrix(XY[, Y.indcs]) ~ data.matrix(XY[, X.indcs]), 
                     data = XY,
                     minbucket= 5,  # minimum required N per terminal node
                     minauto = FALSE,
                     xv = "none", 
                     plot.add=FALSE)
    
    preds <- predict(mvtree, type="matrix") # modif AJ 30.08.18
  }
  
  
  if(binary)
  {
    fit <- R2D2.binom(Y.obs = Y, Y.pred = preds)
  }
  else fit <- R2D2(Y.obs = Y, Y.pred = preds)
  
  return(fit)
}

# R2.MVRegTree(Y = mat.sp, X = mat.env, binary = TRUE, CV=TRUE)


# ------------------------------------------------------------------------------
# UNIVARIATE RANDOM FORESTS from package 'randomForest'
# ------------------------------------------------------------------------------

# A classical univariate randomForest model is fitted to each species separately,
# using the X predictors. 

R2.UniRndForest <- function(Y, X, binary = FALSE)
{
  
  X <- data.frame(X)
  Y <- Y.pred <- data.frame(Y)
  Y.pred[] <- 0 # empty container for predictions
  
  if(binary) 
  {
    Y[Y > 1] <- 1
    Y.orig <- Y
    for(i in 1:ncol(Y)) Y[,i] <- as.factor(Y[,i])
  }
  
  # calculate random forest for each species separately
  for(spec.i in 1:ncol(Y))
  {
    rf.i <- randomForest::randomForest(x = X, 
                                       y = Y[,spec.i], 
                                       ntree = 500,
                                       nodesize = 5)
    if(binary)
    {
      Y.pred[,spec.i] <- predict(rf.i, type = "prob")
    }
    else Y.pred[,spec.i] <- predict(rf.i)
  }
  
  if(binary)
  {
    # extract probabilities for vote oucome = 1 (that is the .1 pattern below)
    Y.pred <- as.matrix(Y.pred)
    Y.pred <- Y.pred[, grep(x = colnames(Y.pred), pattern=".1", fixed=TRUE)]
    fit <- R2D2.binom(Y.obs = Y.orig, Y.pred = Y.pred)
  }
  
  else fit <- R2D2(Y.obs = Y, Y.pred = Y.pred)
  
  return(fit)
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
  
  if(binary) 
  {
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
  
  if(CV) 
  {
    # calculate OOB error for random forests with different node sizes
    # in order to choose the optimal one
    settings <- list(nodesizes = c(2:15))
    
    for(i in 1:length(settings$nodesizes))
    {
      # fit the multivariate random forest model    
      rf <- rfsrc(fmla,
                  nodesize=settings$nodesizes[i],
                  data=XY,
                  tree.err=TRUE,
                  statistics =TRUE)
      # calculate the mean OOB squared error
      Y.pred <- get.mv.predicted(rf, oob=TRUE)
      if(binary)
      {
        Y.pred <- Y.pred[,grep(x = colnames(Y.pred), pattern=".1", fixed=TRUE)]
      }
      
      settings$OOB.err[i] <- mean((Y.orig - Y.pred)^2) 
    }
    # get the nodesize with the lowest OOB error
    best.nodesize <- settings$nodesize[which.min(settings$OOB.err)]
  }
  else best.nodesize = 5
  
  # fit the random forest with the lowest OOB error
  rf.best <- rfsrc(fmla,
                   nodesize=best.nodesize,
                   data=XY)
  
  # predictions
  if(binary)
  {
    Y.pred <- get.mv.predicted(rf.best, oob=FALSE)
    Y.pred <- as.matrix(Y.pred)
    Y.pred <- Y.pred[, grep(x = colnames(Y.pred), pattern=".1", fixed=TRUE)]
    fit <- R2D2.binom(Y.obs = Y.orig, Y.pred = Y.pred)
  }
  else 
  {
    Y.pred <- get.mv.predicted(rf.best, oob=FALSE)
    fit <- R2D2(Y.obs = Y, Y.pred = Y.pred)
    # plot(Y, Y.pred); abline(a=0, b=1, col="red")
  }
  
  return(fit)
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
  if(CV)
  {
    best.n.trees <-  mvbrt$best.trees$best.cv 
    message(best.n.trees)
  }
  else { best.n.trees = 1000  }
  
  # predicted values, using the best number of trees
  if(distr == "bernoulli")
  {
    inverse.logit <- function(x) exp(x)/(1+exp(x))
    preds <- inverse.logit(predict(mvbrt, newdata = X, n.trees = best.n.trees) )
    #hist(preds)
  }
  if(distr == "poisson")
  {
    preds <- exp( predict(mvbrt, newdata = X, n.trees = best.n.trees) )
  }
  if(distr == "gaussian")
  {
    preds <-  predict(mvbrt, newdata = X, n.trees = best.n.trees) 
  }
  
  good.spec <- is.nan(colSums(preds)) == FALSE # some species have NaNs
  # remove the species for which NaNs are produced
  preds <- preds[, good.spec]
  Y <- Y[, good.spec]
  
  # calculate R2s
  if(binary)
  {
    fit <- R2D2.binom(Y.obs = Y, Y.pred = preds)
  }
  else {fit <- R2D2(Y = Y, Y.pred = preds)}
  
  return(fit)
}


# R2.BRT(Y = mat.sp, X = mat.env, binary = TRUE, inter.depth = 1, CV = FALSE)
