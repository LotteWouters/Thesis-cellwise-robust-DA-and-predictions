
cellQDA = function(X, ...) UseMethod("cellQDA")

cellQDA.default = function(X,y){
  
  # This function fits a cellwise robust Quadratic Discriminant
  # Analysis on the given training data.  
  # Cellwise robust discriminant scores are computed based on
  # the cellMCD method.
  #
  # Arguments:
  #   X         : (train) data matrix. Missing values are not allowed.
  #   y         : vector with the given class labels.
  #               Missing y are not allowed, such cases should be
  #               removed first.
  #
  # Returns:
  #   X         : the given training data (converted to a matrix)
  #   yint      : given labels as integer 1, 2, 3, ...
  #   y         : given labels
  #   levels    : levels of y
  #   predint   : predicted class number for train observations. For each case
  #               this is the class with the highest mixture density.
  #   pred      : predicted label (level) of each train observation.
  #   priors    : prior probabilities of each group (computed as the empirical 
  #               probabilities from the train data)
  #   classM    : list with cellMCD centers of each class.
  #   classS    : list with cellMCD covariances of each class.
  #
  X <- as.matrix(X) # in case it is a data frame
  if (nrow(X) == 1) X <- t(X) 
  if (sum(is.na(as.vector(X))) > 0) {
    stop("The data matrix X has NA's.")
  }
  n <- nrow(X)
  d <- ncol(X)
  if (n < 2) stop("The training data should have more than one case.")
  
  # Check whether y and its levels are of the right form:
  checked <- checkLabels(y, n, training = T)
  lab2int <- checked$lab2int 
  indsv   <- checked$indsv 
  levels  <- checked$levels 
  nlab    <- length(levels) 
  yint    <- lab2int(y) 
  yintv   <- yint[indsv]

  classSizes <- rep(0, nlab)
  for (g in seq_len(nlab)) {classSizes[g] <- sum(yintv == g)}
  # classSizes
  MD2s    <- matrix(0, n, nlab) # squared mahalanobis distances
  mixf    <- matrix(0, n, nlab) # component of mixture density
  lmixf   <- matrix(0, n, nlab) # the logarithm of that density
  class.M <- class.S <- list()
  class.prior <- rep(0, nlab)
  
  tinyclasses <- which(classSizes < (d + 1)) 
  if (length(tinyclasses) > 0) {
    stop(paste0(" The following classes have too few members",
                " to compute a covariance matrix:",
                levels[tinyclasses]))
  }
  # Compute the center and covariance matrix of each class:
  iCN <- rep(NA, nlab) # for inverse condition numbers 
  for (g in seq_len(nlab)) {
    clinds       <- indsv[which(yintv == g)] 
    # class indices in 1,...,n
    class.out    <- cellMCD(X[clinds, ]) # compute the cellWise mu and sigma 
    iCN[g]       <- numericalCond(class.out$S) # compute inverse condition number (smallest eigenvalue/largest eigenvalue) of covariance 
    class.M[[g]]   <- class.out$m # store group mean
    class.S[[g]]   <- class.out$S # store group covariance
  }
  tinyCN <- which(iCN < 1e-06)
  if (length(tinyCN) > 0) {
    stop(paste0("\n The following classes have an ill-conditioned",
                " covariance matrix: \n", levels[tinyCN],
                "\n You may want to try LDA or another classifier"))
  } # From here onward we can use the inverse of each S.
  #
  # compute MD2 of each object to all classes:
  for (g in seq_len(nlab)) { # g=1
    MD2s[, g] <- mahalanobis(X, class.M[[g]], class.S[[g]])
  }
  # Compute prior probabilities:
  for (g in seq_len(nlab)) { # g=1
    clinds <- which(yintv == g) # indices in 1, ..., indsv
    class.prior[g] <- length(clinds) / length(indsv) 
  }
  #
  for (g in seq_len(nlab)) { # g=1
    S          <- class.S[[g]]
    prior      <- class.prior[g]
    lmixf[, g] <- -0.5*(log(det(S)) + MD2s[, g]) + log(prior)
    mixf[, g]  <- exp(lmixf[, g])
  }
  
  # Compute predictions for all objects in 1, ..., n
  predint <- rep(NA, n)
  predint <- apply(lmixf[, , drop = FALSE], 1, which.max) 
  cl <- match.call()
  cl[[1L]] <- as.name("cellQDA")
  res <- list(X = X,
              yint = yint,
              y = levels[yint],
              levels = levels,
              predint = predint,
              pred = list(pred.X = levels[predint]),
              priors = class.prior,
              classM = class.M,
              classS = class.S)
  class(res) <- "cellQDA"
  res
}

predict.cellQDA <- function(object, Xnew, ynew = NULL, method = c("QDA","L1","NIMP")){
  
  # A new way of prediction labels for new, unseen observations.
  # Allows for new data with missing values or outlying cells.
  #
  # Arguments:
  # object      : object of class "cellQDA"
  # Xnew        : new, unseen cases to be classified (preferably as matrix)
  # ynew        : true labels of the observations in Xnew
  # method      : determines how the new labels are predicted. With "QDA"
  #               (the default) the labels are predicted based on maximizing the 
  #               cellwise quadratic discriminant scores. For the other methods
  #               each new observation is added to each one of the groups
  #               and an imputed version of the new case is obtained
  #               through the cellHandler algorithm. With the option "L1"
  #               a label is predicted based on the shortest L1-distance 
  #               between xnew and ximp. With "NIMP" a label is given based
  #               on the number of imputed cells in xnew. 
  #
  # Returns:
  # predint     : predicted labels for Xnew (as integers)
  # pred        : predicted labels for Xnew (with the original labels)
  # lPred       : (only for method "QDA") log(mixture density) of Xnew
  #               predicted labels.
  #
  
  if(!inherits(object, "cellQDA")) stop("object not of class \"cellQDA\"")
  method <- match.arg(method) 
  
  Xnew <- as.matrix(Xnew) # in case Xnew is not a matrix 
  X <- object$X
  d <- ncol(X)
  
  # checks
  if (ncol(Xnew) != d){ 
    
    if (ncol(Xnew) == 1 & nrow(Xnew) == d) { Xnew <- t(Xnew) } # user accidently put Xnew in wrong position
    
    else {stop("Xnew and X do not have the same number of variables.")}
  }
  
  n <- nrow(Xnew) # now Xnew is in the right format so we can get size
  levels <- object$levels # same names as in training data
  ngroup <- length(levels) # number of classes
  class.M <- object$classM 
  class.S <- object$classS
  predint <- lPred <- rep(NA, n)
  
  for (i in seq_len(n)){ # consider each observation separately and classify
    
    class.M <- object$classM 
    class.S <- object$classS
    X <- object$X 
    xnew <- Xnew[i, ] # i-th observation
    
    # check for NA's, if so, we continue only with x_observed
    if (any(is.na(xnew))){
      # consider only the non missing variables
      obs.ind <- which(!is.na(xnew)) 
      for (g in seq_len(ngroup)){ # g=1
        class.M[[g]] <- class.M[[g]][obs.ind]
        class.S[[g]] <- class.S[[g]][obs.ind, obs.ind]
      }
      X <- as.matrix(X[,obs.ind])
      xnew <- xnew[obs.ind] # we continue with the observed x
    }
    
    if (method == "QDA"){ 
      
      MD2s <- rep(0, ngroup) # for squared Mahalanobis distances
      mixf <- lmixf <- rep(0, ngroup) # mixture density and its log
      class.prior <- object$priors
      
      for (g in seq_len(ngroup)) { # g=1
        MD2s[g]  <- mahalanobis(xnew, class.M[[g]], class.S[[g]])
        lmixf[g] <- -0.5*(log(det(class.S[[g]])) + MD2s[g]) + log(class.prior[g])
        mixf[g]  <- exp(lmixf[g])
      }
      # get prediction for xnew
      lPred[i]   <- max(lmixf)
      predint[i] <- which.max(lmixf)
    }
    
    else { # method is "L1" or "NIMP" 
      
      yint <- object$yint
      nimputed <- l1.dist <- MD2s <- rep(0,ngroup) 
      xnew.imp <- matrix(0, nrow = ngroup, ncol = ncol(X)) # matrix with the ximp's for each group
      
      for (g in seq_len(ngroup)){ # g=1
        # add xnew to the g-th group
        Xg <- rbind(X[which(yint == g),],xnew)
        # apply cellHandler algorithm
        cellHandler.out <- cellWise::cellHandler(as.matrix(Xg), class.M[[g]], class.S[[g]])
        
        ximp <- cellHandler.out$Ximp[nrow(Xg),]
        xnew.imp[g,] <- as.vector(ximp) # store imputed version of xnew
        
        # store number of imputed cells in xnew, could be zero
        nimputed[g] <- length(intersect(cellHandler.out$indcells, c(nrow(Xg)*1:ncol(Xg))))
        # store L1 distance between xnew and ximp
        l1.dist[g] <- sum(abs(xnew-ximp))
        MD2s[g] <- mahalanobis(ximp, class.M[[g]], class.S[[g]])
      }
      
      if (method == "L1"){ # option 1. assign to class with lowest L1 between xnew and ximp
        predint[i] <- which.min(l1.dist)
        if (sum(l1.dist == l1.dist[predint[i]]) > 1) { # no unique minimum
          tied.groups <- which(l1.dist == l1.dist[predint[i]])
          predint[i] <- tied.groups[which.min(MD2s[tied.groups])]
        }
      }
      
      if (method == "NIMP"){ # option 2. assign to class with lowest nimputed
        predint[i] <- which.min(nimputed)
        if (sum(nimputed == nimputed[predint[i]]) > 1) { # no unique minimum
          tied.groups <- which(nimputed == nimputed[predint[i]])
          predint[i] <- tied.groups[which.min(MD2s[tied.groups])]
        }
      }
    }
  }
  pred <- levels[predint]
  return(list(predint = predint,
              pred = pred,
              lPred = lPred))
}

