
# Extension of generateData() from cellWise package. Function generates 3 types
# of cellwise outliers for data with any center/location. 
generateData2 = function(n, d,
                         mu, Sigma,
                         perout, gamma,
                         outlierType = "casewise", 
                         seed = NULL) {
  
  # outlierType:
  #   casewise: replacement value: gamma*smallest eigenvector where smallest eigenv has md 1
  #   cellwisePlain: cellwise naive
  #   cellwiseStructured: cellwise smallest eigenvector
  #   both
  
  if (!outlierType %in% c("casewise", "cellwisePlain", "cellwiseStructured",
                          "both")) {
    stop("outlierType should be one of \"casewise\", \"cellwisePlain\",
    \"cellwiseStructured\" or \"both\"")
  }
  
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  
  # generating random data
  X = MASS::mvrnorm(n, mu = mu, Sigma = Sigma)
  indcells <- c()
  indrows  <- c()
  W <- array(0, dim(X))
  #adding contamination
  if (perout > 0) {
    if (outlierType == "casewise") {
      replacement <- eigen(Sigma)$vectors[, d] /
        sqrt(mahalanobis(eigen(Sigma)$vectors[, d], mu, Sigma))
      replacement <- gamma * replacement * sqrt(d) * d
      ind         <- seq_len(floor(perout * n)) 
      X[ind, ]    <- matrix(replacement, nrow = length(ind), ncol = d, byrow = TRUE)
      indrows <- ind 
    } else if (outlierType == "cellwisePlain") {
      ind <- replicate(d, sample(seq_len(n), perout * n, replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      X[ind] <- gamma
      indcells <- ind
    } else if (outlierType == "cellwiseStructured") {
      ind <- replicate(d, sample(seq_len(n), perout * n, replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W   <- array(0, dim(X)); W[ind] <- 1
      for (i in seq_len(n)) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)]/
            sqrt(mahalanobis(eigen_out[, length(continds)], rep(0,length(continds)), Sigma[continds, continds]))
          X[i, continds] <- mu[continds] + replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells, i + (continds - 1) * n)
        }
      }
    } else  if (outlierType == "both") {
      replacement <- eigen(Sigma)$vectors[, d] /
        sqrt(mahalanobis(eigen(Sigma)$vectors[, d], mu, Sigma))
      replacement <- gamma * replacement * d * sqrt(d)
      ind         <- seq_len(floor(perout / 2 * n))
      X[ind, ]    <- matrix(replacement, nrow = length(ind), ncol = d, byrow = TRUE)
      indrows <- ind
      
      # Now add cells
      startind <- (floor(perout / 2 * n) + 1) # first row which is not a casewise outliers
      ind <- replicate(d, sample(startind:n, ceiling(perout / 2 * n),
                                 replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W   <- array(0, dim(X)); W[ind] <- 1
      for (i in startind:n) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)] /
            sqrt(mahalanobis(eigen_out[, length(continds)], mu[continds], Sigma[continds, continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells,  i + (continds - 1) * n)
        }
      }
    }  
  }
  return(list(X = X, indcells = indcells, W = W,
              indrows = indrows))
}

