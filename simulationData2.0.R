
# simulationData2.0:
#
# This file contains simulations for testing predictions by cellQDA, L1
# and NIMP on new data with missing cells or/and cellwise contamination. 
# We run 3 simulations with new observations containing (1) only cellwise
# contamination, (2) only missing values, and (3) contamination and missing
# values.


########################
# simulation study 2.0 #
########################

library(caret)
library(classmap)
library(robustHD)
library(cellWise)
library(MASS)
library(ggplot2)
library(reshape2)

### AUXILIARY FUNCTIONS ###

## Function for generating classification data (X,y) with cellwise outliers

cellsimdata <- function(K, d, nk, mu, sigma, perout, gamma){
  
  # K           : number of groups
  # d           : number of variables
  # nk          : vector of length K with number of observations per group
  # mu          : list containing mean of each group
  # sigma       : list containing covariance of each group
  # perout      : K-dim vector containing proportion of contamination per group
  # gamma       : determines outlyingness of the cells (usually between 1 and 10)
  # outlierType : should be either "casewise", "cellwisePlain", "cellwiseStructured" or "both"
  
  if (length(nk) != K){stop("length nk does not match number of groups K")}
  if (length(mu) != K){stop("length of mu does not match K")}
  if (length(sigma) != K){stop("number of covariance matrices does not match K")}
  if (length(perout) != K){"length of perout does not match K"}
  
  # generate X
  X.out = matrix(nrow = 0, ncol = d)
  W.out = matrix(nrow = 0, ncol = d)
  for (i in 1:K){
    g = generateData2(nk[i], d, mu[[i]], sigma[[i]], perout[i], gamma, outlierType = "cellwiseStructured")
    X.out = rbind(X.out,g$X) # add to final data matrix
    # store the contamination matrix of this group
    W.out = rbind(W.out,g$W)
  }
  
  # generate y
  lab = c()
  for (i in 1:K){
    lab = c(lab, rep(i,nk[i])) 
  }
  y.out = as.factor(matrix(lab, ncol = 1))
  
  return(list(X = data.frame(X.out), y = y.out, W = W.out))
}

## Function for getting (columnwise) statistical results from one set of simulation data

get.Colstat = function(X){
  min = apply(X, 2, min)
  max = apply(X, 2, max)
  med = apply(X, 2, median)
  sd = apply(X, 2, sd)
  stats = matrix(rbind(min, max, med, sd), nrow = 4,ncol = ncol(X))
  rownames(stats) = c("min", "max", "med", "sd")
  return(round(stats,digits = 4))
}

## Function for running the simulation: (1) no missing cells

sim.run1 <- function(nsim, K, d, nk, nk.test, gamma, mu, sigma, perout) { 
  
  # Returns list containing 3 objects:
  # $train  : list of length nsim containing train data (X.train,y.train, contamination scheme W)
  # $test   : list of length nsim containing test data, predictions and subsets (X.test, y.test, predictions, subs)
  # $time   : vector containing average execution time of NIMP and L1
  
  # nk : is now the total amount of observations per group (training and testing together)
  # perout : K-dim vector containing percentage of outliers per group  
  
  sim.out.train = list()
  sim.out.test = list()
  methods = c("cellQDA", "NIMP", "L1")
  nimp.pred.time = l1.pred.time = rep(0,nsim)
  
  for (i in 1:nsim){
    
    # Generate classification data
    simdata.out = cellsimdata(K = K, d = d, nk = nk, mu = mu, sigma = sigma, perout = perout, gamma = gamma)
    X = simdata.out$X
    y = simdata.out$y
    W = simdata.out$W
    
    # select at random training and testing data
    test.ind = c() 
    for (g in 1:K){
      rand.ind = sample.int(nk[g], size = nk.test[g]) # sample random indices of the observation to be used per group
      test.ind = c(test.ind, which(y == g)[rand.ind])
    }
    X.test = X[test.ind,]
    y.test = y[test.ind]
    X.train = X[-test.ind,]
    y.train = y[-test.ind]
    W.test = W[test.ind,] 
    
    # train the model
    cellqda = cellQDA(X.train, y.train)
    
    # Get predicted labels for test data from cellQDA, Nimp and L1
    pred.cellqda.test = predict(cellqda, Xnew = X.test, ynew = y.test, method = "QDA")$predint
    start = Sys.time() # we keep track of execution time of NIMP and L1
    pred.nimp.test = predict(cellqda, Xnew = X.test, ynew = y.test, method = "NIMP")$predint
    stop = Sys.time()
    nimp.pred.time[i] = round(stop-start,3)
    start = Sys.time()
    pred.l1.test = predict(cellqda, Xnew = X.test, ynew = y.test, method = "L1")$predint
    stop = Sys.time()
    l1.pred.time[i] = round(stop-start,3)
    
    pred.test = cbind(pred.cellqda.test, pred.nimp.test, pred.l1.test)
    colnames(pred.test) = methods
    
    # extra info: get indices of the observations with 0,1/p,2/p,...100% of outlying cells
    vec = apply(W.test,1,sum)
    subsets = vector(mode = "list", length = p+1)
    
    for (j in 0:p){
      subsets[[j+1]] = as.vector(which(vec == j)) 
    }
    
    # Store training and test data and test predictions for this run
    sim.out.train = c(sim.out.train, list(list(X=X.train,y=y.train, contamination = simdata.out$W)))
    sim.out.test = c(sim.out.test, list(list(X=X.test,y=y.test, predictions = pred.test, subs = subsets)))
  }
  av.l1.pred.time = sum(l1.pred.time) / nsim
  av.nimp.pred.time = sum(nimp.pred.time) / nsim
  return(list(train = sim.out.train, test = sim.out.test, time = c(av.nimp.pred.time, av.l1.pred.time)))
}

## Function for running the simulation: (2) only missing cells

sim.run2 <- function(nsim, K, d, nk, nk.test, mu, sigma, nmissing) { 
  
  # Returns list containing 2 objects (each of length nsim):
  # $train  : train data (X.train,y.train, contamination scheme W)
  # $test   : test data before and after deleting cells and test predictions (X.test, X.test.na, y.test, predictions)

  # nk : is now the total amount of observations per group (training and testing together)
  # perout : K-dim vector containing percentage of outliers per group  
  
  sim.out.train = list()
  sim.out.test = list()
  methods = c("cellQDA", "NIMP", "L1")
  
  for (i in 1:nsim){
    
    # Generate (clean) classification data
    simdata.out = cellsimdata(K = K, d = d, nk = nk, mu = mu, sigma = sigma, perout = rep(0,K), gamma = 5)
    X = simdata.out$X
    y = simdata.out$y
    
    # select random training and testing data
    test.ind = c() 
    for (g in 1:K){
      rand.ind = sample.int(nk[g], size = nk.test[g]) # sample random indices of the observation to be used per group
      test.ind = c(test.ind, which(y == g)[rand.ind])
    }
    X.test = X[test.ind,]
    y.test = y[test.ind]
    X.train = X[-test.ind,]
    y.train = y[-test.ind]
    
    # delete cells from X.test observations
    X.test.na = delete.cells(X.test, nmissing) 
    
    # train the model
    cellqda = cellQDA(X.train, y.train)
    
    # Get predicted labels for test data from cellQDA, NIMP and L1
    pred.cellqda.test = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "QDA")$predint
    pred.nimp.test = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "NIMP")$predint
    pred.l1.test = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "L1")$predint
    
    # Bring predictions together in a matrix
    pred.test = cbind(pred.cellqda.test, pred.nimp.test, pred.l1.test)
    colnames(pred.test) = methods
    
    # Store training and test data and test predictions for this run
    sim.out.train = c(sim.out.train, list(list(X=X.train,y=y.train)))
    sim.out.test = c(sim.out.test, list(list(X=X.test, X.na = X.test.na, y=y.test, predictions = pred.test)))
  }
  return(list(train = sim.out.train, test = sim.out.test))
}

sim.run3 <- function(nsim, K, d, nk, nk.test, gamma, mu, sigma, perout, nmissing) { ##returns list containing 2 objects ($train and $test) each of length nsim containing: ($train) train data (X.train,y.train, W) and ($test) test data and test predictions (X.test, y.test, predictions) per simulation run, predictions = matrix containing predictions from cellQDA, Nimp and L1 of the X.test observations in its columns
  
  # Returns list containing 3 objects:
  # $train  : list of length nsim containing train data (X.train,y.train, contamination scheme W)
  # $test   : list of length nsim containing test data, predictions and subsets (X.test.na, y.test, predictions, subs)
  # $time   : vector containing average execution time of NIMP and L1
  
  # nk : is now the total amount of observations per group (training and testing together)
  # perout : K-dim vector containing percentage of outliers per group 
  
  sim.out.train = list()
  sim.out.test = list()
  methods = c("cellQDA", "NIMP", "L1")
  nimp.pred.time = l1.pred.time = rep(0,nsim)
  
  for (i in 1:nsim){
    
    # Generate classification data with perout% contamination
    simdata.out = cellsimdata(K = K, d = d, nk = nk, mu = mu, sigma = sigma, perout = perout, gamma = gamma)
    X = simdata.out$X
    y = simdata.out$y
    W = simdata.out$W
    
    # select random training and testing data
    test.ind = c() 
    for (g in 1:K){
      rand.ind = sample.int(nk[g], size = nk.test[g]) # sample random indices of the observation to be used per group
      test.ind = c(test.ind, which(y == g)[rand.ind])
    }
    X.test = X[test.ind,]
    y.test = y[test.ind]
    X.train = X[-test.ind,]
    y.train = y[-test.ind]
    W.test = W[test.ind,] 
    
    # delete cells from X.test observations
    X.test.na = delete.cells(X.test, nmissing, W.test)

    # train the model (on possibly contaminated data)
    cellqda = cellQDA(X.train, y.train)
    
    # Get predicted labels for test data from cellQDA, NIMP and L1
    pred.cellqda.test = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "QDA")$predint
    start = Sys.time()
    pred.nimp.test = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "NIMP")$predint
    stop = Sys.time()
    nimp.pred.time[i] = round(stop-start,3)
    start = Sys.time()
    pred.l1.test = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "L1")$predint
    stop = Sys.time()
    l1.pred.time[i] = round(stop-start,3)
    
    pred.test = cbind(pred.cellqda.test, pred.nimp.test, pred.l1.test)
    colnames(pred.test) = methods
    
    # Keep track of subsets in X.test.na
    vec = apply(W.test,1,sum) 
    subsets = vector(mode = "list", length = d) # to store the inds of the test observations with k outlying cells (k=0,..,4)
    for (k in 0:(d-1)){
      subsets[[k+1]] = as.vector(which(vec == k)) # subset[[1]] == indices of test obs with 0 outlying cells (--> 4 clean cells), subset[[5]] == 4 outlying cells 
    }
    
    # Store training and test data and test predictions for this run
    sim.out.train = c(sim.out.train, list(list(X=X.train,y=y.train, contamination = W[-test.ind,])))
    sim.out.test = c(sim.out.test, list(list(X=X.test.na,y=y.test, predictions = pred.test, contamination = W.test, subs = subsets)))
  }
  av.l1.pred.time = sum(l1.pred.time) / nsim
  av.nimp.pred.time = sum(nimp.pred.time) / nsim
  return(list(train = sim.out.train, test = sim.out.test, time = c(av.nimp.pred.time, av.l1.pred.time)))
}


## Function for deleting random cells from a set of observations X

delete.cells <- function(X, nmissing, W = NULL){ 
  
  # Returns: X with nmissing (NA) values in each row chosen at random
  
  p = ncol(X)
  n = nrow(X)
  
  if(is.null(W)){ # in case of no contamination, only missing (sim.run2)
    for (i in 1:n){
      na.ind = sample.int(p, size = nmissing)
      X[i,na.ind] = NA
    }
  } else{ # in case of contamination and missing (sim.run3), we want to avoid 
          # deleting cells that were contaminated, otherwise we create a data
          # set that looks similar to (2) only missing values.
    for (i in 1:n){
      w = W[i,]
      clean.inds = which(w == 0) # get indices of the clean cells
      cont.inds = which(w == 1)
      if (length(clean.inds) >= nmissing){
        na.inds = sample(clean.inds, size = nmissing) # sample only from the clean cells
      } else { # less clean cells than nmissing --> in this case we change all clean.inds to NA's and sample from the contaminated cells the remaining NA indices
        na.inds = clean.inds
        extra = sample(cont.inds, size = nmissing-length(clean.inds))
        na.inds = c(na.inds,extra)
      }
      X[i,na.inds] = NA
    }
  }
  return(X)
}


## Function for evaluating the simulation output

evaluate <- function(sim.run.out){
  
  # Returns list containing 4 objects:
  # $all  : all accuracies of the test predictions for each simulation run  
  # $av   : average accuracy of all nsim for all three methods
  # $stat : basic statistical results of $all
  # $confMs : confusion matrices of all predictions for each simulation and for all methods
  
  # Get train and test simulation output 
  sim.out.test = sim.run.out$test
  methods = colnames(sim.out.test[[1]]$predictions)
  nsim = length(sim.out.test)
  
  # Initialize objects
  all.accuracies = matrix(0,nrow = nsim, ncol = length(methods))
  colnames(all.accuracies) = methods
  conf.list.cellqda = conf.list.nimp = conf.list.l1 = list()
  
  for (i in 1:nsim){
    # Get simulation matrices and predictions
    y.test = sim.out.test[[i]]$y # labels of test data of i-th simulation
    predictions = sim.out.test[[i]]$predictions # n x 3 matrix with columns 'cellqda', 'nimp', 'l1' containing predicted labels for that simulation run.
    
    # Compute accuracy for each method
    for (j in 1:length(methods)){
      all.accuracies[i,j] = length(which(predictions[,j] == y.test)) / sum(nk.test)
    }
    # get confusionMatrix information for cellQDA (,1)
    confusion.cellqda = confusionMatrix(data=factor(predictions[,1], levels = c("1", "2", "3")),reference=y.test)$table
    conf.list.cellqda = append(conf.list.cellqda, list(confusion.cellqda))
    # get confusionMatrix information for nimp (,2)
    confusion.nimp = confusionMatrix(data=factor(predictions[,2], levels = c("1", "2", "3")),reference=y.test)$table
    conf.list.nimp = append(conf.list.nimp, list(confusion.nimp))
    # get confusionMatrix information for l1 (,3)
    confusion.l1 = confusionMatrix(data=factor(predictions[,3], levels = c("1", "2", "3")),reference=y.test)$table
    conf.list.l1 = append(conf.list.l1, list(confusion.l1))
  }
  # Compute average accuracy over all simulations
  average.test.acc = colSums(all.accuracies)/nsim
  # Get basic statistics 
  stats.test = get.Colstat(all.accuracies) # get min, max, med and sd for accuracy scores 
  
  return(list(all = all.accuracies, 
              av = average.test.acc,
              stats = stats.test,
              confMs = list(cellqda = conf.list.cellqda,
                            nimp = conf.list.nimp,
                            l1 = conf.list.l1)))
}

## Function for choosing centers mu1, mu2, mu3

get.3centers = function(p, covars, delta){
  
  mu1 = mu2 = mu3 = rep(0,p)
  
  # get inverse matrices
  inv_covar1 = solve(covars[[1]])
  inv_covar2 = solve(covars[[2]])
  inv_covar3 = solve(covars[[3]])
  
  # get a and b such that mu2 = (-a,0,...0) and mu3 = (0,-b,...,0)
  A = inv_covar1 + inv_covar2
  B = inv_covar1 + inv_covar3
  a = sqrt(delta/A[1,1])
  b = sqrt(delta/B[2,2])
  
  mu2[1] = -a # mu2 = (-a,0,...,0)
  mu3[2] = -b # mu3 = (0,-b,...,0)
  
  # check if the MDs are good
  MD12 = mahalanobis(mu2,mu1,covar1) + mahalanobis(mu1,mu2,covar2)
  MD13 = mahalanobis(mu3,mu1,covar1) + mahalanobis(mu1,mu3,covar3)
  cat("The MD12: ", MD12, "\nThe MD13: ", MD13)
  
  # Extract the solution
  return(mu = list(mu1,mu2,mu3))
}



### 1. COMPARING CELLQDA, NIMPUTED AND L1 (only contamination) ###

nsim = 500
p = 5
K = 3
nk = c(120, 150, 175)
nk.test = c(30, 30, 30)
delta = 81
gamma = 10

covar1 = generateCorMat(p, corrType = "A09", seed = 5) # has rather high correlations
covar2 = generateCorMat(p, corrType = "ALYZ", CN = 2, seed = 6) # normal CN
covar3 = generateCorMat(p, corrType = "ALYZ", CN = 10, seed = 7) # "ill-conditioned", high degree of multicollinearity, but is actually more spread out than covar1!
sigma = list(covar1, covar2, covar3)
mu = get.3centers(p, sigma, delta)

# run the simulation
set.seed(300038)

sim0 = sim.run1(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(0,K))
eval0 = evaluate(sim0)

sim05 = sim.run1(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(0.05,K))
eval05 = evaluate(sim05)

sim10 = sim.run1(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(0.10,K))
eval10 = evaluate(sim10)


# Make data more interpretable/comparable
perout.av.acc.test = matrix(c(eval0$av, eval05$av, eval10$av), 
                            nrow = 3, ncol = 3, byrow = T)
colnames(perout.av.acc.test) = c("cellQDA", "NIMP", "L1")
rownames(perout.av.acc.test) = c("0%", "5%", "10%")

sd.test = matrix(c(eval0$stats[4,], eval05$stats[4,], eval10$stats[4,]),
                 nrow = 3, ncol = 3, byrow = T)
colnames(sd.test) = c("cellQDA", "NIMP", "L1")
rownames(sd.test) = c("0%", "5%", "10%")

perout.time.matrix = matrix(c(as.vector(sim0$time), as.vector(sim05$time),
                              as.vector(sim10$time)),nrow = 3, 
                            ncol = 2, byrow = T)
colnames(perout.time.matrix) = c("av.time.NIMP", "av.time.L1")
rownames(perout.time.matrix) = c("0%", "5%", "10%")


### Boxplot (perout = 0%, 5% and 10%)

# prepare data
bp.data1 = cbind(eval0$all[,1]*100, eval05$all[,1]*100, eval10$all[,1]*100)# = all.testresults of cellqda
bp.data2 = cbind(eval0$all[,2]*100, eval05$all[,2]*100, eval10$all[,2]*100)# = all.testresults of nimp
bp.data3 = cbind(eval0$all[,3]*100, eval05$all[,3]*100, eval10$all[,3]*100)# = all.testresults of l1
bp.data = round(rbind(bp.data1, bp.data2, bp.data3),digits = 4)
bp.data = cbind(bp.data,c(rep("cellQDA",nsim), rep("NIMP",nsim), rep("L1",nsim)))
bp.data = data.frame(bp.data)
colnames(bp.data) = c("0%", "5%", "10%", "method")
bp.data$method = factor(bp.data$method, levels = c("cellQDA","NIMP", "L1"))
melted_data <- reshape2::melt(bp.data, id = "method")
melted_data$value = sapply(melted_data$value, as.numeric) 

# plot (blue palette)
ggplot(melted_data,aes(x=variable,y=value,fill=method))+
  geom_boxplot() +
  theme_minimal() +
  scale_fill_brewer(palette="Blues") +
  labs(x = expression(paste("Percentage of Contamination (", epsilon, ")")), y = "Accuracy (%)")


### Accuracy per subset (perout = 10%)
# initialize objects
methods = c("cellQDA", "NIMP", "L1")
cor.pred.subset = matrix(0, nrow = p+1, ncol = 3) #matrix to keep the number of correct predictions a method made per subset
colnames(cor.pred.subset) = methods
count.subset = rep(0, p+1)

for (i in 1:nsim){
  
  # get information from sim.run
  y.test = sim10$test[[i]]$y
  predictions = sim10$test[[i]]$predictions # n x 3 matrix 
  subsets = sim10$test[[i]]$subs
  
  for (k in 0:p){ # go over each subset
    ind.subset = subsets[[k+1]] # get the observations from this nsim with k cellwise outliers
    # we keep track of how many of each subgroup was present in X.test
    count.subset[k+1] = count.subset[k+1] + length(ind.subset)
    # count number of correct predictions of each method for each subset and add to the sum
    for (j in 1:3){
      cor.pred.subset[k+1,j] = cor.pred.subset[k+1,j] + length(which(predictions[ind.subset,j] == y.test[ind.subset]))
    }
  }
}

# per subgroup, what percentage of this subgroup (if we look at X.test) was correctly predicted? 
per.cor.pred.subset = cor.pred.subset / count.subset * 100

# subgroup plot
x <- seq(0,p)
y <- per.cor.pred.subset

# Combine data into a data frame
data <- data.frame(x = x, cellQDA = y[,1], NIMP = y[,2], L1 = y[,3])

# Reshape the data into long format
data_long <-melt(data, id.vars = "x", variable.name = "method", value.name = "y")

y_labels = paste0(seq(0,100,20), "%")
# Create the plot
ggplot(data_long, aes(x = factor(x), y = y, fill = method)) +
  geom_col(position = "dodge") +  # Add colored bars
  labs(x = "# Contaminated Cells in Observation", y = "Accuracy (%)") +  # Add labels and title
  scale_y_continuous(labels = y_labels, breaks = seq(0, 100, 20)) +  # Set y-axis labels as percentages
  theme_minimal() + # Optional: Use a minimal theme for aesthetics
  scale_fill_brewer(palette="Blues") 




### 2. COMPARING CELLQDA, NIMPUTED AND L1 (only missing values) ###

nsim = 100
p = 5
K = 3
nk = c(120, 150, 175)
nk.test = c(30, 30, 30)
delta = 81

covar1 = generateCorMat(p, corrType = "A09", seed = 5) # has rather high correlations
covar2 = generateCorMat(p, corrType = "ALYZ", CN = 2, seed = 6) # normal CN
covar3 = generateCorMat(p, corrType = "ALYZ", CN = 10, seed = 7) # "ill-conditioned", high degree of multicollinearity, but is actually more spread out than covar1!
sigma = list(covar1, covar2, covar3)
mu = get.3centers(p,sigma, delta)

# run the simulation
set.seed(300038)

sim1 = sim.run2(nsim, K, p, nk, nk.test, mu, sigma, nmissing = 1)
eval1 = evaluate(sim1)

sim2 = sim.run2(nsim, K, p, nk, nk.test, mu, sigma, nmissing = 2)
eval2 = evaluate(sim2)

sim3 = sim.run2(nsim, K, p, nk, nk.test, mu, sigma, nmissing = 3)
eval3 = evaluate(sim3)

# Make data more interpretable/comparable
perout.av.acc.test = matrix(c(eval1$av, eval2$av, eval3$av), 
                            nrow = 3, ncol = 3, byrow = T)
colnames(perout.av.acc.test) = c("cellQDA", "NIMP", "L1")
rownames(perout.av.acc.test) = c("1NA", "2NA", "3NA")

sd.test = matrix(c(eval1$stats[4,], eval2$stats[4,], eval3$stats[4,]),
                 nrow = 3, ncol = 3, byrow = T)
colnames(sd.test) = c("cellQDA", "NIMP", "L1")
rownames(sd.test) = c("1NA", "2NA", "3NA")

# Boxplot

bp.data1 = cbind(eval1$all[,1], eval2$all[,1], eval3$all[,1])# = all.testresults of cellqda
bp.data2 = cbind(eval1$all[,2], eval2$all[,2], eval3$all[,2])# = all.testresults of nimp
bp.data3 = cbind(eval1$all[,3], eval2$all[,3], eval3$all[,3])# = all.testresults of l1
bp.data = round(rbind(bp.data1,bp.data2, bp.data3),digits = 4)
bp.data = cbind(bp.data, c(rep("cellQDA",nsim), rep("NIMP",nsim), rep("L1",nsim)))
bp.data = data.frame(bp.data)
colnames(bp.data) = c("1", "2", "3", "method")
bp.data$method = factor(bp.data$method, levels = c("cellQDA", "NIMP", "L1"))
melted_data <- reshape2::melt(bp.data, id = "method")
melted_data$value = sapply(melted_data$value, as.numeric) 

# blue pallette
ggplot(melted_data,aes(x=variable,y=value,fill=method))+
  geom_boxplot() +
  theme_minimal() +
  scale_fill_brewer(palette="Blues") + 
  labs(x = expression(paste("Number of Missing Cells (", xi, ")")), y = "Accuracy (%)")



### 3. COMPARING CELLQDA, NIMPUTED AND L1 (contamination and missing cells) ###

nsim = 100
p = 5
K = 3
nk = c(120, 150, 175)
nk.test = c(30, 30, 30)
delta = 81
gamma = 10

covar1 = generateCorMat(p, corrType = "A09", seed = 5) # has rather high correlations
covar2 = generateCorMat(p, corrType = "ALYZ", CN = 2, seed = 6) # normal CN
covar3 = generateCorMat(p, corrType = "ALYZ", CN = 10, seed = 7) # "ill-conditioned", high degree of multicollinearity, but is actually more spread out than covar1!
sigma = list(covar1, covar2, covar3)
mu = get.3centers(p,sigma, delta)

# run the simulation
set.seed(300038)

sim1 = sim.run3(nsim, K, p, nk, nk.test, gamma, mu, sigma, 
                perout = rep(0.15,K), nmissing = 1) 
eval1 = evaluate(sim1)

sim2 = sim.run3(nsim, K, p, nk, nk.test, gamma, mu, sigma, 
                perout = rep(0.15,K), nmissing = 2)
eval2 = evaluate(sim2)

sim3 = sim.run3(nsim, K, p, nk, nk.test, gamma, mu, sigma, 
                perout = rep(0.15,K), nmissing = 3)
eval3 = evaluate(sim3)


# interpret the data and figures
perout.av.acc.test = matrix(c(eval1$av, eval2$av, eval3$av), 
                            nrow = 3, ncol = 3, byrow = T)
colnames(perout.av.acc.test) = c("cellQDA", "NIMP", "L1")
rownames(perout.av.acc.test) = c("1NA", "2NA", "3NA")

sd.test = matrix(c(eval1$stats[4,], eval2$stats[4,], eval3$stats[4,]),
                 nrow = 3, ncol = 3, byrow = T)
colnames(sd.test) = c("cellQDA", "NIMP", "L1")
rownames(sd.test) = c("1NA", "2NA", "3NA")

# filter out the four groups: if we make all combinations
# four groups can be formed based on the number of clean cells left
# in the new observation -> 0,1,2,3 clean cells (taking the results of sim1, sim2 and sim3 together)
# Additionally, we also want to check how many of the X.test observations
# had no contamination (group 5)

# nsim1 (nmissing=1)
# options:
# 0 clean cells --> obs had 4 contaminated cells because 1 cell was made NA --> this corresponds with the 5th subset in subs!
# 1 clean cells --> obs had 3 contaminated cells because 1 cell was made NA --> 4th subset
# 2 clean cells --> obs had 2 contaminated cells because 1 cell was made NA --> 3
# 3 clean cells --> obs had 1 contaminated cells because 1 cell was made NA --> 2
# 4 clean cells --> obs had 0 contaminated cells because 1 cell was made NA --> 1
methods = c("cellQDA", "NIMP", "L1")
cor.pred.subset.sim1 = matrix(0, nrow = 5, ncol = 3) # rows reflect the subsets, following the order of the subset --> row1 = test observations from subset[[1]], so with 0 contaminated cells (i.e. 4 clean cells)
colnames(cor.pred.subset.sim1) = methods
rownames(cor.pred.subset.sim1) = c("4clean", "3clean", "2clean", "1clean", "0clean")
count.subset.sim1 = rep(0,5) # to count the number of obs in a particular subset, we also want to count the size of the subset with 4 clean cells, this is the last one
confusion.list = replicate(5, matrix(0, nrow = nsim, ncol = sum(nk.test)), simplify = FALSE)

for (i in 1:nsim){
  
  # get information from sim.run
  y.test = sim1$test[[i]]$y
  predictions = sim1$test[[i]]$predictions # n x 3 matrix (predictions from all 3 methods)
  subsets = sim1$test[[i]]$subs
  
  for (k in 1:5){ # go over ALL subset (4,3,..,0 clean cells)
    ind.subset = subsets[[k]] # get the observations from this nsim that had k-1 contaminated cells
    # we keep track of how many of each subgroup was present in X.test
    count.subset.sim1[k] = count.subset.sim1[k] + length(ind.subset)
    # count number of correct predictions of each method for each subset and add to the sum
    for (j in 1:3){
      cor.pred.subset.sim1[k,j] = cor.pred.subset.sim1[k,j] + length(which(predictions[ind.subset,j] == y.test[ind.subset]))
    }
    # get confusion matrix of cellQDA for this subset
    confusion.list[[k]][i,] = sim1$test[[k]]$predictions[,1]
  }
}


# nsim2 (nmissing=2)
# options:
# 0 clean cells --> obs had 3 contaminated cells because 2 cell were made NA --> 4e subset
# 1 clean cells --> obs had 2 contaminated cells because 2 cell were made NA --> 3e subset
# 2 clean cells --> obs had 1 contaminated cells because 2 cell were made NA --> 2
# 3 clean cells --> obs had 0 contaminated cells because 2 cell were made NA --> 1
# Note: subset[[5]] with 4 contaminated cells in in this case a test observation with 3 contaminated cells and 2 NAs --> add to the count of 3 contaminated 
cor.pred.subset.sim2 = matrix(0, nrow = 4, ncol = 3) # order going down: row1 == obs with 0 contaminated cells (3 clean cells)
colnames(cor.pred.subset.sim2) = methods
rownames(cor.pred.subset.sim2) = c("3clean", "2clean", "1clean", "0clean")
count.subset.sim2 = rep(0,4) 

for (i in 1:nsim){
  
  # get information from sim.run
  y.test = sim2$test[[i]]$y
  predictions = sim2$test[[i]]$predictions 
  subsets = sim2$test[[i]]$subs
  
  for (k in 1:4){ # go over ALL subset (3,2,1,0 clean cells)
    ind.subset = subsets[[k]] # get the observations from this nsim that had k-1 contaminated cells
    # we keep track of how many of each subgroup was present in X.test
    count.subset.sim2[k] = count.subset.sim2[k] + length(ind.subset)
    # count number of correct predictions of each method for each subset and add to the sum
    for (j in 1:3){
      cor.pred.subset.sim2[k,j] = cor.pred.subset.sim2[k,j] + length(which(predictions[ind.subset,j] == y.test[ind.subset]))
    }
  }
  ind.subset = subsets[[5]] 
  count.subset.sim2[4] = count.subset.sim2[4] + length(ind.subset)
  # count number of correct predictions of each method for each subset and add to the sum
  for (j in 1:3){
    cor.pred.subset.sim2[4,j] = cor.pred.subset.sim2[4,j] + length(which(predictions[ind.subset,j] == y.test[ind.subset]))
  }
}

# nsim3 (nmissing=3)
# options:
# 0 clean cells --> obs had 2 contaminated cells because 3 cell were made NA --> 3e subset
# 1 clean cells --> obs had 1 contaminated cells because 3 cell were made NA --> 2e subset
# 2 clean cells --> obs had 0 contaminated cells because 3 cell were made NA --> 1
# Note: subsets[[4]] and subset[[5]], with 3 and 4 contaminated cells become observations with 2 contaminated cells and 3 NAs in this case
cor.pred.subset.sim3 = matrix(0, nrow = 3, ncol = 3) 
colnames(cor.pred.subset.sim3) = methods
rownames(cor.pred.subset.sim3) = c("2clean", "1clean", "0clean")
count.subset.sim3 = rep(0,3) 

for (i in 1:nsim){
  
  # get information from sim.run
  y.test = sim3$test[[i]]$y
  predictions = sim3$test[[i]]$predictions 
  subsets = sim3$test[[i]]$subs
  
  for (k in 1:3){ # go over ALL subset (2,1,0 clean cells)
    ind.subset = subsets[[k]] # get the observations from this nsim that had k-1 contaminated cells
    # we keep track of how many of each subgroup was present in X.test
    count.subset.sim3[k] = count.subset.sim3[k] + length(ind.subset)
    # count number of correct predictions of each method for each subset and add to the sum
    for (j in 1:3){
      cor.pred.subset.sim3[k,j] = cor.pred.subset.sim3[k,j] + length(which(predictions[ind.subset,j] == y.test[ind.subset]))
    }
  }
  ind.subset = subsets[[4]] 
  count.subset.sim3[3] = count.subset.sim3[3] + length(ind.subset)
  # count number of correct predictions of each method for each subset and add to the sum
  for (j in 1:3){
    cor.pred.subset.sim3[3,j] = cor.pred.subset.sim3[3,j] + length(which(predictions[ind.subset,j] == y.test[ind.subset]))
  }
  ind.subset = subsets[[5]] 
  count.subset.sim3[3] = count.subset.sim3[3] + length(ind.subset)
  # count number of correct predictions of each method for each subset and add to the sum
  for (j in 1:3){
    cor.pred.subset.sim3[3,j] = cor.pred.subset.sim3[3,j] + length(which(predictions[ind.subset,j] == y.test[ind.subset]))
  }
}



# plots

# Boxplot (1NA, 2NA, 3NA)

bp.data1 = cbind(eval1$all[,1]*100, eval2$all[,1]*100, eval3$all[,1]*100)# = all.testresults of cellqda
bp.data2 = cbind(eval1$all[,2]*100, eval2$all[,2]*100, eval3$all[,2]*100)# = all.testresults of nimp
bp.data3 = cbind(eval1$all[,3]*100, eval2$all[,3]*100, eval3$all[,3]*100)# = all.testresults of l1
bp.data = round(rbind(bp.data1, bp.data2, bp.data3),digits = 4)
bp.data = cbind(bp.data, c(rep("cellQDA",nsim),rep("NIMP",nsim), rep("L1",nsim)))
bp.data = data.frame(bp.data)
colnames(bp.data) = c("1", "2", "3", "method")
bp.data$method = factor(bp.data$method, levels = c("cellQDA","NIMP", "L1"))
melted_data <- reshape2::melt(bp.data, id = "method")
melted_data$value = sapply(melted_data$value, as.numeric) 


# blue pallette
ggplot(melted_data,aes(x=variable,y=value,fill=method))+
  geom_boxplot() +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_fill_brewer(palette="Blues") + 
  labs(x = expression(paste("Number of Missing Cells (", xi, ")")), y = "Accuracy (%)")


# subgroup plot

# add correct predictions of subgroups together from sim1, sim2 and sim3
cor.pred.subset.total = matrix(0, nrow = 5, ncol = 3)
rownames(cor.pred.subset.total) = c("noCont", "3clean", "2clean", "1clean", "0clean")
colnames(cor.pred.subset.total) = methods
# subset with no contaminated cells, or 4 clean cells (for sim1), or 3 clean cells (for sim2), or 2 clean cells (sim3) --> just all test observations that were clean and only had NAs
cor.pred.subset.total[1,] = cor.pred.subset.sim1[1,] + 
  cor.pred.subset.sim2[1,] + cor.pred.subset.sim3[1,] 
# subset with 3 clean cells
cor.pred.subset.total[2,] = cor.pred.subset.sim1[2,]
# subset with 2 clean cells
cor.pred.subset.total[3,] = cor.pred.subset.sim1[3,] + 
  cor.pred.subset.sim2[2,]
# subset with 1 clean cell
cor.pred.subset.total[4,] = cor.pred.subset.sim1[4,] + 
  cor.pred.subset.sim2[3,] + cor.pred.subset.sim3[2,]
# subset[[5]] with 0 clean cells
cor.pred.subset.total[5,] = cor.pred.subset.sim1[5,] + 
  cor.pred.subset.sim2[4,] + cor.pred.subset.sim3[3,]

# add counts per subgroups together from sim1, sim2, sim3
count.subset.total = rep(0,5)
# subset no contaminated
count.subset.total[1] = count.subset.sim1[1] + count.subset.sim2[1] +
  count.subset.sim3[1]
# subset with 3 clean cells
count.subset.total[2] = count.subset.sim1[2]
# subset with 2 clean cells
count.subset.total[3] = count.subset.sim1[3] + count.subset.sim2[2]
# subset with 1 clean cells
count.subset.total[4] = count.subset.sim1[4] + count.subset.sim2[3] +
  count.subset.sim3[2]
# subset with 0 clean cells
count.subset.total[5] = count.subset.sim1[5] + count.subset.sim2[4] +
  count.subset.sim3[3]

# percentage of correct predictions over sim1, sim2, sim3 per subset
per.cor.pred.subset = (cor.pred.subset.total / count.subset.total) * 100

# Convert the matrix to a data frame
df <- as.data.frame(per.cor.pred.subset)
df$Dataset <- rownames(per.cor.pred.subset)
# Reshape the data frame to long format
scores_long <- melt(df, id.vars = "Dataset", variable.name = "Method", value.name = "Score")
# Set the levels of the Dataset factor explicitly
scores_long$Dataset <- factor(scores_long$Dataset, levels = c("noCont", "3clean", "2clean", "1clean", "0clean"))

ggplot(scores_long, aes(x = Dataset, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Type of Test Observation",
    y = "Accuracy (%)"
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme_minimal() +
  scale_fill_brewer(palette="Blues")