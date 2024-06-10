
# simulationData1.0:
#
# This file contains simulations for testing performance of cellQDA
# when it is trained on datasets with an increasing percentage of cellwise
# contamination. To compare its performance the simulations are also
# run for CQDA and RQDA. 
# The setup is repeated 4 times, for dimensions p=2, 5, 20 and 30. Each 
# section (dimension) should be run and stored separately in the memory.
# In each dimension average execution time is computed for the training
# of cellQDA. The code for the execution times is included at the end of
# this file.


########################
# simulation study 1.0 #
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

cellsimdata <- function(K, d, nk, mu, sigma, perout, gamma, outlierType){
  
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
  if (!outlierType %in% c("casewise", "cellwisePlain", "cellwiseStructured",
                          "both")) {
    stop("jowjow, outlierType should be one of \"casewise\", \"cellwisePlain\",
    \"cellwiseStructured\" or \"both\"")
  }
  
  # generate X
  X.out = matrix(nrow = 0, ncol = d)
  W.out = matrix(nrow = 0, ncol = d)
  for (i in 1:K){
    g = generateData2(nk[i], d, mu[[i]], sigma[[i]], perout[i], gamma, outlierType = outlierType)
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

## Function for running the simulation

sim.run <- function(nsim, K, d, nk, nk.test, gamma, mu, sigma, perout, 
                    outlierType = c("casewise", "cellwisePlain", 
                                    "cellwiseStructured", "both")) { 
  # nsim    : number of simulation runs 
  # perout  : K-dim vector containing percentage of outliers per group 
  
  # Returns: list containing 3 objects (each of length nsim):
  # $train  : train data (X, y, contamination matrix) of that simulation run 
  # $test   : test data and test predictions (X.test, y.test, pred.matrix) of that simulation run
  # $time   : train and prediction times ($train.time, $pred.time) as lists of length nsim
  
  sim.out.train = list()
  sim.out.test = list()
  methods = c("CQDA", "RQDA", "cellQDA")
  runtime.pred = runtime.train = rep(0,nsim)
  
  for (i in 1:nsim){
    
    # Generate training data
    simdata.out = cellsimdata(K = K, d = d, nk = nk, mu = mu, sigma = sigma, perout = perout, gamma = gamma, outlierType = outlierType)
    X = simdata.out$X
    y = simdata.out$y
    
    # Train models
    start = Sys.time() # measure execution time (start)
    cellqda = cellQDA(X,y)
    stop = Sys.time()  
    runtime.train[i] = round(stop - start,3)
    
    cqda = qda(X,y, prior = nk/sum(nk) ,"moment", CV = F)
    rqda = vcr.da.train(X, y, rule = "QDA", estmethod = "DetMCD") 
    
    # Generate test data (without outliers!) from the same populations:
    simdata.out.test = cellsimdata(K = K, d = d, nk = nk.test, mu = mu, sigma = sigma, perout = rep(0,K), gamma = 5, outlierType = "cellwiseStructured")
    X.test = simdata.out.test$X
    y.test = simdata.out.test$y
    
    # Get predicted labels for test data
    start = Sys.time()
    pred.cellqda.test = predict(cellqda, Xnew = X.test, ynew = y.test, method = "QDA")$predint
    stop = Sys.time()
    runtime.pred[i] = round(stop-start,3)
    
    pred.cqda.test = predict(cqda, X.test)$class
    pred.rqda.test = vcr.da.newdata(X.test,y.test, rqda)$predint
    
    pred.test = cbind(pred.cqda.test, pred.rqda.test, pred.cellqda.test)
    colnames(pred.test) = methods
    
    # Store training and test data for this run
    sim.out.train = c(sim.out.train, list(list(X=X,y=y, contamination = simdata.out$W)))
    sim.out.test = c(sim.out.test, list(list(X=X.test,y=y.test, pred.matrix=pred.test)))
  }
  return(list(train = sim.out.train, test = sim.out.test, time = list(train.time = runtime.train, pred.time = runtime.pred)))
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
  methods = colnames(sim.out.test[[1]]$pred.matrix)
  nsim = length(sim.out.test)
  
  # Initialize objects
  all.testresults = matrix(0,nrow = nsim, ncol = length(methods))
  colnames(all.testresults) = methods
  conf.list.cqda = conf.list.rqda = conf.list.cellqda = list()
  
  for (i in 1:nsim){
    # Get simulation matrices and predictions
    X.test = sim.out.test[[i]]$X # test data of the i-th simulation
    y.test = sim.out.test[[i]]$y # labels of test data of i-th simulation
    pred.matrix.test = sim.out.test[[i]]$pred.matrix # n x 3 matrix with columns 'cqda', 'rqda', 'cellqda', containing predicted labels for that simulation run.
    
    # Compute accuracy for each method
    for (j in 1:length(methods)){
      all.testresults[i,j] = length(which(pred.matrix.test[,j] == y.test))/ sum(nk.test)
    }
    # get confusionMatrix information for cQDA (,1)
    confusion.cqda = confusionMatrix(data=factor(pred.matrix.test[,1], levels = c("1", "2", "3")),reference=y.test)$table
    conf.list.cqda = append(conf.list.cqda, list(confusion.cqda))
    # get confusionMatrix information for rQDA (,2)
    confusion.rqda = confusionMatrix(data=factor(pred.matrix.test[,2], levels = c("1", "2", "3")),reference=y.test)$table
    conf.list.rqda = append(conf.list.rqda, list(confusion.rqda))
    # get confusionMatrix information for cellQDA (,3)
    confusion.cellqda = confusionMatrix(data=factor(pred.matrix.test[,3], levels = c("1", "2", "3")),reference=y.test)$table
    conf.list.cellqda = append(conf.list.cellqda, list(confusion.cellqda))
  }
  # Compute average accuracy over all simulations
  average.test.acc = colSums(all.testresults)/nsim
  # Get basic statistics 
  stats.test = get.Colstat(all.testresults) # get min, max, med and sd for accuracy scores 
  # Get average and sd for runtime (pred and train)
  av.train.time = sum(sim.run.out$time$train.time)/nsim
  sd.train.time = sd(sim.run.out$time$train.time)
  av.pred.time = sum(sim.run.out$time$pred.time)/nsim
  sd.pred.time = sd(sim.run.out$time$pred.time)
  time.matrix = matrix(c(av.train.time, sd.train.time, av.pred.time,sd.pred.time),
                       nrow = 2, ncol = 2)
  colnames(time.matrix) = c("train", "pred")
  rownames(time.matrix) = c("av", "sd")
  
  return(list(all = all.testresults, 
              av = average.test.acc,
              stats = stats.test,
              confMs = list(cqda = conf.list.cqda,
                            rqda = conf.list.rqda,
                            cell = conf.list.cellqda),
              time = time.matrix))
}

## Function for choosing centers mu1, mu2, mu3 as described in section 3.1

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


### 1. DIFFERENT PERCENTAGE OF CELLWISE OUTLIERS (p=2) ###

# initializing parameters
nsim = 500
p = 2
K = 3
nk = c(100, 150, 175)
nk.test = c(30, 30, 30)
outliertype = "cellwiseStructured"
delta = 81
gamma = 10

covar1 = generateCorMat(p, corrType = "A09", seed = 1) # has rather high correlations
covar2 = generateCorMat(p, corrType = "ALYZ", CN = 2, seed = 2) # normal CN
covar3 = generateCorMat(p, corrType = "ALYZ", CN = 10, seed = 3) # "ill-conditioned", high degree of multicollinearity, but is actually more spread out than covar1!
sigma = list(covar1, covar2, covar3)
mu = get.3centers(p,sigma,delta)

perout.list = c(0, 0.05, 0.10, 0.15, 0.20) #percentage of outliers to be tested
# we do not go over 20% of cellwise outliers, this scenario is unrealistic, it is also not supported by cellMCD. 

# run and evaluate simulations
set.seed(123)
sim0 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[1],K), outliertype)
eval0 = evaluate(sim0)

sim05 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[2],K), outliertype)
eval05 = evaluate(sim05)

sim10 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[3],K), outliertype)
eval10 = evaluate(sim10)

sim15 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[4],K), outliertype)
eval15 = evaluate(sim15)

sim20 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[5],K), outliertype)
eval20 = evaluate(sim20)

# Make data more interpretable/comparable
perout.av.acc.test = matrix(c(eval0$av, eval05$av, eval10$av,
                              eval15$av, eval20$av), 
                            nrow = length(perout.list), ncol = 3, byrow = T)
colnames(perout.av.acc.test) = c("CQDA", "RQDA", "cellQDA")
rownames(perout.av.acc.test) = c("0%", "5%", "10%", "15%", "20%")

sd.test = matrix(c(eval0$stats[4,], eval05$stats[4,], eval10$stats[4,],
                   eval15$stats[4,], eval20$stats[4,]),
                 nrow = length(perout.list), ncol = 3, byrow = T)
colnames(sd.test) = c("CQDA", "RQDA", "cellQDA")
rownames(sd.test) = c("0%", "5%", "10%", "15%", "20%")

perout.time.matrix = matrix(c(as.vector(eval0$time), as.vector(eval05$time),
                              as.vector(eval10$time), as.vector(eval15$time),
                              as.vector(eval20$time)),nrow = length(perout.list), 
                            ncol = 4, byrow = T)
colnames(perout.time.matrix) = c("Train.time", "sd.train", "Pred.time", "sd.pred")
rownames(perout.time.matrix) = c("0%", "5%", "10%", "15%", "20%")


# Boxplot 

# prepare data 
bp.data1 = cbind(eval0$all[,1], eval05$all[,1], 
                 eval10$all[,1], eval15$all[,1],
                 eval20$all[,1])           # = all.testresults of cqda
bp.data2 = cbind(eval0$all[,2], eval05$all[,2], 
                 eval10$all[,2], eval15$all[,2],
                 eval20$all[,2])           # = all.testresults of rqda
bp.data3 = cbind(eval0$all[,3], eval05$all[,3], 
                 eval10$all[,3], eval15$all[,3],
                 eval20$all[,3])           # = all.testresults of cellqda
bp.data = round(rbind(bp.data1, bp.data2, bp.data3),digits = 4)
bp.data = cbind(bp.data, c(rep("CQDA",nsim),rep("RQDA",nsim), rep("cellQDA",nsim)))
bp.data = data.frame(bp.data)
colnames(bp.data) = c("0%", "5%", "10%", "15%", "20%", "method")
bp.data$method = factor(bp.data$method, levels = c("CQDA", "RQDA", "cellQDA"))
melted_data <- reshape2::melt(bp.data, id = "method")
melted_data$value = sapply(melted_data$value, as.numeric) # ggplot needs numeric input, $value was 'char'...

# blue pallette
ggplot(melted_data,aes(x=variable,y=value,fill=method))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette="Blues") + 
  ggtitle("Simulation results on testset") +
  labs(x = "percentage of outlying cells", y = "accuracy")



### 2. DIFFERENT PERCENTAGE OF CELLWISE OUTLIERS (p=5) ###

# initializing parameters
nsim = 500
p = 5
K = 3
nk = c(100, 150, 175)
nk.test = c(30, 30, 30)
outliertype = "cellwiseStructured"
delta = 81
gamma = 10

covar1 = generateCorMat(p, corrType = "A09", seed = 1) # has rather high correlations
covar2 = generateCorMat(p, corrType = "ALYZ", CN = 2, seed = 2) # normal CN
covar3 = generateCorMat(p, corrType = "ALYZ", CN = 10, seed = 3) # "ill-conditioned", high degree of multicollinearity, but is actually more spread out than covar1!
sigma = list(covar1, covar2, covar3)
mu = get.3centers(p,sigma,delta)

perout.list = c(0, 0.05, 0.10, 0.15, 0.20) 

# run and evaluate simulations
set.seed(123)
sim0 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[1],K), outliertype)
eval0 = evaluate(sim0)

sim05 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[2],K), outliertype)
eval05 = evaluate(sim05)

sim10 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[3],K), outliertype)
eval10 = evaluate(sim10)

sim15 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[4],K), outliertype)
eval15 = evaluate(sim15)

sim20 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[5],K), outliertype)
eval20 = evaluate(sim20)

# Make data more interpretable/comparable
perout.av.acc.test = matrix(c(eval0$av, eval05$av, eval10$av,
                              eval15$av, eval20$av), 
                            nrow = length(perout.list), ncol = 3, byrow = T)
colnames(perout.av.acc.test) = c("CQDA", "RQDA", "cellQDA")
rownames(perout.av.acc.test) = c("0%", "5%", "10%", "15%", "20%")

sd.test = matrix(c(eval0$stats[4,], eval05$stats[4,], eval10$stats[4,],
                   eval15$stats[4,], eval20$stats[4,]),
                 nrow = length(perout.list), ncol = 3, byrow = T)
colnames(sd.test) = c("CQDA", "RQDA", "cellQDA")
rownames(sd.test) = c("0%", "5%", "10%", "15%", "20%")

perout.time.matrix = matrix(c(as.vector(eval0$time), as.vector(eval05$time),
                              as.vector(eval10$time), as.vector(eval15$time),
                              as.vector(eval20$time)),nrow = length(perout.list), 
                            ncol = 4, byrow = T)
colnames(perout.time.matrix) = c("Train.time", "sd.train", "Pred.time", "sd.pred")
rownames(perout.time.matrix) = c("0%", "5%", "10%", "15%", "20%")

# Boxplot

# prepare data 
bp.data1 = cbind(eval0$all[,1], eval05$all[,1], 
                 eval10$all[,1], eval15$all[,1],
                 eval20$all[,1])           # = all.testresults of cqda
bp.data2 = cbind(eval0$all[,2], eval05$all[,2], 
                 eval10$all[,2], eval15$all[,2],
                 eval20$all[,2])           # = all.testresults of rqda
bp.data3 = cbind(eval0$all[,3], eval05$all[,3], 
                 eval10$all[,3], eval15$all[,3],
                 eval20$all[,3])           # = all.testresults of cellqda
bp.data = round(rbind(bp.data1, bp.data2, bp.data3),digits = 4)
bp.data = cbind(bp.data, c(rep("CQDA",nsim),rep("RQDA",nsim), rep("cellQDA",nsim)))
bp.data = data.frame(bp.data)
colnames(bp.data) = c("0%", "5%", "10%", "15%", "20%", "method")
bp.data$method = factor(bp.data$method, levels = c("CQDA", "RQDA", "cellQDA"))
melted_data <- reshape2::melt(bp.data, id = "method")
melted_data$value = sapply(melted_data$value, as.numeric) # ggplot needs numeric input, $value was 'char'...

# blue pallette
ggplot(melted_data,aes(x=variable,y=value,fill=method))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette="Blues") + 
  ggtitle("Simulation results on testset") +
  labs(x = "percentage of outlying cells", y = "accuracy")



### 3. DIFFERENT PERCENTAGE OF CELLWISE OUTLIERS (p=20) ###

# initializing parameters
nsim = 100
p = 20
K = 3
# note that we need to increase the groupsizes since p = 20
nk = c(152, 152*1.5, 152*1.75) 
nk.test = c(30, 30, 30)
outliertype = "cellwiseStructured"
delta = 81
gamma = 10

covar1 = generateCorMat(p, corrType = "A09", seed = 1) # has rather high correlations
covar2 = generateCorMat(p, corrType = "ALYZ", CN = 2, seed = 2) # normal CN
covar3 = generateCorMat(p, corrType = "ALYZ", CN = 10, seed = 3) # "ill-conditioned", high degree of multicollinearity, but is actually more spread out than covar1!
sigma = list(covar1, covar2, covar3)
mu = get.3centers(p,sigma,delta)

perout.list = c(0, 0.05, 0.10, 0.15, 0.20)

# run and evaluate simulations
set.seed(124)
sim0 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[1],K), outliertype)
eval0 = evaluate(sim0)

sim05 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[2],K), outliertype)
eval05 = evaluate(sim05)

sim10 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[3],K), outliertype)
eval10 = evaluate(sim10)

sim15 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[4],K), outliertype)
eval15 = evaluate(sim15)

sim20 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[5],K), outliertype)
eval20 = evaluate(sim20)

# Make data more interpretable/comparable
perout.av.acc.test = matrix(c(eval0$av, eval05$av, eval10$av,
                              eval15$av, eval20$av), 
                            nrow = length(perout.list), ncol = 3, byrow = T)
colnames(perout.av.acc.test) = c("CQDA", "RQDA", "cellQDA")
rownames(perout.av.acc.test) = c("0%", "5%", "10%", "15%", "20%")

sd.test = matrix(c(eval0$stats[4,], eval05$stats[4,], eval10$stats[4,],
                   eval15$stats[4,], eval20$stats[4,]),
                 nrow = length(perout.list), ncol = 3, byrow = T)
colnames(sd.test) = c("CQDA", "RQDA", "cellQDA")
rownames(sd.test) = c("0%", "5%", "10%", "15%", "20%")

perout.time.matrix = matrix(c(as.vector(eval0$time), as.vector(eval05$time),
                              as.vector(eval10$time), as.vector(eval15$time),
                              as.vector(eval20$time)),nrow = length(perout.list), 
                            ncol = 4, byrow = T)
colnames(perout.time.matrix) = c("Train.time", "sd.train", "Pred.time", "sd.pred")
rownames(perout.time.matrix) = c("0%", "5%", "10%", "15%", "20%")

# Boxplot

# prepare data 
bp.data1 = cbind(eval0$all[,1], eval05$all[,1], 
                 eval10$all[,1], eval15$all[,1],
                 eval20$all[,1])           # = all.testresults of cqda
bp.data2 = cbind(eval0$all[,2], eval05$all[,2], 
                 eval10$all[,2], eval15$all[,2],
                 eval20$all[,2])           # = all.testresults of rqda
bp.data3 = cbind(eval0$all[,3], eval05$all[,3], 
                 eval10$all[,3], eval15$all[,3],
                 eval20$all[,3])           # = all.testresults of cellqda
bp.data = round(rbind(bp.data1, bp.data2, bp.data3),digits = 4)
bp.data = cbind(bp.data, c(rep("CQDA",nsim),rep("RQDA",nsim), rep("cellQDA",nsim)))
bp.data = data.frame(bp.data)
colnames(bp.data) = c("0%", "5%", "10%", "15%", "20%", "method")
bp.data$method = factor(bp.data$method, levels = c("CQDA", "RQDA", "cellQDA"))
melted_data <- reshape2::melt(bp.data, id = "method")
melted_data$value = sapply(melted_data$value, as.numeric) # ggplot needs numeric input, $value was 'char'...

# blue pallette
ggplot(melted_data,aes(x=variable,y=value,fill=method))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette="Blues") + 
  ggtitle("Simulation results on testset") +
  labs(x = "percentage of outlying cells", y = "accuracy")



### 4. DIFFERENT PERCENTAGE OF CELLWISE OUTLIERS (p=30) ###

# initializing parameters
nsim = 100
p = 30
K = 3
# increase groupsize again
nk = c(200, 200*1.5, 200*1.75) 
nk.test = c(30, 30, 30)
outliertype = "cellwiseStructured"
delta = 81
gamma = 10

covar1 = generateCorMat(p, corrType = "A09", seed = 1) # has rather high correlations
covar2 = generateCorMat(p, corrType = "ALYZ", CN = 2, seed = 2) # normal CN
covar3 = generateCorMat(p, corrType = "ALYZ", CN = 10, seed = 3) # "ill-conditioned", high degree of multicollinearity, but is actually more spread out than covar1!
sigma = list(covar1, covar2, covar3)
mu = get.3centers(p,sigma,delta)

perout.list = c(0, 0.05, 0.10, 0.15, 0.20) 

# run and evaluate simulations
set.seed(123)
sim0 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[1],K), outliertype)
eval0 = evaluate(sim0)

sim05 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[2],K), outliertype)
eval05 = evaluate(sim05)

sim10 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[3],K), outliertype)
eval10 = evaluate(sim10)

sim15 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[4],K), outliertype)
eval15 = evaluate(sim15)

sim20 = sim.run(nsim, K, p, nk, nk.test, gamma, mu, sigma, rep(perout.list[5],K), outliertype)
eval20 = evaluate(sim20)

# Make data more interpretable/comparable
perout.av.acc.test = matrix(c(eval0$av, eval05$av, eval10$av,
                              eval15$av, eval20$av), 
                            nrow = length(perout.list), ncol = 3, byrow = T)
colnames(perout.av.acc.test) = c("CQDA", "RQDA", "cellQDA")
rownames(perout.av.acc.test) = c("0%", "5%", "10%", "15%", "20%")

sd.test = matrix(c(eval0$stats[4,], eval05$stats[4,], eval10$stats[4,],
                   eval15$stats[4,], eval20$stats[4,]),
                 nrow = length(perout.list), ncol = 3, byrow = T)
colnames(sd.test) = c("CQDA", "RQDA", "cellQDA")
rownames(sd.test) = c("0%", "5%", "10%", "15%", "20%")

perout.time.matrix = matrix(c(as.vector(eval0$time), as.vector(eval05$time),
                              as.vector(eval10$time), as.vector(eval15$time),
                              as.vector(eval20$time)),nrow = length(perout.list), 
                            ncol = 4, byrow = T)
colnames(perout.time.matrix) = c("Train.time", "sd.train", "Pred.time", "sd.pred")
rownames(perout.time.matrix) = c("0%", "5%", "10%", "15%", "20%")


# Boxplot

# prepare data 
bp.data1 = cbind(eval0$all[,1], eval05$all[,1], 
                 eval10$all[,1], eval15$all[,1],
                 eval20$all[,1])           # = all.testresults of cqda
bp.data2 = cbind(eval0$all[,2], eval05$all[,2], 
                 eval10$all[,2], eval15$all[,2],
                 eval20$all[,2])           # = all.testresults of rqda
bp.data3 = cbind(eval0$all[,3], eval05$all[,3], 
                 eval10$all[,3], eval15$all[,3],
                 eval20$all[,3])           # = all.testresults of cellqda
bp.data = round(rbind(bp.data1, bp.data2, bp.data3),digits = 4)
bp.data = cbind(bp.data, c(rep("CQDA",nsim),rep("RQDA",nsim), rep("cellQDA",nsim)))
bp.data = data.frame(bp.data)
colnames(bp.data) = c("0%", "5%", "10%", "15%", "20%", "method")
bp.data$method = factor(bp.data$method, levels = c("CQDA", "RQDA", "cellQDA"))
melted_data <- reshape2::melt(bp.data, id = "method")
melted_data$value = sapply(melted_data$value, as.numeric) # ggplot needs numeric input, $value was 'char'...


# blue pallette
ggplot(melted_data,aes(x=variable,y=value,fill=method))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette="Blues") + 
  ggtitle("Simulation results on testset") +
  labs(x = "percentage of outlying cells", y = "accuracy")


### Execution Time Figures ###

#LOAD DATA p=2:
time.data2 = perout.time.matrix

#LOAD DATA p=5 OVER DATA p=2
time.data5 = perout.time.matrix

#LOAD DATA p=20:
time.data20 = perout.time.matrix

#LOAD DATA p=30:
time.data30 = perout.time.matrix

### prepare data (TRAINING times)
percentage <- c(0, 5, 10, 15, 20)
time_p2 <- time.data2[,1]
time_p5 <- time.data5[,1]
time_p20 <- time.data20[,1]
time_p30 <- time.data30[,1]

data <- data.frame(
  percentage = rep(percentage, 4),
  execution_time = c(time_p2, time_p5, time_p20, time_p30),
  variables = rep(c("p=2", "p=5", "p=20", "p=30"), each = 5)
)

library(ggplot2)
# Ensure 'variables' is a factor with the correct order
data$variables <- factor(data$variables, levels = c("p=2", "p=5", "p=20", "p=30"))

# Define your custom color palette corresponding to the levels
custom_colors <- c("p=2" = "darkolivegreen3", 
                   "p=5" = "firebrick1", 
                   "p=20" = "royalblue1", 
                   "p=30" = "goldenrod1")

# Plot (y=execution, x=percentage)
ggplot(data, aes(x = percentage, y = execution_time, color = variables, group = variables)) +
  geom_line(size = 0.9) +     # Add lines
  geom_point(size = 2.3) +      # Add points
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20), labels = c("0%", "5%", "10%", "15%", "20%")) +
  labs(x = expression(paste("Percentage of Contamination (", epsilon, ")")), y = "Time (s)", color = "Variables") +
  scale_color_manual(values = custom_colors) +  # Use custom colors
  theme_minimal()

# plot (y=execution, x=variables)
# Assuming we have execution times as follows:
execution_times <- list(
  "2" = time_p2,
  "5" = time_p5,
  "20" = time_p20,
  "30" = time_p30
)

# Define the values of p and perout
p_values <- c(2, 5, 20, 30)
perout_values <- c(0, 5, 10, 15, 20)

# Prepare data frame for plotting
data <- data.frame(
  p = rep(p_values, each = length(perout_values)),
  perout = rep(perout_values, times = length(p_values)),
  execution_time = unlist(execution_times)
)

# Define your custom color palette corresponding to the levels
custom_colors <- c("0" = "darkolivegreen3", 
                   "5" = "firebrick1", 
                   "10" = "royalblue1", 
                   "15" = "goldenrod1",
                   "20" = "purple2")


# Create the plot
ggplot(data, aes(x = p, y = execution_time, color = factor(perout))) +
  geom_line(aes(group = perout), size = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = p_values) +
  scale_color_manual(
    values = custom_colors,
    labels = c("0%", "5%", "10%", "15%", "20%")
  ) +
  labs(
    x = "Number of Variables (p)",
    y = "Time (s)",
    color = "Contamination"
  ) +
  theme_minimal()
