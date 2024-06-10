
# WINE DATA

library(corrplot)
library(HDclassif) # for wine data
library(classmap) # for vcr.da.train()
library(cellWise) # for cellQDA()
library(MASS) # for qda()
library(psych) # for pairs
library(gridExtra) # for arranging grid
library(GGally) # for fancy pairs

data(wine)

X = wine[,-1]
y = as.factor(wine[,1])
n=nrow(X)

### PRELIMINARY ANALYSIS

# information about the data --> distribution of the groups
n1 = length(which(y == 1)) 
n2 = length(which(y == 2)) 
n3 = length(which(y == 3)) 

# High degree of multicollinearity, see if we can drop some variables:
corrplot(cov(X), is.corr=FALSE, method = "color", tl.col = "black", addgrid.col = "grey")
# look at correlation to filter out observations with high corr
corrplot(cor(X), is.corr=FALSE, method = "color", tl.col = "black", addgrid.col = "grey")
# compute condition numbers
lambdas = eigen(cor(X))$values
sqrt(lambdas[1]/lambdas)

# reduce original data to p=8 by deleting last 4 features
X.reduced = X[,c(-9,-10,-11,-12,-13)]

# check for multicollinearity again
corrplot(cor(X.reduced), is.corr=FALSE, method = "color", tl.col = "black", addgrid.col = "grey")
# compute condition numbers
lambdas = eigen(cor(X.reduced))$values
sqrt(lambdas[1]/lambdas)
# looks okay
X = X.reduced

# pairs, cor en cov
colors = c('darkgoldenrod1', 'firebrick2', 'royalblue1')
# fancy pairs
df = data.frame(cbind(X,y))
ggpairs(df,
        columns = 1:8, 
        aes(color = y, alpha = 1.2),
        upper = list(continuous = wrap("points", size = 0.5)),
        lower = list(continuous = wrap("points", size = 0.5))) +
  scale_fill_manual(values=colors) +
  scale_colour_manual(values=colors)+
  theme_light()

# check normality assumption per group
X1 = X[which(y==1),]
p = ncol(X)
for (i in 1:p) print(shapiro.test(X1[,i])$p.value)
# select the variables indices for which H_0 is rejected (<0.01)
index = c()
for (i in 1:3) {
  if (shapiro.test(X1[,i])$p.value < 0.01) {
    index = c(index,i)
  }
}

X2 = X[which(y==2),]
for (i in 1:p) print(shapiro.test(X2[,i])$p.value)
# select the variables indices for which H_0 is rejected (<0.01)
index = c()
for (i in 1:p) {
  if (shapiro.test(X2[,i])$p.value < 0.01) {
    index = c(index,i)
  }
}

X3 = X[which(y==3),]
for (i in 1:p) print(shapiro.test(X3[,i])$p.value)
# select the variables indices for which H_0 is rejected (<0.01)
index = c()
for (i in 1:p) {
  if (shapiro.test(X3[,i])$p.value < 0.01) {
    index = c(index,i)
  }
}


# outlier detection in the groups
X1 = X[which(y==1),]
X2 = X[which(y==2),]
X3 = X[which(y==3),]
# cellMCD (cellwise outliers)
cellmcd1 = cellWise::cellMCD(X1)
cellmcd2 = cellWise::cellMCD(X2)
cellmcd3 = cellWise::cellMCD(X3)

# MCD (rowwise outliers)
mcd1 = cov.mcd(X1)
mcd2 = cov.mcd(X2)
mcd3 = cov.mcd(X3)
RMDs1 = mahalanobis(X1,center = mcd1$center, cov = mcd1$cov) # returns squared RMDs
RMDs2 = mahalanobis(X2,center = mcd2$center, cov = mcd2$cov)
RMDs3 = mahalanobis(X3,center = mcd3$center, cov = mcd3$cov)
c = qchisq(0.99, df = ncol(X1))

# cellMaps for X1, X2 and X3  
p1 = cellMap(cellmcd1$Zres, indrows = which(RMDs1>c), mTitle = "X1", drawCircles = T)
p2 = cellMap(cellmcd2$Zres, indrows = which(RMDs2>c), mTitle = "X2", drawCircles = T)
p3 = cellMap(cellmcd3$Zres, indrows = which(RMDs3>c), mTitle = "X3", drawCircles = T)
grid.arrange(grobs = list(p1,p2,p3), ncol = 3)

# get list of observations that contain at least one cellwise outliers.
cellwise_outliers = as.integer(c(which(rowSums(cellmcd1$W)<8),
                                 which(rowSums(cellmcd2$W)<8)+nrow(X1),
                                 which(rowSums(cellmcd3$W)<8)+nrow(X1)+nrow(X2)))
# get list of rowwise outliers
rowwise_outliers = as.integer(c(which(RMDs1>c), 
                                which(RMDs2>c)+nrow(X1), 
                                which(RMDs3>c)+nrow(X1)+nrow(X2)))



### CLASSIFICATON RESULTS ON WINE DATA

# select train and test observations (try cross validation with a (test) fold of 10%)
cv = function(X,y){
  X1 = X[which(y==1),]
  X2 = X[which(y==2),]
  X3 = X[which(y==3),]
  # define folds of 10% for each group
  breakpoints <- quantile(1:nrow(X1), probs = seq(0, 1, by = 0.1))
  # Split the numbers into 10 lists based on the breakpoints
  folds1 <- split(1:nrow(X1), findInterval(1:nrow(X1), breakpoints, all.inside = TRUE))
  
  breakpoints <- quantile(1:nrow(X2), probs = seq(0, 1, by = 0.1))
  folds2 <- split(1:nrow(X2)+nrow(X1), findInterval(1:nrow(X2), breakpoints, all.inside = TRUE))
  
  breakpoints <- quantile(1:nrow(X3), probs = seq(0, 1, by = 0.1))
  folds3 <- split(1:nrow(X3)+nrow(X1)+nrow(X2), findInterval(1:nrow(X3), breakpoints, all.inside = TRUE))
  
  methods = c("CQDA", "RQDA", "cellQDA", "NIMP", "L1")
  score = matrix(0, nrow = 10, ncol = length(methods))
  colnames(score) = methods
  cell.score = matrix(0, nrow = 10, ncol = length(methods))
  colnames(score) = methods
  row.score = matrix(0, nrow = 10, ncol = length(methods))
  colnames(score) = methods
  clean.score = matrix(0, nrow = 10, ncol = length(methods))
  colnames(clean.score) = methods
  
  # train and predict, each time pick a different fold to be X.train
  for (i in 1:10){ 
    test.inds = c(folds1[[i]],folds2[[i]],folds3[[i]])
    X.train = X[-test.inds,]
    y.train = as.factor(y[-test.inds])
    X.test = X[test.inds,]
    y.test = as.factor(y[test.inds])
    
    n1 = length(which(y.train == 1)) 
    n2 = length(which(y.train == 2)) 
    n3 = length(which(y.train == 3))
    n.train = nrow(X.train)
    
    cqda = qda(X.train,y.train, prior = c(n1,n2,n3)/n.train,"moment", CV = F)
    rqda = vcr.da.train(X.train, y.train, rule = "QDA", estmethod = "DetMCD")
    cellqda = cellQDA(X.train,y.train)
    
    pred.cqda = predict(cqda, X.test)$class
    pred.rqda = vcr.da.newdata(X.test,y.test, rqda)$predint
    pred.cellqda = predict(cellqda, Xnew = X.test, ynew = y.test, method = "QDA")$predint
    pred.nimp = predict(cellqda, Xnew = X.test, ynew = y.test, method = "NIMP")$predint
    pred.l1 = predict(cellqda, Xnew = X.test, ynew = y.test, method = "L1")$predint
    
    # get the position of the cell/rowwise outliers that are in y.test
    cell.inds = which(test.inds%in%cellwise_outliers)
    row.inds = which(test.inds%in%rowwise_outliers)
    y.test.cell = y.test[cell.inds]
    y.test.row = y.test[row.inds]
    # percentage of correct predictions, to fully analyse the result we distinguish the scores between cellwise, rowwise or clean observations  
    score[i,] = c(length(which(pred.cqda==y.test)),
                  length(which(pred.rqda==y.test)),
                  length(which(pred.cellqda==y.test)),
                  length(which(pred.nimp==y.test)),
                  length(which(pred.l1==y.test)))/length(y.test) # overall score
    cell.score[i,] = c(length(which(pred.cqda[cell.inds]==y.test.cell)),
                       length(which(pred.rqda[cell.inds]==y.test.cell)),
                       length(which(pred.cellqda[cell.inds]==y.test.cell)),
                       length(which(pred.nimp[cell.inds]==y.test.cell)),
                       length(which(pred.l1[cell.inds]==y.test.cell)))/length(y.test.cell) # score on observations with at least on outlying cell
    row.score[i,] = c(length(which(pred.cqda[row.inds]==y.test.row)),
                      length(which(pred.rqda[row.inds]==y.test.row)),
                      length(which(pred.cellqda[row.inds]==y.test.row)),
                      length(which(pred.nimp[row.inds]==y.test.row)),
                      length(which(pred.l1[row.inds]==y.test.row)))/length(y.test.row) # score on rowwise outliers
    clean.score[i,] =c(length(which(pred.cqda[c(-row.inds,-cell.inds)]==y.test[c(-row.inds,-cell.inds)])),
                       length(which(pred.rqda[c(-row.inds,-cell.inds)]==y.test[c(-row.inds,-cell.inds)])),
                       length(which(pred.cellqda[c(-row.inds,-cell.inds)]==y.test[c(-row.inds,-cell.inds)])),
                       length(which(pred.nimp[c(-row.inds,-cell.inds)]==y.test[c(-row.inds,-cell.inds)])),
                       length(which(pred.l1[c(-row.inds,-cell.inds)]==y.test[c(-row.inds,-cell.inds)])))/length(y.test[c(-row.inds,-cell.inds)]) # score on remaining, clean observations
  }
  return(list(score = score, row = row.score, cell = cell.score, clean = clean.score))
}

score = cv(X,y)
colMeans(score$row) # mean score on rowwise outliers
colMeans(score$cell) # mean score on cellwise outliers
colMeans(score$clean) # mean score on clean observations
average_cv_score = colMeans(score$score) #overall mean score
average_cv_score


# add missing values
delete.cells <- function(X, nmissing){
  # nmissing = list of size nrow(X) with number of cells to be deleted in each observation
  p = ncol(X)
  n = nrow(X)
  for(i in 1:n){ 
    if (nmissing[i]!=0){
      na.ind = sample.int(p, size = nmissing[i])
      X[i,na.ind] = NA
    }
  }
  return(X)
}

cv_missing = function(X,y){
  X1 = X[which(y==1),]
  X2 = X[which(y==2),]
  X3 = X[which(y==3),]
  # define folds of 10% for each group
  breakpoints <- quantile(1:nrow(X1), probs = seq(0, 1, by = 0.1))
  # Split the numbers into 10 lists based on the breakpoints
  folds1 <- split(1:nrow(X1), findInterval(1:nrow(X1), breakpoints, all.inside = TRUE))
  
  breakpoints <- quantile(1:nrow(X2), probs = seq(0, 1, by = 0.1))
  folds2 <- split(1:nrow(X2)+nrow(X1), findInterval(1:nrow(X2), breakpoints, all.inside = TRUE))
  
  breakpoints <- quantile(1:nrow(X3), probs = seq(0, 1, by = 0.1))
  folds3 <- split(1:nrow(X3)+nrow(X1)+nrow(X2), findInterval(1:nrow(X3), breakpoints, all.inside = TRUE))
  
  methods = c("cellQDA", "NIMP", "L1")
  score = matrix(0, nrow = 10, ncol = length(methods))
  colnames(score) = methods
  
  # train and predict, each time pick a different fold to be X.train
  for (i in 1:10){ 
    test.inds = c(folds1[[i]],folds2[[i]],folds3[[i]])
    print(test.inds)
    X.train = X[-test.inds,]
    y.train = as.factor(y[-test.inds])
    X.test = X[test.inds,]
    y.test = as.factor(y[test.inds])
    
    # add missing values to test data in 50% of the observations
    missing.inds = sample(nrow(X.test), nrow(X.test)*0.5)
    # these observations randomly get 1,2 or 3 missing values
    nmissing.list = sample.int(3,length(missing.inds), replace = T)
    nmissing = rep(0,nrow(X.test))
    nmissing[missing.inds] = nmissing.list
    # create test data with missing values
    X.test.na = delete.cells(X.test,nmissing)
    
    # train
    cellqda = cellQDA(X.train,y.train)
    
    # predict
    pred.cellqda = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "QDA")$predint
    pred.nimp = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "NIMP")$predint
    pred.l1 = predict(cellqda, Xnew = X.test.na, ynew = y.test, method = "L1")$predint
    
    # percentage of correct predictions 
    score[i,] = c(length(which(pred.cellqda==y.test)),
                  length(which(pred.nimp==y.test)),
                  length(which(pred.l1==y.test)))/length(y.test)
  }
  return(score)
}
set.seed(4)
score_missing = cv_missing(X,y)
average_cv_score_missing = colMeans(score_missing)
average_cv_score_missing


