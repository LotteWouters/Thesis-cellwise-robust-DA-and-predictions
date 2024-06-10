
# FOREST SOIL DATA

library(classmap) # for vcr.da.train()
library(cellWise) # for cellQDA()
library(MASS) # for qda()
library(gridExtra) # for arranging grids
library(rrcov)
library(corrplot)
library(GGally) # fancy pairs

data(soil)

# (!!) we only work with the observations from 1983

X = data.frame(soil[soil$D == 0,c(4,5,7)]) 
y = soil[soil$D ==0,1]
# information about the data --> distribution of the groups
n1 = length(which(y == 2)) 
n2 = length(which(y == 3))
n3 = length(which(y == 1))
# we work only with two groups: F = 2,3
X = data.frame(X[y != 1,])
y = as.factor(c(rep(1,n1), rep(2,n2)))
n = nrow(X)


### PRELIMINARY ANALYSIS

# corrplot
corrplot(cor(X), is.corr=T, method = "color", tl.col = "black", addgrid.col = "grey")

# fancy pairsplot
colors = c( "firebrick1","royalblue1")
df = data.frame(cbind(X,y))
ggpairs(df,
        columns = 1:3, 
        aes(color = y, alpha = 1.2),
        upper = list(continuous = wrap("points", size = 2)),
        lower = list(continuous = wrap("points", size = 2))) +
  scale_fill_manual(values=colors) +
  scale_colour_manual(values=colors) +
  theme_light()

# check normality assumptions
X1 = X[which(y==1),c(1,2,3)]
rownames(X1) = c(1:23)
for (i in 1:3) print(shapiro.test(X1[,i])$p.value)
# select the variables indices for which H_0 is rejected (<0.01)
index = c()
for (i in 1:3) {
  if (shapiro.test(X1[,i])$p.value < 0.01) {
    index = c(index,i)
  }
}

X2 = X[which(y==2),c(1,2,3)]
rownames(X2) = c(24:47)
for (i in 1:3) print(shapiro.test(X2[,i])$p.value)
# select the variables indices for which H_0 is rejected (<0.01)
index = c()
for (i in 1:3) {
  if (shapiro.test(X2[,i])$p.value < 0.01) {
    index = c(index,i)
  }
}

# 3d scatter plot
# Add a new column with color
mycolors <- c("firebrick1", 'royalblue1')
X$color <- mycolors[ as.numeric(y) ]

# Plot
plot3d(x=X$Ca, y=X$Mg, z=X$Na, col = X$color, 
       type = 's', 
       radius = 1.5,
       xlab="Ca", ylab="Mg", zlab="Na")


# Outlier detection in the groups
# cellMCD
cellmcd1 = cellWise::cellMCD(X1)
cellmcd2 = cellWise::cellMCD(X2)

# rowwise outliers
mcd1 = cov.mcd(X1)
mcd2 = cov.mcd(X2)
RMDs1 = mahalanobis(X1,center = mcd1$center, cov = mcd1$cov) # returns squared RMDs
RMDs2 = mahalanobis(X2,center = mcd2$center, cov = mcd2$cov)
c = qchisq(0.99, df = 3)
which(RMDs1>c) # rowwise outliers X1
which(RMDs2>c) # rowwise outliers X2

# get list of observations that contain at least one cellwise outliers.
cellwise_out <- as.integer(c(which(rowSums(cellmcd1$W)<3),
                             which(rowSums(cellmcd2$W)<3)+nrow(X1)))
# get list of rowwise outliers
rowwise_outliers = as.integer(c(which(RMDs1>c), 
                                which(RMDs2>c)+nrow(X1)))

# cellMaps for X1 and X2  
p1 = cellMap(cellmcd1$Zres, indrows = which(RMDs1>c), mTitle = "X1", drawCircles = T)
p2 = cellMap(cellmcd2$Zres, indrows = which(RMDs2>c), mTitle = "X2", drawCircles = T)
grid.arrange(grobs = list(p1,p2), ncol = 2)



### CLASSIFICATON RESULTS ON FOREST SOIL DATA

X = X[,c(1,2,3)]
y = as.factor(y)
# group sizes
n1 = length(which(y == 1)) #23
n2 = length(which(y == 2)) #24
# train the models
cqda = qda(X,y, prior = c(n1,n2)/n,"moment", CV = F)
rqda = vcr.da.train(X, y, rule = "QDA", estmethod = "DetMCD")
cellqda = cellQDA(X,y)

# predict
pred.cqda = predict(cqda, X)$class
pred.rqda = vcr.da.newdata(X,y, rqda)$predint
pred.cellqda = predict(cellqda, Xnew = X, ynew = y, method = "QDA")$predint
pred.nimp = predict(cellqda, Xnew = X, ynew = y, method = "NIMP")$predint
pred.l1 = predict(cellqda, Xnew = X, ynew = y, method = "L1")$predint

# percentage of correct predictions 
score = data.frame(cqda = length(which(pred.cqda==y)),
                   rqda = length(which(pred.rqda==y)),
                   cellqda = length(which(pred.cellqda==y)),
                   nimp = length(which(pred.nimp==y)),
                   l1 = length(which(pred.l1==y)))

score/n
