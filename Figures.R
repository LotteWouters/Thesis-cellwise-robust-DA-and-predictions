
### FIGURES ###

library(ggforce) # for geom_ellipse
library(ggplot2)
library(emdbook) # (curve3d) for plotting contours
library(cellWise)


### FUNCTIONS ###

# also run : cellsimdata() (from simulationData1.0.R), generateData2.R and
#            cellMCD() (from Auxiliary_Functions_cellQDA.R)


# function for getting data for the tolerance ellipsoid
tol_data = function(sigma){ 
  c = 3.03 # sqrt(chisquare_df=2_alpha=0.99)
  eig = eigen(sigma)
  root_l = sqrt(eig$values) # sqrt l1 and l2
  v1 = eig$vectors[,1]
  angle = atan(v1[2]/v1[1]) # direction of tol.ellips = angle(v1 , pos.x-axis)
  return(list(a = c*root_l[1], b = c*root_l[2], angle = angle))
}

# function for preparing decision boundary
decision_boundary <- function(x1,x2,mu,covar1,covar2){
  
  # mu = list containing two vectors, the estimated mu1_hat, mu2_hat
  # covar1 = estimated covariance for group1
  # covar2 = estimated covariance for group2
  # returns: pi_1*f1(x) - pi_2*f2(x) as a function of x = (x1,x2)
  
  # bivariate normal density function
  phi2 <- function(x1,x2,mu1,mu2,a,b,c,d) {
    # a,b,c,d = elements of cov
    det_inv = 1/(a*d-c*b)
    z1 = x1-mu1
    z2 = x2-mu2
    return(sqrt(det_inv) * exp(-0.5 * det_inv * (d*z1^2 - (c+b)*z1*z2 + a*z2^2)) )
  }
  # group conditional densities
  f1 = phi2(x1,x2,mu[[1]][1],mu[[1]][2],
            covar1[1,1],covar1[1,2],covar1[2,1],covar1[2,2])
  f2 = phi2(x1,x2,mu[[2]][1],mu[[2]][2],
            covar2[1,1],covar2[1,2],covar2[2,1],covar2[2,2])
  return(priors[1]*f1 - priors[2]*f2)
}

# function for adding cellwise outliers to existing data (copied piece from generateData2)
add.cell = function(X, m, s, gamma, perout){
  
  # X := data matrix (2d)
  # mu, sigma := actual mean and covariance of X
  # gamma := outlyingness, between 1 and 10
  # perout := proportion of contamination to be added
  
  indcells = c()
  n = nrow(X)
  ind = replicate(p, sample(seq_len(n), perout * n, replace = FALSE))
  ind = as.vector(t(t(ind) + n * (0:(p - 1))))
  W = array(0, dim(X))
  W[ind] = 1 # the cells that will be made outlying
  for (i in seq_len(n)) {
    continds = which(W[i, ] == 1)
    if (length(continds) > 0) {
      eigen_out <- eigen(s[continds, continds])$vectors
      replacement <- eigen_out[, length(continds)]/
        sqrt(mahalanobis(eigen_out[, length(continds)], rep(0,length(continds)), s[continds, continds]))
      X[i, continds] <- m[continds] + replacement * gamma * sqrt(length(continds))
      indcells <- c(indcells, i + (continds - 1) * n)
    }
  }
  return(list(X = X, W = W))
}


### FIG 1.1: CQDA without outliers ###

# initializing parameters
p = 2
K = 2
nk = c(100, 120)
priors = nk/sum(nk) 

# generate data
mu = list(c(3,3), c(6,5))
covar1 = matrix(c(0.8,-.5, -.5, 0.8), nrow = 2, ncol = 2)
covar2 = diag(c(1.2,0.5))
sigma = list(covar1, covar2)
set.seed(4)
sim.out = cellsimdata(K,p,nk,mu,sigma,perout = c(0,0), gamma = 0, outlierType = "cellwiseStructured")
X = data.frame(sim.out$X[-212,]) 
y = sim.out$y[-212]

# addapt y so that they represent the color of that point in the plot
y = ifelse(y == 1,'firebrick1', 'royalblue1')
rownames(X) = rownames(y) = NULL #reset row numbers
X1 = X[1:nk[1],] # group 1
X2 = X[(nk[1]+1):nrow(X),] # group 2

# get sample mean and cov
s1_cqda = cov(X1)
x_bar1_cqda = colMeans(X1)
s2_cqda = cov(X2)
x_bar2_cqda = colMeans(X2)
tol_data1_cqda = tol_data(s1_cqda)
tol_data2_cqda = tol_data(s2_cqda)

# prepare plot of data (X) and tolerance ellipsoids
p1 = ggplot(X, aes(x=X1, y=X2, color = y)) + 
  geom_point(shape=19, size = 2) +
  geom_ellipse(aes(x0 = x_bar1_cqda[1], y0 = x_bar1_cqda[2], 
                   a = tol_data1_cqda$a, b = tol_data1_cqda$b,  
                   angle = tol_data1_cqda$angle), color = 'firebrick1') +
  geom_ellipse(aes(x0 = x_bar2_cqda[1], y0 = x_bar2_cqda[2], 
                   a = tol_data2_cqda$a, b = tol_data2_cqda$b,  
                   angle = tol_data2_cqda$angle), color = 'royalblue1') +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_identity() + # to overwrite coloring from ggplot with my own
  ggtitle("CQDA") +
  xlab(bquote(X[1])) +
  ylab(bquote(X[2]))
p1

# get contour datapoints for decision boundary (CQDA)
cntr_cqda = curve3d(decision_boundary(x,y,list(x_bar1_cqda, x_bar2_cqda),s1_cqda,s2_cqda),
                    xlim=c(0,11),ylim=c(0,8), sys3d="none")
# reshape data for plotting
dimnames(cntr_cqda$z) = list(cntr_cqda$x,cntr_cqda$y)
mm_cqda = reshape2::melt(cntr_cqda$z)

# final plot 
p1 + geom_contour(data=mm_cqda,aes(x=Var1,y=Var2,z=value),breaks=0,
                  colour="black") 


### FIG 1.3: CQDA with outliers ###

# 5% outliers 

# generate clean data (same as code FIG1)
# initializing parameters
p = 2
K = 2
nk = c(100, 120)
priors = nk/sum(nk) 

# generate data
mu = list(c(3,3), c(6,5))
covar1 = matrix(c(0.8,-.5, -.5, 0.8), nrow = 2, ncol = 2)
covar2 = diag(c(1.2,0.5))
sigma = list(covar1, covar2)
set.seed(4)  
sim.out = cellsimdata(K,p,nk,mu,sigma,perout = c(0,0), gamma = 0, outlierType = "cellwiseStructured")
X = data.frame(sim.out$X[-212,]) 
y = sim.out$y[-212]

# addapt y so that they represent the color of that point in the plot
y = ifelse(y == 1,'firebrick1', 'royalblue1')
rownames(X) = rownames(y) = NULL #reset row numbers
X1 = X[1:nk[1],] # group 1
X2 = X[(nk[1]+1):nrow(X),] # group 2

# add cellwise outliers to X1 and X2 separately 
set.seed(2) # good one
add.cell.out1 = add.cell(X1, mu[[1]], sigma[[1]], gamma = 5, perout = 0.05)
X1.05 = add.cell.out1$X
add.cell.out2 = add.cell(X2, mu[[2]], sigma[[2]], gamma = 5, perout = 0.05)
X2.05 = add.cell.out2$X
X.05 = data.frame(rbind(X1.05,X2.05))
y.05 = y 

# compute location and scale according to CQDR (sample mean and cov)
s1_05 = cov(X1.05)
x_bar1_05 = colMeans(X1.05)
s2_05 = cov(X2.05)
x_bar2_05 = colMeans(X2.05)
tol_data1_05 = tol_data(s1_05)
tol_data2_05 = tol_data(s2_05)

p2 = ggplot(X.05, aes(x=X1, y=X2, color = y.05)) + 
  geom_point(shape=19, size = 2) +
  geom_ellipse(aes(x0 = x_bar1_05[1], y0 = x_bar1_05[2], 
                   a = tol_data1_05$a, b = tol_data1_05$b,  
                   angle = tol_data1_05$angle), color = 'firebrick1') +
  geom_ellipse(aes(x0 = x_bar2_05[1], y0 = x_bar2_05[2], 
                   a = tol_data2_05$a, b = tol_data2_05$b,  
                   angle = tol_data2_05$angle), color = 'royalblue1') +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_identity() + # to overwrite coloring from ggplot with my own
  ggtitle("CQDA (5% cellwise outliers)") +
  xlab(bquote(X[1])) +
  ylab(bquote(X[2]))
p2

# get contour datapoints for decision boundary (CQDA no contamination)
cntr_cqda = curve3d(decision_boundary(x,y,list(x_bar1_cqda, x_bar2_cqda),s1_cqda,s2_cqda),
                    xlim=c(-2,13),ylim=c(-2,10), sys3d="none")
# reshape data for plotting
dimnames(cntr_cqda$z) = list(cntr_cqda$x,cntr_cqda$y)
mm_cqda = reshape2::melt(cntr_cqda$z)

# get contour datapoints for decision boundary (CQDA 5% outliers)
cntr_cqdr.05 = curve3d(decision_boundary(x,y,list(x_bar1_05, x_bar2_05),s1_05,s2_05),
                       xlim=c(-2,13),ylim=c(-2,10), sys3d="none")
# reshape data for plotting
dimnames(cntr_cqdr.05$z) = list(cntr_cqdr.05$x,cntr_cqdr.05$y)
mm_cqdr.05 = reshape2::melt(cntr_cqdr.05$z)

# final plot
p2 + geom_contour(data=mm_cqda,aes(x=Var1,y=Var2,z=value),breaks=0,
                  colour="black", linetype = "dashed") + 
  geom_contour(data=mm_cqdr.05,aes(x=Var1,y=Var2,z=value),breaks=0,
               colour="black")


# 10% outliers

set.seed(4) 
add.cell.out1 = add.cell(X1, mu[[1]], sigma[[1]], gamma = 5, perout = 0.10)
X1.10 = add.cell.out1$X
add.cell.out2 = add.cell(X2, mu[[2]], sigma[[2]], gamma = 5, perout = 0.10)
X2.10 = add.cell.out2$X
X.10 = data.frame(rbind(X1.10,X2.10))
y.10 = y

# compute location and scale according to CQDA (sample mean and cov)
s1_10 = cov(X1.10)
x_bar1_10 = colMeans(X1.10)
s2_10 = cov(X2.10)
x_bar2_10 = colMeans(X2.10)
tol_data1_10 = tol_data(s1_10)
tol_data2_10 = tol_data(s2_10)

p3 = ggplot(X.10, aes(x=X1, y=X2, color = y.10)) + 
  geom_point(shape=19, size = 2) +
  geom_ellipse(aes(x0 = x_bar1_10[1], y0 = x_bar1_10[2], 
                   a = tol_data1_10$a, b = tol_data1_10$b,  
                   angle = tol_data1_10$angle), color = 'firebrick1') +
  geom_ellipse(aes(x0 = x_bar2_10[1], y0 = x_bar2_10[2], 
                   a = tol_data2_10$a, b = tol_data2_10$b,  
                   angle = tol_data2_10$angle), color = 'royalblue1') +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_identity() + # to overwrite coloring from ggplot with my own
  ggtitle("CQDA (10% cellwise outliers)") +
  xlab(bquote(X[1])) +
  ylab(bquote(X[2]))
p3

# get contour datapoints for decision boundary (CQDA no contamination)
cntr_cqda = curve3d(decision_boundary(x,y,list(x_bar1_cqda, x_bar2_cqda),s1_cqda,s2_cqda),
                    xlim=c(-2,13),ylim=c(-2,10), sys3d="none")
# reshape data for plotting
dimnames(cntr_cqda$z) = list(cntr_cqda$x,cntr_cqda$y)
mm_cqda = reshape2::melt(cntr_cqda$z)

# get contour datapoints for decision boundary (CQDA 10% outliers)
cntr_cqdr.10 = curve3d(decision_boundary(x,y,list(x_bar1_10, x_bar2_10),s1_10,s2_10),
                       xlim=c(-2,13),ylim=c(-2,10), sys3d="none")
# reshape data for plotting
dimnames(cntr_cqdr.10$z) = list(cntr_cqdr.10$x,cntr_cqdr.10$y)
mm_cqdr.10 = reshape2::melt(cntr_cqdr.10$z)

# final plot
p3 + geom_contour(data=mm_cqda,aes(x=Var1,y=Var2,z=value),breaks=0,
                  colour="black", linetype = "dashed") + 
  geom_contour(data=mm_cqdr.10,aes(x=Var1,y=Var2,z=value),breaks=0,
               colour="black")



### FIG 2.1: cellQDA with outliers (same dataset as fig 1.3) ###

# 0% cellwise outliers 

# data: X, y, X1, X2
cellmcd1 = cellWise::cellMCD(X1, alpha = 0.98) 
cellmcd2 = cellWise::cellMCD(X2, alpha = 0.94)
s1_cell = cellmcd1$S
x_bar1_cell = cellmcd1$mu
s2_cell = cellmcd2$S
x_bar2_cell = cellmcd2$mu
tol_data1_cell = tol_data(s1_cell)
tol_data2_cell = tol_data(s2_cell)

# plot 
p4 = ggplot(X, aes(x=X1, y=X2, color = y)) + 
  geom_point(shape=19, size = 2) +
  geom_ellipse(aes(x0 = x_bar1_cell[1], y0 = x_bar1_cell[2], 
                   a = tol_data1_cell$a, b = tol_data1_cell$b,  
                   angle = tol_data1_cell$angle), color = 'firebrick1') +
  geom_ellipse(aes(x0 = x_bar2_cell[1], y0 = x_bar2_cell[2], 
                   a = tol_data2_cell$a, b = tol_data2_cell$b,  
                   angle = tol_data2_cell$angle), color = 'royalblue1') +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_identity() + # to overwrite coloring from ggplot with my own
  ggtitle("cellQDA (0% cellwise outliers)") +
  xlab(bquote(X[1])) +
  ylab(bquote(X[2]))
p4

# get contour datapoints for decision boundary (cellQDA)
cntr_cellqdr = curve3d(decision_boundary(x,y,list(x_bar1_cell, x_bar2_cell),s1_cell,s2_cell),
                       xlim=c(0,9),ylim=c(0,8), sys3d="none")
# reshape data for plotting
dimnames(cntr_cellqdr$z) = list(cntr_cellqdr$x,cntr_cellqdr$y)
mm_cellqdr = reshape2::melt(cntr_cellqdr$z)

# final plot 
p4 + geom_contour(data=mm_cellqdr,aes(x=Var1,y=Var2,z=value),breaks=0,
                  colour="black") 


# 5% cellwise outliers

# data: X.05, y, X1.05, X2.05
cellmcd1 = cellWise::cellMCD(X1.05, alpha = 0.90) 
cellmcd2 = cellWise::cellMCD(X2.05, alpha = 0.90)
s1_cell = cellmcd1$S
x_bar1_cell = cellmcd1$mu
s2_cell = cellmcd2$S
x_bar2_cell = cellmcd2$mu
tol_data1_cell = tol_data(s1_cell)
tol_data2_cell = tol_data(s2_cell)

# plot 
p5 = ggplot(X.05, aes(x=X1, y=X2, color = y)) + 
  geom_point(shape=19, size = 2) +
  geom_ellipse(aes(x0 = x_bar1_cell[1], y0 = x_bar1_cell[2], 
                   a = tol_data1_cell$a, b = tol_data1_cell$b,  
                   angle = tol_data1_cell$angle), color = 'firebrick1') +
  geom_ellipse(aes(x0 = x_bar2_cell[1], y0 = x_bar2_cell[2], 
                   a = tol_data2_cell$a, b = tol_data2_cell$b,  
                   angle = tol_data2_cell$angle), color = 'royalblue1') +
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_identity() + # to overwrite coloring from ggplot with my own
  ggtitle("cellQDA (5% cellwise outliers)") +
  xlab(bquote(X[1])) +
  ylab(bquote(X[2]))
p5

# get contour datapoints for decision boundary (cellQDA 5% outliers)
cntr_cellqdr.05 = curve3d(decision_boundary(x,y,list(x_bar1_cell, x_bar2_cell),s1_cell,s2_cell),
                          xlim=c(0,9),ylim=c(0,8), sys3d="none")
# reshape data for plotting
dimnames(cntr_cellqdr.05$z) = list(cntr_cellqdr.05$x,cntr_cellqdr.05$y)
mm_cellqdr.05 = reshape2::melt(cntr_cellqdr.05$z)

# final plot 
p5 + geom_contour(data=mm_cellqdr.05,aes(x=Var1,y=Var2,z=value),breaks=0,
                  colour="black") +
  geom_contour(data=mm_cellqdr,aes(x=Var1,y=Var2,z=value),breaks=0,
               colour="black", linetype = "dashed")




# 10% cellwise outliers

# data: X.10, y, X1.10, X2.10
cellmcd1 = cellWise::cellMCD(X1.10, alpha = 0.85)
cellmcd2 = cellWise::cellMCD(X2.10, alpha = 0.85)

s1_cell = cellmcd1$S
x_bar1_cell = cellmcd1$mu
s2_cell = cellmcd2$S
x_bar2_cell = cellmcd2$mu
tol_data1_cell = tol_data(s1_cell)
tol_data2_cell = tol_data(s2_cell)

# plot 
p6 = ggplot(X.10, aes(x=X1, y=X2, color = y)) + 
  geom_point(shape=19, size = 2) +
  geom_ellipse(aes(x0 = x_bar1_cell[1], y0 = x_bar1_cell[2], 
                   a = tol_data1_cell$a, b = tol_data1_cell$b,  
                   angle = tol_data1_cell$angle), color = 'firebrick1') +
  geom_ellipse(aes(x0 = x_bar2_cell[1], y0 = x_bar2_cell[2], 
                   a = tol_data2_cell$a, b = tol_data2_cell$b,  
                   angle = tol_data2_cell$angle), color = 'royalblue1') +
  theme_minimal() + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_identity() + # to overwrite coloring from ggplot with my own
  ggtitle("cellQDA (10% cellwise outliers)") +
  xlab(bquote(X[1])) +
  ylab(bquote(X[2]))
p6

# get contour datapoints for decision boundary 
cntr_cellqdr.10 = curve3d(decision_boundary(x,y,list(x_bar1_cell, x_bar2_cell),s1_cell,s2_cell),
                          xlim=c(0,9),ylim=c(0,8), sys3d="none")
# reshape data for plotting
dimnames(cntr_cellqdr.10$z) = list(cntr_cellqdr.10$x,cntr_cellqdr.10$y)
mm_cellqdr.10 = reshape2::melt(cntr_cellqdr.10$z)

# final plot 
p6 + geom_contour(data=mm_cellqdr.10,aes(x=Var1,y=Var2,z=value),breaks=0,
                  colour="black") +
  geom_contour(data=mm_cellqdr,aes(x=Var1,y=Var2,z=value),breaks=0,
               colour="black", linetype = "dashed")


