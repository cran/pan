########################################################################
# PAN example command file - works for both the Unix and Windows 
# versions of Splus. Before running this example, make sure that the
# pan() function and object code are loaded into Splus; see the README
# file for details.
########################################################################
# These data, taken from an article by Weil, A.T., Zinberg, N.E. and
# Nelson, J.M. (1968; Clinical and psychological effects of marihuana
# in man; Science, 162, 1234-1242), come from a pilot study to
# investigate the clinical and psychological effects of marijuana use
# in human subjects. Nine subjects subjects each received three
# treatments---low-dose, high-dose, and placebo. Under each treatment,
# changes in heart rate (beats per minute) were measured 15 and 90
# minutes after administration. NA denotes a missing value.
#
#   -----------------------------------------------------------------
#                         15 minutes                 90 minutes
#                   ----------------------     ----------------------
#                   Placebo   Low   High      Placebo   Low   High
#   -----------------------------------------------------------------
#   Subject 1          16     20     16          20     -6     -4       
#           2          12     24     12          -6      4     -8
#           3           8      8     26          -4      4      8
#           4          20      8     NA          NA     20     -4
#           5           8      4     -8          NA     22     -8
#           6          10     20     28         -20     -4     -4
#           7           4     28     24          12      8     18
#           8          -8     20     24          -3      8    -24
#           9          NA     20     24           8     12     NA
#   -----------------------------------------------------------------
#
library(pan)
library(ts)
########################################################################
# Below we show how to multiply impute the missing values using pan(). 
# This example is somewhat atypical because the data consist of a
# single response variable (change in heart rate) measured repeatedly;
# most uses of pan() will involve r > 1 response variables. If we had
# r response variables rather than one, the only difference would be
# that the vector y below would become a matrix with r columns, one
# for each response variable. The dimensions of Sigma (the residual
# covariance matrix for the response) and Psi (the covariance matrix
# for the random effects) would also change to (r x r) and (r*q x r*q),
# respectively, where q is the number of random coefficients in the
# model (in this case q=1 because we have only random intercepts). The
# new dimensions for Sigma and Psi will be reflected in the prior
# distribution, as Dinv and Binv become (r x r) and (r*q x r*q).  
########################################################################
# First we need to enter the response data into a matrix (or, as in
# this case, a vector) with one column for each variable. The first
# six rows contain the data for the first subject; the second six
# contain data for the second subject; and so on. The total number of 
# rows in y is the number of subject-occasions, which in this case is 54.
y <- c(16,20,16,20,-6,-4,
    12,24,12,-6,4,-8,
    8,8,26,-4,4,8,
    20,8,NA,NA,20,-4,
    8,4,-8,NA,22,-8,
    10,20,28,-20,-4,-4,
    4,28,24,12,8,18,
    -8,20,24,-3,8,-24,
    NA,20,24,8,12,NA)
########################################################################
# Next we create the vector subj whose length is equal to the number
# of rows in y. This vector indicates which rows of y belong to the
# first subject, which rows belong to the second subject, etc.
subj <- c(1,1,1,1,1,1,
       2,2,2,2,2,2,
       3,3,3,3,3,3,
       4,4,4,4,4,4,
       5,5,5,5,5,5,
       6,6,6,6,6,6,
       7,7,7,7,7,7,
       8,8,8,8,8,8,
       9,9,9,9,9,9)
########################################################################
# Now we must specify the model to be used for imputation.
# If the six measurements per subject were ordered in time, we might
# consider using a model with time of measurement entered with linear
# (or perhaps higher-order polynomial) effects. But because the
# six measurements are not clearly ordered, let's use a model that has
# an intercept and five dummy codes to allow the population means for
# the six occasions to be estimated freely. We will also allow the
# intercept to randomly vary by subject. For each subject i, the
# covariate matrices are then:
#
#                   1 1 0 0 0 0              1
#                   1 0 1 0 0 0              1
#           Xi =    1 0 0 1 0 0       Zi =   1
#                   1 0 0 0 1 0              1
#                   1 0 0 0 0 1              1
#                   1 0 0 0 0 0              1
#
# When using pan(), these are combined into a single matrix called
# pred. The pred matrix has the same number of rows as y, the number of
# subject-occasions. Each column of Xi and Zi must be represented in
# pred. Because Zi is merely the first column of Xi, we do not need to
# enter that column twice. So pred is simply the matrix Xi, stacked
# upon itself nine times.
#
pred <- cbind(int=rep(1,54),
   dummy1=rep(c(1,0,0,0,0,0),9),
   dummy2=rep(c(0,1,0,0,0,0),9),
   dummy3=rep(c(0,0,1,0,0,0),9),
   dummy4=rep(c(0,0,0,1,0,0),9),
   dummy5=rep(c(0,0,0,0,1,0),9))
#
# Now we must tell pan that all six columns of pred are to be used in
# Xi, but only the first column of pred appears in Zi.
#
xcol <- 1:6
zcol <- 1
########################################################################
# The model specification is now complete. The only task that remains
# is to specify the prior distributions for the covariance matrices
# Sigma and Psi.
#
# Recall that the dimension of Sigma is (r x r) where r
# is the number of response variables (in this case, r=1). The prior
# distribution for Sigma is inverted Wishart with hyperparameters a 
# (scalar) and Binv (r x r), where a is the imaginary degrees of freedom
# and Binv/a is the prior guesstimate of Sigma. The value of a must be
# greater than or equal to r. The "least informative" prior possible
# would have a=r, so here we will take a=1. As a prior guesstimate of 
# Sigma we will use the (r x r) identity matrix, so Binv = 1*1 = 1.
#
# By similar reasoning we choose the prior distribution for Psi. The
# dimension of Psi is (r*q x r*q) where q is the number of random
# effects in the model (i.e. the length of zcol, which in this case is
# one). The hyperparameters for Psi are c and Dinv, where c is the
# imaginary degrees of freedom (which must be greater than or equal to
# r*q) and Cinv/d is the prior guesstimate of Psi. We will take d=1
# and Cinv=1*1 = 1.
#
# The prior is specified as a list with four components named a, Binv,
# c, and Dinv, respectively.
#
prior <- list(a=1,Binv=1,c=1,Dinv=1)
########################################################################
# Now we are ready to run pan(). Let's assume that the pan function
# and the object code have already been loaded into Splus. First we
# do a preliminary run of 1000 iterations. 
#
result <- pan(y,subj,pred,xcol,zcol,prior,seed=13579,iter=1000)
#
# Check the convergence behavior by making time-series plots and acfs
# for the model parameters. Variances will be plotted on a log
# scale. We'll assume that a graphics device has already been opened.
#
plot(1:1000,log(result$sigma[1,1,]),type="l")
acf(log(result$sigma[1,1,]))
plot(1:1000,log(result$psi[1,1,]),type="l")
acf(log(result$psi[1,1,]))
par(mfrow=c(3,2))
for(i in 1:6) plot(1:1000,result$beta[i,1,],type="l")
for(i in 1:6) acf(result$beta[i,1,])
#
# This example appears to converge very rapidly; the only appreciable
# autocorrelations are found in Psi, and even those die down by lag
# 10. With a sample this small we can afford to be cautious, so let's
# impute the missing data m=10 times taking 100 steps between
# imputations. We'll use the current simulated value of y as the first
# imputation, then restart the chain where we left off to produce
# the second thru the tenth.
#
y1 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=9565,iter=100,start=result$last)
y2 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=6047,iter=100,start=result$last)
y3 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=3955,iter=100,start=result$last)
y4 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=4761,iter=100,start=result$last)
y5 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=9188,iter=100,start=result$last)
y6 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=9029,iter=100,start=result$last)
y7 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=4343,iter=100,start=result$last)
y8 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=2372,iter=100,start=result$last)
y9 <- result$y
result <- pan(y,subj,pred,xcol,zcol,prior,seed=7081,iter=100,start=result$last)
y10 <- result$y
########################################################################

