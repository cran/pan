########################################################################
# ECME example command file - works for both the Unix and Windows 
# versions of Splus. Before running this example, make sure that the
# ecme() function and object code are loaded into Splus; see the README
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
########################################################################
# Below we show how to fit a simple linear model to these data using
# ecme(). This will be a traditional repeated-measures style additive
# model with a fixed effect for each column (occasion) and a random 
# intercept for each subject. First we enter the data.
########################################################################
y <- c(16,20,16,20,-6,-4,
    12,24,12,-6,4,-8,
    8,8,26,-4,4,8,
    20,8,20,-4,
    8,4,-8,22,-8,
    10,20,28,-20,-4,-4,
    4,28,24,12,8,18,
    -8,20,24,-3,8,-24,
    20,24,8,12)
occ <- c(1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,5,6,
      1,2,3,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      2,3,4,5)
subj <- c(1,1,1,1,1,1,
       2,2,2,2,2,2,
       3,3,3,3,3,3,
       4,4,4,4,
       5,5,5,5,5,
       6,6,6,6,6,6,
       7,7,7,7,7,7,
       8,8,8,8,8,8,
       9,9,9,9)
########################################################################
# Now we must specify the model. 
# If the six measurements per subject were ordered in time, we might
# consider using a model with time of measurement entered with linear
# (or perhaps higher-order polynomial) effects. But because the
# six measurements are not clearly ordered, let's use a model that has
# an intercept and five dummy codes to allow the population means for
# the six occasions to be estimated freely. We will also allow the
# intercept to randomly vary by subject. For a subject i with no
# missing values, the covariate matrices will be
#
#                   1 1 0 0 0 0              1
#                   1 0 1 0 0 0              1
#           Xi =    1 0 0 1 0 0       Zi =   1
#                   1 0 0 0 1 0              1
#                   1 0 0 0 0 1              1
#                   1 0 0 0 0 0              1
#
# When using ecme(), these are combined into a single matrix called
# pred. The pred matrix has length(y) rows. Each column of Xi and Zi
# must be represented in pred. Because Zi is merely the first column
# of Xi, we do not need to enter that column twice. So pred is simply
# the matrices Xi (i=1,...,9), stacked upon each other.
#
pred <- cbind(int=rep(1,49),dummy1=1*(occ==1),dummy2=1*(occ==2),
    dummy3=1*(occ==3),dummy4=1*(occ==4),dummy5=1*(occ==5))
xcol <- 1:6
zcol <- 1
########################################################################
# Now we can fit the model.
result <- ecme(y,subj,occ,pred,xcol,zcol)
# In this example, ECME converged in 212 steps. The results can be
# viewed by printing out the various components of "result".
########################################################################
