\name{pan}
\alias{pan}
\title{
Imputation of multivariate panel or cluster data
}
\description{
Gibbs sampler for the multivariate linear mixed model with
incomplete data described by Schafer (1997). This function
will typically be used to produce multiple imputations of
missing data values in multivariate panel data or clustered
data. The underlying model is

yi = Xi\%*\%beta + Zi\%*\%bi + ei,    i=1,\dots,m,

where

yi    = (ni x r) matrix of incomplete multivariate
        data for subject or cluster i;

Xi    = (ni x p) matrix of covariates;

Zi    =	(ni x q) matrix of covariates;

beta  = (p x r) matrix of coefficients common to the
        population (fixed effects);

bi    = (q x r) matrix of coefficients specific to
        subject or cluster i (random effects); and

ei    = (ni x r) matrix of residual errors.

The matrix bi, when stacked into a single column, is assumed
to be normally distributed with mean zero and unstructured
covariance matrix psi, and the rows of ei are assumed to be
independently normal with mean zero and unstructured
covariance matrix sigma. Missing values may appear in yi in
any pattern.

In most applications of this model, the first columns of Xi
and Zi will be constant (one) and Zi will contain a subset of
the columns of Xi. 
}

\usage{
pan(y, subj, pred, xcol, zcol, prior, seed, iter=1, start)
}

\arguments{
\item{y}{
matrix of responses. This is simply the individual yi matrices
stacked upon one another. Each column of y corresponds to a
response variable. Each row of y corresponds to a single
subject-occasion, or to a single subject within a cluster.
Missing values (NA) may occur in any pattern.
}

\item{subj}{
vector of length nrow(y) giving the subject (or cluster)
indicators i for the rows of y. For example, suppose 
that y is in fact rbind(y1,y2,y3,y4) where nrow(y1)=2,
nrow(y2)=3, nrow(y3)=2, and nrow(y4)=7. Then subj should
be c(1,1,2,2,2,3,3,4,4,4,4,4,4,4).
}

\item{pred}{
matrix of covariates used to predict y. This should have the
same number of rows as y. The first column will typically be
constant (one), and the remaining columns correspond to other
variables appearing in Xi and Zi.
}

\item{xcol}{
vector of integers indicating which columns of pred will be
used in Xi. That is, pred[,xcol] is the Xi matrices (stacked
upon one another).
}

\item{zcol}{
vector of integers indicating which columns of pred will be
used in Zi. That is, pred[,zcol] is the Zi matrices (stacked
upon one another).
}

\item{prior}{
a list with four components (whose names are a, Binv, c, and
Dinv, respectively) specifying the hyperparameters of the 
prior distributions for psi and sigma. For information on how
to specify and interpret these hyperparameters, see Schafer
(1997) and the example command file "panex.R" distibuted with
this package. Note: This is a slight departure from the
notation in Schafer (1997), where a and Binv were denoted
by "nu1" and "Lambdainv1", and c and Dinv were "nu2" and
"Lambdainv2".
}

\item{seed}{
integer seed for initializing pan()'s internal random number
generator. This argument should be a positive integer. 
}

\item{iter}{
total number of iterations or cycles of the Gibbs sampler
to be carried out.
}

\item{start}{
optional list of quantities to specify the initial state of
the Gibbs sampler. This list has the same form as "last"
(described below), one of the components returned by pan(). 
This argument allows the Gibbs sampler to be restarted from
the final state of a previous run. If "start" is omitted then
pan() chooses its own initial state.
}}

\value{
A list containing the following components. Note that when you
are using pan() to produce multiple imputations, you will
be primarily interested in the component "y" which contains
the imputed data; the arrays "beta", "sigma", and "psi" will
be used primarily for diagnostics (e.g. time-series plots)
to assess the convergence behavior of the Gibbs sampler.

\item{beta}{
array of dimension c(length(xcol),ncol(y),iter) = (p x r x 
number of Gibbs cycles) containing the simulated values of
beta from all cycles. That is, beta[,,T] is the (p x r) matrix of
simulated fixed effects at cycle T.
}

\item{sigma}{
array of dimension c(ncol(y),ncol(y),iter) = (r x r x
number of Gibbs cycles) containing the simulated values of
sigma from all cycles. That is, sigma[,,T] is the simulated version
of the model's sigma at cycle T.
}

\item{psi}{
array of dimension c(length(zcol)*ncol(y), length(zcol)*ncol(y), iter)
= (q*r x q*r x number of Gibbs cycles) containing the simulated values
of psi from all cycles. That is, psi[,,T] is the simulated version of
the model's psi at cycle T.
}

\item{y}{
matrix of imputed data from the final cycle of the Gibbs
sampler. Identical to the input argument y except that the
missing values (NA) have been replaced by imputed values.
If "iter" has been set large enough (which can be determined by
examining time-series plots, etc. of "beta", "sigma", and
"psi") then this is a proper draw from the posterior
predictive distribution of the complete data.
}

\item{last}{
a list of four components characterizing the final state
of the Gibbs sampler. The four components are: "beta", 
"sigma", "psi", and "y", which are the simulated values
of the corresponding model quantities from the final cycle of
Gibbs. This information is already contained in the other
components returned by pan(); we are providing this list merely
as a convenience, to allow the user to start future runs of
the Gibbs sampler at this state.
}}

\details{
The Gibbs sampler algorithm used in pan() is described in
detail by Schafer (1997).
}

\note{
This function assumes that the rows of y (and thus the rows
of subj and pred) have been sorted by subject number. That is,
we assume that subj=sort(subj), y=y[order(subj),], and
pred=pred[order(subj),]. If the matrix y is created by
stacking yi, i=1,\dots,m then this will automatically be the case.
}

\references{
Schafer, J.L. (1997) Imputation of missing covariates under
a multivariate linear mixed model. Technical report, Dept. of
Statistics, The Pennsylvania State University.
}

\examples{
\dontrun{
For a detailed example, see the file "panex.R" distributed
with this function. Here is a simple example of how pan()
might be used to produce three imputations.

# run Gibbs for 1000 cycles
result <- pan(y,subj,pred,xcol,zcol,prior,seed=9565,iter=1000)           
# first imputation
imp1 <- result$y
# another 1000 cycles
result <- pan(y,subj,pred,xcol,zcol,prior,seed=54324,iter=1000,start=result$last)
# second imputation
imp2 <- result$y
# another 1000 cycles
result <- pan(y,subj,pred,xcol,zcol,prior,seed=698212,iter=1000,start=result$last)
# third imputation
imp3 <- result$y
}}

\keyword{models}

% Converted by Sd2Rd version 1.21.
