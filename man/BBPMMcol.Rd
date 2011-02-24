\name{BBPMM}

\alias{BBPMM}

\title{(Multiple) Imputation through Bayesian Bootstrap Predictive Mean
  Matching (BBPMM)}

\description{\sQuote{BBPMM} performs single and multiple imputation
  (MI) of mixed-scale variables using a chained equations approach and
  (Bayesian Bootstrap) Predictive Mean Matching. 
}

\usage{BBPMM(data, M=10, nIter=10, outfile=NULL, ignore=NULL,
vartype=NULL, stepwise=TRUE, maxit=3, maxPerc = 0.98, verbose=TRUE, setSeed,...)}

\arguments{
  \item{data}{A partially incomplete data frame or matrix.}
  \item{M}{Number of multiple imputations. If M=1, no Bayesian
    Bootstrap step is carried out. Default=10.}
  \item{nIter}{Number of iterations of the chained equations algorithm
  before the data set is stored as an 'imputed data set'. Default=10.}
  \item{outfile}{A character string that specifies the path and
    file name for the imputed data sets. If \code{outfile}=NULL (default), no data set is
    stored}
  \item{ignore}{A character or numerical vector that specifies either
    column positions or variable names that are to be excluded from the
    imputation model and process, e.g. an ID variable. If
    \code{ignore}=NULL (default), all variables in \code{data} are used
    in the imputation model.}
  \item{vartype}{A character vector that flags the class of each
    variable in \code{data} (without the variables defined by the
    \code{ignore} argument), with either 'M' for
    metric-scale or 'C' for categorical. The default (NULL) takes over
    the classes of \code{data}. Overruling these classes can sometimes
    make sense: e.g. an ordinal-scale variable is originally classified as 'factor',
    but treating it as metric-scale variable within the imputation
    process might still be a better choice (considering the robust properties
    of PMM to model misspecification).}
  \item{stepwise}{Performs backwise selection for each imputation model based on the 
    Schwarz (Bayes) Information criterion. Default=TRUE.}
  \item{maxit}{Imported argument from the \code{nnet} package that
    specifies the maximum number of iterations for the multinomial logit
    model estimation. Default=3.}
  \item{maxPerc}{The maximum percentage the mode category of a variable is allowed to have in order to 
    try 'regular' imputation. If a variable is approximately Dirac distributed, i.e. if 
    it has (almost) no variance, imputation is carried out by simple hot deck imputation. 
    Default = 0.98.}
  \item{verbose}{The algorithm prints information on imputation and
    iteration numbers. Default=TRUE.}
  \item{setSeed}{Optional argument to fix the pseudo-random number generator in order to
    allow for reproducable results.}
  \item{...}{Further arguments passed to or from other functions.}
}

\value{
  \item{impdata}{The imputed data set, if M=1, or a list containing M
    imputed data sets.}
  \item{misOverview}{The percentage of missing values per incomplete variable.}
  \item{indMatrix}{A matrix with the same dimensions as \code{data}
    minus \code{ignore} containing flags for missing values.}
  \item{M}{Number of (multiple) imputations.}
  \item{nIter}{Number of iterations between two imputations.}
}

\details{\code{BBPMM} is based on a chained equations approach
  that is using a Bayesian Bootstrap approach and Predictive Mean
  Matching variants for metric-scale, binary, and multi-categorical
  variables to generate multiple imputations. In order to emulate a
  monotone missing-data pattern as well as possible, variables are sorted
  by rate of missingness (in ascending order). If no complete variables
  exist, the least incomplete variable is imputed via hot-deck. The
  starting solution then builds the imputation model using the observed values of
  a particular y variable, and the corresponding observed or already
  imputed values of the x variables (i.e. all variables with fewer
  missing values than y).
  Due to the PMM element in the algorithm,
  auto-correlation of subsequent iterations is virtually zero. Therefore, a
  burn-in period is not required, and there is no need to administer
  'high' values (> 20) to \code{nIter} either.
  
  If M=1, no Bayesian Bootstrap step is carried
  out for the chained equations. Note that in this case the algorithm is still unlikely to
  converge to a stable solution, because of the Predictive Mean Matching
  step.
}


\seealso{\code{\link{BBPMM.row}}}

\references{
  
  Koller-Meinfelder, F. (2009) \emph{Analysis of Incomplete Survey Data
    -- Multiple Imputation Via Bayesian Bootstrap Predictive Mean
    Matching}, doctoral thesis.  
}

\keyword{robust}
\keyword{datagen}

\examples{

### sample data set with non-normal variables
n <- 50
x1 <- round(runif(n,0.5,3.5))
x2 <- as.factor(c(rep(1,10),rep(2,25),rep(3,15)))
x3 <- round(rnorm(n,0,3))
y1 <- round(x1-0.25*(x2==2)+0.5*x3+rnorm(n,0,1))
y1 <- ifelse(y1<1,1,y1)
y1 <- as.factor(ifelse(y1>4,5,y1))
y2 <- x1+rnorm(n,0,0.5)
y3 <- round(x3+rnorm(n,0,2))
data1 <- as.data.frame(cbind(x1,x2,x3,y1,y2,y3))
misrow1 <- sample(n,20)
misrow2 <- sample(n,15)
misrow3 <- sample(n,10)
is.na(data1[misrow1, 4]) <- TRUE
is.na(data1[misrow2, 5]) <- TRUE
is.na(data1[misrow2, 6]) <- TRUE

### imputation
imputed.data <- BBPMM(data1, nIter=5, M=5)
}
