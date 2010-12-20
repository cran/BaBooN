\name{BBPMM.row}

\alias{BBPMM.row}

\title{(Multiple) Imputation of variable vectors}

\description{\sQuote{BBPMM.row} performs single and multiple imputation (MI) of metric scale
  variable vectors. For MI, parameter draws from a posterior
  distribution are replaced by a Bayesian Bootstrap step. Imputations
  are generated using Predictive Mean Matching (PMM) as described in
  Little (1988). 
}

\usage{BBPMM.row(misDataPat, blockImp=length(misDataPat$blocks),
 M=10, outfile=NULL, manWeights=NULL, stepwise=TRUE, verbose=TRUE, 
 tol=0.25, setSeed=NULL, ...)}

\arguments{
  \item{misDataPat}{An object created by \code{rowimpPrep} that
    contains information on all identified missing-data patterns.}
  \item{blockImp}{A scalar or vector containing the number(s) of the
    block(s) considered for imputation. Per default only the last
    block is imputed.}
  \item{M}{Number of multiple imputations. If M=1, no Bayesian
    Bootstrap step is carried out.}
  \item{outfile}{A character string that specifies the path and
    file name for the imputed data sets. If \code{outfile}=NULL (default), no data set is
    stored.}
  \item{manWeights}{Optional argument containing manual (non-negative) weights
    for the PMM step. \code{manWeights} can either be a list
    containing a vector for each missingness pattern, or just a
    vector, if only one missingness pattern/block exists. In
    either case, the number of elements in the vector(s) must
    match the number of variables in the corresponding block. Note
    that the higher the weight the higher the importance of a good
    match for the corresponding variable's predictive means. }
  \item{stepwise}{Performs backwise selection for each imputation model based on the 
    Schwarz (Bayes) Information criterion. Default=TRUE.}
  \item{verbose}{The algorithm prints information on weighting
    matrices and imputation numbers. Default=TRUE.}
  \item{tol}{Imported argument from function \code{qr} that specifies the
    tolerance level for linear dependencies among the complete
    variables and defaults to 0.25.}
  \item{setSeed}{Optional argument to fix the pseudo-random number generator in order to
    allow for reproducable results.}
  \item{...}{Further arguments passed to or from other functions.}
}

\value{
  \item{impdata}{A list containing M completed data sets.}
  \item{weightMatrix}{A list containing weight matrices for all imputations and blocks.}
  \item{model}{A list containing the lm-objects for all imputations and blocks.}
  \item{pairlist}{A list containing the donor/recipient pairlist data frames for all imputations and blocks.}
  \item{indMatrix}{A matrix with the same dimensions as \code{data}
    containing flags for missing values.}
}

\details{The simultaneous
  imputation of several variables is useful for missing-by-design
  patterns, such as data fusion or split questionnaire designs. The
  predictive means of the imputation variables are weighted by the
  inverse of the covariance matrix of the residuals from the regression
  of these variables on the complete variables. The intuitive idea
  behind is that distances between predictive means should be punished
  more severely, if the particular variable can be explained well by the
  (completely observed) imputation model variables.
  Through partialization and subsequent usage of the residuals the
  weight matrix is transformed into a diagonal matrix.
  The calculated weights can be adjusted by manual weights. Since the weight matrix
  is a Mahalanobis type of distance matrix, the weights are in the denominator and
  therefore the lower the weight, the higher the influence. As this is somewhat 
  counterintuitive, the reciprocal of the manual weights is taken. Therefore, the higher 
  the manual weight the higher in the influence of the corresponding variable's 
  predictor on the overall distance.
  The donor/recipient ID pairlist for each imputation and identified pattern
  ('block') is stored. In general, \code{weightMatrix}, \code{model} and \code{pairlist}
  are list objects named \sQuote{M1} to \sQuote{M<M>}, and each in return is a list object
  named \sQuote{block1} to \sQuote{block<length(blockImp)>}. \code{model} contains another list object
  with lm-objects for all variables in a particular block.
  Unlike \code{BBPMM} this algorithm is not based on sequential
  regression. Therefore, imputed variables are conditionally independent
  given the completely observed variables (of which at least one must exist).}


\seealso{\code{\link{rowimpPrep}}, \code{\link{BBPMM}}}

\references{
  Little, R.J.A. (1988) \emph{Missing-Data Adjustments in Large
    Surveys}, Journal of Business and Economic Statistics, Vol.6,
  No.3,pp. 287-296.

  Koller-Meinfelder, F. (2009) \emph{Analysis of Incomplete Survey Data
    -- Multiple Imputation Via Bayesian Bootstrap Predictive Mean
    Matching}, doctoral thesis.  
}

\keyword{robust}
\keyword{datagen}

\examples{

### sample data set with non-normal variables and a single
### missingness pattern
set.seed(1000)
n <- 50
x1 <- round(runif(n,0.5,3.5))
x2 <- as.factor(c(rep(1,10),rep(2,25),rep(3,15)))
x3 <- round(rnorm(n,0,3))
y1 <- round(x1-0.25*(x2==2)+0.5*x3+rnorm(n,0,1))
y1 <- ifelse(y1<1,1,y1)
y1 <- ifelse(y1>4,5,y1)
y2 <- y1+rnorm(n,0,0.5)
y3 <- round(x3+rnorm(n,0,2))
data <- as.data.frame(cbind(x1,x2,x3,y1,y2,y3))
misrow1 <- sample(n,20)
data[misrow1, c(4:6)] <- NA

### preparation step
impblock <- rowimpPrep(data)

### imputation
imputed.data <- BBPMM.row(impblock, M=5)
}
