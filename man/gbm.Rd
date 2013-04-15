\name{gbm}
\alias{gbm}
\alias{gbm.more}
\alias{gbm.fit}
\title{Generalized Boosted Regression Modeling}
\description{Fits generalized boosted regression models.}
\usage{
gbm(formula = formula(data),
    distribution = "bernoulli",
    data = list(),
    weights,
    var.monotone = NULL,
    n.trees = 100,
    interaction.depth = 1,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    bag.fraction = 0.5,
    train.fraction = 1.0,
    cv.folds=0,
    keep.data = TRUE,
    verbose = TRUE)

gbm.fit(x,y,
        offset = NULL,
        misc = NULL,
        distribution = "bernoulli",
        w = NULL,
        var.monotone = NULL,
        n.trees = 100,
        interaction.depth = 1,
        n.minobsinnode = 10,
        shrinkage = 0.001,
        bag.fraction = 0.5,
        train.fraction = 1.0,
        keep.data = TRUE,
        verbose = TRUE,
        var.names = NULL,
        response.name = NULL)

gbm.more(object,
         n.new.trees = 100,
         data = NULL,
         weights = NULL,
         offset = NULL,
         verbose = NULL)
}
\arguments{
\item{formula}{a symbolic description of the model to be fit. The formula may 
   include an offset term (e.g. y~offset(n)+x). If \code{keep.data=FALSE} in 
   the initial call to \code{gbm} then it is the user's responsibility to 
   resupply the offset to \code{\link{gbm.more}}.}
\item{distribution}{a character string specifying the name of the distribution 
   to use or a list with a component \code{name} specifying the distribution 
   and any additional parameters needed. 
   Currently available options are "gaussian" (squared error), "laplace"
   (absolute loss), "bernoulli" (logistic regression for 0-1 outcomes), 
   "adaboost" (the AdaBoost exponential loss for 0-1 outcomes), "poisson" 
   (count outcomes), and "coxph" (right censored observations). Quantile 
   regression is also available and \code{distribution} must a list of the form 
   \code{list(name="quantile",alpha=0.25)} where \code{alpha} is the quantile 
   to estimate. The current version's Laplace and quantile regression methods do 
   not handle non-constant weights and will stop.}
\item{data}{an optional data frame containing the variables in the model. By
   default the variables are taken from \code{environment(formula)}, typically 
   the environment from which \code{gbm} is called. If \code{keep.data=TRUE} in 
   the initial call to \code{gbm} then \code{gbm} stores a copy with the 
   object. If \code{keep.data=FALSE} then subsequent calls to 
   \code{\link{gbm.more}} must resupply the same dataset. It becomes the user's 
   responsibility to resupply the same data at this point.}
\item{weights}{an optional vector of weights to be used in the fitting process. 
   Must be positive but do not need to be normalized. If \code{keep.data=FALSE} 
   in the initial call to \code{gbm} then it is the user's responsibility to 
   resupply the weights to \code{\link{gbm.more}}.}
\item{var.monotone}{an optional vector, the same length as the number of
predictors, indicating which variables have a monotone increasing (+1),
decreasing (-1), or arbitrary (0) relationship with the outcome.}
\item{n.trees}{the total number of trees to fit. This is equivalent to the
number of iterations and the number of basis functions in the additive
expansion.}
\item{cv.folds}{Number of cross-validation folds to perform. If \code{cv.folds}>1 then
\code{gbm}, in addition to the usual fit, will perform a cross-validation, calculate
an estimate of generalization error returned in \code{cv.error}.}
\item{interaction.depth}{The maximum depth of variable interactions. 1 implies
an additive model, 2 implies a model with up to 2-way interactions, etc.}
\item{n.minobsinnode}{minimum number of observations in the trees terminal
nodes. Note that this is the actual number of observations not the total
weight.}
\item{shrinkage}{a shrinkage parameter applied to each tree in the expansion.
Also known as the learning rate or step-size reduction.}
\item{bag.fraction}{the fraction of the training set observations randomly
selected to propose the next tree in the expansion. This introduces randomnesses
into the model fit. If \code{bag.fraction}<1 then running the same model twice
will result in similar but different fits. \code{gbm} uses the R random number
generator so \code{set.seed} can ensure that the model can be
reconstructed. Preferably, the user can save the returned
\code{\link{gbm.object}} using \code{\link{save}}.}
\item{train.fraction}{The first \code{train.fraction * nrows(data)}
observations are used to fit the \code{gbm} and the remainder are used for
computing out-of-sample estimates of the loss function.}
\item{keep.data}{a logical variable indicating whether to keep the data and
an index of the data stored with the object. Keeping the data and index makes
subsequent calls to \code{\link{gbm.more}} faster at the cost of storing an
extra copy of the dataset.}
\item{object}{a \code{gbm} object created from an initial call to
\code{\link{gbm}}.}
\item{n.new.trees}{the number of additional trees to add to \code{object}.}
\item{verbose}{If TRUE, gbm will print out progress and performance indicators.
If this option is left unspecified for gbm.more then it uses \code{verbose} from
\code{object}.}

\item{x, y}{For \code{gbm.fit}: \code{x} is a data frame or data matrix containing the
predictor variables and \code{y} is the vector of outcomes. The number of rows
in \code{x} must be the same as the length of \code{y}.}
\item{offset}{a vector of values for the offset}
\item{misc}{For \code{gbm.fit}: \code{misc} is an R object that is simply passed on to
the gbm engine. It can be used for additional data for the specific distribution.
Currently it is only used for passing the censoring indicator for the Cox
proportional hazards model.}
\item{w}{For \code{gbm.fit}: \code{w} is a vector of weights of the same
length as the \code{y}.}
\item{var.names}{For \code{gbm.fit}: A vector of strings of length equal to the
number of columns of \code{x} containing the names of the predictor variables.}
\item{response.name}{For \code{gbm.fit}: A character string label for the response
variable.}

}
\details{

See \code{vignette("gbm")} for technical details of the package. Also available 
at \url{../doc/gbm.pdf} (if you are using HTML help).

This package implements the generalized boosted modeling framework.
Boosting is the process of iteratively adding basis functions in a greedy
fashion so that each additional basis function further reduces the selected
loss function. This implementation closely follows Friedman's Gradient
Boosting Machine (Friedman, 2001).

In addition to many of the features documented in the Gradient Boosting Machine,
\code{gbm} offers additional features including the out-of-bag estimator for
the optimal number of iterations, the ability to store and manipulate the
resulting \code{gbm} object, and a variety of other loss functions that had not
previously had associated boosting algorithms, including the Cox partial
likelihood for censored data, the poisson likelihood for count outcomes, and a
gradient boosting implementation to minimize the AdaBoost exponential loss
function.

\code{gbm.fit} provides the link between R and the C++ gbm engine. \code{gbm}
is a front-end to \code{gbm.fit} that uses the familiar R modeling formulas.
However, \code{\link[stats]{model.frame}} is very slow if there are many
predictor variables. For power-users with many variables use \code{gbm.fit}.
For general practice \code{gbm} is preferable.

}
\value{
\code{gbm}, \code{gbm.fit}, and \code{gbm.more} return a
\code{\link{gbm.object}}.
}
\references{
Y. Freund and R.E. Schapire (1997) \dQuote{A decision-theoretic generalization of
on-line learning and an application to boosting,} \emph{Journal of Computer and
System Sciences,} 55(1):119-139.

G. Ridgeway (1999). \dQuote{The state of boosting,} \emph{Computing Science and
Statistics} 31:172-181.

J.H. Friedman, T. Hastie, R. Tibshirani (2000). \dQuote{Additive Logistic Regression:
a Statistical View of Boosting,} \emph{Annals of Statistics} 28(2):337-374.

J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient Boosting
Machine,} \emph{Annals of Statistics} 29(5):1189-1232.

J.H. Friedman (2002). \dQuote{Stochastic Gradient Boosting,} \emph{Computational Statistics
and Data Analysis} 38(4):367-378.

B. Kriegler (2007). \emph{Cost-Sensitive Stochastic Gradient Boosting Within a 
Quantitative Regression Framework}. PhD dissertation, UCLA Statistics. 
\url{http://theses.stat.ucla.edu/57/KrieglerDissertation.pdf}

\url{http://www.i-pensieri.com/gregr/gbm.shtml}

\url{http://www-stat.stanford.edu/~jhf/R-MART.html}
}

\author{Greg Ridgeway \email{gregr@rand.org}

Quantile regression code developed by Brian Kriegler \email{bk@stat.ucla.edu}}

\seealso{
\code{\link{gbm.object}},
\code{\link{gbm.perf}},
\code{\link{plot.gbm}},
\code{\link{predict.gbm}},
\code{\link{summary.gbm}},
\code{\link{pretty.gbm.tree}}.
}
\examples{
# A least squares regression example
# create some data

N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- ordered(sample(letters[1:4],N,replace=TRUE),levels=letters[4:1])
X4 <- factor(sample(letters[1:6],N,replace=TRUE))
X5 <- factor(sample(letters[1:3],N,replace=TRUE))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]

SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + mu
sigma <- sqrt(var(Y)/SNR)
Y <- Y + rnorm(N,0,sigma)

# introduce some missing values
X1[sample(1:N,size=500)] <- NA
X4[sample(1:N,size=300)] <- NA

data <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

# fit initial model
gbm1 <- gbm(Y~X1+X2+X3+X4+X5+X6,         # formula
    data=data,                   # dataset
    var.monotone=c(0,0,0,0,0,0), # -1: monotone decrease,
                                 # +1: monotone increase,
                                 #  0: no monotone restrictions
    distribution="gaussian",     # bernoulli, adaboost, gaussian,
                                 # poisson, coxph, and quantile available
    n.trees=3000,                # number of trees
    shrinkage=0.005,             # shrinkage or learning rate,
                                 # 0.001 to 0.1 usually work
    interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
    bag.fraction = 0.5,          # subsampling fraction, 0.5 is probably best
    train.fraction = 0.5,        # fraction of data for training,
                                 # first train.fraction*N used for training
    n.minobsinnode = 10,         # minimum total weight needed in each node
    cv.folds = 5,                # do 5-fold cross-validation
    keep.data=TRUE,              # keep a copy of the dataset with the object
    verbose=TRUE)                # print out progress

# check performance using an out-of-bag estimator
# OOB underestimates the optimal number of iterations
best.iter <- gbm.perf(gbm1,method="OOB")
print(best.iter)

# check performance using a 50\% heldout test set
best.iter <- gbm.perf(gbm1,method="test")
print(best.iter)

# check performance using 5-fold cross-validation
best.iter <- gbm.perf(gbm1,method="cv")
print(best.iter)

# plot the performance
# plot variable influence
summary(gbm1,n.trees=1)         # based on the first tree
summary(gbm1,n.trees=best.iter) # based on the estimated best number of trees

# compactly print the first and last trees for curiosity
print(pretty.gbm.tree(gbm1,1))
print(pretty.gbm.tree(gbm1,gbm1$n.trees))

# make some new data
N <- 1000
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- ordered(sample(letters[1:4],N,replace=TRUE))
X4 <- factor(sample(letters[1:6],N,replace=TRUE))
X5 <- factor(sample(letters[1:3],N,replace=TRUE))
X6 <- 3*runif(N)
mu <- c(-1,0,1,2)[as.numeric(X3)]

Y <- X1**1.5 + 2 * (X2**.5) + mu + rnorm(N,0,sigma)

data2 <- data.frame(Y=Y,X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6)

# predict on the new data using "best" number of trees
# f.predict generally will be on the canonical scale (logit,log,etc.)
f.predict <- predict.gbm(gbm1,data2,best.iter)

# least squares error
print(sum((data2$Y-f.predict)^2))

# create marginal plots
# plot variable X1,X2,X3 after "best" iterations
par(mfrow=c(1,3))
plot.gbm(gbm1,1,best.iter)
plot.gbm(gbm1,2,best.iter)
plot.gbm(gbm1,3,best.iter)
par(mfrow=c(1,1))
# contour plot of variables 1 and 2 after "best" iterations
plot.gbm(gbm1,1:2,best.iter)
# lattice plot of variables 2 and 3
plot.gbm(gbm1,2:3,best.iter)
# lattice plot of variables 3 and 4
plot.gbm(gbm1,3:4,best.iter)

# 3-way plots
plot.gbm(gbm1,c(1,2,6),best.iter,cont=20)
plot.gbm(gbm1,1:3,best.iter)
plot.gbm(gbm1,2:4,best.iter)
plot.gbm(gbm1,3:5,best.iter)

# do another 100 iterations
gbm2 <- gbm.more(gbm1,100,
                 verbose=FALSE) # stop printing detailed progress
}
\keyword{models}
\keyword{nonlinear}
\keyword{survival}
\keyword{nonparametric}
\keyword{tree}
