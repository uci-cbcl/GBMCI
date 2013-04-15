\name{gbm-package}
\alias{gbm-package}
\docType{package}
\title{
Generalized Boosted Regression Models
}
\description{
This package implements extensions to Freund and
Schapire's AdaBoost algorithm and J. Friedman's gradient
boosting machine. Includes regression methods for least
squares, absolute loss, logistic, Poisson, Cox proportional
hazards partial likelihood, and AdaBoost exponential loss.
}
\details{
\tabular{ll}{
Package: \tab gbm\cr
Version: \tab 1.5-6\cr
Date: \tab 2006-1-20\cr
Depends: \tab R (>= 2.1.0), survival, lattice, mgcv\cr
License: \tab GPL (version 2 or newer)\cr
URL: \tab http://www.i-pensieri.com/gregr/gbm.shtml\cr
Built: \tab R 2.2.1; i386-pc-mingw32; 2006-02-24 18:09:42; windows\cr
}

Index:
\preformatted{
basehaz.gbm             Baseline hazard function
calibrate.plot          Calibration plot
gbm                     Generalized Boosted Regression Modeling
gbm.object              Generalized Boosted Regression Model Object
gbm.perf                GBM performance
plot.gbm                Marginal plots of fitted gbm objects
predict.gbm             Predict method for GBM Model Fits
pretty.gbm.tree         Print gbm tree components
quantile.rug            Quantile rug plot
relative.influence      Methods for estimating relative influence
shrink.gbm              L1 shrinkage of the predictor variables in a
                        GBM
shrink.gbm.pred         Predictions from a shrunked GBM
summary.gbm             Summary of a gbm object
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{gbm} \tab Generalized Boosted Models: A guide to the gbm package (source, pdf)\cr
}
}
\author{
Greg Ridgeway \email{gregr@rand.org}
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

\url{http://www.i-pensieri.com/gregr/gbm.shtml}

\url{http://www-stat.stanford.edu/~jhf/R-MART.html}}
\keyword{ package }
