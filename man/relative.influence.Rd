\name{relative.influence}
\alias{relative.influence}
\alias{permutation.test.gbm}
\alias{gbm.loss}
\title{ Methods for estimating relative influence }
\description{
Helper functions for computing the relative influence of each variable in the gbm object.
}
\usage{
relative.influence(object, n.trees)
permutation.test.gbm(object, n.trees)
gbm.loss(y,f,w,offset,dist,baseline)
}
\arguments{
\item{object}{a \code{gbm} object created from an initial call to \code{\link{gbm}}.}
\item{n.trees}{ the number of trees to use for computations.}
\item{y,f,w,offset,dist,baseline}{For \code{gbm.loss}: These components are the
outcome, predicted value, observation weight, offset, distribution, and comparison
loss function, respectively.}
}
\details{
This is not intended for end-user use. These functions offer the different
methods for computing the relative influence in \code{\link{summary.gbm}}.
\code{gbm.loss} is a helper function for \code{permutation.test.gbm}.
}
\value{
Returns an unprocessed vector of estimated relative influences.
}
\references{
J.H. Friedman (2001). "Greedy Function Approximation: A Gradient Boosting
Machine," Annals of Statistics 29(5):1189-1232.

L. Breiman (2001). "Random Forests," Available at \url{ftp://ftp.stat.berkeley.edu/pub/users/breiman/randomforest2001.pdf}.

}
\author{Greg Ridgeway \email{gregr@rand.org}}

\seealso{ \code{\link{summary.gbm}} }

\keyword{ hplot }
