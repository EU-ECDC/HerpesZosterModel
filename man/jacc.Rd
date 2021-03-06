% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures.R
\name{jacc}
\alias{jacc}
\title{Jaccard index}
\usage{
jacc(tp, tn, fp, fn, ...)
}
\arguments{
\item{tp}{True positives}

\item{tn}{True negatives}

\item{fp}{False positives}

\item{fn}{False negatives}
}
\description{
The Jaccard index (also known as Tanimoto similarity) is given by
\deqn{Jaccard = TP / (TP + FP + FN)}
where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
\eqn{FP} denotes false positives, and \eqn{FN} denotes false negatives.
}
\examples{
jacc(tp = 45, fp = 15, fn = 25, tn = 15)
}
\references{
\insertRef{Jaccard1908}{zostmod}
\insertRef{Tversky1977}{zostmod}
}
\seealso{
Other classification scores: \code{\link{auc}},
  \code{\link{bcr}}, \code{\link{brier}},
  \code{\link{dor}}, \code{\link{dp}}, \code{\link{err}},
  \code{\link{et}}, \code{\link{f1}}, \code{\link{fai}},
  \code{\link{fm}}, \code{\link{gain}}, \code{\link{gl}},
  \code{\link{gm}}, \code{\link{ignr}},
  \code{\link{kappa}}, \code{\link{lr}}, \code{\link{mcc}},
  \code{\link{op}}, \code{\link{rand}}, \code{\link{rt}},
  \code{\link{rus}}, \code{\link{sm}}, \code{\link{ss}},
  \code{\link{tss}}
}
\concept{classification scores}
\keyword{classif}
