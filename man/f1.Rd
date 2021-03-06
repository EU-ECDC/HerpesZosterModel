% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measures.R
\name{f1}
\alias{f1}
\title{F score}
\usage{
f1(tp, tn, fp, fn, beta = 1, ...)
}
\arguments{
\item{tp}{True positives}

\item{tn}{True negatives}

\item{fp}{False positives}

\item{fn}{False negatives}

\item{beta}{The weighting parameter, defaults to 1}
}
\description{
The F score is given by
\deqn{F = ((1 + \beta^2) \cdot TP) / ((1 + \beta^2) \cdot TP + \beta^2 \cdot FP + FN)}
where \eqn{TP} denotes true positives, \eqn{TN} denotes true negatives,
\eqn{FP} denotes false positives, \eqn{FN} denotes false negatives, and
\eqn{\beta} is a weighting parameter. In our implementation \eqn{\beta}
has a default value of \eqn{1}, meaning our F score defaults to the F1 score
(hence the function name).
}
\examples{
f1(tp = 45, fp = 15, fn = 25, tn = 15)
}
\references{
\insertRef{Soerensen1948}{zostmod}
\insertRef{Tversky1977}{zostmod}
}
\seealso{
\code{\link[=mcc]{mcc()}}

Other classification scores: \code{\link{auc}},
  \code{\link{bcr}}, \code{\link{brier}},
  \code{\link{dor}}, \code{\link{dp}}, \code{\link{err}},
  \code{\link{et}}, \code{\link{fai}}, \code{\link{fm}},
  \code{\link{gain}}, \code{\link{gl}}, \code{\link{gm}},
  \code{\link{ignr}}, \code{\link{jacc}},
  \code{\link{kappa}}, \code{\link{lr}}, \code{\link{mcc}},
  \code{\link{op}}, \code{\link{rand}}, \code{\link{rt}},
  \code{\link{rus}}, \code{\link{sm}}, \code{\link{ss}},
  \code{\link{tss}}
}
\concept{classification scores}
\keyword{classif}
